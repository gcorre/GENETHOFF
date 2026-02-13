# ========================================================= #
# Guillaume CORRE @ GENETHON 2025
# This script is used to generate HTML tables and figures for GUIDE-seq data analysis.
# It is used in the snakemake pipeline for GUIDE-seq analysis.
# ========================================================= #

# Set options for R environment
options(tidyverse.quiet = TRUE,warn = -1,verbose = F,warn = -10,conflicts.policy = list(warn = FALSE))
options(knitr.kable.NA = '')


# Load necessary libraries
library(tidyverse,quietly = T,warn.conflicts = F,verbose = F)
library(readxl,quietly = T,warn.conflicts = F,verbose = F)
library(kableExtra,quietly = T,warn.conflicts = F,verbose = F)
library(GenomicRanges,quietly = T,warn.conflicts = F,verbose = F)
library(yaml,quietly = T,warn.conflicts = F,verbose = F)
library(rmdformats,quietly = T,warn.conflicts = F,verbose = F)
library(ggrepel,quietly = T,warn.conflicts = F,verbose = F)
library(data.table,quietly = T,warn.conflicts = F,verbose = F)
library(DT, quietly = T,warn.conflicts = F,verbose = F)



# debug

args <-  c("05-Report//Lib_1_CD34_ODN_neg_summary_grna.csv 05-Report/Lib_1_K562_iODN_neg_summary_grna.csv 05-Report/Lib_1_K562_ODN_neg_summary_grna.csv 05-Report//Lib_2_CD34_ODN_pos_summary_grna.csv 05-Report//Lib_3_CD34_ODN_pos_summary_grna.csv 05-Report//Lib_3_K562_iODN_pos_summary_grna.csv 05-Report//Lib_3_K562_ODN_pos_summary_grna.csv",
           "sample_info.xlsx",
           "configuration.yml",
           1, 
           100,
           "results/test.xlsx")




args <- commandArgs(trailingOnly = T)

# 1 : csv files from all libraries
# 2 : sample info file
# 3 : configuration file
# 4 : min number of UMI for a cluster to be reported in html files
# 5 : Number of 



summary_files = args[1] 


if(str_ends(args[2],pattern = "xlsx")){
  sampleInfo <- read_xlsx(args[2])
} else {
  sampleInfo <- read.delim(args[2], sep=";")
}

config <- read_yaml(args[3])
minUMI_alignments_figure <- as.numeric(args[4])
max_clusters <- as.numeric(args[5])



#----------------------------------------------------------#
# 01 - get statistics of libraries ----
#----------------------------------------------------------#

files <- list.files("05-Report/", pattern = "stat$", full.names = T)
names(files) <- str_remove(basename(files),"\\.stat")



stats <- lapply(files, read.delim)
stats <- stats %>% 
  bind_rows() %>% 
  distinct() %>% 
  select(file,num_seqs) %>% 
  mutate(library = str_match(file,"/(.+)_R1")[,2]) %>% 
  mutate(step = case_when(str_detect(file, "_R1.fastq.gz")~"Demultiplexed",
                          str_ends(file, "R1.ODN.UMI.fastq.gz")~"ODN_match",
                          str_ends(file, "R1.UMI.ODN.trimmed.filtered.fastq.gz")~"Filtered",
                          TRUE ~ NA)) %>% 
  filter(!is.na(step)) %>% 
  select(-file) %>% 
  group_by(library) %>% 
  mutate(prop = round(num_seqs/dplyr::first(num_seqs)*100,digits = 2),
         value = paste(format(num_seqs,big.mark=",")," (",prop,"%)",sep="")) %>% 
  pivot_wider(id_cols = "library",names_from = "step", values_from = "value")


#----------------------------------------------------------#
# 02 - Load summary files ----
#----------------------------------------------------------#

summary_files <- unlist(str_split(summary_files," "))

summary <- lapply(summary_files, function(x){
  
  if(file.size(x)>0){
    
   read.delim(x, sep="\t", header = T)

  }
})


names(summary) <- str_remove(basename(summary_files),"\\_summary_grna.csv")

libraries_count <- length(summary)


# --> convert columns type

summary <- lapply(summary, function(x){
  
  x %>% mutate(., across(everything(), as.character))
  
})





#----------------------------------------------------------#
# 03 - Aggregate libraries ----
#----------------------------------------------------------#

summary_df <- bind_rows(summary, .id = "library") %>% 
  type.convert(as.is = TRUE)


#----------------------------------------------------------#
# 03 - Flag the most probable ON-target ----
#----------------------------------------------------------#
# based on smallest number EDITS in crRNA & PAM

# get cluster with best gRNa alignment
OT_candidate <- summary_df %>% 
  filter(!is.na(Alignment)) %>% 
  arrange(desc(N_UMI_cluster)) %>% 
  group_by(library) %>% 
  mutate(Rank = dense_rank(x = -N_UMI_cluster)) %>% 
  slice_min(n = 1,with_ties = T, order_by = N_edits) %>% 
  slice_min(n = 1,with_ties = T, order_by = PAM_indel_count) %>% 
  select(clusterID,Rank) %>% 
  mutate(best="TRUE")

# calculate relative abundance of gRNA matches clusters

OT_candidate <- summary_df %>% 
  filter(!is.na(Alignment)) %>% 
  left_join(OT_candidate, by = c("library","clusterID")) %>% 
  group_by(library) %>% 
  mutate(Relative_abundance = round(N_UMI_cluster / sum(N_UMI_cluster) *100,digits = 2)) %>% 
  select(clusterID, best,Relative_abundance, Rank )

# add relative abundance and flag to table
summary_df <- summary_df %>% 
  left_join(OT_candidate, by = c("library","clusterID"))


## get best match information by library:
best_aligns <- summary_df %>%  
  filter(best=="TRUE") %>%
  arrange(library) %>% 
  select(library,
         clusterID,
         chromosome ,
         cut_modal_position,
         cut_gRNa_alignment,
         Edits_gRNA = N_edits,
         Edits_PAM=PAM_indel_count,
         N_UMI_cluster,
         Relative_abundance, Rank)

# In case one library has no "best" alignment, add empty informations
best_aligns <- best_aligns %>% 
  full_join(data.frame(library = names(summary)))



#----------------------------------------------------------#
# 04 - Get stats per library ----
#----------------------------------------------------------#

stats_summary <- summary_df %>% 
  group_by(library) %>% 
  summarise(
    Reads = sum(N_reads_cluster),
    UMIs = sum(N_UMI_cluster),
    Insertions = sum(N_IS_cluster),
    Clusters = n(),
    "With gRNA match .."=  length(which(!is.na(Alignment))),
    ".. 2 PCR orientations" = length(which(!is.na(Alignment) & N_orientations_PCR==2)),
    ".. 2 ODN orientations" = length(which(!is.na(Alignment) & N_orientations_cluster==2)),
    ".. In Oncogene" = length(which(!is.na(Alignment) & str_detect(Onco_annotation,"[A-Za-z]") & (Onco_annotation != "Not_Done")))) %>% 
  pivot_longer(cols =starts_with(".."),names_to = " And ..", values_to = "count") 


#----------------------------------------------------------#
# 05 - Make figures ----
#----------------------------------------------------------#

## rank-abundance curve

RankAbundance_data <- summary_df %>% 
  filter(N_UMI_cluster>=minUMI_alignments_figure, !is.na(Alignment)) %>% 
  group_by(library) %>% 
  mutate(rank_desc = row_number(dplyr::desc(N_UMI_cluster))) %>% 
  ungroup() %>% 
  arrange(rank_desc)



fig_RankAbundance <- ggplot(RankAbundance_data, aes(rank_desc,N_UMI_cluster)) +
  geom_step(col = "black", direction = "hv") +
  geom_point(pch=19)+
  facet_wrap(~library, ncol = 3, scales = "free") + 
  theme_bw(base_size = 12)+
  scale_y_log10() + 
  scale_x_log10() + 
  ggrepel::geom_text_repel(data = . %>% 
                             group_by(library) %>% 
                             slice_min(n=3,order_by = rank_desc), 
                           aes(rank_desc,N_UMI_cluster,label=paste(chromosome,cut_gRNa_alignment,sep=":")),
                           inherit.aes = F,col = "purple",cex=3,nudge_x = 1,nudge_y = 0.5,force = 5,direction= "x")+
  geom_point(data = . %>% 
               filter(N_edits==0,PAM_indel_count==0), 
             aes(rank_desc,N_UMI_cluster),col = "red",cex= 2,inherit.aes = F)+
  labs(x = "Ranked clusters by decreasing abundance",
       y = "UMI count per cluster",
       caption = paste("Only clusters with >=",minUMI_alignments_figure, "UMIs",sep=" "))+
  scale_color_manual(values = c("black","green3"))






## Chromosome distribution

ChromDistr_data <-  summary_df %>% 
  filter(!is.na(Alignment) & (str_starts(chromosome,"[0-9XYMT]") | str_starts(chromosome,"chr"))) %>% 
  group_by(library,chromosome,Predicted=predicted) %>% 
  summarise(clusters = n(), reads = sum(N_reads_cluster), UMI = sum(N_UMI_cluster)) %>% 
  mutate(chromosome = recode(chromosome, "MT" = "M")) 

if(!any(str_starts(ChromDistr_data$chromosome,"chr"))){
  ChromDistr_data <- ChromDistr_data %>%
    mutate(chromosome = paste("chr",chromosome,sep=""))
} 

ChromDistr_data <- ChromDistr_data %>% 
  mutate(chromosome = factor(chromosome, levels = paste("chr",c(1:22,"X","Y","MT"),sep="")),
         Predicted = factor(Predicted, levels = c("yes","no")))


# figures
fig_ChromDistr_clusters <- ggplot(ChromDistr_data, aes( chromosome, clusters,fill = Predicted)) +
  geom_col(col = "black",na.rm = F) +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  coord_polar() +
  theme_bw(base_size = 12)+
  scale_fill_manual(values = c("green4","grey"),drop = F) +
  scale_x_discrete(drop=F)+
  theme(panel.grid.major.x = element_line(linetype = 2,colour = "grey40"),axis.text =  element_text(color = "black"))


fig_ChromDistr_UMI <- ggplot(ChromDistr_data, aes( chromosome, UMI,fill = Predicted)) +
  geom_col(col = "black") +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  coord_polar() +
  theme_bw(base_size = 12)+
  scale_fill_manual(values = c("green4","grey"),drop = F) +
  scale_x_discrete(drop=F)+
  theme(panel.grid.major.x = element_line(linetype = 2,colour = "grey40"),axis.text =  element_text(color = "black"))



# Distribution of UMI around cut best cut site

files <- list.files("04-IScalling/", pattern="UMIs_per_IS_in_Cluster.bed", full.names = T)
names(files) <- str_remove(basename(files),".UMIs_per_IS_in_Cluster.bed")

IS <- lapply(files, read.delim,header=F)
IS <- lapply(IS, function(x){
  x %>% mutate(V1 = as.character(V1))
})


IS <- IS %>% bind_rows(.id="library")

IS <- IS %>% inner_join(best_aligns, by = c("library","V12"="clusterID"))

IS <- IS %>% mutate(rel_dist_gRNA = V2 - cut_gRNa_alignment,
                    rel_dist_mod = V2 - as.numeric(cut_modal_position))

IS_stat <- IS %>%
  group_by(library,V6,chromosome,cut_gRNa_alignment ,rel_dist_gRNA,clusterID=V12 ) %>%
  summarise(count=sum(V8),
            UMI = log10(sum(V8)+1)) %>%
  mutate(UMI = case_when(V6 == "+"~UMI,
                         TRUE ~ -UMI))



fig_distrAroundBestSite <- ggplot(IS_stat, aes(rel_dist_gRNA,UMI,fill = factor(sign(UMI)))) +
  geom_vline(xintercept = 0, lty = 2,col = "black") +
  geom_col(show.legend = F) + 
  facet_wrap(~library+paste(chromosome,cut_gRNa_alignment,sep=":"),scales= "free_y",ncol = 3) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  labs(x = "Distance to cut site (bp)", y = "log10(UMI)")+
  scale_fill_manual(values = c("blue3","red3"))





##  Define function

compare_strings <- function(DNA,RNA){
  
  DNA <- str_split_1(DNA,"")
  RNA <- str_split_1(RNA,"")
  
  DNA[DNA==RNA]<- "."
  return(paste(DNA,collapse=""))
  
}



#----------------------------------------------------------#
# 06 - HTML tables for off targets  ----
#----------------------------------------------------------#

summary_pred_bulge <- split(summary_df, summary_df$library)

# --> save libraries as an excel file with one tab per library

writexl::write_xlsx(summary_pred_bulge, path = args[6] ,col_names = T, format_headers = T)




tables_off <- lapply(seq_along(summary_pred_bulge),function(x){
  
  
  
  dt <- summary_pred_bulge[[x]] %>% 
    filter(!is.na(Alignment)) 
  
  if(nrow(dt)>0){
    
    cat("Generating html table for sample : ",names(summary_pred_bulge)[x],"\n")
  
    
    dt <- dt %>% 
      mutate("UMIs (%)" = round(N_UMI_cluster / sum(N_UMI_cluster)*100,digits = 1),.after=N_UMI_cluster) %>% 
      filter(N_UMI_cluster > minUMI_alignments_figure) %>% 
      slice_max(n = max_clusters,order_by = N_UMI_cluster,with_ties = F) %>% 
      as.data.table()
    
    
    dt  <- dt %>% 
      rowwise() %>% 
      mutate(seq_gDNA_plot = compare_strings(seq_gDNA,seq_gRNA),
             PAM_plot = compare_strings(pam_gDNA,pam_gRNA)) %>%
      ungroup %>%
      data.table()
    
    
    
    ### FORMAT alignments to HTML to add some colors
    ## use a monospace font to keep equal character width
    
    dt[, seq_gRNA_html := paste(text_spec(background_as_tile = FALSE, monospace = TRUE,
                                          strsplit(seq_gRNA, split = "")[[1]],
                                          background = recode(strsplit(seq_gRNA, split = "")[[1]],
                                                              A = "#129749", T = "#d62839", C = "#255c99", G = "#f7b32b")),
                                collapse = ""), by = 1:nrow(dt)]
    
    dt[, seq_gDNA_html := paste(text_spec(background_as_tile = FALSE, monospace = TRUE,
                                          strsplit(seq_gDNA_plot, split = "")[[1]],
                                          background = recode(strsplit(seq_gDNA_plot, split = "")[[1]],
                                                              A = "#129749", T = "#d62839", C = "#255c99", G = "#f7b32b", "-" = "grey70")),
                                collapse = ""), by = 1:nrow(dt)]
    
    dt[, pam_gRNA_html := paste(text_spec(background_as_tile = FALSE, monospace = TRUE,
                                          strsplit(pam_gRNA, split = "")[[1]],
                                          background = recode(strsplit(pam_gRNA, split = "")[[1]],
                                                              N = "grey", A = "#129749", T = "#d62839", C = "#255c99", G = "#f7b32b")),
                                collapse = ""), by = 1:nrow(dt)]
    
    dt[, pam_gDNA_html := paste(text_spec(background_as_tile = FALSE, monospace = TRUE,
                                          strsplit(PAM_plot, split = "")[[1]],
                                          background = recode(strsplit(PAM_plot, split = "")[[1]],
                                                              A = "#129749", T = "#d62839", C = "#255c99", G = "#f7b32b")),
                                collapse = ""), by = 1:nrow(dt)]
    
    
    if(unique(dt$pam_side) == "3"){
      dt[, alignment_html := paste("gRNA: ", seq_gRNA_html, " ", pam_gRNA_html, " <br>gDNA: ", seq_gDNA_html, " ", pam_gDNA_html, sep = "")]
    } else {
      dt[, alignment_html := paste("gRNA: ",pam_gRNA_html," ", seq_gRNA_html, " <br>gDNA: ", pam_gDNA_html," ",seq_gDNA_html,   sep = "")]
    }
    
    convert_to_text_spec <- function(gene_symbols) {
      if(!is.na(gene_symbols)){
        base_url <- "http://www.genecards.org/cgi-bin/carddisp.pl?gene="
        symbols <- unlist(strsplit(gene_symbols, ","))
        text_spec_list <- lapply(symbols, function(symbol) {
          text_spec(symbol, link = paste0(base_url, str_trim(symbol,side = "both")),new_tab = TRUE)
        })
        # do.call(c, text_spec_list)
        toString(text_spec_list)} 
      else {
        ""
      }
      
    }
    
    dt[,gene_links := sapply(Symbol, convert_to_text_spec, simplify = TRUE)]
    dt[, Symbol_html := str_replace_all(gene_links, ", ", " <br>")]
    
    dt[, position_html := str_replace_all(position, ", ", " <br>")]
    dt[, Oncogene_html := str_replace_all(Onco_annotation, ",", " <br>")]
    
    
    dt$predicted_alignment_html <- cell_spec(dt$predicted, background = ifelse(dt$predicted == "yes", "#129749", "white"))
    
    
    
    
    ## which columns to add to the alignment html file
    dt <- dt %>% 
      mutate("cut offset" = relative_distance,
             PAM_indel_count = as.numeric(PAM_indel_count),
             PCR = case_when(N_UMI_positive>0 & N_UMI_negative > 0 ~ "+/-",
                             N_UMI_positive>0 ~"+",
                             N_UMI_negative>0 ~"-",
                             TRUE ~ NA)) %>% 
      select(chromosome,
             cut_gRNa_alignment, 
             `cut offset`,
             Alignment=alignment_html,
             UMIs=N_UMI_cluster,
             `UMIs (%)`,
             Reads=N_reads_cluster,
             "Edits crRNA"=N_edits, 
             "Edits pam"= PAM_indel_count,
             Symbol=Symbol_html,
             Feature=position_html,
             Predicted= predicted_alignment_html,
             PCR,
             Oncogene = Oncogene_html) %>% 
      unite(col = "Position",chromosome,cut_gRNa_alignment,sep = ":") 
    
    
    
    ## Make kable static table
    
    kb <- kbl(dt,
              escape = F,
              align=c(rep("r",6),rep('c', 7))) %>%
      kable_classic_2(full_width = F,html_font = "helvetica") %>%
      kable_styling(bootstrap_options = c("condensed","hover","stripped"),
                    font_size = 12,
                    fixed_thead = T) %>%
      column_spec(1:(ncol(dt)),extra_css = "vertical-align:middle;")
    
    save_kable(kb,file = paste("results/report-files/",names(summary_pred_bulge)[x],"_offtargets.html",sep=""),self_contained=T)
    
    
    # Make DT dynamic table
    
    dt2 <- DT::datatable(dt,
                         escape = F, 
                         extensions = 'Buttons',
                         rownames = FALSE,
                         filter = 'top',
                         class = 'nowrap display compact',
                         options = list(
                           pageLength = 20,
                           autoWidth = F,
                           lengthMenu = list(c(20, 50, -1), c('20', '50', 'All')),
                           fixedHeader = TRUE,
                           dom = 'Blcfrtip',
                           buttons = c('copy', 'csv', 'excel'),
                           columnDefs = list(list(className = 'dt-body-center', targets = c(6,7,8,9,10)),
                                             list(className = 'dt-body-right', targets = c(0,1,2,3,4,5)),
                                             list(className = 'dt-head-center', targets = seq(1:(ncol(dt)-1))),
                                             list(width = "2%", targets = 5)),
                           initComplete = DT::JS(
                             "function(settings, json) {",
                             "$('body').css({'font-family': 'Calibri', 'font-size': '12px'});",
                             "$('table').css({'width': '100%'});
                            }"
                           )
                         )
    ) %>% 
      formatStyle(
        'UMIs',valueColumns = "UMIs",
        background = styleColorBar(data = dt$UMIs, 'steelblue'),
        backgroundSize = '95% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )%>% 
      formatStyle(
        'UMIs (%)',
        background = styleColorBar(data = dt$`UMIs (%)`,'steelblue'),
        backgroundSize = '95% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )%>%
      formatStyle(
        'Reads',
        background = styleColorBar(dt$Reads, 'orange'),
        backgroundSize = '95% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
    
    
    
    htmlwidgets::saveWidget(widget = dt2, file = paste("results/report-files/",names(summary_pred_bulge)[x],"_offtargets_dynamic.html",sep=""),selfcontained = T)
    
    return(paste("results/report-files/",names(summary_pred_bulge)[x],"_offtargets_dynamic.html",sep="")) # name the output in list
  } else {  
    return("")
  }
}
)


names(tables_off) <- names(summary_pred_bulge)



# Save report data ------------------------------------------------------------------

save(list = c("tables_off", "summary_pred_bulge","sampleInfo", "stats", "stats_summary","best_aligns", grep(x=ls(),"fig",value = T)), file = "results/report.rdata")







