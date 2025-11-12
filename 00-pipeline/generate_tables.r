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
  #library(UpSetR, quietly = T,warn.conflicts = F,verbose = F)  
  


# get command line arguments 
  args <- commandArgs(trailingOnly = T)


# Assign command-line arguments to variables
  summary_files = args[1] 
  
  if(str_ends(args[2],pattern = "xlsx")){
    sampleInfo <- read_xlsx(args[2])
  } else {
    sampleInfo <- read.delim(args[2], sep=";")
  }

  config <- read_yaml(args[3])
  predicted_files = args[4] 
  max_clusters <- as.numeric(args[5])
  minUMI_alignments_figure <- as.numeric(args[6])
  min_predicted_distance <- as.numeric(args[7])

## debug
# 
#   summary_files= "results/GUIDEseq_Sp_NT_summary.xlsx results/GUIDEseq_SpRY_gRNA_3_HBG_summary.xlsx results/GUIDEseq_SpRY_gRNA_21_HBG_summary.xlsx results/GUIDEseq_Sp_gRNA_12_HBG_summary.xlsx results/GUIDEseq_SpRY_NT_summary.xlsx"
# 
#  sampleInfo <- read_xlsx("sampleInfo.xlsx")
#  
# #= read.delim("sampleInfo.csv",sep=";")
# config = read_yaml("guideSeq_GNT.yml")
# predicted_files ="06-offPredict/GRCh38_CCTTCCCCACACTATCTCAA_NNN_3.csv 06-offPredict/GRCh38_GCCCCTTCCCCACACTATCT_NNN_3.csv 06-offPredict/GRCh38_GTGGGGAAGGGGCCCCCAAG_NGG_3.csv"
# max_clusters=100
# minUMI_alignments_figure=1
# min_predicted_distance=100

  

# Load demux stat -------------------------------------------------------------------------

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



# Load summary excel file -------------------------------------------------------
  summary_files <- unlist(str_split(summary_files," "))
  summary <- lapply(summary_files,read_excel)
  
  summary <- lapply(summary, function(x){

    x %>%
      mutate(chromosome = case_when(!str_starts(chromosome ,"chr") ~ paste("chr",chromosome ,sep=""), 
                           TRUE ~ chromosome))
  })
  
  names(summary) <- str_remove(basename(summary_files),"\\_summary.xlsx")
  
  libraries_count <- length(summary)


# Get probable ON-targets -----------------------------------------------------------
# based on smallest EDITS in crRNA & PAM
  
  summary <- lapply(summary,function(x){
    
    best_df <- x %>%
      filter(!is.na(Alignment)) %>% ## remove alignments without gRNA match
      arrange(desc(N_UMI_cluster)) %>% 
      mutate(Rank = dense_rank(x = -N_UMI_cluster)) %>% 
      slice_min(n = 1,with_ties = T, order_by = N_edits) %>% 
      slice_min(n = 1,with_ties = T, order_by = PAM_indel_count) %>% 
      #filter(N_edits==0) %>% 
      select(clusterID,Rank) %>% 
      mutate(best=TRUE)
    
    abundance <- x  %>% 
      filter(!is.na(Alignment)) %>% ## remove alignments without gRNA match
      left_join(best_df,by = c("clusterID")) %>%
      mutate(Relative_abundance = round(N_UMI_cluster / sum(N_UMI_cluster) *100,digits = 2)) %>% 
      select(clusterID, best,Relative_abundance, Rank )
    
    x <-x %>% left_join(abundance) 
    return(x)
    
  })

## get best match information by library:
  best_aligns <- summary %>%  
    bind_rows(.id="library") %>% 
    filter(best==T) %>%
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







# Get stats -------------------------------------------------------------------------

stats_summary <- lapply(summary, function(x){
  
    x %>% summarise(
    Reads = sum(N_reads_cluster),
    UMIs = sum(N_UMI_cluster),
    Insertions = sum(N_IS_cluster),
    Clusters = n(),
    "With gRNA match .."= length(which(!is.na(Alignment))),
    ".. 2 PCR orientations" = length(which(!is.na(Alignment) & N_orientations_PCR==2)),
    ".. 2 ODN orientations" = length(which(!is.na(Alignment) & N_orientations_cluster==2)),
    #"& with >3 sites" = length(which(!is.na(Alignment) & N_IS_cluster>3)),
    ".. In Oncogene" = length(which(!is.na(Alignment) & str_detect(Onco_annotation,"[A-Za-z]")))) %>% 
    pivot_longer(cols =starts_with(".."),names_to = " And ..", values_to = "count") 
}
)

stats_summary <- stats_summary %>% bind_rows(.id="library")



# Annotate predictions ------------------------------------------------------------------
cat("Reading predicted cutting sites\n")
predicted_files <- unlist(str_split(predicted_files," "))
predict_gRNA <- lapply(predicted_files,data.table::fread  )
names(predict_gRNA) <- str_remove(basename(predicted_files),"\\.csv")



# Perform analysis for each library

summary_pred_bulge <- lapply(seq_along(summary), function(x){
  lib_name <- names(summary)[x]
  cat("processing",lib_name,"\n")
  
  
  ## take summary excel file
  df <- summary[[x]] 
  
  # remove clusters without gRNA match
  df_grna <- df %>% 
    dplyr::filter(!is.na(Alignment))
  
  if(nrow(df_grna)>0){
    
    ### load prediction for corresponding gRNA/genome/PAM
    if(nrow(df_grna)>0){
      grna_name <- sampleInfo %>% 
        filter(sampleName == names(summary)[x]) %>% 
        mutate("grna_name" = paste(Genome,"_",gRNA_sequence,"_",PAM_sequence,"_",PAM_side,sep="")) %>%
        distinct(grna_name) %>% 
        pull(grna_name)
      
      if(length(grna_name)==1){
        cat("Using .....", grna_name,"\n")
        
        if(nrow(predict_gRNA[[grna_name]])>0){
          
          predict <- predict_gRNA[[grna_name]]  %>% 
            mutate(Chromosome = case_when(!str_starts(Chromosome,"chr") ~ paste("chr",str_remove_all(Chromosome,"_.+$"),sep=""), 
                                          TRUE ~ Chromosome))
          
          
          
          ## intersect cluster range and prediction ranges to find predicted clusters
          gr_predict <- makeGRangesFromDataFrame(predict,
                                                 start.field = "Position",
                                                 end.field = "Position",
                                                 seqnames.field = "chromosome",
                                                 keep.extra.columns = F)
          
          gr_data <- makeGRangesFromDataFrame(df_grna, 
                                              start.field = "cut_gRNa_alignment",
                                              end.field = "cut_gRNa_alignment",
                                              seqnames.field = "chromosome",
                                              keep.extra.columns = F)
          
          # distance to nearest predicted cut site
          near_df <- distanceToNearest(gr_data,gr_predict)
          
          # filter predicted that are within window defined in configuration file
          near_df <- near_df[near_df@elementMetadata$distance <= min_predicted_distance,]
          
          # clear
          rm(gr_predict);rm(gr_data)
          
          # get identified clusters with gRNA
          x3 <- df_grna %>% 
            unite("DNA", seq_gDNA,pam_gDNA,sep="",remove = F ) %>% 
            unite("crRNA",seq_gRNA,pam_gRNA,sep="",remove=F) %>% 
            distinct(chromosome, clusterID,DNA,crRNA)
          
          
          # keep those predicted
          x3 <- data.frame(x3[near_df@from,])
          
          # Add predicted alignments
          pred <- predict[near_df@to,] %>% data.frame()
          
          pred <- pred %>%
            select(chromosome_pred=Chromosome,DNA_pred=AlignedText,crRNA_pred=AlignedTarget)
          
          # merge identified clusters and the predicted alignments
          x3 <- x3 %>% bind_cols(pred)
          
          # annotate clusters if predicted by position (in window) or exact alignment pattern
          x3 <- x3 %>% 
            mutate(predicted_position = "yes",
                   predicted_pattern = case_when(DNA == DNA_pred & crRNA == crRNA_pred ~ "yes",
                                                 TRUE ~"no"))
          
          
          
          ## add prediction status to table of identified OT with gRNA match
          x3 <- df_grna %>% 
            left_join(x3, by = "clusterID") %>% 
            replace_na(list(predicted_position="no",predicted_pattern="no"))
          
          
          ## define if there are bulges in the alignments --> presence of indels
          x3 <- x3 %>% mutate(bulge = case_when(str_detect(seq_gDNA,"-") & str_detect(seq_gRNA,"-") ~ "both",
                                                str_detect(seq_gRNA,"-") ~ "gDNA",
                                                str_detect(seq_gDNA,"-") ~ "gRNA",
                                                TRUE ~ "none")) %>%
            select(clusterID,predicted=predicted_position,bulge)
          
          
          # annotate table of clusters with prediction / bulges
          df <- df %>% left_join(x3, by = "clusterID")
          
          # update xls file with prediction annotation
          df_list <- list(df)
          names(df_list) <- lib_name
          writexl::write_xlsx(df_list,path = summary_files[x],col_names = T,format_headers = T)
          
          return(df)
          
          
        }
      }
    }
    
  }
  
}
)

names(summary_pred_bulge) <- names(summary)

# Make figures ----------------------------------------------------------------------

## rank-abundance curve ----

RankAbundance_data <- lapply(summary_pred_bulge, function(x){
  
  x %>%
    filter(N_UMI_cluster>=minUMI_alignments_figure, !is.na(Alignment))
  }
)

RankAbundance_data <- RankAbundance_data %>%
  bind_rows(.id="library") %>% 
  group_by(library) %>% 
  mutate(rank_desc = row_number(dplyr::desc(N_UMI_cluster)))



fig_RankAbundance <- ggplot(RankAbundance_data, aes(rank_desc,N_UMI_cluster)) +
  geom_step(col = "black", direction = "hv") +
  geom_point(pch=19)+
  facet_wrap(~library, ncol = 3, scales = "free") + 
  theme_bw(base_size = 12)+
  scale_y_log10() + 
  scale_x_log10() + 
  ggrepel::geom_text_repel(data = . %>% 
                             group_by(library) %>% 
                             slice_head(n=3), 
                           aes(rank_desc,N_UMI_cluster,label=paste(chromosome,cut_gRNa_alignment,sep=":")),
                           inherit.aes = F,col = "purple",cex=3,nudge_x = 1,nudge_y = 0.5,force = 5,direction= "x")+
  geom_point(data = . %>% 
               filter(N_edits==0,PAM_indel_count==0), 
             aes(rank_desc,N_UMI_cluster),col = "red",cex= 2,inherit.aes = F)+
  labs(x = "Ranked clusters by decreasing abundance",
       y = "UMI count per cluster",
       caption = paste("Only clusters with >=",minUMI_alignments_figure, "UMIs",sep=" "))+
  scale_color_manual(values = c("black","green3"))






## Chromosome distribution --------------------------------------------------------------------------

ChromDistr_data <- lapply(summary_pred_bulge, function(x){
  
  x <- x %>% 
    filter(!is.na(Alignment) & (str_starts(chromosome,"[0-9XYM]") | str_starts(chromosome,"chr"))) %>% 
    group_by(chromosome,Predicted=predicted) %>% 
    summarise(clusters = n(), reads = sum(N_reads_cluster), UMI = sum(N_UMI_cluster)) %>% 
    mutate(chromosome = recode(chromosome, "MT" = "M")) 
  
  if(!any(str_starts(x$chromosome,"chr"))){
    x <- x %>%
      mutate(chromosome = paste("chr",chromosome,sep=""))
  } else{
    x
  }
}
)

ChromDistr_data <- ChromDistr_data %>% bind_rows(.id="library")

ChromDistr_data <- ChromDistr_data %>% 
  mutate(chromosome = factor(chromosome, levels = paste("chr",c(1:30,"X","Y","M"),sep="")),
         Predicted = factor(Predicted, levels = c("yes","no")))


# figures
fig_ChromDistr_clusters <- ggplot(ChromDistr_data, aes( chromosome, clusters,fill = Predicted)) +
  geom_col(col = "black") +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  coord_polar() +
  theme_bw(base_size = 12)+
  scale_fill_manual(values = c("green4","grey")) +
  theme(panel.grid.major.x = element_line(linetype = 2,colour = "grey40"),axis.text =  element_text(color = "black"))


fig_ChromDistr_UMI <- ggplot(ChromDistr_data, aes( chromosome, UMI,fill = Predicted)) +
  geom_col(col = "black") +
  facet_wrap(~library, ncol = 3, scales = "free_y") + 
  coord_polar() +
  theme_bw(base_size = 12)+
  scale_fill_manual(values = c("green4","grey")) +
  theme(panel.grid.major.x = element_line(linetype = 2,colour = "grey40"),axis.text =  element_text(color = "black"))



# Distribution of UMI around cut best cut site ------------------------------

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



# HTML tables for off targets ----------------------------------------------------------------


tables_off <- lapply(seq_along(summary_pred_bulge),function(x){
  
  cat("Generating html table for sample : ",names(summary_pred_bulge)[x],"\n")
  
  df <- summary_pred_bulge[[x]] %>% filter(!is.na(Alignment))  ## keep only clusters with gRNA match
  
  ## subset top n clusters by UMI count (defined in the configuration file)
  dt <-  df  %>% ungroup %>% 
    mutate("UMIs (%)" = round(N_UMI_cluster / sum(N_UMI_cluster)*100,digits = 1),.after=N_UMI_cluster) %>% 
    filter(N_UMI_cluster > minUMI_alignments_figure) %>% 
    slice_max(n = max_clusters,order_by = N_UMI_cluster,with_ties = F) %>% 
    as.data.table(dt)
  
  
  compare_strings <- function(DNA,RNA){
    
    DNA <- str_split_1(DNA,"")
    RNA <- str_split_1(RNA,"")
    
    DNA[DNA==RNA]<- "."
    return(paste(DNA,collapse=""))
    
  }
  
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
    mutate("cut offset" = as.numeric(cut_modal_position) - as.numeric(cut_gRNa_alignment),
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
           #N_mismatches,
           #n_indels,
           #soft_trim,
           "Edits pam"= PAM_indel_count,
           Symbol=Symbol_html,
           Feature=position_html,
           Predicted= predicted_alignment_html,
           #bulge,
           PCR,
           #multiHit,
           #n_UMI_multiSum,
           Oncogene = Oncogene_html) %>% 
    unite(col = "Position",chromosome,cut_gRNa_alignment,sep = ":") 
    #unite("Mismatch_indels",N_mismatches, n_indels,sep = "_",remove = T) %>% 
    
  
  
  
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
  
})

names(tables_off) <- names(summary_pred_bulge)



# Save report data ------------------------------------------------------------------

save(list = c("tables_off", "summary_pred_bulge","sampleInfo", "stats", "stats_summary","best_aligns", grep(x=ls(),"fig",value = T)), file = "results/report.rdata")







