# CORRE Guillaume @ GENETHON
# 2026-08-14
# Benchmark analysis of GENETHOFF, GUIDE-seq (v1.0.2 / v2) and iGUIDE-seq (v1.2.0)


#---------------------------------------------##
# Load libraries ----
#---------------------------------------------##

library(tidyverse)
library(GenomicRanges)
library(patchwork)
library(ggprism)

slop_window = 10 # use to cluster cleavage sites across pipelines

#---------------------------------------------##
# Run time analysis----
#---------------------------------------------##

runtime <- read.delim("runtime_benchmark.csv", sep =";", check.names = F)



runtime_raw <- ggplot(runtime, aes(total_threads, as.numeric(real_min), fill = factor(tools))) +
  geom_point(size = 3, pch = 21, col = "black")+
  ggprism::theme_prism(base_size = 12)+
  stat_summary(
    fun = mean,
    geom = "line",
    aes(group = tools, col = tools),
    linetype = 2,
    linewidth = 1, show.legend = F
  )+
  labs(x = "Provided threads", y = "Runtime (min)")+
  scale_x_continuous(breaks = c(6,12,24))+
  scale_y_continuous()

scalability <- ggplot(runtime, aes(total_threads,`user+sys`/as.numeric(real_min), fill = factor(tools))) +
  geom_point(size = 3, pch = 21, col = "black")+
  geom_abline(lty = 2, col = "grey") + 
  ggprism::theme_prism(base_size = 12)+
  stat_summary(
    fun = mean,
    geom = "line",
    aes(group = tools, col = tools),
    linetype = 2,
    linewidth = 1, show.legend = F
  )+
  labs(x = "Provided threads", y = "Estimated threads")+
  scale_x_continuous(breaks = c(6,12,24))+
  scale_y_continuous(breaks = c(0,6,12,18,24,30))
  

#---------------------------------------------##
# Load each pipeline results table ----
#---------------------------------------------##

  #---------------------------------------------##
  ## GUIDESEQ v1----
  #---------------------------------------------##
  
  # list files of interest 

  files <- list.files("GUIDE-seq_v1.0.2/identified/", full.names = T)
  names(files) <- str_replace(str_remove(basename(files),"_identifiedOfftargets.txt"),"-","_")
  files <- grep(files, pattern = "control", value = T, invert = T)
  
  # read them all 
  guideseq <- lapply(files, read.delim)
  
  guideseq_v1 <- lapply(guideseq, function(x){
    
    x %>% 
      mutate(Chromosome = as.character(Chromosome))%>% 
      select(Chromosome, Position, bi.sum.mi) 
  })
  
  rm(files);
  

  # merge replicates for each gRNA
  targets <- names(guideseq_v1)
  targets_unique <- str_extract(targets,"^[A-Za-z0-9]+") %>% unique()
  
  guideseq_v1_aggregated <- list()
  
  for(i in targets_unique){
    cat("processing",i,"\n")
    
    guideseq_v1_aggregated[[i]] <- bind_rows(guideseq_v1[str_starts(names(guideseq_v1),i)])%>% 
      group_by(Chromosome, Position) %>% 
      summarise(counts = sum(bi.sum.mi))
  }
  

  guideseq_v1 <- bind_rows(guideseq_v1_aggregated,.id = "gRNA")
  
  
  # clean objects
  rm(targets); rm(guideseq); rm(guideseq_v1_aggregated)
  
  #---------------------------------------------##
  ## GUIDESEQ v2 ----
  #---------------------------------------------##
  
  # list files of interest 
    
  files <- list.files("GUIDE-seq_v2/identified/", full.names = T, pattern = "_identifiedOfftargets.txt$")
  files <- grep(files,pattern = "Control",value = T, invert = T)
  names(files) <- str_replace(str_remove(basename(files),"_identifiedOfftargets.txt"),"-","_")

  # read them all 
  guideseq <- lapply(files, read.delim)

  guideseq_v2 <- lapply(guideseq, function(x){
    
    x %>% 
      mutate(Chromosome = as.character(X.BED_Chromosome))%>% 
      select(Chromosome, Position, bi.sum.mi)  
    
  })

  rm(files); rm(guideseq)

  
  # merge replicates for each gRNA
  targets <- names(guideseq_v2)
  targets_unique <- str_extract(targets,"^[A-Za-z0-9]+") %>% unique()
  
  guideseq_v2_aggregated <- list()
  
  for(i in targets_unique){
    
    guideseq_v2_aggregated[[i]] <- bind_rows(guideseq_v2[str_starts(names(guideseq_v2),i)])%>%
      group_by(Chromosome, Position) %>% 
      summarise(counts = sum(bi.sum.mi))
  }
  
  guideseq_v2 <- bind_rows(guideseq_v2_aggregated,.id = "gRNA")
  
  
  # clean objects
  rm(targets); rm(guideseq); rm(guideseq_v2_aggregated)
  
  


  
  #---------------------------------------------##
  ## GENETHOFF ----
  #---------------------------------------------##
  path <- "GENETHOFF/results/paper_review.xlsx"
  
  genethoff <- path %>% 
    readxl::excel_sheets() %>% 
    set_names() %>% 
    map(readxl::read_excel, path = path)
  
  rm(path)
  
  
  # remove sites without gRNA match
  
  genethoff_sub <- lapply(genethoff, function(x){
    x %>%
      filter(!is.na(Alignment),
             abs(relative_distance) < 15) %>%
      select(chromosome, N_UMI_cluster,cut_modal_position)  
    
  }
  )
  
  
  # merge libraries for each gRNA
  targets <- names(genethoff_sub)
  targets_unique <- str_extract(targets,"^[A-Za-z0-9]+") %>% unique()
  targets_unique <- grep(targets_unique,pattern = "Mock",ignore.case = T,invert = T, value = T)
  
  genethoff_aggregated <- list()
  
  for(i in targets_unique){
  
  genethoff_aggregated[[i]] <- bind_rows(genethoff_sub[str_starts(names(genethoff_sub),i)])%>% 
    group_by(Chromosome = chromosome, Position = cut_modal_position) %>% 
    summarise(counts = sum(N_UMI_cluster))
  }
  
  
  genethoff <- bind_rows(,.id = "gRNA")

  rm(targets); rm(genethoff_aggregated); rm(genethoff_sub)

  
  #---------------------------------------------##
  ## iGUIDE-seq ----
  #---------------------------------------------##
  path <- "iGUIDE-seq_v1.2.0/manual_data_extraction.xlsx"
  
  iguide <- path %>% 
    readxl::excel_sheets() %>% 
    set_names() %>% 
    map(readxl::read_excel, path = path)
  
  rm(path)
  
  
  iguide_aggregated <- lapply(iguide, function(x){
    
    x %>% select("Edit Site", counts="Abund.") %>%
      separate(convert = T,"Edit Site", into = c("Chromosome","strand","Position"), sep =":") %>% 
      select(-strand)
  })
  
  iguideseq <- bind_rows(iguide_aggregated,.id = "gRNA")

  rm(path); rm(iguide_aggregated);rm(iguide)
  
  
#---------------------------------------------##
# Aggregate the 4 pipelines results ----
#---------------------------------------------##
  
  all_gRNAs_pipelines <- bind_rows("GUIDE-seq_v1" = guideseq_v1,
                   "GUIDE-seq_v2"= guideseq_v2,
                   "iGUIDE-seq" = iguideseq,
                   "GENETHOFF"= genethoff, .id = "workflow")
  
  

#---------------------------------------------##
# Cluster cleavage sites across pipelines ----
#---------------------------------------------##
  
    #---------------------------------------------##
    ## Reformat columns ----
    #---------------------------------------------##
  
  all_gRNAs_pipelines <- all_gRNAs_pipelines %>%
    mutate(start = Position ,
           end = Position , 
           Chromosome = str_remove_all(Chromosome,"chr"),
           gRNA = case_when(gRNA=="TRAC5"~ "TRAC",
                            TRUE ~ gRNA),
           workflow = factor(workflow))
  
  
    #---------------------------------------------##
    ## Convert to gRange object for clustering ----
    #---------------------------------------------##
  all_gRNAs_pipelines_grange <- makeGRangesFromDataFrame(df = all_gRNAs_pipelines, 
                                                         keep.extra.columns = T,
                                                         seqnames.field = "Chromosome",
                                                         start.field = "start",
                                                         end.field = "end")

  
    #---------------------------------------------##
    ## Cluster all cleavage sites ----
    #---------------------------------------------##
      
    # find overlapping regions in a window
    window = 10
    
    hits <- findOverlaps(all_gRNAs_pipelines_grange, reduce(all_gRNAs_pipelines_grange, min.gapwidth = window))
    cluster_id <- subjectHits(hits)[order(queryHits(hits))]

    # annotate cleavage sites with cluster ID
    all_gRNAs_pipelines$cluster <- cluster_id


    
    #---------------------------------------------##
    # Get relative abundance of cleavage sites per gRNA and cluster ----
    #---------------------------------------------##
    
    # calculate relative abundance of cleavage site per gRNA and method
    all_gRNAs_pipelines <- all_gRNAs_pipelines %>% 
      group_by(workflow,gRNA, Chromosome, cluster) %>% 
      summarise(start = min(start),
                end = max(end),
                counts = sum(counts)) %>% 
      group_by(workflow,gRNA) %>% 
      mutate(prop = counts / sum(counts) * 100) %>% 
      unite(col = "Position",Chromosome,start,end,remove = F,sep = "_")
    
    
    
    
    #---------------------------------------------##
    # Annotate cluster if they are the on target expected site ----
    #---------------------------------------------##
    
    all_gRNAs_pipelines <- all_gRNAs_pipelines %>% 
      mutate(OT = case_when(gRNA=="B2M" & Position %in% c("15_44711567_44711568" ,"15_44711569_44711569","15_44711568_44711568") ~ T,
                                     gRNA=="TRAC" & Position %in% c("14_22547664_22547664","14_22547663_22547664") ~ T,
                                     gRNA=="VEGFAs2" & Position %in% c("6_43770825_43770825","6_43770824_43770824") ~ T,
                                     gRNA=="VEGFAs3" & Position %in%c("6_43769732_43769732", "6_43769733_43769733") ~ T,
                                     TRUE ~ FALSE)
             )

    
    #---------------------------------------------##
    # Pivot table to get pipeline in columns ----
    #---------------------------------------------##
    
    ## add method OT positions
    all_gRNAs_pipelines_wide <- all_gRNAs_pipelines %>% 
      pivot_wider(names_from = "workflow", values_from = c(Position,counts, prop),names_glue ="{workflow}_{.value}", id_cols = c("gRNA","cluster") )
    

#---------------------------------------------##
# Save tables  ----
#---------------------------------------------##

write.table(all_gRNAs_pipelines, paste("complete_OT_table_all_methods_",window,"bp.csv",sep=""), sep=";",row.names = F, quote = F)

write.table(all_gRNAs_pipelines_wide, paste("complete_OT_table_all_methods_wide_",window,"bp.csv",sep=""), sep=";",row.names = F, quote = F)





#---------------------------------------------##
# Make some plots ----
#---------------------------------------------##
library(UpSetR)
library(tidyverse)

#reload tables if necessary 
all_gRNAs_pipelines <- read.delim("complete_OT_table_all_methods_10bp.csv", header = T, sep =";")
all_gRNAs_pipelines_wide <- read.delim("complete_OT_table_all_methods_wide_10bp.csv", header = T, sep =";")



  #---------------------------------------------##
  ## Rank-abundance plot ----
  #---------------------------------------------##

rank_ab <- ggplot(all_gRNAs_pipelines %>%
                    filter(counts >0) %>% 
         group_by(gRNA,workflow) %>% 
         mutate(rank = row_number(-counts)), 
       aes(rank,counts, col = workflow)) + 
  geom_line(show.legend = F) +
  geom_point(show.legend = F)+
  facet_wrap(~gRNA, nrow = 2, scale = "free_x") +
  scale_y_log10(breaks = c(1,10,100,1000,10000,100000))+
  scale_x_log10(breaks = c(1,10,100,1000), limits = c(1,NA))+
  ggprism::theme_prism(base_size = 12, border = T) +
  geom_point(data = . %>% filter(OT==T),col = "black", shape = 21, size = 1,stroke = 2,show.legend = F)+
  labs(x = "Rank (descending abundance)", y = "Abundance" ) 

rank_ab


  #---------------------------------------------##
  ## Rank of top cleavage sites ----
  #---------------------------------------------##
  
  # keep site that are in the 4 pipelines
  keep <- all_gRNAs_pipelines %>%
  ungroup %>% 
    count(gRNA,cluster) %>%
    filter(n>=4)

  ggplot(data = all_gRNAs_pipelines %>% 
           group_by(gRNA,workflow) %>% 
           mutate(rank = row_number(-counts)) %>% # calculate rank per gRNA and pipeline
           semi_join(keep) %>%                    # keep sites present in 4 pipelines
           filter(#rank<= 25,
                  gRNA %in% c("VEGFAs2","VEGFAs3")), 
         aes(workflow,rank, col = factor(cluster))) + 
    geom_point(show.legend = c(fill=F),pch = 21,col = "black", aes(size = prop, fill = factor(cluster)), alpha = 0.6) + 
    geom_line(aes(group = cluster), show.legend = F)+
    facet_wrap(~gRNA, scale="free")+
    ggprism::theme_prism(base_size = 12,border = T,axis_text_angle = 45) +
    labs(x = NULL, y = "Cleavage site Rank")
  


  
  #---------------------------------------------##
  # Venn diagrams and upset plots ----
  #---------------------------------------------##
  library(UpSetR)
  library(ComplexUpset)
  library(ggvenn)

# set parameters

margins <- c(0,0,0,0)
text_size = 4


  venn_plot <- list()
  upset_plot <- list()

  for(i in c("B2M","TRAC","VEGFAs2", "VEGFAs3")){
    
    cat("processing",i,"\n")
    
    df <- all_gRNAs_pipelines %>% filter(gRNA == i )
    
    # make a list
    lst <- split(df$cluster,f=df$workflow)
    
    # make the venn
    venn_plot[[i]] <- ggvenn::ggvenn(lst,
                               set_name_color = "white",
                               fill_alpha = 0.8,
                               fill_color = scale_fill_hue()$palette(4),
                               set_name_size =3,
                               text_size = text_size,
                               show_percentage = FALSE) +
      ggtitle(i) + 
      ggprism::theme_prism(base_size = 12)+
      theme(plot.margin = unit(margins,units = "cm"), 
            axis.line = element_blank(),
            axis.title = element_blank(), 
            axis.text = element_blank(),
            axis.ticks = element_blank())
    
    
    
    # make the upset plot
    df <- fromList(lst)
    
    upset_plot[[i]]  <- ComplexUpset::upset(
      df,
      intersect = colnames(df),
      min_size=1,
      width_ratio=0.2,name = NULL,height_ratio = 0.6, 
      stripes='white',
      base_annotations=list(
        'Cleavage sites'= intersection_size(
          counts=TRUE,
          fill = "black")
      ),
      set_sizes = upset_set_size(
        geom = ggplot2::geom_bar(fill = "black")
      )+
        labs(y="Cleavage sites")
    ) 
    
    
  }



  #---------------------------------------------##
  ## Make the final montage ----
  #---------------------------------------------##

  venn_plots <- cowplot::as_grob((venn_plot$B2M + venn_plot$TRAC) / (venn_plot$VEGFAs2 + venn_plot$VEGFAs3))

x11();
(runtime_raw / scalability |plot_spacer()| rank_ab | plot_spacer()| venn_plots) +  
  plot_annotation(tag_levels = 'A')+
  plot_layout(guides = "collect",widths = c(1,0.25, 2,0.25,2))&
  theme(legend.position = "bottom")


#---------------------------------------------##
# Calculate Kendall W concordance score on ranks ----
#---------------------------------------------##

library(irr)

for(grna in c("B2M","TRAC","VEGFAs2","VEGFAs3")){
  x_df <- all_gRNAs_pipelines_wide %>% ungroup %>% 
    filter(gRNA == grna) %>% 
    select(cluster,ends_with("prop")) %>%
    column_to_rownames("cluster") %>%
    mutate(across(everything(), ~replace_na(.x, 0))) %>% 
    filter(if_any(everything(), ~ . >= 1)) %>%             # keep site that repsent more than x % of total abundance
    mutate(across(everything(), ~row_number(-.x)))         # calculte the rank
  if(nrow(x_df)>1){
    x = kendall(x_df,correct = T)
  }else {
    x = 1
  }
  cat("####################\n")
  print(x)
  cat("####################\n")
  }


