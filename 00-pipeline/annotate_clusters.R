# ========================================================= #
# Guillaume CORRE @ GENETHON 2025
# This script is used to annotate cutting sites from GUIDE-seq data with gene and oncogenes
# It is used in the snakemake pipeline for GUIDE-seq analysis.
# ========================================================= #


## the main issue is to have common chromosom names between all databases.

# Set R options
options(tidyverse.quiet = TRUE,warn = -1,verbose = F,warn = -10,conflicts.policy = list(warn = FALSE))


# load necessary libraries
library(tidyverse,quietly = T, verbose = F,warn.conflicts = F)
library(GenomicRanges,quietly = T, verbose = F,warn.conflicts = F)
library(annotatr,quietly = T, verbose = F,warn.conflicts = F)
library(writexl,quietly = T, verbose = F,warn.conflicts = F)
library(rtracklayer,quietly = T, verbose = F,warn.conflicts = F)
library(pwalign,quietly = T, verbose = F,warn.conflicts = F)


# get argument from command line 

#debug
# args <- c("04-IScalling/SaKo8_ODN.cluster_slop.bed",
#           "04-IScalling/SaKo8_ODN.UMIs_per_IS_in_Cluster.bed",
#           "GRCm39.rds",
#           "None",
#           "06-offPredict/GRCm39_TGTCCAGGGCTAGCTTAACG_NNGRRT_6.csv",
#           50,
#           "results/SaKo8_ODN_summary.csv"
#           )


args <- commandArgs(trailingOnly = T)

# 1 : Slop.bed
# 2 : UMI_perIS_per_cluster
# 3 : annotation rds file
# 4 : onco_list  -> converted to "None" if empty in configuration file
# 5 : prediction file
# 6 : min_predicted distance --> disatnce for a cut to be predicted by SWOffinder
# 7 : output file

min_predicted_distance <- as.numeric(args[6])
#----------------------------------------------------------#
# 01 - Get clusters informations ----
#----------------------------------------------------------#

bed <- args[1]


clusters_df <- read.delim(bed,
                          header=F, 
                          col.names = c("chromosome",
                                        "start_cluster",
                                        "end_cluster",
                                        "clusterID",
                                        "medianMAPQ_cluster",
                                        "N_IS_cluster",
                                        "N_orientations_cluster",
                                        "N_orientations_PCR",
                                        "N_UMI_cluster",
                                        "N_reads_cluster"))


#----------------------------------------------------------#
# 02 - Extract cluster information by PCR orientation ----
#----------------------------------------------------------#
bed <- args[2]

## Load data 

bed_df <- read.delim(bed,
                     header=F, 
                     col.names = c("chromosome",
                                   "start_IS",
                                   "end_IS",
                                   "IS_ID",
                                   "MedianMAPQ_IS",
                                   "strand_IS",
                                   "PCR_orientation",
                                   "N_UMI_IS",
                                   "Nreads_IS",
                                   "UMI_list",
                                   "ReadPerUMI",
                                   "clusterID"))


## get the chromosome format from bed files
bed_chromosome_name_format <- ifelse(all(str_starts(bed_df$chromosome, "chr")),yes = "long",no = "short")



rm(bed)

## Breakdown positive and negative PCRs
# PCR orientation is a factor with the 2 levels. If one level is missing, it is still reported with value 0.


clusters_split_orientation <- bed_df %>% 
  group_by(clusterID,PCR_orientation,.drop = F) %>% 
  summarise(N_IS = n_distinct(IS_ID),
            N_UMI = sum(N_UMI_IS),
            N_Reads = sum(Nreads_IS),
            N_Orientation = n_distinct(strand_IS),
            MapQ = median(MedianMAPQ_IS)
  ) %>%
  pivot_wider(names_from = PCR_orientation, values_from = c(N_IS:MapQ), values_fill = 0)



# get the position with most abundant UMIs in a cluster
modal_cut_position <- bed_df %>% 
  group_by(clusterID, start_IS) %>%
  summarise(count_UMI = sum(N_UMI_IS)) %>% 
  mutate(UMI_proportion = round(count_UMI / sum(count_UMI)*100,digits = 1)) %>% 
  group_by(clusterID) %>% 
  slice_max(n = 1, order_by = UMI_proportion,with_ties = F) %>%    ### here pick the first if multiple modes (especialty when low number of UMI and number of IS)
  rename("start_IS"="cut_modal_position")



## annotate clusters with identified gRNA sequence
clusters_split_orientation <- clusters_split_orientation %>%
  left_join(modal_cut_position, by = c("clusterID")) 


rm(modal_cut_position)
rm(bed_df)


#----------------------------------------------------------#
# 03 - Merge total and PCR orientation specific cluster informations ----
#----------------------------------------------------------#

clusters_df <- clusters_df %>% 
  left_join(clusters_split_orientation, by = "clusterID")

rm(clusters_split_orientation)




#----------------------------------------------------------#
# 04 - Perform gene feature annotation ----
#----------------------------------------------------------#



if(nrow(clusters_df)>0){
  
  # load pre-prepared annotation GRange object
  gtf <- readRDS(args[3])
  
  
  # check that all sources have the same chromosome format (chrXX or just XX)
  gtf_chromosome_name_format <- ifelse(all(str_starts(levels(gtf@seqnames), "chr")),yes = "long",no = "short")
  
  if(bed_chromosome_name_format != gtf_chromosome_name_format){
    # convert all sources to gtf chromosome format
    
    clusters_df <- clusters_df %>% 
      mutate(chromosome = case_when(gtf_chromosome_name_format == "long" ~ paste("chr",chromosome, sep= ""), 
                                    TRUE ~ str_match(chromosome,"chr([0-9XYMT]+)")[,2]))
  }
  

  
  
  
  # convert OT to gRanges
  results_granges <- makeGRangesFromDataFrame(clusters_df, 
                                              ignore.strand = T,
                                              keep.extra.columns = T,
                                              start.field = "start_cluster",
                                              end.field = "end_cluster",
                                              seqnames.field = "chromosome",
                                              na.rm = T)
  
  
  
  # annotate gRanges
  
  results_granges_annot <- annotate_regions(regions = results_granges,
                                            annotations = gtf,
                                            ignore.strand = T,
                                            minoverlap = 1)
  
  
  
  # convert to data.frame
  results_granges_df <- data.frame(results_granges_annot) %>%
    dplyr::select(clusterID,starts_with('anno') ) 
  
  rm(gtf)
  rm(results_granges)
  rm(results_granges_annot)
  
  genes_annot <- results_granges_df %>%
    filter(annot.type == "gene") %>%
    select(-annot.type,-annot.transcript_id)%>% ungroup()  %>%
    distinct()
  
  
  features_annot <- results_granges_df %>% 
    filter(annot.type != "gene") %>% 
    select(clusterID,annot.seqnames,annot.gene_id,annot.type)%>% 
    distinct() %>% 
    ungroup() %>% 
    mutate(annot.type = as.character(annot.type))
  

  
    #----------------------------------------------------------#
  # 05 - Perform oncogene annotation ----
  #----------------------------------------------------------#
  
  
  # human list examples : https://bioinfo.uth.edu/TSGene/download.cgi & https://bioinfo-minzhao.org/ongene/download.html 
  
  if(args[4]!="None"){
    
    onco_list_df <- read.delim(args[4],sep=";") 
    ## this file must contain columns : "transcriptID", "Onco_annotation"
    if(all(c("transcriptID","Onco_annotation") %in% colnames(onco_list_df)) & nrow(onco_list_df)>0){
      onco_list_df <- onco_list_df %>%
        distinct(transcriptID,Onco_annotation)%>% 
        filter(!is.na(Onco_annotation))
      
      ## this file must contain columns : "transcriptID", "Onco_annotation"
      # convert to gene Id
      onco_annotations <- results_granges_df %>% 
        ungroup %>% 
        filter(!annot.transcript_id=="") %>% 
        distinct(annot.gene_id,annot.transcript_id) %>% 
        left_join(onco_list_df, by = c("annot.transcript_id" = "transcriptID") ) %>% 
        select(-annot.transcript_id) %>% 
        distinct() %>% 
        filter(!is.na(Onco_annotation))
      
      rm(onco_list_df)
    } else{
      
      errorCondition("Oncogene list is empty or is not correctly formatted. 
                     \nColnames must be : 'ensembl.transcriptID','Is.Oncogene','Is.Tumor.Suppressor.Gene'")
      
    }
  } else {
    warning("Oncogene list was empty. Annotating with NAs")
    onco_annotations <- results_granges_df %>% 
      distinct(annot.gene_id) %>%  
      mutate(Onco_annotation = "Not_Done")
  }
  
  
  #----------------------------------------------------------#
  # 06 - aggregate annotations ----
  #----------------------------------------------------------#
  
  results_granges_df_annot <- genes_annot %>% 
    left_join(features_annot, by = c("clusterID","annot.seqnames", "annot.gene_id")) %>% 
    left_join(onco_annotations, by = c("annot.gene_id")) %>% 
    replace_na(list(annot.type = "intron",onco_annotations = "")) %>% 
    replace_na(list("Onco_annotation" = ""))
  
  rm(genes_annot)
  rm(features_annot)
  rm(onco_annotations)
  rm(results_granges_df)
  
  #----------------------------------------------------------#
  # 07 - collapse annotations per cluster ----
  #----------------------------------------------------------#
  
  results_granges_df_annot <- results_granges_df_annot %>%
    distinct(clusterID,annot.type,annot.gene_id,annot.gene_name,annot.gene_biotype,Onco_annotation)%>%
    group_by(clusterID) %>% 
    summarise(gene_ensemblID = toString(annot.gene_id),
              Symbol = toString(annot.gene_name),
              gene_biotype = toString(annot.gene_biotype),
              position = toString(annot.type),
              Onco_annotation = toString(Onco_annotation))
  
 
 
  #----------------------------------------------------------#
  # 08 - Annotate clusters ----
  #----------------------------------------------------------#
  
  clusters_df <- clusters_df %>%
    left_join(results_granges_df_annot, by = "clusterID")
  
  rm(results_granges_df_annot)
  
  
  #----------------------------------------------------------#
  # 09 - Perform cluster/prediction annotation ----
  #----------------------------------------------------------#
  
  prediction_gRNA <- data.table::fread(args[5], sep =",") %>% 
    filter(str_starts(Chromosome, "(chr)?[0-9XYMT]+"))
  prediction_gRNA <- prediction_gRNA %>% 
    mutate(Chromosome = str_match(Chromosome,"^((chr)?[0-9XYMT]+)")[,2]) %>% 
    filter(!is.na(Chromosome))
  
  prediction_chromosome_name_format <- ifelse(all(str_starts(prediction_gRNA$Chromosome, "chr")),yes = "long",no = "short")
  
  if(bed_chromosome_name_format != prediction_chromosome_name_format){
    # convert all sources to prediction chromosome format
    
    clusters_df <- clusters_df %>% 
      mutate(chromosome = case_when(prediction_chromosome_name_format == "long" ~ paste("chr",str_match(chromosome,"^([0-9XYMT]+)")[,2], sep= ""), 
                                    TRUE ~ str_match(chromosome,"chr([0-9XYMT]+)")[,2]))
  }
  
  
  gr_predict <- makeGRangesFromDataFrame(prediction_gRNA,
                                         start.field = "EndPosition",
                                         end.field = "EndPosition",
                                         seqnames.field = "chromosome",
                                         keep.extra.columns = F)
  
  
  gr_data <- makeGRangesFromDataFrame(clusters_df, 
                                      start.field = "cut_modal_position",
                                      end.field = "cut_modal_position",
                                      seqnames.field = "chromosome",
                                      keep.extra.columns = F)
  
  
  
  
  # distance to nearest predicted cut site
  near_df <- distanceToNearest(gr_data,gr_predict) %>% data.frame
  
  near_df <- near_df %>% 
    select("clusterID"=queryHits,"closest_predicted_distance" = "distance") %>% 
    mutate(predicted = case_when(closest_predicted_distance < min_predicted_distance ~ "yes",
                                 TRUE ~"no"))
  
  
  # clear 
  rm(gr_predict);rm(gr_data)
  
  
  
  ## annotate clusters with prediction (by position, NOT exact match)

  clusters_df <- clusters_df %>% left_join(near_df, by = "clusterID")
  
  clusters_df <- clusters_df %>% 
    replace_na(list("Onco_annotation" = ""))
  
  
  
  
  ## convert chromosome names to long format

  clusters_df_chromosome_name_format <- ifelse(all(str_detect(clusters_df$chromosome, "chr")),yes = "long",no = "short")
  
  if(clusters_df_chromosome_name_format != "long"){
    clusters_df <- clusters_df %>% 
      mutate(chromosome = paste("chr",chromosome, sep= ""))
  }
  
  
  # format oncogene annotation -> remove empty cells
  
  clusters_df <- clusters_df %>% 
    mutate(Onco_annotation = case_when(!str_detect(Onco_annotation,"[A-Za-z0-9]+")~ "",
                                       TRUE ~ Onco_annotation))
}
  #----------------------------------------------------------#
  # 10 - save summary file as xlsx ----
  #----------------------------------------------------------#
  
 
  write.table(clusters_df,file= args[7],sep = "\t", row.names = F, quote = F)
  

