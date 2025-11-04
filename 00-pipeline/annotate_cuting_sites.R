# ========================================================= #
# Guillaume CORRE @ GENETHON 2025
# This script is used to annotate cutting sites from GUIDE-seq data with gene and oncogenes
# It is used in the snakemake pipeline for GUIDE-seq analysis.
# ========================================================= #

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
args <- commandArgs(trailingOnly = T)

# Parse arguments to objects
annotation <- args[1]
file <- args[2]
output<- args[3]
onco_list <- args[4]


###########################################################
## debug
 # file = "05-Report/EBS_hDMD.rdata"
 # annotation <- "GRCh38"
 # onco_list <- "../02-ressources/OncoList_OncoKB_GRCh38_2025-07-04.tsv"
 # output <- "results/EBS_hDMD_summary.xlsx"
###########################################################

# load sample rdata file
library <- str_remove(basename(file),pattern = ".rdata")
load(file)


## if file is not empty ...
if(exists("clusters_grna_match") && nrow(clusters_grna_match)>0){
  
  ## sort & convert some features type
  results_df <- clusters_grna_match %>% 
    arrange(desc(N_UMI_cluster)) %>% 
    mutate(gRNA = as.character(grna),
           gRNA_name = grna@metadata$name) 
  rm(clusters_grna_match);
  
  
  # load pre-prepared annotation GRange object
  gtf <- readRDS(annotation)
  
  # convert OT to gRanges -----------------------------------------------------
  results_granges <- makeGRangesFromDataFrame(results_df, 
                                              ignore.strand = T,
                                              keep.extra.columns = T,
                                              strand.field = "strand_guide",
                                              start.field = "start_cluster",
                                              end.field = "end_cluster",
                                              seqnames.field = "chromosome",
                                              na.rm = T)
  
  
  
  
  # annotate gRanges -------------------------------------------------------------------
  
  results_granges_annot <- annotate_regions(regions = results_granges,
                                            annotations = gtf,
                                            ignore.strand = T,
                                            minoverlap = 1)
  
  
  # convert to data.frame
  results_granges_df <- data.frame(results_granges_annot) %>%
    dplyr::select(clusterID,starts_with('anno') ) 
  
  
# Add oncogene / tumors suppressing gene annotation ---------------
  
  # human list examples : https://bioinfo.uth.edu/TSGene/download.cgi & https://bioinfo-minzhao.org/ongene/download.html 
  
  if(onco_list!=""){
    
    onco_list_df <- read.delim(onco_list,sep=";") 
    
    if(all(c("ensembl.transcriptID","Onco_annotation") %in% colnames(onco_list_df)) & nrow(onco_list_df)>0){
      
      ## this file must contain columns : "ensembl.transcriptID", "Onco_annotation"
      
      onco_list_df <- onco_list_df %>%
        distinct(ensembl.transcriptID,Onco_annotation)
      
      results_granges_df <- results_granges_df %>%
        left_join(onco_list_df, by = c("annot.transcript_id" = "ensembl.transcriptID") ) %>% 
        select(-annot.transcript_id) %>%  
        distinct() %>% 
        replace_na(replace = list(Onco_annotation=""))
    } else{
      
      errorCondition("Oncogene list is empty or is not correctly formatted. 
                     \nColnames must be : 'ensembl.transcriptID','Is.Oncogene','Is.Tumor.Suppressor.Gene'")
      
    }
  } else {
    warning("Oncogene list was empty. Annotating with NAs")
    results_granges_df <- results_granges_df %>% mutate(Onco_annotation = "")
  }
  

  
  # collapse annotation per cluster -----------------------------------------------------
  results_granges_df_annot <- results_granges_df %>%
    distinct(clusterID,annot.type,annot.gene_id,annot.gene_name,annot.gene_biotype,Onco_annotation)%>%
    group_by(clusterID,annot.gene_id,annot.gene_name,annot.gene_biotype,Onco_annotation) %>%
    summarise(pos =toString(annot.type)) %>% 
    mutate(position=case_when(str_detect(pos,"exon") ~ "exon", TRUE ~ "intron")) %>%
    dplyr::select(-pos) %>% 
    group_by(clusterID) %>% 
    summarise(gene_ensemblID = toString(annot.gene_id),
              Symbol = toString(annot.gene_name),
              gene_biotype = toString(annot.gene_biotype),
              position = toString(position),
              Onco_annotation = toString(Onco_annotation))
  
  

  # Annotate clusters with oncogenes------------------------------------------------------------------
  
  results_granges_df_annot <- results_df %>%
    left_join(results_granges_df_annot, by = "clusterID")
  
  out_list <- list(results_granges_df_annot )
  names(out_list) <- library
  
  # save to excel file -----------------------------------------------------
  
  writexl::write_xlsx(out_list,path = output,col_names = T,format_headers = T)
  
}
