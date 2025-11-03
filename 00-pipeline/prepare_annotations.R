# ========================================================= #
# Guillaume CORRE @ GENETHON 2025
# This script is used to prepare gene annotations from a GTF file.
# It is used in the snakemake pipeline for GUIDE-seq analysis.
# ========================================================= #

options(tidyverse.quiet = TRUE,warn = -1,verbose = F,warn = -10,conflicts.policy = list(warn = FALSE))


library(rtracklayer,quietly = T, verbose = F,warn.conflicts = F)
library(tidyverse,quietly = T, verbose = F,warn.conflicts = F)
## convert gtf file to R object

args <- commandArgs(trailingOnly = T)


species <- args[1]
file <- args[2]

# build annotation objects ----------------------------------------------------------


  gtf <- rtracklayer::import(file)
  gtf <- gtf[gtf$type %in% c("gene","exon")]
  
  # replace rename gene_type gene_biotype if necessary
  if ("gene_type" %in% colnames(GenomicRanges::mcols(gtf))) {
    gtf$gene_biotype <- gtf$gene_type
  } else {
    cat("Error: Missing column. Please make sure your GTF file includes a 'gene_biotype' or 'gene_type' column \n", file = stderr())
    sys.exit(1)
  }
  
  # select only columns of interest
  mcols(gtf) <- gtf %>%
    data.frame %>% 
    dplyr::select(gene_id,gene_name,transcript_id,gene_biotype,type) %>% 
    mutate(gene_name = case_when(is.na(gene_name)~ gene_id,TRUE ~ gene_name))
  
  saveRDS(object = gtf, file = paste("../02-ressources/",species,".rds",sep = "" ))