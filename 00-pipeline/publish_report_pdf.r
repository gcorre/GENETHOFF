# ========================================================= #
# Guillaume CORRE @ GENETHON 2025
# This script is used to publish the GUIDE-seq report in HTML format.
# It is used in the snakemake pipeline for GUIDE-seq analysis.
# ========================================================= #

# Set global options
options(tidyverse.quiet = TRUE,warn = -1,verbose = F,conflicts.policy = list(warn = FALSE))

# Load necessary libraries
library(rmdformats,quietly = T,warn.conflicts = F,verbose = F)

# get arguments from command line
args <- commandArgs(trailingOnly = T)

# debug
# args <- c("M:\\DB_Unite de Genomique\\Bioinformatique\\guideSeq\\PRJNA752948_bcl11_mario\\results\\report.rdata", 
#           "M:\\DB_Unite de Genomique\\Bioinformatique\\guideSeq\\PRJNA752948_bcl11_mario\\guideSeq_GNT.yml",
#           "./",
#           "test.pdf")


# parse arguments to objects
load(args[1])
config_path <- args[2]
template <- paste(args[4],"/genethOFF_report_pdf.rmd",sep="") 


# Generate the html report

rmarkdown::render(input = template, 
                  output_format = "pdf_document", 
                  output_dir = dirname(args[3]),
                  output_file = basename(args[3]),
                  clean = T)