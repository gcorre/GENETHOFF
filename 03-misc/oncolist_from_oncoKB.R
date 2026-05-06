###  Format onco gene list from oncoKG

# The output must be a ";" separated table with 2 columns only and named as below: 
# - transcriptID
# - Onco_annotation


library(tidyverse)

#update <- format.Date(x = Sys.Date(),format = "%Y-%m-%d")
update <- "2026-01-29"  # get from onkodb website 
ids <- "refseq"  # or ensembl

## Example for GRCh38 using oncoKB database
# onco_list <- read.delim("https://www.oncokb.org/api/v1/utils/cancerGeneList.txt")
onco_list <- read.delim("../03-misc/cancerGeneList.tsv")


onco_list_output <- onco_list %>% 
  filter(GRCh38.RefSeq!="",X..of.occurrence.within.resources..Column.J.P.>1) %>%
  rename("GRCh38.RefSeq"="transcriptID") %>% 
  select(transcriptID,Gene.Type) %>% 
  group_by(transcriptID) %>% 
  summarise(Onco_annotation = paste(unique(Gene.Type),collapse = "|"))

write.table(onco_list_output, file = paste("02-ressources/OncoList_OncoKB_",ids,"_",update,".tsv",sep = ""), sep=";")

