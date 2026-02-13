# ========================================================= #
# Guillaume CORRE @ GENETHON 2025
# This script is used to align clusters of sequences to a gRNA sequence and report the best matches.
# It is used in the snakemake pipeline for GUIDE-seq analysis.
# ========================================================= #


# Set global options
options(tidyverse.quiet = TRUE,warn = -1,verbose = F,warn = -10,conflicts.policy = list(warn = FALSE))


#R.Version()$version.string == "R version 4.3.2 (2023-10-31 ucrt)"






#----------------------------------------------------------#
# 01 - Load necessary libraries ----
#----------------------------------------------------------#

library(tidyverse,quietly = T, verbose = F,warn.conflicts = F)
library(Biostrings,quietly = T,verbose = F, warn.conflicts = F)
library(DECIPHER,quietly = T,warn.conflicts = F,verbose = F)
library(pwalign,quietly = T, verbose = F,warn.conflicts = F)



#----------------------------------------------------------#
# 02- Parse arguments ----
#----------------------------------------------------------#

# debug:

args <- c("04-IScalling/Lib_3_K562_ODN_pos.cluster_slop.fa",
          "05-Report/Lib_3_K562_ODN_pos_summary.csv",
          "CTAACAGTTGCTTTTATCAC",
          "gRNA3",
          "NNN",
          -4,
          6,
          TRUE,
          3,
          "05-Report//Lib_3_K562_ODN_pos_summary_grna.csv")

# Get arguments from command line 
args <- commandArgs(trailingOnly = T)


fasta <- readDNAStringSet(args[1],
                          use.names = T)


grna <- DNAString(args[3])
gRNA_name <- args[4]
grna@metadata$name <- gRNA_name


pam <- DNAString(args[5])
pam_length <- nchar(pam)

offset <- as.numeric(args[6])

max_edits <- as.numeric(args[7])

bulges <- as.logical(str_to_upper(args[8]))

pam_side <- args[9]




## check that offset is in the good direction
## if PAM is in the 3', offset should be negative and opposite
offset <- if((pam_side == "3" & sign(offset) < 0) |(pam_side == "5" & sign(offset) > 0) ){
  offset <- offset
} else {
  offset <- offset * -1 # change sign
  
}




# Change gap penalty if bulges are tolerated or not

gap_open_penalty <- 10 ## default (this tolerate gap (bulges))
if(bulges == FALSE){
  gap_open_penalty <- gap_open_penalty *100  # this forces alignment without gap
}







#----------------------------------------------------------#
# 03 - Define matching function with IUPAC letters ----
#----------------------------------------------------------#

# This function is for the detection of mismatches between 2 DNAstring objects using IUAPC code

count_iupac_mismatches <- function(sequence, motif) {
  # Convert inputs to DNAString objects if they aren't already
  seq <- DNAString(sequence)
  mot <- DNAString(motif)
  
  # Get lengths
  seq_len <- length(seq)
  mot_len <- length(mot)
  
  if (seq_len != mot_len) {
    stop("Sequence and motif must be the same length")
  }
  
  # Convert to character vectors for comparison
  seq_chars <- strsplit(as.character(seq), "")[[1]]
  mot_chars <- strsplit(as.character(mot), "")[[1]]
  
  # Count mismatches considering IUPAC codes
  mismatches <- 0
  mm_position <- NULL
  for (i in 1:seq_len) {
    # Get the IUPAC alternatives for both positions
    seq_bases <- IUPAC_CODE_MAP[seq_chars[i]]
    mot_bases <- IUPAC_CODE_MAP[mot_chars[i]]
    
    # If there's no overlap between the possible bases, it's a mismatch
    if (length(intersect(unlist(strsplit(seq_bases, "")), 
                         unlist(strsplit(mot_bases, "")))) == 0) {
      mismatches <- mismatches + 1
      mm_position <- c(mm_position,i)
    }
  }
  
  return(paste(mismatches,paste(mm_position,collapse = ","),sep = "_"))
}




#----------------------------------------------------------#
# 04 - Pairwise alignments gRNA/gDNA sequences ----
#----------------------------------------------------------#

watson = pairwiseAlignment(pattern = fasta,subject = grna,type = "local-global",gapOpening=gap_open_penalty)

crick = pairwiseAlignment(pattern = reverseComplement(fasta),subject = grna,type = "local-global",gapOpening=gap_open_penalty)







# aggregate results to dataframe
align_stat <- data.frame(position = names(fasta),
                         watson_score = score(watson),
                         crick_score =  score(crick),
                         watson_edits = nedit(watson),
                         crick_edits = nedit(crick),
                         watson_mismatches = nmismatch(watson),
                         crick_mismatches = nmismatch(crick),
                         watson_matches = nmatch(watson),
                         crick_matches = nmatch(crick),
                         watson_pid = pid(watson),
                         crick_pid = pid(crick) ) %>% 
  separate("position", into = c("clusterID","cluster"),sep = "::",convert = T)



#----------------------------------------------------------#
# 05 - Process watson and crick orientations ----
#----------------------------------------------------------#

## if result table is not empty ...

if(nrow(align_stat)>0){
  
  # find orientation of gRNA with best score
  align_stat <- align_stat %>%
    mutate(grna_orientation = case_when(watson_score > crick_score ~ "watson",
                                        watson_score < crick_score  ~ "crick",
                                        TRUE ~ NA))
  
  
  
  
  
  
  ## WATSON strand processing
  cat("Analyzing WATSON strand")
  
  
  # get clusters with watson alignment of gRNA
  
  watson_best <- align_stat %>% 
    filter(grna_orientation == "watson", watson_edits <= max_edits) %>% ### keep clusters with less than n edits
    dplyr::select(clusterID,cluster,grna_orientation,starts_with("watson")) %>%
    rename_all(~str_remove(.,"watson_")) 
  
  
  # if there are some clusters with gRNA alignment in watson orientation...
  if(nrow(watson_best)>0){
    
    # extract alignments where gRNA matches in watson orientation
    watson_sub <- watson[watson_best$clusterID]
    
    #get alignment width (should match gRNA length +/- indels)
    width <- width(alignedPattern(watson_sub))
    
    
    
    
    
    # get INDELS counts and total width (INS + DEL separated)
    indels <- nindel(watson_sub)
    
    ins <- indels@insertion %>% data.frame
    names(ins) <- c('Insertion_Events','Insertion_length')
    
    del <- indels@deletion %>% data.frame
    names(del) <- c('Deletion_Events','Deletion_length')
    
    
    
    
    ## We need to correct the edits as unaligned sequenced in extremities are not counted as gaps.
    watson_best <- watson_best %>% 
      mutate(alignment.width = width,
             soft_trim = width - edits -matches) %>% 
      dplyr::select(-alignment.width) %>% 
      bind_cols(ins) %>% 
      bind_cols(del) %>% 
      mutate(indels = Insertion_length+ Deletion_length) %>% 
      mutate(edits = indels + soft_trim + mismatches) %>% 
      filter(edits <= max_edits)
    
    # update watson alignments  with updated number of edits
    watson_sub <- watson[watson_best$clusterID]
    
    
    
    # subset fasta file with watson clusters
    fasta_temp <- fasta[watson_best$clusterID]
    
    
    
    # get aligned gRNA and gDNA sequences with INDELS
    seq_gDNA <- alignedPattern(watson_sub)   # this step can be long
    seq_gRNA <- alignedSubject(watson_sub)   # this step can be long
    # add sequences to table
    watson_best$seq_gDNA <- seq_gDNA %>% as.character
    watson_best$seq_gRNA <- seq_gRNA %>% as.character
    
    
    
    
    # get hybridization energy
    energy <- CalculateEfficiencyArray(seq_gDNA,seq_gRNA,temp = 37,FA = 0)
    colnames(energy) <- c("Gibbs_Hybridization_efficacy","Gibbs_dG_0")
    # add to table
    watson_best <- watson_best %>%
      bind_cols(energy) 
    
    
    
    
    # format alignment
    watson_best <- watson_best %>% 
      rowwise() %>%
      mutate("Alignment" = paste(paste("gDNA :",seq_gDNA), paste("gRNA :", seq_gRNA),sep = "_"))
    
    
    
    
    
    watson_best$GC_content = letterFrequency(seq_gDNA, letters = "GC",as.prob = T)
    
    watson_best$mismatches_position_gRNA <- lapply(watson_sub@subject@mismatch,toString) %>% unlist 
    
    # add gRNA match start and end positions in cluster sequence
    watson_best <- watson_best %>% 
      bind_cols(data.frame(watson_sub@pattern@range))
    
    
    # add cluster sequence to table
    watson_best <- watson_best %>% 
      left_join(fasta_temp %>%
                  data.frame %>% 
                  select(sequence_window = ".") %>%
                  rownames_to_column('id') %>% 
                  separate("id", into = c("clusterID","cluster"),sep = "::",convert = T))
    
    
    
    # identify PAM sequence
    
    pams <- lapply(seq_along(watson_sub), function(x) {
      # if PAM is in 3' if gRNA, the PAM starts 1 nt after gRNA alignment end
      if(pam_side == "3"){
        start = watson_best$end[x]+1
        end = watson_best$end[x] + pam_length
        
        # get cluster fasta sequence length
        length = seqlengths(fasta[watson_best$clusterID[x]]) 
        # if PAM end position is not in fasta sequence add "..." as PAM
        if((end + pam_length)> length){
          DNAStringSet(paste(rep(".",pam_length),collapse = ""))
        } else {
          subseq(fasta[watson_best$clusterID[x]],start  , end ) # else extract PAM motif from coordinates
        }
      } else {
        # if PAM is 5' of gRNA, do the same stuff but the other side
        end = watson_best$start[x]-1
        start = watson_best$start[x] - pam_length
        length = seqlengths(fasta[watson_best$clusterID[x]])
        if((start - pam_length)< 1){
          DNAStringSet(paste(rep(".",pam_length),collapse = ""))
        } else {
          subseq(fasta[watson_best$clusterID[x]],start  , end )
          
        }
      }
    }
    )
    
    pams <- do.call(c,pams)
    
    watson_best$pam_gDNA <- pams %>% as.character
    watson_best$pam_gRNA <- pam %>%  as.character
    watson_best$pam_side <- pam_side
    watson_best <- watson_best %>%
      ungroup %>% 
      mutate(rank=row_number())
    
    
    ## get indels coordinates in alignment between gRNA and gDNA
    
    indels_list <- indel(watson_sub)
    
    indels_table <- bind_rows(
      
      
      bind_rows(
        "insertions" = insertion(indels_list)%>%
          data.frame,
        "deletions"= deletion(indels_list) %>%
          data.frame, .id = "indel") %>% 
        
        
        bind_rows(data.frame(
          indel = c("deletions","insertions"),
          group=-1,
          start=-1,
          end=-1,
          width = 1)))%>% # this line is here in case there are no indels to avoid missing columns later
      group_by(indel,group) %>% 
      summarise(start = toString(start),
                end = toString(end),
                width = toString(width)) %>% 
      pivot_wider(names_from = indel, values_from = c(start,end,width)) %>% 
      select(group,ends_with("deletions"),ends_with("insertions")) %>% 
      filter(group>0) %>%  
      arrange(group)
    
    
    watson_best <- watson_best %>% 
      left_join(indels_table, by = c("rank"="group"))
    
  } 
  
  
  ## CRICK strand analysis
  cat("Analyzing CRICK strand")
  
  crick_best <- align_stat %>% 
    filter(grna_orientation == "crick", crick_edits <= max_edits) %>% 
    dplyr::select(clusterID,cluster,grna_orientation,starts_with("crick")) %>%
    rename_all(~str_remove(.,"crick_"))
  
  if(nrow(crick_best)>0) {
    fasta_temp <- reverseComplement(fasta[crick_best$clusterID])
    
    crick_sub <- crick[crick_best$clusterID]
    
    width <- width(alignedPattern(crick_sub))
    
    indels <- nindel(crick_sub)
    ins <- indels@insertion %>% data.frame
    names(ins) <- c('Insertion_Events','Insertion_length')
    del <- indels@deletion %>% data.frame
    names(del) <- c('Deletion_Events','Deletion_length')
    
    ## We need to correct the edits as unaligned sequenced in extremities are not counted as gaps.
    crick_best <- crick_best %>% 
      mutate(alignment.width = width,
             soft_trim = width - edits -matches) %>% 
      dplyr::select(-alignment.width) %>% 
      bind_cols(ins) %>% 
      bind_cols(del) %>% 
      mutate(indels = Insertion_length+ Deletion_length) %>% 
      mutate(edits = indels + soft_trim + mismatches) %>% 
      filter(edits <= max_edits)
    
    crick_sub <- crick[crick_best$clusterID]
    
    seq_gDNA <- alignedPattern(crick_sub)
    seq_gRNA <- alignedSubject(crick_sub)
    
    energy <- CalculateEfficiencyArray(seq_gDNA,seq_gRNA,temp = 37,FA = 0)
    colnames(energy) <- c("Gibbs_Hybridization_efficacy","Gibbs_dG_0")
    
    crick_best$seq_gDNA <- seq_gDNA %>% as.character
    crick_best$seq_gRNA <- seq_gRNA %>% as.character
    crick_best <- crick_best %>% 
      rowwise() %>%
      mutate("Alignment" = paste(paste("gDNA :",seq_gDNA), paste("gRNA :", seq_gRNA),sep = "_"))
    
    crick_best <- crick_best %>%
      bind_cols(energy) 
    
    crick_best$GC_content = letterFrequency(seq_gDNA, letters = "GC",as.prob = T)
    
    crick_best$mismatches_position_gRNA <- lapply(crick_sub@subject@mismatch,toString) %>% unlist 
    
    crick_best <- crick_best %>% bind_cols(data.frame(crick_sub@pattern@range))
    
    crick_best <- crick_best %>% 
      left_join(fasta_temp %>% data.frame %>% select(sequence_window = ".") %>% rownames_to_column('id') %>% separate("id", into = c("clusterID","cluster"),sep = "::",convert = T))
    
    
    pams <- lapply(seq_along(crick_sub), function(x) {
      #print(x)
      if(pam_side == "3"){
        start = crick_best$end[x]+1
        end = crick_best$end[x] + pam_length
        length = seqlengths(fasta[crick_best$clusterID[x]])
        if((end + pam_length)> length){
          DNAStringSet(paste(rep(".",pam_length),collapse = ""))
        } else {
          subseq(reverseComplement(fasta[crick_best$clusterID[x]]),start  , end )
          
        }
      } else {
        end = crick_best$start[x]-1
        start = crick_best$start[x] - pam_length
        length = seqlengths(fasta[crick_best$clusterID[x]])
        if((start - pam_length)< 1){
          DNAStringSet(paste(rep(".",pam_length),collapse = ""))
        } else {
          subseq(reverseComplement(fasta[crick_best$clusterID[x]]),start,end )
          
        }
      }
    }
    )
    
    pams <- do.call(c,pams)
    
    
    crick_best$pam_gDNA <- pams %>% as.character
    crick_best$pam_gRNA <- pam %>%  as.character
    crick_best$pam_side <- pam_side
    crick_best <-crick_best %>% ungroup %>%  mutate(rank=row_number())
    
    indels_list <- indel(crick_sub)
    
    indels_table <- bind_rows(
      bind_rows(
        "insertions" = insertion(indels_list)%>%
          data.frame,
        "deletions"= deletion(indels_list) %>%
          data.frame, .id = "indel") %>% 
        bind_rows(data.frame(indel = c("deletions","insertions"),group=-1,start=-1,end=-1,width = 1))) %>% 
        group_by(indel,group = group) %>% 
        summarise(start = toString(start),
                  end = toString(end),
                  width = toString(width)) %>% 
      pivot_wider(names_from = indel, values_from = c(start,end,width)) %>% 
      select(group,ends_with("deletions"),ends_with("insertions")) %>% 
      filter(group>0) %>%  
      arrange(group) 
    
    
    crick_best <- crick_best %>% 
      left_join(indels_table, by = c("rank"="group"))
  }
  
  
  
  
  
  ## bind results from watson and crick alignments
  
  if(nrow(watson_best) > 0 & nrow(crick_best) > 0){
    best <- bind_rows(watson_best,crick_best)
    
  } else {
    if(nrow(watson_best)>0){
      best <- watson_best
    } else { 
      if(nrow(crick_best)>0){
        best <- crick_best
      } else { # if both datatables are empty :
        col.idx <- c("start" ,"end", "width", "clusterID", "chromosome", "start_chr",  "end_chr","sequence", "mismatches", "strand_guide","cut_gRNa_alignment","seq_gDNA","seq_gRNA")
        best <- data.frame(matrix(nrow = 0, ncol = length(col.idx) ))
        names(best) <- col.idx
      }
      
    }
  }
  
  rm(crick_best)
  rm(watson_best)
  
  
  
  if(nrow(best)>0){
    # reorder columns : 
    best <- best %>% 
      dplyr::select(clusterID, 
                    cluster, 
                    sequence_window,
                    grna_orientation,
                    seq_gDNA,seq_gRNA,
                    Alignment,
                    starts_with("Gibbs"),
                    GC_content,
                    alignment_start_gDNA = start,
                    alignment_end_gDNA = end,
                    alignment_width_gDNA = width,
                    alignment_score=score,
                    Identity_pct = pid,
                    N_edits = edits, 
                    N_mismatches = mismatches,
                    mismatches_position_gRNA,
                    soft_trim, 
                    n_indels=indels,
                    Insertion_length,
                    start_insertions,
                    end_insertions,
                    width_insertions,
                    Deletion_length,
                    start_deletions,
                    end_deletions,
                    width_deletions,
                    pam_gDNA,
                    pam_gRNA,
                    pam_side )
    
    
    best$pam_iupac <- sapply(best$pam_gDNA, function(x) {count_iupac_mismatches(x,pam)})
    best <- best %>% separate(pam_iupac, into = c("PAM_indel_count","PAM_indel_pos"),sep = "_",remove = T)
    
    
    
    
    ## get cutting site position based on CAS offset and strand orientation
    
    
    best <- best %>% 
      separate(cluster, into = c("chromosome","start_chr","end_chr"), sep = "[:]+|-",convert = T) %>% 
      mutate(chromosome = as.character(chromosome))
    
    if(pam_side == "3"){
      
      best <- best %>%
        mutate(cut_gRNa_alignment = case_when(grna_orientation == "watson" ~ start_chr + alignment_end_gDNA + offset + 1 ,
                                              TRUE ~ end_chr - alignment_end_gDNA - offset ))
    } else {
      best <- best %>%
        mutate(cut_gRNa_alignment = case_when(grna_orientation == "watson" ~ start_chr + alignment_start_gDNA + offset +1,
                                              TRUE ~ end_chr - alignment_start_gDNA - offset ))
    }
  } 
  
  


  
  

  

} else {
  # if none if the clusters has a gRNA match ...
  col.idx <- c("start" ,"end", "width", "clusterID", "chromosome", "start_chr",  "end_chr","sequence", "mismatches", "strand_guide","cut_gRNa_alignment","seq_gDNA","seq_gRNA")
  best <- data.frame(matrix(nrow = 0, ncol = length(col.idx) ))
  names(best) <- col.idx

}



rm(align_stat)
rm(crick)
rm(crick_sub)
rm(fasta_temp)
rm(watson)
rm(pam);rm(pams);rm(seq_gDNA); rm(seq_gRNA)


#----------------------------------------------------------#
# 06 - Aggregate cluster information and gRNA match ----
#----------------------------------------------------------#
clusters_df <- read.delim(args[2], sep="\t", header = T) %>% 
  mutate(chromosome = case_when(str_starts(chromosome, "chr") ~ as.character(chromosome),
                                TRUE ~ paste("chr",chromosome,sep="")))



best <- best %>% mutate(chromosome = case_when(str_starts(chromosome, "chr") ~ as.character(chromosome),
                                               TRUE ~ paste("chr",chromosome,sep="")))


clusters_df <- clusters_df %>%
  left_join(best, by = c("clusterID","chromosome"))



#----------------------------------------------------------#
# 07 - Add additional informations ----
#----------------------------------------------------------#


# calculate distance between cut sites and cluster gRNA matched theoretical cutting site

clusters_df <- clusters_df %>%
    mutate(relative_distance = cut_modal_position - cut_gRNa_alignment)  # distance between IS and predicted gRNA cut site

  
# add sequence of the cluster

seqs <- fasta %>% 
  data.frame  %>% 
  select(sequence_cluster = ".") %>%
  rownames_to_column('id')%>%
  separate("id", into = c("clusterID","cluster"),sep = "::",convert = T) %>% 
  select(clusterID, sequence_cluster)

clusters_df <- clusters_df %>% 
  left_join(seqs, by = "clusterID")


# annotate bulges
clusters_df <- clusters_df %>% 
  mutate(bulge = case_when(str_detect(seq_gDNA,"-") & str_detect(seq_gRNA,"-") ~ "both",
                           str_detect(seq_gRNA,"-") ~ "gDNA",
                           str_detect(seq_gDNA,"-") ~ "gRNA",
                           TRUE ~ "none"))


write.table(clusters_df,file = args[10], sep="\t", row.names = F, quote = F)

