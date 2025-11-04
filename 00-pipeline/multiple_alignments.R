# ========================================================= #
# Guillaume CORRE @ GENETHON 2025
# This script is used to align clusters of sequences to a gRNA sequence and report the best matches.
# It is used in the snakemake pipeline for GUIDE-seq analysis.
# ========================================================= #


# Set global options
options(tidyverse.quiet = TRUE,warn = -1,verbose = F,warn = -10,conflicts.policy = list(warn = FALSE))


#R.Version()$version.string == "R version 4.3.2 (2023-10-31 ucrt)"

# Get arguments fro mcommand line 
args <- commandArgs(trailingOnly = T)


# Load necessary libraries --------------------------------------------------------------------

library(tidyverse,quietly = T, verbose = F,warn.conflicts = F)
library(Biostrings,quietly = T,verbose = F, warn.conflicts = F)
library(DECIPHER,quietly = T,warn.conflicts = F,verbose = F)
library(pwalign,quietly = T, verbose = F,warn.conflicts = F)


# for debugging
 # args <- c( "04-IScalling/Cpf1.cluster_slop.fa", "04-IScalling/Cpf1.cluster_slop.bed", "04-IScalling/Cpf1.UMIs_per_IS_in_Cluster.bed", "AAAAGAGAGAGAGAGGGGGGAAA", "SIPRa-S1",  "TTTN",  -4, 6, TRUE, "5", "05-Report/Cpf1.rdata")



# Parse arguments to objects ---------------------------------------------------------------


fasta <- readDNAStringSet(args[1],
                          use.names = T)

clusters <- read.delim(args[2],
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

bed_collapsedUMI <- read.delim(args[3],
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


bed_collapsedUMI <- bed_collapsedUMI %>% mutate(chromosome = as.character(chromosome),
                                                PCR_orientation = factor(PCR_orientation,levels = c("negative","positive")))




grna <- DNAString(args[4])
gRNA_name <- args[5]

grna@metadata$name <- gRNA_name

pam <- DNAString(args[6])
pam_length <- nchar(pam)


max_edits <- as.numeric(args[8])

bulges <- toupper(args[9])

pam_side <- args[10]

offset <- as.numeric(args[7])

## check that offset is in the good direction
## if PAM is in the 3', offset should be negative and opposite
offset <- if((pam_side == "3" & sign(offset) < 0) |(pam_side == "5" & sign(offset) > 0) ){
  offset <- offset
} else {
  offset <- offset * -1 # change sign
  
}

# Change gap penalty if bulges are tolerated or not

gap_open_penalty <- 10 ## default (this tolerate gap (bulges))
if(bulges == "FALSE"){
  gap_open_penalty <- gap_open_penalty *100  # this forces alignment without gap
}




# Define matching function with IUPAC letters ---------------------------------------
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




## Calculate pairwise alignments between cluster sequence and gRNA sequence in forward and reverse orientation
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



if(nrow(align_stat)>0){
  
  # find strand with best score
  align_stat <- align_stat %>%
    mutate(grna_orientation = case_when(watson_score > crick_score ~ "watson",
                                        watson_score < crick_score  ~ "crick",
                                        TRUE ~ NA))
  
  
  
  
  
  ## WATSON strand analysis
  cat("Analyzing WATSON strand")
  
  watson_best <- align_stat %>% 
    filter(grna_orientation == "watson", watson_edits <= max_edits) %>% ### keep clusters with less than n edits
    dplyr::select(clusterID,cluster,grna_orientation,starts_with("watson")) %>%
    rename_all(~str_remove(.,"watson_")) 
  
  
  
  
  
  
  if(nrow(watson_best)>0){
    
    watson_sub <- watson[watson_best$clusterID]
    
    width <- width(alignedPattern(watson_sub))
    
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
    
    watson_sub <- watson[watson_best$clusterID]
    
    
    
    
    
    fasta_temp <- fasta[watson_best$clusterID]
    seq_gDNA <- alignedPattern(watson_sub)   # this step can be long
    seq_gRNA <- alignedSubject(watson_sub)   # this step can be long
    
    energy <- CalculateEfficiencyArray(seq_gDNA,seq_gRNA,temp = 37,FA = 0)
    colnames(energy) <- c("Gibbs_Hybridization_efficacy","Gibbs_dG_0")
    
    
    watson_best$seq_gDNA <- seq_gDNA %>% as.character
    watson_best$seq_gRNA <- seq_gRNA %>% as.character
    watson_best <- watson_best %>% 
      rowwise() %>%
      mutate("Alignment" = paste(paste("gDNA :",seq_gDNA), paste("gRNA :", seq_gRNA),sep = "\r\n"))
    watson_best <- watson_best %>%
      bind_cols(energy) 
    watson_best$GC_content = letterFrequency(seq_gDNA, letters = "GC",as.prob = T)
    watson_best$mismatches_position_gRNA <- lapply(watson_sub@subject@mismatch,toString) %>% unlist 
    watson_best <- watson_best %>% bind_cols(data.frame(watson_sub@pattern@range))
    watson_best <- watson_best %>% 
      left_join(fasta_temp %>%
                  data.frame %>% select(sequence_window = ".") %>%
                  rownames_to_column('id') %>% 
                  separate("id", into = c("clusterID","cluster"),sep = "::",convert = T))
    
    
    
    pams <- lapply(seq_along(watson_sub), function(x) {
      #print(x)
      if(pam_side == "3"){
        start = watson_best$end[x]+1
        end = watson_best$end[x] + pam_length
        length = seqlengths(fasta[watson_best$clusterID[x]])
        if((end + pam_length)> length){
          DNAStringSet(paste(rep(".",pam_length),collapse = ""))
        } else {
          subseq(fasta[watson_best$clusterID[x]],start  , end )
          
        }
      } else {
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
    watson_best <- watson_best %>% ungroup %>%  mutate(rank=row_number())
    
    
    indels_list <- indel(watson_sub)
    
    indels_table <- bind_rows(
      bind_rows(
        "insertions" = insertion(indels_list)%>%
          data.frame,
        "deletions"= deletion(indels_list) %>%
          data.frame, .id = "indel") %>% 
        bind_rows(data.frame(indel = c("deletions","insertions"),group=-1,start=-1,end=-1,width = 1)) %>% # this line is here in case there are no indels to avoid missing columns later
        group_by(indel,group) %>% 
        summarise(start = toString(start),
                  end = toString(end),
                  width = toString(width))) %>% 
      pivot_wider(names_from = indel, values_from = c(start,end,width)) %>% 
      select(group,ends_with("deletions"),ends_with("insertions")) %>% 
      filter(group>0) %>%  
      arrange(group)
    
    
    watson_best <- watson_best %>% left_join(indels_table, by = c("rank"="group"))
    
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
      mutate("Alignment" = paste(paste("gDNA :",seq_gDNA), paste("gRNA :", seq_gRNA),sep = "\r\n"))
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
        bind_rows(data.frame(indel = c("deletions","insertions"),group=-1,start=-1,end=-1,width = 1)) %>% 
        group_by(indel,group = group) %>% 
        summarise(start = toString(start),
                  end = toString(end),
                  width = toString(width))) %>% 
      pivot_wider(names_from = indel, values_from = c(start,end,width)) %>% 
      select(group,ends_with("deletions"),ends_with("insertions")) %>% 
      filter(group>0) %>%  
      arrange(group) 
    
    
    crick_best <- crick_best %>% left_join(indels_table, by = c("rank"="group"))
  }
  
  
  
  if(nrow(watson_best)>0 & nrow(crick_best)>0){
    # collapse watson and crick
    best <- bind_rows(watson_best,crick_best)
  } else {
    if(nrow(watson_best)>0){
      best <- watson_best
    } else {
      best <- crick_best
    }
  }
  
  
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
                  pam_gDNA,pam_gRNA,pam_side )
  
  
  best$pam_iupac <- sapply(best$pam_gDNA, function(x) {count_iupac_mismatches(x,pam)})
  best <- best %>% separate(pam_iupac, into = c("PAM_indel_count","PAM_indel_pos"),sep = "_",remove = T)
  
  
  
  ## get cutting site position based on CAS offset and strand orientation
  
  if(nrow(best)>0){
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
  } else{
    
    col.idx <- c("start" ,"end", "width", "clusterID", "chromosome", "start_chr",  "end_chr","sequence", "mismatches", "strand_guide","cut_gRNa_alignment")
    best <- data.frame(matrix(nrow = 0, ncol = length(col.idx) ))
    names(best) <- col.idx 
  }
  
  
  
  
  ## Breakdown positive and negative PCRs
  # PCR orientation is a factor with the 2 levels. If one level is missing, it is still reported with value 0.
  
  
  cluster_PCR <- bed_collapsedUMI %>% 
    mutate(chromosome = as.character(chromosome)) %>% 
    select(-c(start_IS:end_IS,UMI_list:ReadPerUMI)) %>%
    group_by(clusterID,chromosome,PCR_orientation,.drop = F) %>% 
    summarise(N_IS = n_distinct(IS_ID),
              N_UMI = sum(N_UMI_IS),
              N_Reads = sum(Nreads_IS),
              N_Orientation = n_distinct(strand_IS),
              MapQ = median(MedianMAPQ_IS)
    ) %>%
    pivot_wider(names_from = PCR_orientation, values_from = c(N_IS:MapQ), values_fill = 0)
  
  
  
  ## annotate clusters with identified gRNA sequence
  cluster_annotated <- clusters %>% mutate(chromosome = as.character(chromosome)) %>% 
    left_join(cluster_PCR, by = c("clusterID","chromosome")) %>% 
    left_join(best %>% 
                select(-start_chr,-end_chr,-alignment_width_gDNA,-alignment_start_gDNA,-alignment_end_gDNA) ,by = c("clusterID","chromosome"))
  
  
  
  # calculate distance between cut sites and cluster gRNA matched theoretical cutting site 
  IS_in_clusters = bed_collapsedUMI %>% 
    left_join(cluster_annotated %>%  
                select(clusterID,chromosome,cut_gRNa_alignment), 
              by = c("clusterID","chromosome")) 
  
  IS_in_clusters = IS_in_clusters %>% 
    mutate(relative_distance = start_IS - cut_gRNa_alignment)  # distance between IS and predicted gRNA cut site
  
  
  # get the position with most abundant UMIs in a cluster
  modal_cut_position <- IS_in_clusters %>% 
    group_by(clusterID,chromosome,start_IS) %>%
    summarise(count_UMI = sum(N_UMI_IS)) %>% 
    mutate(UMI_proportion = round(count_UMI / sum(count_UMI)*100,digits = 1)) %>% 
    group_by(clusterID,chromosome) %>% 
    slice_max(n = 1, order_by = UMI_proportion,with_ties = F) %>%    ### here pick the first if multiple modes (especialty when low number of UMI and number of IS)
    mutate(cut_modal_position = start_IS)
  
  
  
  clusters_grna_match <- cluster_annotated %>%
    left_join(modal_cut_position, by = c("clusterID","chromosome"))
  
  save(clusters_grna_match, grna,  file = args[11])
} else {
  save(grna, file = args[11])
  
}


