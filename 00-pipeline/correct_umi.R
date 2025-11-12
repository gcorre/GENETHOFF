# ========================================================= #
# Guillaume CORRE @ GENETHON 2025
# This script is used to correct UMIs based on Hamming distance and filter UMIs based on a motif.
# It is used in the snakemake pipeline for GUIDE-seq analysis.
# ========================================================= #


options(tidyverse.quiet = TRUE,warn = -1,verbose = F,warn = -10,conflicts.policy = list(warn = FALSE))

library(stringdist,quietly = T, verbose = F,warn.conflicts = F)
library(igraph,quietly = T, verbose = F,warn.conflicts = F)
library(Biostrings,quietly = T, verbose = F,warn.conflicts = F)
library(tidyverse,quietly = T, verbose = F,warn.conflicts = F)


#R.Version()$version.string == "R version 4.3.2 (2023-10-31 ucrt)"


args <- commandArgs(trailingOnly = T)
file <- args[1]
bed <- read.delim(file, header = F) %>% mutate(V1 = as.character(V1))
motif <- args[2]
filt.umi <- as.logical(str_to_upper(args[3]))
hamming_threshold <- as.numeric(args[4])
method <- str_to_lower(args[5])
output <- args[6]
rescueR2 <- as.logical(str_to_upper(args[7]))

## debug
# 
# file <- "C:\\Users\\gcorre\\OneDrive - Institut\\Bureau\\VEGFA_s1_K562_neg.reads_per_UMI_per_IS.bed"
# bed <- read.delim(file, header = F) %>% mutate(V1 = as.character(V1))
# motif <- "NNWNNWNN"
# filt.umi <- as.logical("TRUE")
# hamming_threshold <- as.numeric(1)
# method <- str_to_lower("Adjacency")
# output <- "C:\\Users\\gcorre\\OneDrive - Institut\\Bureau\\VEGFA_s1_K562_neg.reads_per_UMI_per_IS_corrected2.bed"
# rescueR2 <- TRUE



## If R2 reads were rescued from unmapped or too short pairs, merge them with PE read alignments
if(rescueR2==TRUE){
  rescuedR2_file <- str_replace(file,".reads_per_UMI_per_IS.bed","_R2rescued.reads_per_UMI_per_IS.bed")
  bed_rescuedR2 <- read.delim(rescuedR2_file,header=F) %>% 
    mutate(V1 = as.character(V1))
  
  bed <- bind_rows(bed,bed_rescuedR2)
}

names(bed) <- c('chr','start','end','orientation','UMI','strand', 'count','Qual')



# Select sites to correct UMIs -----------------------------------------------------------

bed <- bed %>%  
  unite("site", chr,start,end,strand,remove = F)

to_process <- bed %>% dplyr::count(site) %>% filter(n>1) %>% distinct(site) # get sites with multiple UMIs

bed_sub <- bed %>% semi_join(to_process)

sp <- split(bed_sub,f = paste(bed_sub$site, bed_sub$orientation,sep=";"))





# Define function to perform UMI correction -----------------------------------------
group_umis <- function(umis, counts, hamming_threshold = 1,graph=F,type = "Adjacency") {
  # Calculate Hamming distance matrix
  dist_matrix <- stringdistmatrix(umis, umis, method = "hamming")
  
  # Create adjacency matrix based on threshold
  adj_matrix <- dist_matrix <= hamming_threshold
  diag(adj_matrix) <- FALSE  # Remove self-loops
  
  # Create graph
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  if(graph){
    plot(g,vertex.color="black",vertex.label.cex=0.8,vertex.size = 1,edge.width = 5	)
  }
  # Find connected components
  components <- igraph::components(g)
  
  if(type=="cluster"){
    # Initialize results list
    results <- list()
    
    # Process each connected component as a single read group
    for (comp_id in 1:components$no) {
      # Get nodes in this component
      comp_nodes <- which(components$membership == comp_id)
      comp_umis <- umis[comp_nodes]
      comp_counts <- counts[comp_nodes]
      
      # Find UMI with highest count as representative
      rep_idx <- which.max(comp_counts)
      
      # Create group
      results[[length(results) + 1]] <- list(
        representative = comp_umis[rep_idx],
        members = comp_umis,
        member_counts = comp_counts,
        member_indices = comp_nodes,
        total_count = sum(comp_counts)
      )
    }
    #return(results)
    result_df <- lapply(results,function(x){
      
      data.frame(UMI_node=x$representative,UMI=x$members,row.names = NULL)
    }) %>% bind_rows()
    
    return(result_df)
  }
  
  if(type == "adjacency"){
    
    # Initialize results list
    results <- list()
    
    # Process each connected component
    for (comp_id in 1:components$no) {
      # Get nodes in this component
      comp_nodes <- which(components$membership == comp_id)
      comp_umis <- umis[comp_nodes]
      comp_counts <- counts[comp_nodes]
      
      # Skip if component has only one node
      if (length(comp_nodes) == 1) {
        results[[length(results) + 1]] <- list(
          representative = comp_umis,
          members = comp_umis,
          member_counts = comp_counts,
          member_indices = comp_nodes
        )
        next
      }
      
      # Create subgraph for this component
      subg <- induced_subgraph(g, comp_nodes)
      
      # Initialize unassigned nodes
      unassigned <- seq_along(comp_nodes)
      groups <- list()
      
      while (length(unassigned) > 0) {
        # Get unassigned nodes and their counts
        unassigned_counts <- comp_counts[unassigned]
        
        # Select top nodes by count
        n_select <- min(2, length(unassigned))
        top_nodes <- order(unassigned_counts, decreasing = TRUE)[1:n_select]
        selected_nodes <- unassigned[top_nodes]
        
        # Process each selected node
        for (node_idx in seq_along(selected_nodes)) {
          node <- selected_nodes[node_idx]
          if (!(node %in% unassigned)) next
          
          # Find all unassigned nodes one edge away
          neighbors <- neighbors(subg, which(comp_nodes == comp_nodes[node]))
          neighbor_indices <- comp_nodes[as.numeric(neighbors)]
          available_neighbors <- neighbor_indices[neighbor_indices %in% comp_nodes[unassigned]]
          
          # Create group
          group_nodes <- c(comp_nodes[node], available_neighbors)
          group_umis <- umis[group_nodes]
          group_counts <- counts[group_nodes]
          
          # Add to results if group is not empty
          if (length(group_nodes) > 0) {
            groups[[length(groups) + 1]] <- list(
              representative = umis[comp_nodes[node]],
              members = group_umis,
              member_counts = group_counts
            )
          }
          
          # Remove assigned nodes from unassigned
          unassigned <- unassigned[!(comp_nodes[unassigned] %in% group_nodes)]
        }
      }
      
      # Merge groups for this component into results
      results <- c(results, groups)
    }
    #
    result_df <- lapply(results,function(x){
      
      data.frame(UMI_node=x$representative,UMI=x$members,row.names = NULL)
    }) %>% bind_rows()
    
    return(result_df)
  }
}






# Perform correction ----------------------------------------------------------------

sp2 <- lapply(sp, function(x){
  
  umis <- x$UMI
  counts <- x$count
  
  group_umis(umis, counts, hamming_threshold = hamming_threshold, graph = F, type = method)
}
)

# collapse all
corrected <- sp2 %>% 
  bind_rows(.id="site") %>% 
  distinct() %>% 
  separate(site, into =c("site","orientation"),sep = ";")

# annotate UMI with their UMI node sequence and aggregate
w <- bed_sub %>% 
  left_join(corrected , by = c("site","orientation","UMI")) %>% 
  group_by(chr,start,end,UMI=UMI_node,strand,orientation) %>% 
  summarise(
    avg_qual = sum(count*Qual)/sum(count),
    UMIs = toString(UMI),
    UMIs_count = toString(count),
    count = sum(count),
    nUMIs = n()) %>% 
  ungroup


# Aggregate with single UMI sites ---------------------------------------------------

bed_corrected <- bed %>% 
  anti_join(to_process) %>% 
  bind_rows(w) %>% 
  select(-site) %>% 
  replace_na(list(nUMIs = 1)) %>% 
  mutate(UMIs  = case_when(is.na(UMIs )~ UMI, TRUE ~ UMIs ),
         UMIs_count  = case_when(is.na(UMIs_count )~ as.character(count), TRUE ~ UMIs_count ),
         avg_qual =case_when(is.na(avg_qual )~ Qual, TRUE ~ avg_qual ) ) %>% 
  select(-Qual)


# Check if UMI pattern is correct ---------------------------------------------------

has_motif <- vcountPattern(motif, DNAStringSet(bed_corrected$UMI), fixed = FALSE)

bed_corrected <- bed_corrected %>% 
  ungroup %>% 
  mutate("has_motif"= as.logical(has_motif))


if(filt.umi==TRUE){
  bed_corrected <- bed_corrected %>% filter(has_motif == T)
}

bed_corrected <- bed_corrected %>% 
  select(chr,start,end,UMI,orientation,strand,count,avg_qual,UMIs,UMIs_count,nUMIs,has_motif)

write.table(bed_corrected, file = output, sep="\t", quote = F, row.names = F)



