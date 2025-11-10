resolve_smallblock_bridge <- function(block_annotation, filtered_dt) {
  # Ensure data.table semantics if available
  if (!data.table::is.data.table(block_annotation)) block_annotation <- data.table::as.data.table(block_annotation)
  
  # Define the roles that count as events
  event_roles <- c("inversion", "inversion+translocation", "translocation")
  
  # If fewer than 3 rows, nothing to do
  if (nrow(block_annotation) < 3) return(block_annotation)
  
  # Build a small lookup from gene id to position (V1 -> V9)
  # Keep only the two columns for speed
  pos_map <- filtered_dt[, .(gene = V1, pos = V9)]
  
  # Slide a window of size 3 over rows 1..n
  for (i in 1:(nrow(block_annotation) - 2)) {
    
    # Row indices of the triplet
    i1 <- i
    i2 <- i + 1
    i3 <- i + 2
    
    # Extract roles
    r1 <- block_annotation$role[i1]
    r2 <- block_annotation$role[i2]
    r3 <- block_annotation$role[i3]
    
    # Quick gate: event, small_block with n_genes < 3, same event again
    if (!(r1 %in% event_roles)) next
    if (!(r2 == "small_block" && block_annotation$n_genes[i2] < 3)) next
    if (!(r3 == r1)) next
    
    #####################################
    # I) Fix end breakpoint first block # 
    #####################################
    # Get direction of the first event ("positive" or "negative")
    dir1 <- block_annotation$direction[i1]
    
    # Identify the last gene of the first event block
    end_gene_block1 <- block_annotation$end_gene[i1]
    
    # Identify the first gene of the small_block
    start_gene_small <- block_annotation$start_gene[i2]
    
    # Identify the last gene of the small_block
    end_gene_small <- block_annotation$end_gene[i2]
    
    # Map gene ids to positions using filtered_dt V9
    pos_end_block1 <- pos_map[gene == end_gene_block1, pos][1]
    pos_start_small <- pos_map[gene == start_gene_small, pos][1]

    # Check trend continuity:
    # positive: small_block first gene position must be greater than end position of first event
    # negative: small_block first gene position must be less than end position of first event
    follows_trend <- (dir1 == "positive" && pos_start_small > pos_end_block1) ||
      (dir1 == "negative" && pos_start_small < pos_end_block1)
    
    # If trend is followed, set the END breakpoint of the first event
    if (follows_trend) {
      block_annotation$end_breakpoint[i1] <- paste(start_gene_small, end_gene_small, sep = "_")
      block_annotation$end_gene <- start_gene_small
    }else{
      block_annotation$end_breakpoint[i1] <- paste(end_gene_block1, start_gene_small, sep = "_")
    }
    
    #########################################
    # II) Fix start breakpoint second block # 
    #########################################
    # Get direction of the first event ("positive" or "negative")
    dir2 <- block_annotation$direction[i3]
    
    # Identify the last gene of the first event block
    start_gene_block3 <- block_annotation$start_gene[i3]
    
    # Identify the first gene of the small_block
    end_gene_small <- block_annotation$end_gene[i2]
    
    # Identify the last gene of the small_block
    start_gene_small <- block_annotation$start_gene[i2]
    
    # Map gene ids to positions using filtered_dt V9
    pos_end_block3 <- pos_map[gene == start_gene_block3, pos][1]
    pos_end_small <- pos_map[gene == start_gene_small, pos][1]
    
    # Check trend continuity:
    follows_trend <- (dir2 == "positive" && start_gene_block3 > end_gene_small) ||
      (dir1 == "negative" && start_gene_block3 < end_gene_small)
    
    # If trend is followed, set the start breakpoint of the last event
    if (follows_trend) {
      block_annotation$start_breakpoint[i3] <- paste(start_gene_small, end_gene_small, sep = "_")
      block_annotation$start_gene[i3] <- start_gene_small
    }else{
      block_annotation$start_breakpoint[i3] <- paste(end_gene_small, start_gene_block3, sep = "_")
    }
    
  }
  
  # Return the modified table
  return(block_annotation)
}
