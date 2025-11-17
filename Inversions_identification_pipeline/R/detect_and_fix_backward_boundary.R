detect_and_fix_backward_boundary <- function(dt, pos_col_index, gene_col_index, block_data) {
  ann <- copy(block_data)                                            # copy block annotation so we do not modify in place
  
  for (i in seq_len(nrow(ann) - 1)) {                                # loop over each pair of consecutive blocks
    
    g1 <- unlist(ann$genes[[i]])                                     # extract gene vector for block i
    g2 <- unlist(ann$genes[[i + 1]])                                 # extract gene vector for block i+1
    
    if (length(g1) < 2 || length(g2) < 2) next                       # skip if any block has fewer than 2 genes (cannot move boundary)
    
    # compute interval for block 1
    g1_pos <- dt[[pos_col_index]][match(g1, dt[[gene_col_index]])]
    g1_low  <- median(g1_pos[1:min(3, length(g1_pos))])
    g1_high <- median(g1_pos[max(1, length(g1_pos)-2):length(g1_pos)])
    #
    ## compute interval for block 2
    g2_pos <- dt[[pos_col_index]][match(g2, dt[[gene_col_index]])]
    g2_low  <- median(g2_pos[1:min(3, length(g2_pos))])
    g2_high <- median(g2_pos[max(1, length(g2_pos)-2):length(g2_pos)])
    #
    ## skip if intervals do NOT overlap
    if (g1_high < g2_low || g2_high < g1_low) next
    
    all_genes <- c(g1, g2)                                           # concatenate genes of both blocks in current block order
    dt_pair <- dt[dt[[gene_col_index]] %in% all_genes]               # subset original table to these genes
    dt_pair <- dt_pair[match(all_genes, dt_pair[[gene_col_index]])]  # reorder subset so that rows follow all_genes order exactly
    
    pos_pair <- dt_pair[[pos_col_index]]                             # extract position vector for this pair
    genes_pair <- dt_pair[[gene_col_index]]                          # extract gene id vector for this pair
    n_pair <- length(genes_pair)                                     # number of genes in the pair
    
    orig_boundary <- length(g1)                                      # index of the current breakpoint (last gene of block i) in genes_pair
    
    score_vec <- rep(Inf, 5)                                         # initialize score vector for up to 5 shifts to the left
    
    for (shift in 0:4) {                                             # loop over candidate shifts 0..4 to the left
      
      idx_bef <- orig_boundary - shift                               # index of gene just before candidate breakpoint
      idx_aft <- idx_bef + 1                                         # index of gene just after candidate breakpoint
      
      if (idx_bef - 5 < 1) next                                      # skip if we cannot take 5 left neighbors
      if (idx_aft + 5 > n_pair) next                                 # skip if we cannot take 5 right neighbors
      
      pos_bef <- pos_pair[idx_bef]                                   # position of left gene at candidate breakpoint
      pos_aft <- pos_pair[idx_aft]                                   # position of right gene at candidate breakpoint
      
      left_neighbors  <- pos_pair[(idx_bef - 5):(idx_bef - 1)]       # positions of 5 genes immediately left of pos_bef
      right_neighbors <- pos_pair[(idx_aft + 1):(idx_aft + 5)]       # positions of 5 genes immediately right of pos_aft
      
      left_mean  <- mean(abs(left_neighbors - pos_bef))              # mean absolute distance between left neighbors and pos_bef
      right_mean <- mean(abs(right_neighbors - pos_aft))             # mean absolute distance between right neighbors and pos_aft
      
      score_vec[shift + 1] <- left_mean + right_mean                 # store total score for this shift (lower is better)
    }
    
    best_idx <- which.min(score_vec)                                 # index of minimal score in score_vec (1..5)
    if (!is.finite(score_vec[best_idx])) next                        # skip if no finite candidate score was found
    
    best_shift <- best_idx - 1                                       # convert back to shift value (0..4)
    
    cut1 <- orig_boundary - best_shift                               # new cut position in genes_pair (after this index is breakpoint)
    if (cut1 < 1 || cut1 >= n_pair) next                             # safety check to avoid empty blocks
    
    g1_new <- genes_pair[1:cut1]                                     # new gene set for block i (left side)
    g2_new <- genes_pair[(cut1 + 1):n_pair]                          # new gene set for block i+1 (right side)
    
    ann$genes[[i]]     <- g1_new                                     # update genes for block i
    ann$genes[[i + 1]] <- g2_new                                     # update genes for block i+1
    
    ann$n_genes[i]     <- length(g1_new)                             # update gene count for block i
    ann$n_genes[i + 1] <- length(g2_new)                             # update gene count for block i+1
    
    ann$end_gene[i]           <- g1_new[length(g1_new)]              # update end gene for block i
    ann$start_gene[i + 1]     <- g2_new[1]                           # update start gene for block i+1
    
    ann$median_pos[i]     <- median(dt[dt[[gene_col_index]] %in% g1_new][[pos_col_index]])  # update median position for block i
    ann$median_pos[i + 1] <- median(dt[dt[[gene_col_index]] %in% g2_new][[pos_col_index]])  # update median position for block i+1
    
    boundary_left_gene  <- g1_new[length(g1_new)]                    # gene id at left of corrected breakpoint
    boundary_right_gene <- g2_new[1]                                 # gene id at right of corrected breakpoint
    
    ann$end_breakpoint[i]     <- paste0(boundary_left_gene, "_", boundary_right_gene)  # update end breakpoint of block i
    ann$start_breakpoint[i+1] <- paste0(boundary_left_gene, "_", boundary_right_gene)  # update start breakpoint of block i+1
  }
  
  return(ann)                                                        # return updated block annotation
}
