remove_too_fuzzy_breakpoints <- function(block_annotation, gene_order, max_genes_dist = 3) {
  block_annotation <- data.table::as.data.table(block_annotation)  # ensure data.table
  
  if (is.null(block_annotation) || nrow(block_annotation) == 0) {  # if no rows, return as is
    return(block_annotation)                                       # return input unchanged
  }
  
  if (missing(gene_order) || is.null(gene_order) || length(gene_order) == 0) {  # require gene_order
    stop("gene_order must be provided and must contain at least one gene.")     # error if missing
  }
  
  if (is.character(gene_order) && length(gene_order) == 1 && file.exists(gene_order)) {  # if gene_order is a file path
    gene_order <- data.table::fread(gene_order, header = FALSE)[[1]]                    # read file as vector
  }
  
  gene_order <- as.character(gene_order)                                # coerce gene_order to character vector
  gene_pos <- setNames(seq_along(gene_order), gene_order)               # map gene name to index in gene_order
  
  event_roles <- c("inversion", "inversion+translocation", "translocation")  # roles for which breakpoints matter
  
  keep_row <- rep(TRUE, nrow(block_annotation))                          # logical index of rows to keep
  
  for (i in seq_len(nrow(block_annotation))) {                           # loop over each row of block_annotation
    if (!block_annotation$role[i] %in% event_roles) {                    # skip non event roles
      next                                                               
    }
    
    bp_vec <- c(block_annotation$start_breakpoint[i],                    # collect start breakpoint
                block_annotation$end_breakpoint[i])                      # and end breakpoint
    
    for (bp in bp_vec) {                                                 # loop over the two breakpoints
      if (is.na(bp) || bp %in% c("Undefined", "Not_applicable")) {       # ignore undefined or not applicable
        next
      }
      
      genes <- strsplit(bp, "_", fixed = TRUE)[[1]]                      # split breakpoint string into two genes
      if (length(genes) != 2) {                                          # skip malformed breakpoint strings
        next
      }
      
      g1 <- genes[1]                                                     # first gene at breakpoint
      g2 <- genes[2]                                                     # second gene at breakpoint
      
      if (!(g1 %in% names(gene_pos)) || !(g2 %in% names(gene_pos))) {    # if any gene not in gene_order, skip
        next
      }
      
      idx1 <- gene_pos[[g1]]                                             # index of first gene in gene_order
      idx2 <- gene_pos[[g2]]                                             # index of second gene in gene_order
      
      dist_genes <- abs(idx1 - idx2) - 1L                                # number of genes between them
                                                                        # (difference minus one)
      if (!is.na(dist_genes) && dist_genes > max_genes_dist) {           # if distance exceeds user threshold
        keep_row[i] <- FALSE                                             # mark this block row for removal
        break                                                            # no need to inspect other breakpoint
      }
    }
  }
  
  return(block_annotation[keep_row])                                     # return filtered block_annotation
}
