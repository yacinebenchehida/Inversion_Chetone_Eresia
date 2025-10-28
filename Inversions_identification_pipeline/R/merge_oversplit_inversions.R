merge_oversplit_inversions <- function(inv_table, dt, pos_col_index = 9, gene_col_index = 1, direction = "positive") {
  
  # Load data.table package for data.table operations
  library(data.table)
  
  # Ensure the inversion table is a data.table
  inv_table <- as.data.table(inv_table)
  
  # Ensure the full dataset is a data.table
  dt <- as.data.table(dt)
  
  # If the inversion table is NULL or has zero rows, return it immediately
  if (is.null(inv_table) || nrow(inv_table) == 0) return(inv_table)
  
  # Helper function to get the index of a gene in the dataset dt
  gene_idx <- function(gene) {
    if (is.na(gene) || length(gene) == 0) return(NA_integer_)  # Return NA if gene is missing
    i <- which(dt[[gene_col_index]] == gene)                    # Find positions of the gene in dt
    if (length(i) == 0) return(NA_integer_) else return(i[1])  # Return first match or NA if not found
  }
  
  # Get a sorted list of unique inversion IDs
  inv_ids <- sort(unique(inv_table$inversion_id))
  
  # Initialize output list to store merged inversions
  out <- list()
  
  # Initialize new inversion ID counter
  new_id <- 1
  
  # Initialize loop counter
  i <- 1
  
  # Loop over all inversion IDs
  while (i <= length(inv_ids)) {
    
    # Current inversion ID
    id1 <- inv_ids[i]
    
    # Subset the rows corresponding to the current inversion
    rows1 <- inv_table[inversion_id == id1]
    
    # Extract roles present in the current inversion
    roles1 <- rows1$role
    
    # Extract start gene if present
    start1 <- if ("start" %in% roles1) rows1$gene_name[rows1$role == "start"][1] else NA_character_
    
    # Extract end gene if present
    end1   <- if ("end" %in% roles1)   rows1$gene_name[rows1$role == "end"][1]   else NA_character_
    
    # Extract undefined gene if present
    undef1 <- if ("undefined" %in% roles1) rows1$gene_name[rows1$role == "undefined"][1] else NA_character_
    
    # Extract type of the inversion (should be single value)
    type1 <- unique(rows1$type)[1]
    
    # Flag to indicate if a merge occurred in this iteration
    merged_this_round <- FALSE
    
    # Check if a next inversion exists
    if (i < length(inv_ids)) {
      
      # Next inversion ID
      id2 <- inv_ids[i + 1]
      
      # Subset rows for the next inversion
      rows2 <- inv_table[inversion_id == id2]
      
      # Extract roles for next inversion
      roles2 <- rows2$role
      
      # Extract start gene for next inversion
      start2 <- if ("start" %in% roles2) rows2$gene_name[rows2$role == "start"][1] else NA_character_
      
      # Extract end gene for next inversion
      end2   <- if ("end" %in% roles2)   rows2$gene_name[rows2$role == "end"][1]   else NA_character_
      
      # Extract undefined gene for next inversion
      undef2 <- if ("undefined" %in% roles2) rows2$gene_name[rows2$role == "undefined"][1] else NA_character_
      
      # Extract type for next inversion
      type2 <- unique(rows2$type)[1]
      
      # Only consider merging if both inversions are of the same type
      if (type1 == type2) {
        
        # Determine position for orientation check
        pos1 <- if (!is.na(undef1)) dt[[pos_col_index]][gene_idx(undef1)] else dt[[pos_col_index]][gene_idx(end1)]
        
        # Determine position of the second inversion breakpoint
        pos2 <- if (!is.na(end2)) dt[[pos_col_index]][gene_idx(end2)] else dt[[pos_col_index]][gene_idx(start2)]
        
        # Compute indices in dt for gene gap calculation
        idx1_end <- if (!is.na(end1)) gene_idx(end1) else gene_idx(undef1)
        idx2_start <- if (!is.na(start2)) gene_idx(start2) else gene_idx(undef2)
        
        # Compute number of genes between the inversions
        gap <- idx2_start - idx1_end - 1L
        
        # Initialize merge condition to FALSE
        merge_condition <- FALSE
        
        # Check if positions are available and gap <= 2 to allow merging
        if (!is.na(pos1) && !is.na(pos2) && gap <= 2) {
          
          # Orientation check depending on direction
          if (direction == "positive") merge_condition <- pos2 < pos1   # Positive trend → inversion goes down
          if (direction == "negative") merge_condition <- pos2 > pos1   # Negative trend → inversion goes up
        }
        
        # If merge condition is met, determine merge type
        if (merge_condition) {
          
          # CASE 1: undefined + start/end → merge as undefined, breakpoint = end of second inversion
          if (!is.na(undef1) && (!is.na(start2) || !is.na(end2))) {
            chosen_bp <- if (!is.na(end2)) end2 else start2
            out[[length(out) + 1]] <- data.table(
              inversion_id = new_id,
              gene_name = chosen_bp,
              role = "undefined",
              type = type1
            )
            merged_this_round <- TRUE
          }
          
          # CASE 2: start/end + start/end → merge start of first with end of second
          else if (!is.na(start1) && !is.na(end1) && !is.na(start2) && !is.na(end2)) {
            out[[length(out) + 1]] <- data.table(
              inversion_id = new_id,
              gene_name = c(start1, end2),
              role = c("start", "end"),
              type = type1
            )
            merged_this_round <- TRUE
          }
          
          # CASE 3: start/end + undefined → merge as undefined, breakpoint = start of first inversion
          else if ((!is.na(start1) || !is.na(end1)) && !is.na(undef2)) {
            chosen_bp <- if (!is.na(start1)) start1 else end1
            out[[length(out) + 1]] <- data.table(
              inversion_id = new_id,
              gene_name = chosen_bp,
              role = "undefined",
              type = type1
            )
            merged_this_round <- TRUE
          }
        }
      }
    }
    
    # If merge occurred, increment new ID and skip next inversion
    if (merged_this_round) {
      new_id <- new_id + 1
      i <- i + 2L
    } else {
      # Otherwise, keep current inversion as is
      rows1$inversion_id <- new_id
      out[[length(out) + 1]] <- rows1[, .(inversion_id, gene_name, role, type)]
      new_id <- new_id + 1
      i <- i + 1L
    }
  }
  
  # Combine all merged/unmerged inversions into a single data.table
  result <- rbindlist(out, use.names = TRUE, fill = TRUE)
  
  # Sort output by inversion ID
  setorder(result, inversion_id)
  
  # Return final merged inversion table
  return(result[])
}
