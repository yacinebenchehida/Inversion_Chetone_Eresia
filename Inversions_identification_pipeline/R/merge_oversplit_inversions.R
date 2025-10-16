merge_oversplit_inversions <- function(inv_table, dt, pos_col_index = 9, gene_col_index = 1, direction) {  # Define function with inputs: inversion table, full dataset, column indices, and inversion direction
  if (is.null(inv_table) || nrow(inv_table) == 0) return(inv_table)  # If inversion table is NULL or empty, return it immediately
  
  merged <- list()  # Initialize empty list to store merged inversions
  inversion_counter <- 1  # Initialize counter for new inversion IDs
  skip_next <- FALSE  # Initialize flag to skip next inversion if it was merged
  
  inv_ids <- unique(inv_table$inversion_id)  # Extract unique inversion IDs from the input table
  
  for (i in seq_along(inv_ids)) {  # Loop over each inversion ID
    if (skip_next) {  # Check if previous iteration merged the next inversion
      skip_next <- FALSE  # Reset skip flag
      next  # Skip current iteration
    }
    
    id1 <- inv_ids[i]  # Current inversion ID
    start_gene1 <- inv_table$gene_name[inv_table$inversion_id == id1 & inv_table$role == "start"]  # Get start gene of current inversion
    end_gene1 <- inv_table$gene_name[inv_table$inversion_id == id1 & inv_table$role == "end"]  # Get end gene of current inversion
    type1 <- inv_table$type[inv_table$inversion_id == id1 & inv_table$role == "start"]  # Get type of current inversion
    
    if (i < length(inv_ids)) {  # Check if a next inversion exists
      id2 <- inv_ids[i + 1]  # Next inversion ID
      start_gene2 <- inv_table$gene_name[inv_table$inversion_id == id2 & inv_table$role == "start"]  # Get start gene of next inversion
      end_gene2 <- inv_table$gene_name[inv_table$inversion_id == id2 & inv_table$role == "end"]  # Get end gene of next inversion
      type2 <- inv_table$type[inv_table$inversion_id == id2 & inv_table$role == "start"]  # Get type of next inversion
      
      if (type1 == type2) {  # Only attempt merge if both inversions have the same type
        idx1_end <- which(dt[[gene_col_index]] == end_gene1)  # Find index of end gene of first inversion in dataset
        idx2_start <- which(dt[[gene_col_index]] == start_gene2)  # Find index of start gene of second inversion in dataset
        gap_genes <- idx2_start - idx1_end - 1  # Compute number of genes between the inversions
        
        if (gap_genes <= 2) {  # Check if gap is small enough to consider merging
          pos1 <- dt[[pos_col_index]][idx1_end]  # Position of end gene of first inversion
          pos2 <- dt[[pos_col_index]][idx2_start]  # Position of start gene of second inversion
          
          merge_condition <- if (direction == "positive") pos2 < pos1 else pos2 > pos1  # Determine if positions indicate a merged inversion
          
          if (merge_condition) {  # If condition met, merge the two inversions
            merged[[inversion_counter]] <- data.table(  # Append merged inversion to list
              inversion_id = inversion_counter,  # Assign new inversion ID
              gene_name = c(start_gene1, end_gene2),  # Start gene is first inversion start, end gene is second inversion end
              role = c("start", "end"),  # Roles of genes
              type = type1  # Type remains the same
            )
            inversion_counter <- inversion_counter + 1  # Increment new inversion counter
            skip_next <- TRUE  # Set flag to skip the next inversion in loop
            next  # Skip rest of loop for this iteration
          }
        }
      }
    }
    
    merged[[inversion_counter]] <- data.table(  # If no merge, keep current inversion
      inversion_id = inversion_counter,  # Assign new inversion ID
      gene_name = c(start_gene1, end_gene1),  # Start and end genes of current inversion
      role = c("start", "end"),  # Roles of genes
      type = type1  # Type remains the same
    )
    inversion_counter <- inversion_counter + 1  # Increment inversion counter
  }
  
  return(data.table::rbindlist(merged))  # Combine list of merged inversions into single data.table and return
}
