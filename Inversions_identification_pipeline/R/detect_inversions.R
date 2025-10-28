#------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------
# Function to detect inversion based on the direction of successive points
#-------------------------------------------------------------------------
detect_inversions <- function(dt_segment, pos_col_index = 9, gene_col_index = 1, min_consecutive = 3, direction) {
  
  dt_segment <- as.data.table(dt_segment)  # Ensure input is a data.table
  orig_genes <- as.character(dt_segment[[gene_col_index]])  # Store original gene names
  dt_segment[[gene_col_index]] <- orig_genes  # Keep original gene names
  
  positions <- dt_segment[[pos_col_index]]  # Extract positions
  expected_sign <- if (direction == "positive") 1 else -1  # Determine expected direction
  adj_diff <- diff(positions)  # Compute differences between consecutive positions
  reverse_steps <- sign(adj_diff) != expected_sign  # Identify steps opposite to expected direction
  
  results_inversion_table <- data.table(inversion_id = integer(), gene_name = character(), role = character())  # Initialize results table
  run_start <- NULL  # Initialize start index for current run
  run_length <- 0  # Initialize run length counter
  inversion_id <- 1  # Initialize inversion counter
  
  for (i in seq_along(reverse_steps)) {  # Loop over reverse steps
    if (reverse_steps[i]) {  # If step is reversed
      if (is.null(run_start)) run_start <- i  # Start new run if none active
      run_length <- run_length + 1  # Increment run length
    } else {  # If step returns to expected direction
      if (!is.null(run_start) && run_length >= min_consecutive) {  # Check if run qualifies as inversion
        inversion_start_idx <- run_start  # Default start index of inversion
        inversion_end_idx <- run_start + run_length  # Default end index of inversion
        
        # Determine which side touches edge
        if (inversion_start_idx == 1) {  # Starts at first row → ignore start, mark end as undefined
          results_inversion_table <- rbind(results_inversion_table,
                                           data.table(
                                             inversion_id = inversion_id,
                                             gene_name = dt_segment[[gene_col_index]][inversion_end_idx],  # Only end row
                                             role = "undefined"  # Mark as undefined
                                           ),
                                           use.names = TRUE)
        } else if (inversion_end_idx == nrow(dt_segment)) {  # Ends at last row → ignore end, mark start as undefined
          results_inversion_table <- rbind(results_inversion_table,
                                           data.table(
                                             inversion_id = inversion_id,
                                             gene_name = dt_segment[[gene_col_index]][inversion_start_idx],  # Only start row
                                             role = "undefined"  # Mark as undefined
                                           ),
                                           use.names = TRUE)
        } else {  # Normal case, neither side touches edge
          results_inversion_table <- rbind(results_inversion_table,
                                           data.table(
                                             inversion_id = inversion_id,
                                             gene_name = c(dt_segment[[gene_col_index]][inversion_start_idx],
                                                           dt_segment[[gene_col_index]][inversion_end_idx]),
                                             role = c("start", "end")  # Standard roles
                                           ),
                                           use.names = TRUE)
        }
        
        inversion_id <- inversion_id + 1  # Increment inversion counter
      }
      run_start <- NULL  # Reset run
      run_length <- 0  # Reset length
    }
  }
  
  # Final check if last run reaches end
  if (!is.null(run_start) && run_length >= min_consecutive) {
    inversion_start_idx <- run_start  # Default start index
    inversion_end_idx <- run_start + run_length  # Default end index
    
    # Edge handling
    if (inversion_start_idx == 1) {  # Starts at first row → ignore start, mark end as undefined
      results_inversion_table <- rbind(results_inversion_table,
                                       data.table(
                                         inversion_id = inversion_id,
                                         gene_name = dt_segment[[gene_col_index]][inversion_end_idx],  # Only end row
                                         role = "undefined"  # Mark as undefined
                                       ),
                                       use.names = TRUE)
    } else if (inversion_end_idx == nrow(dt_segment)) {  # Ends at last row → ignore end, mark start as undefined
      results_inversion_table <- rbind(results_inversion_table,
                                       data.table(
                                         inversion_id = inversion_id,
                                         gene_name = dt_segment[[gene_col_index]][inversion_start_idx],  # Only start row
                                         role = "undefined"  # Mark as undefined
                                       ),
                                       use.names = TRUE)
    } else {  # Normal case
      results_inversion_table <- rbind(results_inversion_table,
                                       data.table(
                                         inversion_id = inversion_id,
                                         gene_name = c(dt_segment[[gene_col_index]][inversion_start_idx],
                                                       dt_segment[[gene_col_index]][inversion_end_idx]),
                                         role = c("start", "end")
                                       ),
                                       use.names = TRUE)
    }
  }
  
  return(results_inversion_table)  # Return final inversion table
}
