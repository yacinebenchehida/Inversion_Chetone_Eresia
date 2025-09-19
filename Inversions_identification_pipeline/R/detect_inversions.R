detect_inversions <- function(dt_segment, pos_col_index = 9, gene_col_index = 1, min_consecutive = 3, direction) {
  dt_segment <- data.table::as.data.table(dt_segment)  # Ensure input is a data.table
  orig_genes <- as.character(dt_segment[[gene_col_index]])  # Store original gene names
  dt_segment[[gene_col_index]] <- orig_genes  # Keep original gene names
  
  positions <- dt_segment[[pos_col_index]]  # Extract positions
  adj_diff <- diff(positions)  # Compute differences between consecutive positions
  expected_sign <- if (direction == "positive") 1 else -1  # Determine expected sign
  reverse_steps <- sign(adj_diff) != expected_sign  # Identify steps opposite to expected direction
  
  results_inversion_table <- data.table::data.table(  # Initialize results table
    inversion_id = integer(), 
    gene_name = character(),
    role = character()
  )
  
  run_start <- NULL  # Initialize start index for a run
  run_length <- 0  # Initialize length counter for a run
  inversion_id <- 1  # Initialize inversion ID counter
  
  for (i in seq_along(reverse_steps)) {  # Loop through reverse steps
    if (reverse_steps[i]) {  # If current step is reversed
      if (is.null(run_start)) run_start <- i  # Start new run if none active
      run_length <- run_length + 1  # Increment run length
    } else {  # If current step follows expected direction
      if (!is.null(run_start) && run_length >= min_consecutive) {  # Check if run qualifies as inversion
        start_name <- dt_segment[[gene_col_index]][run_start]  # Get start gene of inversion
        end_name <- dt_segment[[gene_col_index]][run_start + run_length]  # Get end gene of inversion
        results_inversion_table <- rbind(  # Append inversion to results table
          results_inversion_table,
          data.table::data.table(
            inversion_id = inversion_id, 
            gene_name = c(start_name, end_name), 
            role = c("start", "end")
          ),
          use.names = TRUE
        )
        inversion_id <- inversion_id + 1  # Increment inversion ID
      }
      run_start <- NULL  # Reset run start
      run_length <- 0  # Reset run length
    }
  }
  
  if (!is.null(run_start) && run_length >= min_consecutive) {  # Final check for last run
    start_name <- dt_segment[[gene_col_index]][run_start]  # Start gene
    end_name <- dt_segment[[gene_col_index]][run_start + run_length]  # End gene
    results_inversion_table <- rbind(  # Append final inversion
      results_inversion_table,
      data.table::data.table(
        inversion_id = inversion_id,
        gene_name = c(start_name, end_name),
        role = c("start", "end")
      ),
      use.names = TRUE
    )
  }
  
  return(results_inversion_table)  # Return inversion table
} # End of detect_inversions_simple
