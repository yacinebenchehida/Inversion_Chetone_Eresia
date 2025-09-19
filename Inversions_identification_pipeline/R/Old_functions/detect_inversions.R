detect_inversions <- function(dt, pos_col_index = 9, gene_col_index = 1, min_consecutive = 3) {
  # Convert input to a data.table to allow fast indexing and operations
  dt <- as.data.table(dt)
  
  # Modify the gene names by prefixing them with their row number to ensure uniqueness
  # Example: "HMEL000013-RA" becomes "18_HMEL000013-RA"
  dt[[gene_col_index]] <- paste0(seq_len(nrow(dt)), "_", dt[[gene_col_index]])
  
  # Step 1: Assess the overall direction (positive or negative) using an external helper function
  # This function must return either "positive", "negative", or "non-monotonous"
  direction <- assess_direction(dt, pos_col_index = pos_col_index)
  
  # Step 2: If the window is non-monotonous, stop the analysis and return NULL
  if (direction == "non-monotonous") {
    message("Window is non-monotonous. No inversion detection performed.")  # Inform user
    return(NULL)  # Exit early if no clear direction
  } # End of check for non-monotonous direction
  
  # Step 3: Extract positions (column specified by pos_col_index)
  positions <- dt[[pos_col_index]]  # Store positions in a vector
  
  # Compute the difference between adjacent positions to detect changes
  adj_diff <- diff(positions)  # Difference between each consecutive position
  
  # Determine the expected sign based on overall direction:
  # - If direction is "positive", positions should increase, so expected sign is +1
  # - If direction is "negative", positions should decrease, so expected sign is -1
  expected_sign <- if (direction == "positive") 1 else -1
  
  # Step 4: Flag reverse steps
  # A reverse step occurs when the actual sign does NOT match the expected sign
  reverse_steps <- sign(adj_diff) != expected_sign
  
  # Step 5: Prepare an empty results table to store inversion boundaries
  results_inversion_table <- data.table(
    inversion_id = integer(),     # Unique ID for each inversion detected
    gene_position = integer(),    # Position index of the boundary gene in the input table
    gene_name = character(),      # Gene name at the boundary
    role = character()            # "start" or "end" of inversion
  )
  
  # Initialize variables to track runs of reverse steps
  run_start <- NULL    # Starting index of a reverse run
  run_length <- 0      # Length of the current reverse run
  inversion_id <- 1    # Counter to uniquely identify each inversion
  
  # Loop through each step to detect consecutive reverse runs
  for (i in seq_along(reverse_steps)) {
    # If this step is a reverse step, either start a new run or extend the current run
    if (reverse_steps[i]) {
      if (is.null(run_start)) run_start <- i  # Mark start of run if not already started
      run_length <- run_length + 1            # Increment run length
    } else {
      # If we reach a non-reverse step, check if the run length is long enough to count as inversion
      if (!is.null(run_start) && run_length >= min_consecutive) {
        # Add start and end genes of this inversion to the results table
        results_inversion_table <- rbind(
          results_inversion_table,
          data.table(
            inversion_id = inversion_id,
            gene_position = c(run_start, run_start + run_length),
            gene_name = c(dt[[gene_col_index]][run_start],
                          dt[[gene_col_index]][run_start + run_length]),
            role = c("start", "end")
          )
        )
        inversion_id <- inversion_id + 1  # Increment inversion ID for next detection
      }
      # Reset the run tracking variables
      run_start <- NULL
      run_length <- 0
    } # End of handling non-reverse step
  } # End of loop over reverse_steps
  
  # Final check at the end of the vector in case the last steps were a valid inversion run
  if (!is.null(run_start) && run_length >= min_consecutive) {
    results_inversion_table <- rbind(
      results_inversion_table,
      data.table(
        inversion_id = inversion_id,
        gene_position = c(run_start, run_start + run_length),
        gene_name = c(dt[[gene_col_index]][run_start],
                      dt[[gene_col_index]][run_start + run_length]),
        role = c("start", "end")
      )
    )
  } # End of final run check
  
  # Return the full table of inversion boundaries
  return(results_inversion_table)
} # End of detect_inversions function
