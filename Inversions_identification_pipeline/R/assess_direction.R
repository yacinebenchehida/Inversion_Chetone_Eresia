assess_direction <- function(dt, pos_col_index = 9, keep_global = TRUE) {
  
  # Convert input to data.table for consistent indexing
  dt <- as.data.table(dt)
  
  # Extract the positions column based on provided index
  positions <- dt[[pos_col_index]]
  
  # Count number of positions
  n <- length(positions)
  
  # Stop if there are fewer than 6 points
  if (n < 6) stop("Not enough points to assess direction (need at least 6 rows).")
  
  # Extract first 3 positions
  first3 <- positions[1:3]
  # Compute median of first 3 positions for global trend comparison
  med_first3 <- median(first3)
  
  # Extract last 3 positions
  last3 <- positions[(n-2):n]
  # Compute median of last 3 positions for global trend comparison
  med_last3 <- median(last3)
  
  # Determine global direction based on median comparison
  if (med_first3 < med_last3) global_direction <- "positive"
  else if (med_first3 > med_last3) global_direction <- "negative"
  else global_direction <- "non-monotonous"
  
  # Check local monotony in first 3 points
  first_increasing <- all(diff(first3) > 0)
  first_decreasing <- all(diff(first3) < 0)
  
  # Check local monotony in last 3 points
  last_increasing <- all(diff(last3) > 0)
  last_decreasing <- all(diff(last3) < 0)
  
  # Decide overall direction combining global trend, local monotony, and keep_global option
  # 1) Both start and end monotonic in same direction
  if ((first_increasing && last_increasing) || (first_decreasing && last_decreasing)) {
    # a) Local trend matches global trend → keep global trend
    if ((first_increasing && global_direction == "positive") || (first_decreasing && global_direction == "negative")) {
      direction <- global_direction
    } else { # b) Local trend conflicts with global trend
      if (keep_global) {
        warning("Global and local trends are inconsistent; using global trend.")
        direction <- global_direction
      } else {
        direction <- "non-monotonous"
      }
    }
    
  } else if ((first_increasing || first_decreasing) && !(last_increasing || last_decreasing)) { # 2) Only start is monotonic, end is non-monotonic
    # a) Global trend matches start → keep global trend
    if ((first_increasing && global_direction == "positive") || (first_decreasing && global_direction == "negative")) {
      warning("Only the start of the window is monotonous; proceeding with caution.")
      direction <- global_direction
    } else { # b) Conflict between start and global trend
      if (keep_global) {
        warning("Global and local trends are inconsistent; using global trend.")
        direction <- global_direction
      } else {
        direction <- "non-monotonous"
      }
    }
    
  } else if ((first_increasing && last_decreasing) || (first_decreasing && last_increasing)) { # 3) Start and end monotonic but in opposite directions → use global trend
    
    direction <-  global_direction
    
  } else {
    # Both start and end non-monotonic → always non-monotonous
    direction <- global_direction
  }
  
  # Handle rare case where global trend is exactly flat (medians equal)
  if (global_direction == "non-monotonous") {
    warning("Global trend is flat; cannot determine direction.")
    direction <- "non-monotonous"
  }
  
  # Return the assessed direction
  return(direction)
}
