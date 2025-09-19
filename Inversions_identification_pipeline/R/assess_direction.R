assess_direction <- function(dt, pos_col_index = 9) {
  
  dt <- as.data.table(dt)
  
  # Extract the positions column
  positions <- dt[[pos_col_index]]
  n <- length(positions)
  
  if (n < 6) {
    stop("Not enough points to assess direction (need at least 6 rows).")
  }
  
  # First 3 points
  first3 <- positions[1:3]
  
  # Last 3 points
  last3 <- positions[(n-2):n]
  
  # Check monotony
  first_increasing <- all(diff(first3) > 0)
  first_decreasing <- all(diff(first3) < 0)
  last_increasing  <- all(diff(last3) > 0)
  last_decreasing  <- all(diff(last3) < 0)
  
  # Decide overall direction
  if ((first_increasing & last_increasing) | (first_decreasing & last_decreasing)) {
    direction <- ifelse(first_increasing, "positive", "negative")
  } else if(first_increasing | first_decreasing) {
    warning("Only the start of the window is monotonous; proceeding with caution.")
    direction <- ifelse(first_increasing, "positive", "negative")
  } else {
    direction <- "non-monotonous"
  }
  
  return(direction)
}
