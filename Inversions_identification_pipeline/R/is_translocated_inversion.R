is_translocated_inversion <- function(direction, inv_pos, first_positions, last_positions) {
  # Ensure inputs are numeric
  inv_pos <- as.numeric(inv_pos)
  first_positions <- as.numeric(first_positions)
  last_positions <- as.numeric(last_positions)
  
  if(length(first_positions) == 0 | length(last_positions) == 0) {
    stop("first_positions and last_positions must have at least one value each")
  }
  
  # Check for translocation depending on direction
  if(direction == "positive") {
    # Positive direction: translocated if inv_pos is smaller than the first positions
    # or larger than the last positions
    if(inv_pos < min(first_positions) | inv_pos > max(last_positions)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
    
  } else if(direction == "negative") {
    # Negative direction: translocated if inv_pos is larger than the first positions
    # or smaller than the last positions
    if(inv_pos > max(first_positions) | inv_pos < min(last_positions)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
    
  } else {
    stop("direction must be either 'positive' or 'negative'")
  }
}
