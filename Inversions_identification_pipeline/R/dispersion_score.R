dispersion_score <- function(dt, pos_col_index = 9, chaos_excess_score = 1.5) {
  
  pos <- dt[[pos_col_index]]                   # extract position vector
  
  # Linear fit and r^2
  r2 <- summary(lm(pos ~ seq_along(pos)))$r.squared
  
  # Fraction of sign changes in second differences
  sgn_changes_frac <- mean(diff(sign(diff(pos))) != 0, na.rm = TRUE)
  
  # Longest monotonic run (normalized)
  d <- diff(pos)
  runs <- rle(sign(d))
  longest_run <- max(runs$lengths) + 1
  longest_run_ratio <- longest_run / length(pos)
  
  # Combined score
  score <- (1 - abs(r2)) + sgn_changes_frac + (1 - longest_run_ratio)
  
  # Boolean output
  high_chaos <- score > chaos_excess_score
  return(high_chaos)
}
