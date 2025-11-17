
is_big_enough <- function(dt, pos_col_index = 9, threshold = 10000000) {
  
  dt <- as.data.table(dt)
  
  start_val <- dt[[pos_col_index]][1L]
  end_val   <- dt[[pos_col_index]][nrow(dt)]
  
  estimated_genome_size <- end_val - start_val
  print(paste("The genome size is around", estimated_genome_size))
  
  if (estimated_genome_size > threshold) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
