fix_small_edge_blocks <- function(block_annotation) {
  
  # Check first row
  if (block_annotation$role[1] %in% c("inversion", "inversion+translocation", "translocation") &&
      block_annotation$n_genes[1] < 4) {
    block_annotation$role[1] <- "unclear"
  }
  
  # Check last row
  last <- nrow(block_annotation)
  if (block_annotation$role[last] %in% c("inversion", "inversion+translocation", "translocation") &&
      block_annotation$n_genes[last] < 4) {
    block_annotation$role[last] <- "unclear"
  }
  
  return(block_annotation)
}
