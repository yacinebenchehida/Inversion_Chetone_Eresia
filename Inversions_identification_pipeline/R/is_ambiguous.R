is_ambiguous <- function(block_annotation){
  
  # Remove rows where the role is "small_block"
  block_annotation <- block_annotation[role != "small_block"]
  
  # Default output is FALSE (not ambiguous)
  is_ambiguous = FALSE
  
  # Continue only if exactly two blocks remain
  if(nrow(block_annotation) == 2){
    
    # Check that the pair contains one inversion-type role
    # and one colinear/translocated-type role
    if (
      any(block_annotation$role %in% c("inversion", "inversion+translocation")) &&
      any(block_annotation$role %in% c("colinear", "translocated"))
    ) {
      
      # Check if at least one block has an undefined breakpoint
      if (any(block_annotation$start_breakpoint == "Undefined" |
              block_annotation$end_breakpoint == "Undefined")) {
        
        # Mark as ambiguous
        is_ambiguous = TRUE
      }
    }
  }
  
  # Return the logical result
  return(is_ambiguous)
}
