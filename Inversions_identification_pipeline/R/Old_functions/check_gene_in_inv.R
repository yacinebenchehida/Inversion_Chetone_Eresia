check_gene_in_inv <- function(inv, gene_order, gene_name = NULL) {
  # If inv is NULL or empty, just return FALSE
  if (is.null(inv) || nrow(inv) == 0) {
    return(FALSE)
  }
  
  # Convert to data.table safely
  inv <- data.table::as.data.table(inv)
  
  # If required columns are missing, skip quietly
  if (!all(c("inversion_id", "gene_name", "role") %in% names(inv))) {
    return(FALSE)
  }
  
  # Create lookup: gene name → numeric order index
  gene_pos <- setNames(seq_along(gene_order), gene_order)
  
  # Keep only start/end roles and reshape wide
  inv_pairs <- inv[role %in% c("start", "end"), .(gene_name, role, inversion_id)]
  inv_pairs <- data.table::dcast(inv_pairs, inversion_id ~ role, value.var = "gene_name")
  
  # Skip if no valid start–end pairs exist
  if (nrow(inv_pairs) == 0) {
    return(FALSE)
  }
  
  # Add numeric positions for start and end
  inv_pairs[, start_idx := gene_pos[start]]
  inv_pairs[, end_idx   := gene_pos[end]]
  
  # Use midpoint gene if none provided
  if (is.null(gene_name)) {
    mid_idx <- ceiling(length(gene_order) / 2)
    gene_name <- gene_order[mid_idx]
    message(sprintf("No gene_name provided → using midpoint gene: %s", gene_name))
  }
  
  # If gene not in gene_order, return FALSE
  if (!gene_name %in% names(gene_pos)) {
    return(FALSE)
  }
  
  # Index of target gene
  target_idx <- gene_pos[gene_name]
  
  # Check if target lies between start and end (any orientation)
  inside <- inv_pairs[
    (target_idx >= pmin(start_idx, end_idx)) &
      (target_idx <= pmax(start_idx, end_idx))
  ]
  
  # Return TRUE if gene inside any inversion
  return(nrow(inside) > 0)
}
