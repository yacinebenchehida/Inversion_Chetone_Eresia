check_gene_in_block <- function(block_data,
                                gene_order,
                                gene_name) {
  
  # block_data missing or empty -> return FALSE
  if (is.null(block_data) || nrow(block_data) == 0) {
    return(FALSE)
  }
  
  # require gene_order
  if (missing(gene_order) || is.null(gene_order) || length(gene_order) == 0) {
    stop("gene_order must be provided and must contain at least one gene.")
  }
  
  # If gene_order is a character vector of length 1 that matches a file --> read it
  if (is.character(gene_order) && length(gene_order) == 1 && file.exists(gene_order)) {
    # Read file, one gene per line
    gene_order <- data.table::fread(gene_order, header = FALSE)
    # transform gene_order to a vector
    gene_order <- gene_order[[1]]
  }
  
  # require gene_name
  if (missing(gene_name) || is.null(gene_name)) {
    stop("A gene_name must be provided.")
  }
  
  # convert to data.table
  block_data <- data.table::as.data.table(block_data)
  
  # required columns in block_data
  required_cols <- c("start_gene","end_gene","role")
  if (!all(required_cols %in% names(block_data))) {
    stop("block_data is missing required columns: start_gene, end_gene, role.")
  }
  
  # lookup from gene name to index in canonical order
  gene_pos <- setNames(seq_along(gene_order), gene_order)
  
  # check gene_name exists in gene_order
  if (!gene_name %in% names(gene_pos)) {
    stop(sprintf("The gene '%s' is not found in gene_order.", gene_name))
  }
  
  # target gene index
  target_idx <- gene_pos[[gene_name]]
  
  # roles that count as inversions
  allowed_roles <- c("inversion","translocation","inversion+translocation")
  inv_rows <- block_data[role %in% allowed_roles]
  
  # no inversion-like blocks -> FALSE
  if (nrow(inv_rows) == 0) {
    return(FALSE)
  }
  
  # map terminal genes to canonical indices
  inv_rows[, start_idx := gene_pos[start_gene]]
  inv_rows[, end_idx   := gene_pos[end_gene]]
  
  # drop intervals missing endpoints
  inv_rows <- inv_rows[!is.na(start_idx) & !is.na(end_idx)]
  if (nrow(inv_rows) == 0) {
    return(FALSE)
  }
  
  # inclusive interval bounds
  inv_rows[, lo := pmin(start_idx, end_idx)]
  inv_rows[, hi := pmax(start_idx, end_idx)]
  
  # check whether target lies inside any block
  inside_any <- inv_rows[(target_idx >= lo) & (target_idx <= hi)]
  
  return(nrow(inside_any) > 0)
}
