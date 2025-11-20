check_start_end_melpo_genes_present <- function(dt, genes_set1, genes_set2, gene_col_index) {

  # extract gene names from the specified column
  gene_col <- dt[[gene_col_index]]

  # check presence
  found1 <- any(gene_col %in% genes_set1)
  found2 <- any(gene_col %in% genes_set2)

  # condition: at least one from each set
  if (found1 && found2) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
