pct_genes_in_scaffold <- function(filtered_dt,gene_order){
  gene_order <- read.table(gene_order)
  return(nrow(filtered_dt)/nrow(gene_order))
}
