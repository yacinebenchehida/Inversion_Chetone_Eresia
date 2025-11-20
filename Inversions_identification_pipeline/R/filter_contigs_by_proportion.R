#-------------------------------
# Function: get the right contig
#-------------------------------
filter_contigs_by_proportion <- function(blast_dt, contig_col_index = 2, pos_col_index = 9,
                                         bit_score_index = 12, threshold = 0.5, min_bit_score = 70, outlier = TRUE, gene_order) {
  
  # Ensure input is a data.table
  dt <- as.data.table(blast_dt)
  
  # Filter rows with bit score below threshold
  dt <- dt[dt[[bit_score_index]] >= min_bit_score]
  
  # Remove outliers if requested
  if (outlier) {
    mean_pos <- mean(dt[[pos_col_index]], na.rm = TRUE)   # compute mean of positions
    sd_pos <- sd(dt[[pos_col_index]], na.rm = TRUE)       # compute standard deviation of positions
    lower_margin <- mean_pos -   2 * sd_pos
    upper_margin <- mean_pos +   2 * sd_pos
    dt <- dt[dt[[pos_col_index]] >= lower_margin & dt[[pos_col_index]] <= upper_margin]  # keep values within bounds
  }
  
  # Handle duplicate positions: keep the row with highest bit score
  dt <- dt[!duplicated(dt[[pos_col_index]]), ]
  
  # Count hits per contig
  contig_counts <- dt[, .N, by = dt[[contig_col_index]]]
  setnames(contig_counts, "dt", "contig")                     # rename for clarity
  
  # Compute proportion of total queries
  contig_counts[, proportion := N / sum(N)]
  
  # Identify contigs meeting threshold
  selected_contigs <- contig_counts[proportion >= threshold, contig]
  
  # Check if any contigs pass threshold
  if (length(selected_contigs) == 0) {
    max_prop <- max(contig_counts$proportion)
    max_contig <- contig_counts[proportion == max_prop, contig]
    message(sprintf("No contigs have frequency above the indicated threshold. Most abundant contig(s): %s with proportion %.1f%%", 
                    paste(max_contig, collapse = ", "), max_prop * 100))
    stop("Dropping this genome")
  }
  
  # Filter original table for those contigs
  filtered_dt <- dt[dt[[contig_col_index]] %in% selected_contigs]
  
  # Call the Chaos_noise_score function on the filtered data
  high_chaos <- dispersion_score(filtered_dt, pos_col_index = pos_col_index, chaos_excess_score = 1.4)
  
  # If chaos is too high, print message and stop
  if (high_chaos) {
    message("Detected excessive noise in blast hits. Dispersion score above threshold. Dropping this genome.")
    filtered_dt <- NULL
  }
  
  # Check if there is enough genomics data (scaffold >= 1Mb)
  if (is_big_enough(filtered_dt,threshold = 1000000)==FALSE){
    message("The scaffold is smaller than 1Mb. Dropping this genome.")
    filtered_dt <- NULL
  }


  # Check if 80% of the melpomene chromosomes 15 genes are there
  # Read gene order and transform it into a vector
    if (is.character(gene_order) && length(gene_order) == 1 && file.exists(gene_order)) {
    gene_order <- data.table::fread(gene_order, header = FALSE)[[1]]
  }

  pct_genes_present <- pct_genes_in_scaffold(filtered_dt=filtered_dt, gene_order="genes_order_busco.txt")
  if(pct_genes_present < 0.80){
    filtered_dt <- NULL
  }

  # Make sure at least one of the first 10 and 10 last melpo genes are present
  #if(check_start_end_melpo_genes_present(filtered_dt, 
   # genes_set1 = gene_order[1:10],
  #genes_set2 = gene_order[(length(gene_order) - 9):length(gene_order)],
  #  gene_col_index=1)==FALSE){
  #  filtered_dt <- NULL
  #}

  # Return filtered table
  return(filtered_dt)
}
