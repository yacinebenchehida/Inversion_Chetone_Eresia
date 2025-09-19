filter_contigs_by_proportion <- function(blast_dt, contig_col_index = 2, threshold = 0.7) {
  
  # Ensure input is a data.table
  dt <- as.data.table(blast_dt)
  
  # Count hits per contig
  contig_counts <- dt[, .N, by = dt[[contig_col_index]]]
  setnames(contig_counts, "dt", "contig")  # rename for clarity
  
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
    return(NULL)
  }
  
  # Filter original table for those contigs
  filtered_dt <- dt[dt[[contig_col_index]] %in% selected_contigs]
  
  # Return filtered table
  return(filtered_dt)
}
