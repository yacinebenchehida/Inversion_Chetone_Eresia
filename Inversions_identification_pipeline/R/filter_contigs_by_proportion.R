#-------------------------------
# Function: get the right contig
#-------------------------------
filter_contigs_by_proportion <- function(blast_dt, contig_col_index = 2, pos_col_index = 9,
                                         bit_score_index = 12, threshold = 0.7, min_bit_score = 80, outlier = TRUE) {
  
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
  high_chaos <- dispersion_score(filtered_dt, pos_col_index = pos_col_index, chaos_excess_score = 1.7)
  
  # If chaos is too high, print message and stop
  if (high_chaos) {
    message("Detected excessive noise in blast hits. Dispersion score above threshold. Dropping this genome.")
    filtered_dt <- NULL
  }
  
  # Return filtered table
  return(filtered_dt)
}
