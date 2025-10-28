wrapper_detect_inversions <- function(dt, pos_col_index = 9, gene_col_index = 1,
                                      min_consecutive = 3, boundary_genes = NULL,
                                      plot_window = FALSE) {
  dt <- data.table::as.data.table(dt)  # Convert input to data.table
  orig_genes <- as.character(dt[[gene_col_index]])  # Store original gene names
  dt[[gene_col_index]] <- orig_genes  # Keep original names
  
  direction <- assess_direction(dt, pos_col_index = pos_col_index, keep_global = TRUE)  # Determine direction on full dataset
  if (direction == "non-monotonous") {  # Check if direction is valid
    message("Full data non-monotonous. No inversion detection performed.")  
    return(NULL)  # Exit if non-monotonous
  }
  
  count_dt <- dt  # Copy full dataset for counting inversions
  if (!is.null(boundary_genes)) {  # If boundaries provided
    if (length(boundary_genes) != 2) stop("boundary_genes must contain exactly two gene names.")  # Enforce 2 genes
    positions_idxs <- which(orig_genes %in% boundary_genes)  # Find indices of boundary genes
    if (length(positions_idxs) != 2) stop("One or both boundary genes not found in the data.")  # Validate presence
    start_idx <- min(positions_idxs)  # Determine left boundary
    end_idx <- max(positions_idxs)  # Determine right boundary
    count_dt <- dt[start_idx:end_idx, ]  # Subset data for counting inversions
  }
  
  results_inversion_table <- detect_inversions(  # Call simple inversion function
    dt_segment = count_dt,
    pos_col_index = pos_col_index,
    gene_col_index = gene_col_index,
    min_consecutive = min_consecutive,
    direction = direction
  )
  
  if (!is.null(results_inversion_table) && nrow(results_inversion_table) > 0) {  # Only proceed if table exists
    # Handle "undefined" inversions separately and assign type to start/end inversions
    results_inversion_table[, type := {
      if (all(role == "undefined")) {  # All undefined → simple inversion type
        rep("simple inversion", .N)
      } else {  # Start/end inversions
        start_gene <- gene_name[role == "start"]
        if (length(start_gene) == 0) {  # No start gene → nothing to compute
          rep(NA_character_, .N)
        } else {
          inv_start_pos <- count_dt[[pos_col_index]][match(start_gene, count_dt[[gene_col_index]])]
          first3_positions <- head(count_dt[[pos_col_index]], 3)
          last3_positions <- tail(count_dt[[pos_col_index]], 3)
          translocated <- is_translocated_inversion(direction, inv_start_pos, first3_positions, last3_positions)
          ifelse(translocated, "inversion+translocation", "simple inversion")
        }
      }
    }, by = inversion_id]
    
    # Merge oversplit/adjacent inversions for all types
    merged <- merge_oversplit_inversions(
      inv_table = results_inversion_table,
      dt = count_dt,
      pos_col_index = pos_col_index,
      gene_col_index = gene_col_index,
      direction = direction
    )
    if (!is.null(merged) && nrow(merged) > 0) results_inversion_table <- merged  # Keep merged table if not empty
    
    # Fix short interruptions in inversions (syntheny change or local blast errors)
    results_inversion_table <- Fix_fake_interruptions(inv_table = results_inversion_table, dt= count_dt, max_gap = 2,direction = direction)
  }
  
  if (plot_window) {  # Plot only if requested
    library(ggplot2)  # Load ggplot2
    plot_dt <- dt  # Use full dataset for plotting
    plot_dt[[gene_col_index]] <- factor(plot_dt[[gene_col_index]], levels = plot_dt[[gene_col_index]])  # Preserve gene order
    
    xmin_val <- -Inf  # Default left boundary for plot
    xmax_val <- Inf  # Default right boundary for plot
    if (!is.null(boundary_genes)) {  # Adjust if boundaries provided
      xmin_val <- which(plot_dt[[gene_col_index]] == boundary_genes[1])  # Left boundary
      xmax_val <- which(plot_dt[[gene_col_index]] == boundary_genes[2])  # Right boundary
    }
    
    p <- ggplot(data = plot_dt, aes_string(x = names(plot_dt)[gene_col_index], y = names(plot_dt)[pos_col_index])) +  # Map x and y
      theme_bw() +  # Clean theme
      geom_point() +  # Scatter plot
      ggtitle("Species") +  # Plot title
      {if(!is.null(boundary_genes)) annotate("rect", xmin = xmin_val, xmax = xmax_val, ymin = -Inf, ymax = Inf, alpha = 0.3, fill = "red")} +  # Highlight region
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) +  # Rotate x labels
      xlab("gene order") +  # X axis label
      ylab("Gene start position")  # Y axis label
    
    print(p)  # Display plot
  }
  
  # Count total inversions including undefined
  inv_numb <- length(unique(results_inversion_table$inversion_id))
  message(paste("Youhouuuu! There are", inv_numb, "inversion(s) in the specified window"))
  return(results_inversion_table)  # Return inversion table
}
