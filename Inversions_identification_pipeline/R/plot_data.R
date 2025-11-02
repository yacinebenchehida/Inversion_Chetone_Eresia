plot_data <- function(dt, block_annotation = NULL, pos_col_index = 9, gene_col_index = 1,
                       plot_title = "Gene order with blocks", annotate_blocks = TRUE) {
  # Convert dt to data.table
  dt <- data.table::as.data.table(dt)
  
  # Convert block_annotation to data.table if provided
  if (!is.null(block_annotation)) {
    block_annotation <- data.table::as.data.table(block_annotation)
  }
  
  # Factor genes in their order
  gene_names <- as.character(dt[[gene_col_index]])
  dt[[gene_col_index]] <- factor(gene_names, levels = gene_names)
  gene_to_index <- setNames(seq_along(gene_names), gene_names)
  
  library(ggplot2)
  
  # Base scatter plot
  p <- ggplot(dt, aes_string(x = names(dt)[gene_col_index], y = names(dt)[pos_col_index])) +
    geom_point(size = 1) +
    theme_bw() +
    ggtitle(plot_title) +
    xlab("Gene order") + ylab("Gene start position") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))
  
  # Overlay structural variant blocks if requested
  if (annotate_blocks && !is.null(block_annotation)) {
    highlight_roles <- c("translocation", "inversion", "inversion+translocation")
    highlight_blocks <- block_annotation[role %in% highlight_roles]
    
    if (nrow(highlight_blocks) > 0) {
      # Compute rectangle edges
      rects_df <- data.frame(
        xmin = sapply(highlight_blocks$start_gene, function(g) gene_to_index[[g]]) - 0.5,
        xmax = sapply(highlight_blocks$end_gene, function(g) gene_to_index[[g]]) + 0.5,
        role = highlight_blocks$role,
        start_breakpoint = highlight_blocks$start_breakpoint,
        end_breakpoint = highlight_blocks$end_breakpoint
      )
      
      # Plot rectangles
      p <- p + geom_rect(data = rects_df, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = role),
                         alpha = 0.3, inherit.aes = FALSE)
      
      # Compute dotted line positions: only edges without Undefined
      dotted_x <- c()
      for (i in seq_len(nrow(rects_df))) {
        if (!grepl("Undefined", rects_df$start_breakpoint[i])) dotted_x <- c(dotted_x, rects_df$xmin[i])
        if (!grepl("Undefined", rects_df$end_breakpoint[i])) dotted_x <- c(dotted_x, rects_df$xmax[i])
      }
      dotted_x <- unique(dotted_x)  # remove duplicates
      
      if (length(dotted_x) > 0) {
        p <- p + geom_vline(xintercept = dotted_x, linetype = "dotted", color = "black")
      }
      
      # Legend
      p <- p + scale_fill_manual(name = "Structural Variant",
                                 values = c("translocation" = "skyblue",
                                            "inversion" = "lightcoral",
                                            "inversion+translocation" = "palegreen"))
    }
  }
  
  return(p)
}
