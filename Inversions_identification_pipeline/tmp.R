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
    lower_margin <- mean_pos -   sd_pos
    upper_margin <- mean_pos +   sd_pos
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
    stop("Please decrease threshold")
  }
  
  # Filter original table for those contigs
  filtered_dt <- dt[dt[[contig_col_index]] %in% selected_contigs]
  
  # Return filtered table
  return(filtered_dt)
}

#------------------------------------------------------------------------------------------
#----------------------------------------------------------------------
# Function to assess overall direction based on first and last 3 points
#----------------------------------------------------------------------
assess_direction <- function(dt, pos_col_index = 9, keep_global = TRUE) {
  
  # Convert input to data.table for consistent indexing
  dt <- as.data.table(dt)
  
  # Extract the positions column based on provided index
  positions <- dt[[pos_col_index]]
  
  # Count number of positions
  n <- length(positions)
  
  # Stop if there are fewer than 6 points
  if (n < 6) stop("Not enough points to assess direction (need at least 6 rows).")
  
  # Extract first 3 positions
  first3 <- positions[1:3]
  # Compute median of first 3 positions for global trend comparison
  med_first3 <- median(first3)
  
  # Extract last 3 positions
  last3 <- positions[(n-2):n]
  # Compute median of last 3 positions for global trend comparison
  med_last3 <- median(last3)
  
  # Determine global direction based on median comparison
  if (med_first3 < med_last3) global_direction <- "positive"
  else if (med_first3 > med_last3) global_direction <- "negative"
  else global_direction <- "non-monotonous"
  
  # Check local monotony in first 3 points
  first_increasing <- all(diff(first3) > 0)
  first_decreasing <- all(diff(first3) < 0)
  
  # Check local monotony in last 3 points
  last_increasing <- all(diff(last3) > 0)
  last_decreasing <- all(diff(last3) < 0)
  
  # Decide overall direction combining global trend, local monotony, and keep_global option
  if ((first_increasing && last_increasing) || (first_decreasing && last_decreasing)) {
    # Both start and end monotonic in same direction
    if ((first_increasing && global_direction == "positive") || (first_decreasing && global_direction == "negative")) {
      # Local trend matches global trend → keep global trend
      direction <- global_direction
    } else {
      # Local trend conflicts with global trend
      if (keep_global) {
        warning("Global and local trends are inconsistent; using global trend.")
        direction <- global_direction
      } else {
        direction <- "non-monotonous"
      }
    }
    
  } else if ((first_increasing || first_decreasing) && !(last_increasing || last_decreasing)) {
    # Only start is monotonic, end is non-monotonic
    if ((first_increasing && global_direction == "positive") || (first_decreasing && global_direction == "negative")) {
      # Global trend matches start → keep global trend, with partial warning
      warning("Only the start of the window is monotonous; proceeding with caution.")
      direction <- global_direction
    } else {
      # Conflict between start and global trend
      if (keep_global) {
        warning("Global and local trends are inconsistent; using global trend.")
        direction <- global_direction
      } else {
        direction <- "non-monotonous"
      }
    }
    
  } else if ((first_increasing && last_decreasing) || (first_decreasing && last_increasing)) {
    # Start and end monotonic but in opposite directions → always non-monotonous
    direction <- "non-monotonous"
    
  } else {
    # Both start and end non-monotonic → always non-monotonous
    direction <- "non-monotonous"
  }
  
  # Handle rare case where global trend is exactly flat (medians equal)
  if (global_direction == "non-monotonous") {
    warning("Global trend is flat; cannot determine direction.")
    direction <- "non-monotonous"
  }
  
  # Return the assessed direction
  return(direction)
}



#------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------
# Function to detect inversion based on the direction of successive points
#-------------------------------------------------------------------------
detect_inversions <- function(dt_segment, pos_col_index = 9, gene_col_index = 1, min_consecutive = 3, direction = "positive") {
  dt_segment <- as.data.table(dt_segment)  # Ensure input is a data.table
  orig_genes <- as.character(dt_segment[[gene_col_index]])  # Store original gene names
  dt_segment[[gene_col_index]] <- orig_genes  # Keep original gene names
  
  positions <- dt_segment[[pos_col_index]]  # Extract positions
  expected_sign <- if (direction == "positive") 1 else -1  # Determine expected direction
  adj_diff <- diff(positions)  # Compute differences between consecutive positions
  reverse_steps <- sign(adj_diff) != expected_sign  # Identify steps opposite to expected direction
  
  results_inversion_table <- data.table(inversion_id = integer(), gene_name = character(), role = character())  # Initialize results table
  run_start <- NULL  # Initialize start index for current run
  run_length <- 0  # Initialize run length counter
  inversion_id <- 1  # Initialize inversion counter
  
  for (i in seq_along(reverse_steps)) {  # Loop over reverse steps
    if (reverse_steps[i]) {  # If step is reversed
      if (is.null(run_start)) run_start <- i  # Start new run if none active
      run_length <- run_length + 1  # Increment run length
    } else {  # If step returns to expected direction
      if (!is.null(run_start) && run_length >= min_consecutive) {  # Check if run qualifies as inversion
        # Adjust start index using i-1 and i+2 rule
        prev_pos <- ifelse(run_start > 1, positions[run_start - 1], positions[run_start])
        next_pos <- ifelse(run_start + 2 <= length(positions), positions[run_start + 2], positions[run_start + 1])
        if (direction == "positive") {
          inversion_start_idx <- ifelse(next_pos > prev_pos, run_start, run_start + 1)
        } else {
          inversion_start_idx <- ifelse(next_pos < prev_pos, run_start, run_start + 1)
        }
        inversion_end_idx <- run_start + run_length  # End index of inversion
        
        results_inversion_table <- rbind(results_inversion_table,  # Append inversion
                                         data.table(
                                           inversion_id = inversion_id,
                                           gene_name = c(dt_segment[[gene_col_index]][inversion_start_idx],
                                                         dt_segment[[gene_col_index]][inversion_end_idx]),
                                           role = c("start", "end")
                                         ),
                                         use.names = TRUE)
        inversion_id <- inversion_id + 1  # Increment inversion counter
      }
      run_start <- NULL  # Reset run
      run_length <- 0  # Reset length
    }
  }
  
  # Final check if last run reaches end
  if (!is.null(run_start) && run_length >= min_consecutive) {
    prev_pos <- ifelse(run_start > 1, positions[run_start - 1], positions[run_start])
    next_pos <- ifelse(run_start + 2 <= length(positions), positions[run_start + 2], positions[run_start + 1])
    inversion_start_idx <- ifelse(next_pos > prev_pos, run_start, run_start + 1)
    inversion_end_idx <- run_start + run_length
    
    results_inversion_table <- rbind(results_inversion_table,
                                     data.table(
                                       inversion_id = inversion_id,
                                       gene_name = c(dt_segment[[gene_col_index]][inversion_start_idx],
                                                     dt_segment[[gene_col_index]][inversion_end_idx]),
                                       role = c("start", "end")
                                     ),
                                     use.names = TRUE)
  }
  
  return(results_inversion_table)  # Return final inversion table
}
# End of detect_inversions_simple
#------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------

# Wrapper function for the above function that defines an interval-and make a plot of the region of interest with the interval if requested.
#--------------------------------------------------------------------------------------------------------------------------------------------------
wrapper_detect_inversions <- function(dt, pos_col_index = 9, gene_col_index = 1,
                                      min_consecutive = 3, boundary_genes = NULL,
                                      plot_window = FALSE) {
  dt <- data.table::as.data.table(dt)  # Convert input to data.table
  orig_genes <- as.character(dt[[gene_col_index]])  # Store original gene names
  dt[[gene_col_index]] <- orig_genes  # Keep original names
  
  direction <- assess_direction(dt, pos_col_index = pos_col_index, keep_global = FALSE)  # Determine direction on full dataset
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
  
  if (!is.null(results_inversion_table) && nrow(results_inversion_table) > 0) {
    results_inversion_table[, type := {
      start_gene <- gene_name[role == "start"]
      inv_start_pos <- count_dt[[pos_col_index]][match(start_gene, count_dt[[gene_col_index]])]
      first3_positions <- head(count_dt[[pos_col_index]], 3)
      last3_positions <- tail(count_dt[[pos_col_index]], 3)
      translocated <- is_translocated_inversion(direction, inv_start_pos, first3_positions, last3_positions)
      if (translocated) "inversion+translocation" else "simple inversion"
    }, by = inversion_id]
    
    # Merge oversplit/adjacent inversions (often happens when syntheny is modified within the inversion)
    results_inversion_table <- merge_oversplit_inversions(                                                         # Call merging function after type set
      inv_table = results_inversion_table,
      dt = count_dt,
      pos_col_index = pos_col_index,
      gene_col_index = gene_col_index,
      direction = direction
    )
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
  
  if (nrow(results_inversion_table) == 0) {
    message("There are no inversion in the specified window")
    return(NULL)
  } else {
    inv_numb <- nrow(results_inversion_table)/2
    message(paste("Youhouuuu! There are ",inv_numb," inversion(s) in the specified window"))
    return(results_inversion_table)  # Return inversion table
  }
  
} # End of detect_inversions wrapper
#------------------------------------------------------------------------------------------
is_translocated_inversion <- function(direction, inv_pos, first_positions, last_positions) {
  # Ensure inputs are numeric
  inv_pos <- as.numeric(inv_pos)
  first_positions <- as.numeric(first_positions)
  last_positions <- as.numeric(last_positions)
  
  if(length(first_positions) == 0 | length(last_positions) == 0) {
    stop("first_positions and last_positions must have at least one value each")
  }
  
  # Check for translocation depending on direction
  if(direction == "positive") {
    # Positive direction: translocated if inv_pos is smaller than the first positions
    # or larger than the last positions
    if(inv_pos < min(first_positions) | inv_pos > max(last_positions)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
    
  } else if(direction == "negative") {
    # Negative direction: translocated if inv_pos is larger than the first positions
    # or smaller than the last positions
    if(inv_pos > max(first_positions) | inv_pos < min(last_positions)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
    
  } else {
    stop("direction must be either 'positive' or 'negative'")
  }
}



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


merge_oversplit_inversions <- function(inv_table, dt, pos_col_index = 9, gene_col_index = 1, direction) {  # Define function with inputs: inversion table, full dataset, column indices, and inversion direction
  if (is.null(inv_table) || nrow(inv_table) == 0) return(inv_table)  # If inversion table is NULL or empty, return it immediately
  
  merged <- list()  # Initialize empty list to store merged inversions
  inversion_counter <- 1  # Initialize counter for new inversion IDs
  skip_next <- FALSE  # Initialize flag to skip next inversion if it was merged
  
  inv_ids <- unique(inv_table$inversion_id)  # Extract unique inversion IDs from the input table
  
  for (i in seq_along(inv_ids)) {  # Loop over each inversion ID
    if (skip_next) {  # Check if previous iteration merged the next inversion
      skip_next <- FALSE  # Reset skip flag
      next  # Skip current iteration
    }
    
    id1 <- inv_ids[i]  # Current inversion ID
    start_gene1 <- inv_table$gene_name[inv_table$inversion_id == id1 & inv_table$role == "start"]  # Get start gene of current inversion
    end_gene1 <- inv_table$gene_name[inv_table$inversion_id == id1 & inv_table$role == "end"]  # Get end gene of current inversion
    type1 <- inv_table$type[inv_table$inversion_id == id1 & inv_table$role == "start"]  # Get type of current inversion
    
    if (i < length(inv_ids)) {  # Check if a next inversion exists
      id2 <- inv_ids[i + 1]  # Next inversion ID
      start_gene2 <- inv_table$gene_name[inv_table$inversion_id == id2 & inv_table$role == "start"]  # Get start gene of next inversion
      end_gene2 <- inv_table$gene_name[inv_table$inversion_id == id2 & inv_table$role == "end"]  # Get end gene of next inversion
      type2 <- inv_table$type[inv_table$inversion_id == id2 & inv_table$role == "start"]  # Get type of next inversion
      
      if (type1 == type2) {  # Only attempt merge if both inversions have the same type
        idx1_end <- which(dt[[gene_col_index]] == end_gene1)  # Find index of end gene of first inversion in dataset
        idx2_start <- which(dt[[gene_col_index]] == start_gene2)  # Find index of start gene of second inversion in dataset
        gap_genes <- idx2_start - idx1_end - 1  # Compute number of genes between the inversions
        
        if (gap_genes <= 2) {  # Check if gap is small enough to consider merging
          pos1 <- dt[[pos_col_index]][idx1_end]  # Position of end gene of first inversion
          pos2 <- dt[[pos_col_index]][idx2_start]  # Position of start gene of second inversion
          
          merge_condition <- if (direction == "positive") pos2 < pos1 else pos2 > pos1  # Determine if positions indicate a merged inversion
          
          if (merge_condition) {  # If condition met, merge the two inversions
            merged[[inversion_counter]] <- data.table(  # Append merged inversion to list
              inversion_id = inversion_counter,  # Assign new inversion ID
              gene_name = c(start_gene1, end_gene2),  # Start gene is first inversion start, end gene is second inversion end
              role = c("start", "end"),  # Roles of genes
              type = type1  # Type remains the same
            )
            inversion_counter <- inversion_counter + 1  # Increment new inversion counter
            skip_next <- TRUE  # Set flag to skip the next inversion in loop
            next  # Skip rest of loop for this iteration
          }
        }
      }
    }
    
    merged[[inversion_counter]] <- data.table(  # If no merge, keep current inversion
      inversion_id = inversion_counter,  # Assign new inversion ID
      gene_name = c(start_gene1, end_gene1),  # Start and end genes of current inversion
      role = c("start", "end"),  # Roles of genes
      type = type1  # Type remains the same
    )
    inversion_counter <- inversion_counter + 1  # Increment inversion counter
  }
  
  return(data.table::rbindlist(merged))  # Combine list of merged inversions into single data.table and return
}



# args <- commandArgs(trailingOnly=TRUE)
folder_path <- "./mmseq2"                                 # Folder with GCA files
#gca_files <- list.files(folder_path, pattern = "^GCA.*", full.names = TRUE)  # List GCA files
gca_files <- list.files(folder_path, pattern = "*.tsv", full.names = TRUE)  # List GCA files

counter_ivory <- 0
gca_with_gene <- list()

# === Process each file and collect results ===
inversion_results <- lapply(seq_along(gca_files), function(i) {
  file <- gca_files[i]                                   # Current file path
  gca_name <- tools::file_path_sans_ext(basename(file))  # GCA name without extension
  message(sprintf("[%d/%d] Processing %s...", i, length(gca_files), gca_name))  # Progress message
  
  dt <- fread(file, header = FALSE)                      # Read BLAST output (no header)
  
  result <- tryCatch({                                   # Wrap per-file processing in tryCatch
    suppressMessages(suppressWarnings({
      filtered_dt <- filter_contigs_by_proportion(dt, contig_col_index = 2, threshold = 0.7,outlier = TRUE)  # Filter contigs
      if (!is.null(filtered_dt) && nrow(filtered_dt) > 0) {
        gene_order <- dt$V1    # Save gene order (as character) from this file
        
        # pdf(paste(gca_name,".pdf",sep=""))
        inv <- wrapper_detect_inversions(                # Detect inversions (existing function)
          filtered_dt,
          min_consecutive = 2,
          plot_window = FALSE)
        # Initialize list before loop
        
        
        # Inside the loop, replace your print line with:
        if ("HMEL000020g1.t1" %in% inv$gene_name) {
          print(gca_name)
          print(inv[inv$gene_name == "HMEL000020g1.t1", ])
          gca_with_gene[[length(gca_with_gene) + 1]] <<- gca_name
        }
          # boundary_genes = c("HMEL000004-RA", "HMEL000055g1.t1"))
        # dev.off()
        is_ivory_in_inv <- check_gene_in_inv(inv, gene_order, gene_name = "HMEL000025-RB")
        if(is_ivory_in_inv==TRUE){
          counter_ivory <<- counter_ivory + 1
        }
        
        if (!is.null(inv) && nrow(inv) > 0) {
          return(list(name = gca_name, data = inv, gene_order = gene_order))  # Return when inversions found
        }
      }
      NULL                                                # Return NULL if nothing found
    }))
  }, error = function(e) {
    message(sprintf("Skipping %s: %s", basename(file), e$message))  # On error, report and skip
    NULL
  })
  return(result)                                         # Return per-file result (or NULL)
})


# === Number of times ivory is found in an inversion ===
message(sprintf("Ivory is present in %d inversion(s)", counter_ivory))

# === Remove NULL results ===
inversion_results <- Filter(Negate(is.null), inversion_results)  # Remove files with no inversions

# === If no results, stop gracefully ===
if (length(inversion_results) == 0) stop("No inversion results found in any GCA files.")

# === Use gene order from the first successful file as reference ===
gene_order <- dt$V1       # Character vector of gene names in genome order

# === Combine all inversion tables into one data.table ===
inversion_results_df <- rbindlist(
  lapply(inversion_results, function(x) {
    df <- copy(x$data)                                 # Copy results table
    df[, GCA := x$name]                                # Add GCA column
    setcolorder(df, c("GCA", setdiff(names(df), "GCA"))) # Put GCA first
    return(df)                                         # Return table
  }),
  fill = TRUE                                          # Fill missing columns if any
)

# === Ensure gene_name is character (prevent factor-induced NAs) ===
inversion_results_df[, gene_name := as.character(gene_name)]  # Force character

# === Create unique breakpoint identifier BEFORE factor conversion ===
inversion_results_df[, unique_breakpoints := paste(gene_name, role, sep = "_")]  # Create unique id

# === Determine ordered levels for unique_breakpoints using gene_order and role order ===
# Build baseline breakpoint order from gene_order: for each gene, we want "<gene>_start" then "<gene>_end" (or existing role)
bp_levels <- c()                                       # Initialize vector
for (g in gene_order) {                                # Loop genome-ordered genes
  # Add possible roles in a consistent order to preserve visual grouping
  bp_levels <- c(bp_levels, paste0(g, "_start"), paste0(g, "_end"))
}
# Keep only levels that actually appear in your data (avoid NA/unused levels)
bp_levels <- intersect(bp_levels, unique(inversion_results_df$unique_breakpoints))  # Preserve order by genome

# === Convert unique_breakpoints to factor with genome-based ordering ===
inversion_results_df[, unique_breakpoints := factor(unique_breakpoints, levels = bp_levels)]

# === Count frequencies by unique_breakpoints and type (data.table grouped count) ===
df_type <- inversion_results_df[, .N, by = .(unique_breakpoints, type)]  # Count occurrences
setnames(df_type, "N", "Frequency")                     # Rename count column

# === Ensure that for any breakpoint missing a type we include a zero row (no NA rows created) ===
# Build complete set of (breakpoint, type) combinations present in data
all_breaks <- unique(df_type$unique_breakpoints)        # Breakpoints present
all_types  <- unique(df_type$type)                      # Types present (e.g., simple inversion, inversion+translocation)
grid <- CJ(unique_breakpoints = all_breaks, type = all_types, unique = TRUE)  # Cartesian join
setkey(grid, unique_breakpoints, type)
setkey(df_type, unique_breakpoints, type)
df_type <- merge(grid, df_type, by = c("unique_breakpoints", "type"), all.x = TRUE)  # Left join to ensure zeros
df_type[is.na(Frequency), Frequency := 0]               # Replace NA counts with 0

# === Force factor ordering on df_type to preserve plotting order ===
df_type[, unique_breakpoints := factor(as.character(unique_breakpoints), levels = levels(inversion_results_df$unique_breakpoints))]

# === Plot: stacked bar (absolute counts) with genome-ordered breakpoints ===
ggplot(df_type, aes(x = Frequency, y = unique_breakpoints, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +      # Stacked counts
  scale_fill_manual(values = c("simple inversion" = "steelblue", "inversion+translocation" = "darkred")) +
  theme_minimal() +
  labs(title = "Frequencies of inversion types per breakpoint",
       x = "Count",
       y = "Breakpoint",
       fill = "Type") +
  theme(axis.text.y = element_text(size = 7))            # Adjust y-axis label size

library(ggplot2)
library(data.table)

setwd("/Users/yacinebenchehida/Desktop/Convergent_evolution/Chetone_histrio/R_package")
# Read BLAST output
dt <- fread("Inputs/GCA_917862395.2_iHelSar1.2.txt",header = FALSE)
filtered_dt <- filter_contigs_by_proportion(dt,2, threshold =0.6,outlier = TRUE)
my_dir <- assess_direction(filtered_dt)
print(wrapper_detect_inversions(filtered_dt, min_consecutive=3, plot_window=TRUE))
invers_table <- wrapper_detect_inversions(filtered_dt, min_consecutive=3, plot_window=TRUE)

merge_oversplit_inversions(invers_table, dt, pos_col_index = 9, gene_col_index = 1, direction)

detect_translocations(filtered_dt, pos_col_index = 9, gene_col_index = 1, min_consecutive = 3, direction=my_dir)
