detect_inversions_translo_using_blocks <- function(dt,
                                                   pos_col_index = 9,
                                                   gene_col_index = 1,
                                                   min_consecutive = 3,
                                                   direction,
                                                   block_data) {
  # --------------# 
  #  Prepare data # 
  # --------------# 
  # Guard: need enough rows to compute first/last medians and detect runs
  if (nrow(dt) < 6) stop("Not enough points to detect translocations (need ≥ 6).")
  
  # Ensure input is a data.table
  dt <- as.data.table(dt)
  
  # Extract and store gene names
  dt[[gene_col_index]] <- as.character(dt[[gene_col_index]])
  
  # Extract positions
  positions <- dt[[pos_col_index]]
  
  # Pull convenient vectors (names unchanged)
  genes  <- as.character(dt[[gene_col_index]])   # gene identifiers (ordered)
  
  # Define in-bounds using medians of the first 3 and the last 3 positions
  med_first3 <- median(positions[1:3])                                         # median of first 3
  med_last3  <- median(positions[(length(positions) - 2):length(positions)])   # median of last 3
  lower_bound <- min(med_first3, med_last3)                                    # inclusive lower bound
  upper_bound <- max(med_first3, med_last3)                                    # inclusive upper bound
  
  
  # Initialize role and breakpoint columns
  block_data <- block_data[, -c("mean_step"), with = FALSE]
  block_data[, role := NA_character_]           # create empty role column
  block_data[, start_breakpoint := NA_character_]  # create empty start breakpoint column
  block_data[, end_breakpoint := NA_character_]    # create empty end breakpoint column
  
  # -------------------# 
  #  Mark small blocks # 
  # -------------------# 
  # Assign role "small block" to any block with <= min_consecutive genes
  block_data[n_genes <= min_consecutive, role := "small_block"]
  block_data[n_genes <= min_consecutive, start_breakpoint := "Not_applicable"]
  block_data[n_genes <= min_consecutive, end_breakpoint := "Not_applicable"]
  
  # -------------------# 
  #  Detect inversions # 
  # -------------------# 
  # Loop over rows and assign role and breakpoints
  for (i in seq_len(nrow(block_data))) {
    
    # Skip already-marked small blocks
    if (!is.na(block_data$role[i]) && block_data$role[i] == "small_block") next
    
    if (block_data$direction[i] != direction) {
      # Block is an inversion
      block_data$role[i] <- "inversion"
      
      # Determine start breakpoint
      if (i == 1 || block_data$role[i-1] == "small_block") {
        block_data$start_breakpoint[i] <- "Undefined"
      } else {
        block_data$start_breakpoint[i] <- paste(block_data$end_gene[i-1], block_data$start_gene[i], sep="_")
      }
      
      # Determine end breakpoint
      if (i < nrow(block_data) && !is.na(block_data$role[i+1]) && block_data$role[i+1] == "small_block") {
        block_data$end_breakpoint[i] <- "Undefined"
      } else if (i == nrow(block_data)) {
        block_data$end_breakpoint[i] <- "Undefined"
      } else {
        block_data$end_breakpoint[i] <- paste(block_data$end_gene[i], block_data$start_gene[i+1], sep="_")
      }
      
    } else {
      # Block is colinear
      block_data$role[i] <- "colinear"
    }
  }
  
  # ---------------------------------# 
  #  Detect translocation+inversions # 
  # ---------------------------------#
  for (i in seq_len(nrow(block_data))) {
    if (block_data$role[i] == "inversion") {
      if(i > 1 && i < nrow(block_data)) {
        inversion_median <- block_data$median_pos[i]
        if (block_data$median_pos[i] < lower_bound || block_data$median_pos[i] > upper_bound) {
          block_data$role[i] <- "inversion+translocation"
        }
      }
    }
  }
  
  # ----------------------# 
  #  Detect translocation # 
  # ----------------------#
  for (i in seq_len(nrow(block_data))) {
    if (block_data$role[i] == "colinear") {
      block_data$start_breakpoint[i] <- "Not_applicable"
      block_data$end_breakpoint[i] <- "Not_applicable"
      if(i > 1 && i < nrow(block_data)) {
        inversion_median <- block_data$median_pos[i]
        if (block_data$median_pos[i] < lower_bound || block_data$median_pos[i] > upper_bound) {
          block_data$role[i] <- "translocation"
          # Determine start breakpoint
          if (i == 1 || block_data$role[i-1] == "small_block") {
            block_data$start_breakpoint[i] <- "Undefined"
          } else {
            block_data$start_breakpoint[i] <- paste(block_data$end_gene[i-1], block_data$start_gene[i], sep="_")
          }
          
          # Determine end breakpoint
          if (i < nrow(block_data) && !is.na(block_data$role[i+1]) && block_data$role[i+1] == "small_block") {
            block_data$end_breakpoint[i] <- "Undefined"
          } else if (i == nrow(block_data)) {
            block_data$end_breakpoint[i] <- "Undefined"
          } else {
            block_data$end_breakpoint[i] <- paste(block_data$end_gene[i], block_data$start_gene[i+1], sep="_")
          }
          
        }
      }
    }
  }
  
  # ---------------------------------------------------# 
  #  Correct for translos wrongly assign to inversions # 
  # ---------------------------------------------------#
  cleaner_block_data <- block_data[role != "small_block"]
  
  if (all(cleaner_block_data$role != "colinear")) {
    if (all(cleaner_block_data$role %in% c("inversion","inversion+translocation"))) {
      if (all(cleaner_block_data$direction != direction)) {
        set(block_data, i = which(block_data$role %in% c("inversion","inversion+translocation")),
            j = "role", value = "translocation")
      }
    }
  }
  
  # -----------------------------#
  #  Adjust defined breakpoints  #
  # -----------------------------#
  for(i in seq_len(nrow(block_data))) {
    # Only for inversion/translocation/inversion+translocation
    if(block_data$role[i] %in% c("inversion","translocation","inversion+translocation")) {
      
      # Fix start breakpoint if defined
      if(!is.na(block_data$start_breakpoint[i]) && block_data$start_breakpoint[i] != "Undefined") {
        if(i > 1) {
          prev_median <- block_data$median_pos[i-1]
          curr_median <- block_data$median_pos[i]
          
          # find the median of the first min_consecutive point of the block i and of the last min_consecutive of the block i - 1
          debut_block <- block_data$start_gene[i]
          fin_block <- block_data$end_gene[i]
          start_row <- which(dt[[gene_col_index]] == debut_block)
          end_row   <- which(dt[[gene_col_index]] == fin_block)
          block_values <- dt[seq(from = start_row, to = end_row), ][[pos_col_index]]
          curr_median <- median(block_values[1:(min_consecutive + 1)])
            
          
          debut_block <- block_data$start_gene[i-1]
          fin_block <- block_data$end_gene[i-1]
          start_row <- which(dt[[gene_col_index]] == debut_block)
          end_row   <- which(dt[[gene_col_index]] == fin_block)
          block_values <- dt[seq(from = start_row, to = end_row), ][[pos_col_index]]
          prev_median <- median(block_values[(length(block_values) - min_consecutive):length(block_values)])
          
          
          # Compare previous block to current, adjust start breakpoint
          if((block_data$direction[i-1] == "positive" && prev_median < curr_median) ||
             (block_data$direction[i-1] == "negative" && prev_median > curr_median)) {
            curr_gene_start <- block_data$start_gene[i]
            block_data$start_breakpoint[i] <- paste(genes[which(genes==curr_gene_start)-2], genes[which(genes==curr_gene_start)-1] , sep="_")
            block_data$start_gene[i] <- genes[which(genes==curr_gene_start)-1]
            
          }
        }
      }
      
      # Fix end breakpoint if defined
      if(!is.na(block_data$end_breakpoint[i]) && block_data$end_breakpoint[i] != "Undefined") {
        if (i < nrow(block_data) && !(block_data$role[i+1] %in% c("inversion","translocation","inversion+translocation"))) {

          # find the median of the last min_consecutive point of the block i and of the first min_consecutive of the block i + 1
          debut_block <- block_data$start_gene[i]
          fin_block <- block_data$end_gene[i]
          start_row <- which(dt[[gene_col_index]] == debut_block)
          end_row   <- which(dt[[gene_col_index]] == fin_block)
          block_values <- dt[seq(from = start_row, to = end_row), ][[pos_col_index]]
          curr_median <- median(block_values[(length(block_values) - min_consecutive):length(block_values)])

          debut_block <- block_data$start_gene[i+1]
          fin_block <- block_data$end_gene[i+1]
          start_row <- which(dt[[gene_col_index]] == debut_block)
          end_row   <- which(dt[[gene_col_index]] == fin_block)
          block_values <- dt[seq(from = start_row, to = end_row), ][[pos_col_index]]
          next_median <- median(block_values[1:(min_consecutive + 1)])
          

          # Compare next block to current, adjust end breakpoint
          if((block_data$direction[i] == "positive" && curr_median < next_median) ||
             (block_data$direction[i] == "negative" && curr_median > next_median)) {
            curr_gene_end <- block_data$end_gene[i]
            
            block_data$end_breakpoint[i] <- paste(genes[which(genes==curr_gene_end)-1], genes[which(genes==curr_gene_end)] , sep="_")
            block_data$end_gene[i] <- genes[which(genes==curr_gene_end)-1]
            block_data$start_gene[i+1] <- genes[which(genes==curr_gene_end)]
            #block_data$n_genes[i] <- block_data$n_genes[i] - 1
            #}
          }
         }
      }
    }
  }
  # ---------------------------------------------------------#
  #  Classified instable direction blocks to unclear_blocks  #
  # ---------------------------------------------------------#
  for(i in seq_len(nrow(block_data))) {
    # Only for inversion/translocation/inversion+translocation
    if(block_data$role[i] %in% c("inversion","translocation","inversion+translocation")) {
      debut_block <- block_data$start_gene[i]
      fin_block <- block_data$end_gene[i]
      # Find the corresponding row numbers
      start_row <- which(dt[[gene_col_index]] == debut_block)
      end_row   <- which(dt[[gene_col_index]] == fin_block)
      
      # Extract all rows between them (inclusive)
      block_values <- dt[seq(from = start_row, to = end_row), ][[pos_col_index]]
      consecutive_sign_diff <- sign(diff(block_values))
      tbl <- table(consecutive_sign_diff)
      if(length(tbl) < 2) next
  
      ratio_sign_change_in_block <- max(tbl) / sum(tbl)
      if(ratio_sign_change_in_block < 0.75){
        block_data$role[i] <- "unclear"
        block_data$start_breakpoint[i] <- "Not_applicable"
        block_data$end_breakpoint[i] <- "Not_applicable"
      }
    }
  }
  
  # -----------------------------------------------------#
  #  Fix breakpoints assigned to different breakpoints   #
  # -----------------------------------------------------#
  for(i in seq_len(nrow(block_data))) {
    # Only for inversion/translocation/inversion+translocation
    if(block_data$role[i] %in% c("inversion","translocation","inversion+translocation")) {
      if(i<nrow(block_data)){
        if(!is.na(block_data$end_breakpoint[i]) && !is.na(block_data$start_breakpoint[i+1])) {
          block1_end <- unlist(strsplit(block_data$end_breakpoint[i], "_"))[1]
          block2_start <- unlist(strsplit(block_data$start_breakpoint[i+1], "_"))[2]
            if(block1_end==block2_start){
              block_data$end_breakpoint[i]=block_data$start_breakpoint[i+1]
              block_data$end_gene[i] <- unlist(strsplit(block_data$start_breakpoint[i+1], "_"))[1]
            }
        }
      }
    }
  }

  # -----------------------#
  #  Return the final data #
  # -----------------------#
  return(block_data)
}
