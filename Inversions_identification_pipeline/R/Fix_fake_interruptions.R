Fix_fake_interruptions <- function(inv_table, dt, pos_col_index = 9, gene_col_index = 1, max_gap = 2, direction) {
    
    # Ensure inv_table is a data.table
    inv_table <- as.data.table(inv_table)  # Convert inversion table to data.table
    
    # Ensure dt is a data.table
    dt <- as.data.table(dt)  # Convert full dataset to data.table
    
    # If inversion table is NULL or empty, return it immediately
    if (is.null(inv_table) || nrow(inv_table) == 0) return(inv_table)  # Nothing to process
    
    # Initialize list to store extended inversions
    extended <- list()  # Empty list for storing extended inversions
    
    # Initialize inversion counter for assigning new inversion IDs
    inversion_counter <- 1  # Counter for new inversions
    
    # Extract unique inversion IDs in order
    inv_ids <- unique(inv_table$inversion_id)  # Get all inversion IDs in the table
    
    # Loop over each inversion ID
    for (i in seq_along(inv_ids)) {
      
      # Subset rows corresponding to current inversion
      rows <- inv_table[inversion_id == inv_ids[i]]  # Current inversion rows
      
      # Identify roles present in current inversion
      roles <- unique(rows$role)  # Can be "start", "end", or "undefined"
      
      # Get start gene if present
      start_gene <- if ("start" %in% roles) rows$gene_name[rows$role == "start"] else NA_character_  # Start gene
      
      # Get end gene if present
      end_gene <- if ("end" %in% roles) rows$gene_name[rows$role == "end"] else NA_character_  # End gene
      
      # Get inversion type
      type <- unique(rows$type)  # Type of inversion
      
      # Determine start index in full dataset
      start_idx <- if (!is.na(start_gene)) which(dt[[gene_col_index]] == start_gene) else 1  # Default to 1 if NA
      
      # Determine end index in full dataset
      end_idx <- if (!is.na(end_gene)) which(dt[[gene_col_index]] == end_gene) else nrow(dt)  # Default to last row if NA
      
      # Determine next inversion start index to avoid overlapping
      if (i < length(inv_ids)) {
        next_start_gene <- inv_table[inversion_id == inv_ids[i + 1] & role == "start"]$gene_name  # Next inversion start gene
        next_start_idx <- if (length(next_start_gene) > 0) which(dt[[gene_col_index]] == next_start_gene) else nrow(dt) + 1  # Index in dt
      } else {
        next_start_idx <- nrow(dt) + 1  # No next inversion, use beyond dataset
      }
      
      # Determine expected inversion direction based on global direction
      expected_sign <- if (direction == "positive") -1 else 1  # Positive → inversion goes down, negative → goes up
      
      # Initialize extended end index as original end
      extended_end_idx <- end_idx  # Will potentially extend this
      
      # Initialize current index for scanning
      current_idx <- end_idx + 1  # Start scanning after current end
      
      # Initialize gap counter to track interruptions
      gap_count <- 0  # Count consecutive genes opposite to inversion trend
      
      # Loop forward to extend inversion while respecting max_gap and next inversion
      while (current_idx < next_start_idx && current_idx <= nrow(dt)) {  # Stop at next inversion or dataset end
        
        # Compute step between consecutive positions
        step <- dt[[pos_col_index]][current_idx] - dt[[pos_col_index]][current_idx - 1]  # Difference in positions
        
        # Compute sign of step (+1, -1, or 0)
        step_sign <- sign(step)  # Step direction
        
        # Check if step follows inversion trend
        if (step_sign == expected_sign) {  # If in inversion direction
          
          # Extend breakpoint to current gene
          extended_end_idx <- current_idx  # Update extended end
          
          # Reset gap counter
          gap_count <- 0  # No interruption
          
        } else {  # Step opposite to inversion trend
          
          # Increment gap counter
          gap_count <- gap_count + 1  # Count interruption
          
          # Stop extending if gap exceeds allowed max_gap
          if (gap_count > max_gap) break  # Terminate extension
          
        }
        
        # Move to next gene
        current_idx <- current_idx + 1  # Advance scanning index
      }
      
      # Determine new end gene after extension
      new_end_gene <- dt[[gene_col_index]][extended_end_idx]  # Gene at extended end
      
      # Check if current inversion is undefined
      if ("undefined" %in% roles) {
        
        # Update gene_name to extended gene; role remains undefined
        rows$gene_name[rows$role == "undefined"] <- new_end_gene  # Extend undefined
        # Type remains unchanged
        rows$inversion_id <- inversion_counter  # Assign new inversion ID
        extended[[inversion_counter]] <- rows  # Append to result
        inversion_counter <- inversion_counter + 1  # Increment counter
        
      } else if ("start" %in% roles && "end" %in% roles) {  # Start/End inversion
        
        # Check if extended end hits last gene in dataset
        if (extended_end_idx == nrow(dt)) {  # Edge case: end reached dataset end
          
          # Replace original start with undefined
          new_rows <- data.table(
            inversion_id = inversion_counter,  # New inversion ID
            gene_name = start_gene,  # Original start becomes undefined breakpoint
            role = "undefined",  # Role is undefined
            type = type  # Type remains
          )
          
          # Append to result
          extended[[inversion_counter]] <- new_rows
          inversion_counter <- inversion_counter + 1  # Increment counter
          
        } else {  # Normal start/end inversion
          
          # Update end gene to extended end; start remains unchanged
          rows$gene_name[rows$role == "end"] <- new_end_gene  # Update end
          rows$inversion_id <- inversion_counter  # Assign new inversion ID
          extended[[inversion_counter]] <- rows  # Append
          inversion_counter <- inversion_counter + 1  # Increment counter
          
        }
      }
      
    }  # End loop over inversions
    
    # Combine all extended inversions into one data.table
    return(rbindlist(extended))  # Return fully extended inversion table
  }
