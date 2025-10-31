Fix_fake_interruptions <- function(dt, inv_table,
                                   pos_col_index = 9, gene_col_index = 1,
                                   direction = "positive", max_gap = 2) {
  # Convert dt to data.table (ensure correct type)
  dt <- as.data.table(dt)  # convert input sequence to data.table
  
  # Convert inv_table to data.table (ensure correct type)
  inv_table <- as.data.table(inv_table)  # convert input inversions table
  
  # If there are no inversions to process, return the table immediately
  if (is.null(inv_table) || nrow(inv_table) == 0) return(inv_table)  # nothing to do
  
  # Pull convenience vectors for quick lookups
  positions <- dt[[pos_col_index]]  # numeric vector of genomic coordinates in compared species
  genes <- dt[[gene_col_index]]     # character vector of gene names in melpomene order
  
  # Define the sign we expect for inversion continuation:
  # If direction == "positive" then non-inverted steps are positive,
  # so inversion steps are negative (-1). Conversely if direction == "negative",
  # inversion steps are +1. This value is compared to sign(position[i] - position[i-1]).
  expected_inv_sign <- if (direction == "positive") -1L else 1L  # expected sign for inversion steps
  
  # Prepare output list to collect extended inversion rows
  out_list <- list()  # will collect data.table rows for each processed inversion
  
  # Get unique inversion IDs in original order (preserve input order)
  inv_ids <- unique(inv_table$inversion_id)  # vector of inversion ids in order of appearance
  
  # Loop through each inversion id and process one at a time
  new_id_counter <- 1L  # new sequential id for output (keeps results compact and reindexed)
  for (k in seq_along(inv_ids)) {
    
    # Current inversion id
    cur_id <- inv_ids[k]  # numeric / integer id from input
    
    # Subset all rows belonging to this inversion
    cur_rows <- inv_table[inversion_id == cur_id]  # rows for this inversion (1 or 2 rows)
    
    # Determine roles present (could be "start" & "end" or single "undefined")
    roles <- unique(cur_rows$role)  # roles in this inversion
    
    # Retrieve inversion type (should be identical for rows of same inversion)
    inv_type <- unique(cur_rows$type)[1]  # type string (use first if multiple)
    
    # If the inversion has both start and end, identify them; else identify undefined
    start_gene <- if ("start" %in% roles) cur_rows$gene_name[cur_rows$role == "start"][1] else NA_character_  # start gene if present
    end_gene   <- if ("end"   %in% roles) cur_rows$gene_name[cur_rows$role == "end"][1]   else NA_character_  # end gene if present
    undef_gene <- if ("undefined" %in% roles) cur_rows$gene_name[cur_rows$role == "undefined"][1] else NA_character_  # undefined gene if present
    
    # Compute index positions in dt for start/end/undefined genes
    # If gene not found, which(...) returns integer(0); we convert to NA_integer_
    idx_start <- if (!is.na(start_gene) && length(which(genes == start_gene)) > 0) which(genes == start_gene)[1] else NA_integer_  # index of start
    idx_end   <- if (!is.na(end_gene)   && length(which(genes == end_gene))   > 0) which(genes == end_gene)[1]   else NA_integer_  # index of end
    idx_undef <- if (!is.na(undef_gene) && length(which(genes == undef_gene)) > 0) which(genes == undef_gene)[1] else NA_integer_  # index of undefined
    
    # Compute index of next inversion's start gene to avoid overlapping extension
    if (k < length(inv_ids)) {
      # Find start gene of next inversion if present, otherwise set beyond dataset boundary
      next_rows <- inv_table[inversion_id == inv_ids[k + 1]]
      next_start_gene <- if ("start" %in% unique(next_rows$role)) next_rows$gene_name[next_rows$role == "start"][1] else NA_character_
      next_start_idx <- if (!is.na(next_start_gene) && length(which(genes == next_start_gene)) > 0) which(genes == next_start_gene)[1] else (nrow(dt) + 1L)
    } else {
      # No next inversion: set next_start_idx beyond dataset end
      next_start_idx <- nrow(dt) + 1L
    }
    
    # Initialize extended_end_index variable (this will hold the index of the extended breakpoint)
    # For undefined inversions we will move from the undefined breakpoint outward in the inversion direction,
    # for start/end inversions we will extend the end only.
    extended_end_idx <- NA_integer_  # placeholder
    
    # ---------- CASE A: undefined inversion (single row) ----------
    if (!is.na(idx_undef)) {
      # Start scanning from the undefined breakpoint outward in the inversion direction.
      # We scan *forward* in gene order (increasing index).
      # When moving from i to i+1 we compute sign(positions[i+1] - positions[i]).
      # A value equal to expected_inv_sign indicates continuation of the inversion.
      
      # Set initial breakpoint index to the undefined gene index
      extended_end_idx <- idx_undef  # start from the undefined bp
      
      # Set scanning start index (we scan forward from the breakpoint)
      scan_idx <- idx_undef + 1L  # begin immediately after the undefined bp
      
      # Initialize interruption counter
      gap_counter <- 0L  # counts consecutive non-inversion steps
      
      # Loop while we haven't reached the start of the next inversion and within dataset
      while (scan_idx < next_start_idx && scan_idx <= nrow(dt)) {
        # Compute sign of the step between scan_idx-1 and scan_idx
        # (we use positions[scan_idx] - positions[scan_idx - 1])
        step_sign <- sign(positions[scan_idx] - positions[scan_idx - 1L])  # +1, 0, or -1
        
        # If step matches expected inversion sign, update breakpoint and reset gap counter
        if (!is.na(step_sign) && step_sign == expected_inv_sign) {
          extended_end_idx <- scan_idx  # extend bp to this gene
          gap_counter <- 0L  # reset interruption counter
        } else {
          # Step does not follow inversion sign → increment gap counter
          gap_counter <- gap_counter + 1L  # count interruption
          # If interruptions exceed allowed max, stop scanning
          if (gap_counter > max_gap) break  # exit loop
        }
        # Advance scanning index by one gene (always forward scanning for undefined)
        scan_idx <- scan_idx + 1L  # increment scan index
      } # end while for undefined
      
      # After scanning, set new_bp gene based on extended_end_idx
      new_bp_gene <- genes[extended_end_idx]  # gene id at extended breakpoint index
      
      # Build output row: undefined remains undefined and gene updated to new_bp_gene
      out_list[[length(out_list) + 1L]] <- data.table(inversion_id = new_id_counter,
                                                      gene_name = new_bp_gene,
                                                      role = "undefined",
                                                      type = inv_type)
      # increment output id counter
      new_id_counter <- new_id_counter + 1L
      
      # proceed to next inversion
      next  # continue for-loop
      
    } # end CASE A (undefined)
    
    # ---------- CASE B: start & end both present ----------
    if (!is.na(idx_start) && !is.na(idx_end)) {
      # We will attempt to extend the end only (never change start), scanning forward from end index+1.
      extended_end_idx <- idx_end  # initialize extended end to original end index
      scan_idx <- idx_end + 1L  # start scanning immediately after original end
      gap_counter <- 0L  # reset gap counter
      
      # Loop while within dataset and not overlapping next inversion
      while (scan_idx < next_start_idx && scan_idx <= nrow(dt)) {
        # Compute step_sign between scan_idx-1 and scan_idx
        step_sign <- sign(positions[scan_idx] - positions[scan_idx - 1L])  # +1,0,-1
        
        # If this step continues inversion (matches expected_inv_sign), extend end
        if (!is.na(step_sign) && step_sign == expected_inv_sign) {
          extended_end_idx <- scan_idx  # extend end to this gene index
          gap_counter <- 0L  # reset gap counter
        } else {
          # This step is opposite to inversion trend; count it as interruption
          gap_counter <- gap_counter + 1L  # increase interruption count
          # If interruptions exceed allowed max, we must stop extending
          if (gap_counter > max_gap) break  # stop scanning extension
        }
        # Move scanning index forward by one gene
        scan_idx <- scan_idx + 1L  # increment scan index
      } # end while scanning extension for start/end
      
      # At this point extended_end_idx is the furthest index we can extend to without overlapping next inversion
      
      # If we did not extend (extended_end_idx equals original idx_end), then keep original rows unchanged
      if (extended_end_idx == idx_end) {
        # Keep original start and end rows; preserve type
        out_list[[length(out_list) + 1L]] <- data.table(inversion_id = new_id_counter,
                                                        gene_name = c(start_gene, end_gene),
                                                        role = c("start", "end"),
                                                        type = c(inv_type, inv_type))
        new_id_counter <- new_id_counter + 1L  # increment output id
        # Move to next inversion
        next
      }
      
      # If we extended and the extended end hit the dataset last row, apply Case 3 behavior:
      if (extended_end_idx == nrow(dt)) {
        # RULE: when end extends to the dataset end, the original start becomes undefined
        # and the start/end pair disappears; we create a single undefined row with gene_name = start_gene.
        out_list[[length(out_list) + 1L]] <- data.table(inversion_id = new_id_counter,
                                                        gene_name = start_gene,
                                                        role = "undefined",
                                                        type = inv_type)
        new_id_counter <- new_id_counter + 1L  # increment output id
        # Continue to next inversion
        next
      }
      
      # Otherwise we extended but did NOT hit dataset end: update end only (start unchanged).
      new_end_gene <- genes[extended_end_idx]  # gene name at extended end index
      out_list[[length(out_list) + 1L]] <- data.table(inversion_id = new_id_counter,
                                                      gene_name = c(start_gene, new_end_gene),
                                                      role = c("start", "end"),
                                                      type = c(inv_type, inv_type))
      new_id_counter <- new_id_counter + 1L  # increment output id
      # Continue to next inversion
      next
    } # end CASE B
    
    # If we reach here something unexpected occurred (neither undefined nor start/end)
    # To be safe, append the original rows unchanged
    out_list[[length(out_list) + 1L]] <- copy(cur_rows)[, inversion_id := new_id_counter]  # append original rows reindexed
    new_id_counter <- new_id_counter + 1L  # increment id
  } # end for loop over inversions
  
  # Combine all output pieces into a single data.table
  result <- rbindlist(out_list, use.names = TRUE, fill = TRUE)  # combine results
  
  # Ensure result columns are ordered nicely: inversion_id, gene_name, role, type
  setcolorder(result, c("inversion_id", "gene_name", "role", "type"))  # reorder columns
  
  # Return the final extended inversion table
  return(result)  # return result to caller
}
