detect_and_split_nested_translocations <- function(block_data,
                                                   dt,
                                                   direction,
                                                   pos_col_index = 9,
                                                   gene_col_index = 1,
                                                   min_consecutive = 3,
                                                   max_gap = 1) {
  # convert block_data to data.table
  block_data <- data.table::as.data.table(block_data)
  # convert dt to data.table
  dt <- data.table::as.data.table(dt)
  # extract numeric positions from dt using pos_col_index
  positions <- dt[[pos_col_index]]
  # extract gene names from dt using gene_col_index
  genes_dt <- as.character(dt[[gene_col_index]])
  
  # helper function to compute the sign of a step
  step_sign <- function(dx) {
    # return 0 if dx is NA
    if (is.na(dx)) return(0)
    # return +1 if dx greater than zero
    if (dx > 0) return(1)
    # return -1 if dx less than zero
    if (dx < 0) return(-1)
    # return 0 if dx equal zero
    return(0)
  }
  
  # list to store all final blocks (original or split)
  final_blocks <- list()
  # index to fill final_blocks list
  fb_index <- 1
  
  # loop over each block in block_data
  for (i in seq_len(nrow(block_data))) {
    
    # extract role of current block
    role_i <- block_data$role[i]
    # extract direction of current block
    dir_i <- block_data$direction[i]
    # extract start index of current block in dt
    start_i <- block_data$start_index[i]
    # extract end index of current block in dt
    end_i <- block_data$end_index[i]
    # compute length of current block in number of genes
    block_len <- end_i - start_i + 1
    
    # if current block is not colinear or inversion, keep it unchanged
    if (!(role_i %in% c("colinear", "inversion"))) {
      # append original block row to final_blocks
      final_blocks[[fb_index]] <- block_data[i]
      # increment final_blocks index
      fb_index <- fb_index + 1
      # continue to next block
      next
    }
    
    # if block is too short to hold start, end, and internal nested region
    if (block_len < (3 + 3 + min_consecutive)) {
      # append original block row unchanged
      final_blocks[[fb_index]] <- block_data[i]
      # increment final_blocks index
      fb_index <- fb_index + 1
      # continue to next block
      next
    }
    
    # extract positions inside this block
    block_vals <- positions[start_i:end_i]
    # extract gene names inside this block
    block_genes <- genes_dt[start_i:end_i]
    
    # compute reference median for first three positions of block
    block_start_ref <- median(block_vals[1:3])
    # compute reference median for last three positions of block
    block_end_ref <- median(block_vals[(block_len - 2):block_len])
    
    # compute lower bound of reference interval
    ref_lower <- min(block_start_ref, block_end_ref)
    # compute upper bound of reference interval
    ref_upper <- max(block_start_ref, block_end_ref)
    
    # create vector to classify each local index in block as outside high, outside low, or inside
    outside_state <- integer(block_len)
    # loop over each index in block_vals
    for (k in seq_len(block_len)) {
      # get value at position k
      val <- block_vals[k]
      # if value greater than ref_upper, mark as +1 (outside high)
      if (val > ref_upper) {
        outside_state[k] <- 1L
      } else if (val < ref_lower) {
        # if value less than ref_lower, mark as -1 (outside low)
        outside_state[k] <- -1L
      } else {
        # otherwise mark as 0 (inside)
        outside_state[k] <- 0L
      }
    }
    
    # define start of interior region local index (exclude first three)
    interior_start <- 4L
    # define end of interior region local index (exclude last three)
    interior_end <- block_len - 3L
    
    # if interior region is empty, no nested detection possible
    if (interior_start > interior_end) {
      # append original block unchanged
      final_blocks[[fb_index]] <- block_data[i]
      # increment final_blocks index
      fb_index <- fb_index + 1
      # continue to next block
      next
    }
    
    # list to store nested segments found in this block
    nested_segments <- list()
    # index to fill nested_segments list
    ns_index <- 1
    
    # initialize index to walk through interior region
    idx <- interior_start
    
    # loop while idx stays inside interior region
    while (idx <= interior_end) {
      # get state at current index (1, 0, or -1)
      state_here <- outside_state[idx]
      
      # if current state is inside (0), move to next index
      if (state_here == 0L) {
        # increment idx
        idx <- idx + 1L
        # continue to next iteration
        next
      }
      
      # if current state is outside, start a new run
      run_side <- state_here
      # local index where run begins
      run_start <- idx
      # current last index inside run
      run_end <- idx
      # count of outside points in this run
      count_outside <- 1L
      # count of inside points in this run
      count_inside <- 0L
      # local index of first outside point in run
      first_out_local <- idx
      # local index of last outside point in run
      last_out_local <- idx
      
      # move j to next index after idx
      j <- idx + 1L
      
      # extend the run while j stays inside interior region
      while (j <= interior_end) {
        # read state at index j
        s <- outside_state[j]
        
        # if state matches run_side, it is an outside point of same side
        if (s == run_side) {
          # update run_end
          run_end <- j
          # increment outside count
          count_outside <- count_outside + 1L
          # update last outside local index
          last_out_local <- j
          # reset inside count since this is outside
          # but we do not actually need to reset because we only care about total inside
          # move j forward
          j <- j + 1L
          # continue loop
          next
        }
        
        # if state is 0, it is an inside point
        if (s == 0L) {
          # increment inside count
          count_inside <- count_inside + 1L
          # if inside count exceeds max_gap, stop extending run
          if (count_inside > max_gap) {
            # break out of while loop
            break
          } else {
            # still within allowable gap, extend run to j
            run_end <- j
            # move j forward
            j <- j + 1L
            # continue loop
            next
          }
        }
        
        # if we reach here, state is opposite side (-run_side)
        # break because run cannot include opposite side
        break
      }
      
      # after extending run, compute how many indices are in run
      run_indices <- run_start:run_end
      # extract values for this run
      run_vals <- block_vals[run_indices]
      # build outside mask for this run using run_side
      outside_mask <- if (run_side == 1L) run_vals > ref_upper else run_vals < ref_lower
      # build inside mask for this run
      inside_mask <- !outside_mask
      # compute count of outside points in run from mask
      n_outside <- sum(outside_mask)
      # compute count of inside points in run from mask
      n_inside <- sum(inside_mask)
      
      # check if run has enough outside points and limited inside points
      if (n_outside >= min_consecutive && n_inside <= max_gap) {
        # indices of outside points in run_vals (1 based within run)
        out_idx_run <- which(outside_mask)
        # index in run of first outside point in outside-only region
        first_out_run <- out_idx_run[1]
        # index in run of last outside point in outside-only region
        last_out_run <- out_idx_run[length(out_idx_run)]
        # compute local block index of first outside point
        seg_start_local <- run_start + first_out_run - 1L
        # compute local block index of last outside point
        seg_end_local <- run_start + last_out_run - 1L
        
        # compute median of first min_consecutive outside points in run
        first_med_indices <- out_idx_run[seq_len(min_consecutive)]
        # compute first median over outside values
        first_med <- median(run_vals[first_med_indices])
        # compute median of last min_consecutive outside points in run
        last_med_indices <- tail(out_idx_run, min_consecutive)
        # compute last median over outside values
        last_med <- median(run_vals[last_med_indices])
        
        # check high side condition if run_side is +1
        if (run_side == 1L &&
            first_med > ref_upper &&
            last_med > ref_upper) {
          # record high side nested segment
          nested_segments[[ns_index]] <- list(
            seg_start = seg_start_local,
            seg_end = seg_end_local,
            side = "high"
          )
          # increment nested segment index
          ns_index <- ns_index + 1L
        }
        
        # check low side condition if run_side is -1
        if (run_side == -1L &&
            first_med < ref_lower &&
            last_med < ref_lower) {
          # record low side nested segment
          nested_segments[[ns_index]] <- list(
            seg_start = seg_start_local,
            seg_end = seg_end_local,
            side = "low"
          )
          # increment nested segment index
          ns_index <- ns_index + 1L
        }
      }
      
      # move idx to first index after this run
      idx <- run_end + 1L
    }
    
    # if no nested segments were detected, keep block unchanged
    if (length(nested_segments) == 0) {
      # append original block row to final_blocks
      final_blocks[[fb_index]] <- block_data[i]
      # increment final_blocks index
      fb_index <- fb_index + 1
      # continue to next block
      next
    }
    
    # sort nested_segments by seg_start to ensure increasing order
    nested_segments <- nested_segments[order(vapply(nested_segments,
                                                    function(z) z$seg_start,
                                                    integer(1)))]
    
    # list to store subsegments (prefix, nested, suffix, etc.)
    subsegments <- list()
    # index to fill subsegments list
    sub_index <- 1
    # local current position in block (1 based)
    current_local <- 1L
    
    # loop over each nested segment in order
    for (seg in nested_segments) {
      # local start index for nested segment
      ns_start <- seg$seg_start
      # local end index for nested segment
      ns_end <- seg$seg_end
      
      # if there is a normal segment before this nested segment
      if (current_local < ns_start) {
        # record normal subsegment from current_local to ns_start - 1
        subsegments[[sub_index]] <- list(
          type = "normal",
          local_start = current_local,
          local_end = ns_start - 1L
        )
        # increment subsegment index
        sub_index <- sub_index + 1L
      }
      
      # record nested subsegment itself
      subsegments[[sub_index]] <- list(
        type = "nested",
        local_start = ns_start,
        local_end = ns_end
      )
      # increment subsegment index
      sub_index <- sub_index + 1L
      
      # move current_local to position after nested segment
      current_local <- ns_end + 1L
    }
    
    # if there is remaining tail after last nested segment
    if (current_local <= block_len) {
      # record trailing normal subsegment
      subsegments[[sub_index]] <- list(
        type = "normal",
        local_start = current_local,
        local_end = block_len
      )
      # increment subsegment index
      sub_index <- sub_index + 1L
    }
    
    # list to store new block rows for this original block
    new_blocks <- list()
    # index to fill new_blocks list
    nb_index <- 1
    
    # loop over all subsegments to build new block rows
    for (sub in subsegments) {
      # local start index inside block for this subsegment
      local_start <- sub$local_start
      # local end index inside block for this subsegment
      local_end <- sub$local_end
      # absolute start index in dt
      abs_start <- start_i + local_start - 1L
      # absolute end index in dt
      abs_end <- start_i + local_end - 1L
      
      # extract gene names for this subsegment
      gsub <- genes_dt[abs_start:abs_end]
      # extract positions for this subsegment
      psub <- positions[abs_start:abs_end]
      # compute number of genes in this subsegment
      seg_len <- length(gsub)
      
      # default sub_role equals original block role
      sub_role <- role_i
      
      # if this subsegment is nested, assign translocation type
      if (sub$type == "nested") {
        # if original role is colinear, use "translocation"
        if (role_i == "colinear") sub_role <- "translocation"
        # if original role is inversion, use "inversion+translocation"
        if (role_i == "inversion") sub_role <- "inversion+translocation"
      } else {
        # for normal subsegment, if too short, mark as small_block
        if (seg_len < min_consecutive) {
          sub_role <- "small_block"
        }
      }
      
      # initialize start_breakpoint for this subsegment as "Not_applicable"
      sub_start_bp <- "Not_applicable"
      # initialize end_breakpoint for this subsegment as "Not_applicable"
      sub_end_bp <- "Not_applicable"
      
      # if this subsegment is nested, compute true gene-level breakpoints
      if (sub$type == "nested") {
        # if there is a previous gene before this nested region in dt
        if (abs_start > 1L) {
          # define start breakpoint using previous gene and first gene of nested region
          sub_start_bp <- paste(genes_dt[abs_start - 1L],
                                genes_dt[abs_start],
                                sep = "_")
        } else {
          # if nested region starts at first gene, breakpoint undefined
          sub_start_bp <- "Undefined"
        }
        
        # if there is a next gene after this nested region in dt
        if (abs_end < length(genes_dt)) {
          # define end breakpoint using last gene of nested region and next gene after it
          sub_end_bp <- paste(genes_dt[abs_end],
                              genes_dt[abs_end + 1L],
                              sep = "_")
        } else {
          # if nested region ends at last gene, breakpoint undefined
          sub_end_bp <- "Undefined"
        }
      }
      
      # create new block row for this subsegment
      new_blocks[[nb_index]] <- data.table(
        block_id = NA_real_,
        start_index = abs_start,
        end_index = abs_end,
        start_gene = gsub[1],
        end_gene = gsub[length(gsub)],
        direction = dir_i,
        n_genes = seg_len,
        median_pos = median(psub),
        genes = list(gsub),
        role = sub_role,
        start_breakpoint = sub_start_bp,
        end_breakpoint = sub_end_bp
      )
      # increment new_blocks index
      nb_index <- nb_index + 1L
    }
    
    # bind all new_blocks for this block into a data.table
    new_dt <- data.table::rbindlist(new_blocks, use.names = TRUE)
    
    # append rows of new_dt to final_blocks
    for (r in seq_len(nrow(new_dt))) {
      # store each row into final_blocks
      final_blocks[[fb_index]] <- new_dt[r]
      # increment final_blocks index
      fb_index <- fb_index + 1L
    }
  }
  
  # bind all final_blocks rows into a single data.table
  out <- data.table::rbindlist(final_blocks, use.names = TRUE)
  # renumber block_id sequentially from 1 to number of rows
  out[, block_id := seq_len(.N)]
  # ensure columns are in expected order
  out <- out[, .(block_id,
                 start_index,
                 end_index,
                 start_gene,
                 end_gene,
                 direction,
                 n_genes,
                 median_pos,
                 genes,
                 role,
                 start_breakpoint,
                 end_breakpoint)]
  # return final updated block_data table
  return(out)
}
