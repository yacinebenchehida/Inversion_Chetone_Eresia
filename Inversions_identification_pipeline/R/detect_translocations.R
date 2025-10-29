detect_translocations <- function(dt_segment, inv_table,
                                  pos_col_index = 9, gene_col_index = 1,
                                  min_consecutive = 3, direction) {
  # Ensure inputs are data.table for fast, consistent ops
  dt_segment <- as.data.table(dt_segment)   # convert synteny segment to data.table
  inv_table  <- as.data.table(inv_table)    # convert inversion table to data.table
  
  # Pull convenient vectors (names unchanged)
  genes     <- as.character(dt_segment[[gene_col_index]])   # gene identifiers (ordered)
  positions <- as.numeric(dt_segment[[pos_col_index]])      # genomic coordinates (numeric)
  
  # Guard: need enough rows to compute first/last medians and detect runs
  if (nrow(dt_segment) < 6) stop("Not enough points to detect translocations (need ≥ 6).")
  
  # Define in-bounds using medians of the first 3 and the last 3 positions
  med_first3 <- median(positions[1:3])                                         # median of first 3
  med_last3  <- median(positions[(length(positions) - 2):length(positions)])   # median of last 3
  lower_bound <- min(med_first3, med_last3)                                    # inclusive lower bound
  upper_bound <- max(med_first3, med_last3)                                    # inclusive upper bound
  
  # Precompute per-point "side":
  #   -1 = below in-bounds (positions < lower_bound)
  #    0 = inside in-bounds (lower_bound ≤ pos ≤ upper_bound)
  #   +1 = above in-bounds (positions > upper_bound)
  side <- integer(length(positions))                              # allocate integer vector
  side[positions <  lower_bound] <- -1L                           # mark below
  side[positions >  upper_bound] <-  1L                           # mark above
  # points already 0L are inside
  
  # Expected overall trend per your global direction:
  # positive → net increase; negative → net decrease
  expected_sign <- if (direction == "positive")  1L else -1L      # +1 or -1
  
  # Prepare output table (schema fixed)
  results_translocation_table <- data.table(
    translocation_id = integer(),  # sequential id per detected translocation
    gene_name        = character(),# start/end gene names
    role             = character(),# "start"/"end"
    type             = character() # "simple translocation"
  )
  
  # Helper: finalize a candidate run [s:e] after all validations
  finalize_run <- function(s, e) {
    # Ensure length criterion
    if ((e - s + 1L) < min_consecutive) return(invisible(NULL))   # too short → ignore
    
    # Ensure all points in the run stay on the **same side** (all -1 or all +1)
    run_side <- unique(side[s:e])                                  # unique side codes in run
    run_side <- run_side[run_side != 0L]                           # drop any 0 (should not exist by construction)
    if (length(run_side) != 1L) return(invisible(NULL))            # mixed sides → split earlier; discard here if it happens
    
    # Ensure **overall direction** matches: net delta sign must match expected_sign
    net_sign <- sign(positions[e] - positions[s])                  # sign of start→end position change
    if (is.na(net_sign) || net_sign == 0L || net_sign != expected_sign) {
      return(invisible(NULL))                                      # no movement or wrong overall direction → ignore
    }
    
    # Build gene block and check overlap with inversions
    block_genes <- genes[s:e]                                      # genes in this candidate block
    inv_genes   <- inv_table$gene_name                             # genes referenced by inversions (start/end/undefined)
    overlap     <- intersect(block_genes, inv_genes)               # any overlap?
    
    if (length(overlap) > 0L) {
      # Overlap is allowed **only** if it is confined to the block's edges (start/end genes)
      edge_genes <- c(genes[s], genes[e])                          # allowed overlap set
      if (!all(overlap %in% edge_genes)) {
        return(invisible(NULL))                                    # interior overlap with inversion → discard
      }
      # If we are here, overlap (if any) touches edges only → allowed
    }
    
    # Record start/end for this translocation (type fixed string per your spec)
    results_translocation_table <<- rbind(
      results_translocation_table,
      data.table(
        translocation_id = trans_id,                               # current id
        gene_name        = c(genes[s], genes[e]),                  # start then end
        role             = c("start", "end"),                      # roles
        type             = c("simple translocation", "simple translocation")  # same type on both rows
      ),
      use.names = TRUE
    )
    
    # Increment global counter in parent frame
    trans_id <<- trans_id + 1L                                     # next id
  }
  
  # Scan through points, splitting runs whenever we hit in-bounds (0) or flip side (-1 ↔ +1)
  trans_id   <- 1L                         # initialize translocation id counter
  run_start  <- NA_integer_                # current run start index (NA = no run open)
  run_side   <- NA_integer_                # current run side (-1 or +1)
  
  # Iterate all indices
  for (i in seq_along(side)) {
    si <- side[i]                          # side at position i
    
    if (si == 0L) {                        # in-bounds → close any open run
      if (!is.na(run_start)) finalize_run(run_start, i - 1L)  # finalize previous run [run_start .. i-1]
      run_start <- NA_integer_             # reset run start
      run_side  <- NA_integer_             # reset run side
      next                                      # continue scanning
    }
    
    # si is out-of-bounds (-1 or +1)
    if (is.na(run_start)) {
      # start a new run
      run_start <- i                       # mark run start
      run_side  <- si                      # fix the side for this run
    } else {
      # run already open → check side consistency
      if (si != run_side) {
        # side flipped (below ↔ above) → finalize old run, start a new one at i
        finalize_run(run_start, i - 1L)    # finalize up to previous index
        run_start <- i                     # new run starts here
        run_side  <- si                    # new run side
      }
      # else: same side → continue the run
    }
  }
  
  # After loop: close any trailing run that reached the end
  if (!is.na(run_start)) finalize_run(run_start, length(side))     # finalize last open run
  
  # Return results (could be empty if nothing valid passes)
  return(results_translocation_table)
}
