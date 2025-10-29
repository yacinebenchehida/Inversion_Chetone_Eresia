fix_blocks <- function(dt_segment, inv_table,
                             pos_col_index = 9, gene_col_index = 1,
                             direction, max_gap = 2) {
  # Ensure synteny data is a data.table
  dt_segment <- data.table::as.data.table(dt_segment)     # convert input to data.table
  # Ensure inversion table is a data.table (or NULL allowed)
  inv_table  <- if (is.null(inv_table)) NULL else data.table::as.data.table(inv_table)  # normalize inv_table
  
  # Extract gene names in order (index = gene order within the window)
  genes <- as.character(dt_segment[[gene_col_index]])     # vector of gene identifiers
  # Extract numeric positions (genomic coordinates)
  positions <- as.numeric(dt_segment[[pos_col_index]])    # vector of coordinates
  
  # Sanity check: need enough data to define blocks robustly
  if (length(positions) < 6) stop("Not enough points to reconcile blocks (need at least 6 rows).")
  
  # Define the expected global step sign from the assessed global direction
  expected_sign <- if (direction == "positive") 1L else -1L  # +1 for increasing, -1 for decreasing
  
  # Precompute first differences and their signs for stepwise logic
  diffs <- diff(positions)                                   # step-wise delta (length n-1)
  step_sign <- sign(diffs)                                   # sign of each step: -1, 0, +1
  
  # ------------------------------
  # BLOCK DETECTION (with max_gap)
  # ------------------------------
  # We build blocks that follow the global direction, but we allow brief counter-steps
  # up to 'max_gap' in length. If a counter-run is short (<= max_gap), we normally keep
  # the same block. However, if after such a jump the trend resumes in the global
  # direction *on the other side of the starting level of the block* (your rule), we
  # split and start a new block at the resume point.
  
  # Initialize container for blocks (start/end indices in gene index space)
  blocks <- list()                                           # will collect (start_idx, end_idx)
  # Current block start (gene index)
  block_start <- 1L                                          # first gene starts block 1
  # Track current interruption against the global direction
  in_counter <- FALSE                                        # are we inside a counter-direction run?
  counter_len <- 0L                                          # length of the current counter run
  
  # Helper to finalize a block up to 'end_idx' and start a new one at 'new_start'
  finalize_start_new <- function(end_idx, new_start) {
    blocks[[length(blocks) + 1L]] <<- list(start_idx = block_start, end_idx = end_idx)  # save current block
    block_start <<- new_start                                                             # begin new block
    in_counter <<- FALSE                                                                  # reset counter state
    counter_len <<- 0L                                                                    # reset length
  }
  
  # Iterate across steps (between gene i and i+1)
  for (i in seq_along(step_sign)) {
    # Current step's sign
    s <- step_sign[i]                                  # -1, 0, or +1
    
    # Treat zero-step (identical positions) as aligned with global direction (neutral)
    aligned <- (s == expected_sign) || (s == 0L)       # TRUE if not opposing the global direction
    
    if (aligned) {
      # If we were in a counter-run, it ends here
      if (in_counter) {
        # We have just *resumed* the global trend at index i+1 (next gene)
        resume_idx <- i + 1L                           # candidate resume index
        # Apply your "offset" rule to decide split despite short interruption:
        # - If direction is negative: new coordinates after resume are "higher" than the FIRST point of the block → start new block
        # - If direction is positive: new coordinates after resume are "lower"  than the FIRST point of the block → start new block
        if (resume_idx <= length(positions)) {
          if (expected_sign == -1L && positions[resume_idx] > positions[block_start]) {
            # End old block just before resume, start new block at resume
            finalize_start_new(end_idx = resume_idx - 1L, new_start = resume_idx)
          } else if (expected_sign == 1L && positions[resume_idx] < positions[block_start]) {
            # End old block just before resume, start new block at resume
            finalize_start_new(end_idx = resume_idx - 1L, new_start = resume_idx)
          } else {
            # Offset rule not triggered → continue same block
            in_counter <- FALSE
            counter_len <- 0L
          }
        } else {
          # Defensive: resume beyond last index → close block at end
          in_counter <- FALSE
          counter_len <- 0L
        }
      }
      # Nothing else to do; the block continues
    } else {
      # Opposes the global direction → inside/extend a counter-run
      if (!in_counter) {
        in_counter <- TRUE                              # we just entered a counter run
        counter_len <- 1L                               # initialize its length
      } else {
        counter_len <- counter_len + 1L                 # extend counter-run
      }
      # If the counter-run exceeds the maximum allowed gap → force a split
      if (counter_len > max_gap) {
        # We end the current block before the counter-run started:
        end_of_block <- i - counter_len                 # last aligned index before counter
        # Next block starts where the counter-run began:
        new_block_start <- i - counter_len + 1L
        # Finalize and start new block at the first counter index
        finalize_start_new(end_idx = end_of_block, new_start = new_block_start)
        # Note: we stay in_counter=FALSE now; the loop will continue accumulating regular logic from here
      }
    }
  }
  # After the loop: close the last block at the final gene
  blocks[[length(blocks) + 1L]] <- list(start_idx = block_start, end_idx = length(positions))
  
  # Materialize blocks as a data.table
  blocks_dt <- data.table::rbindlist(blocks)                 # combine list to data.table
  # Compute each block's net direction from start→end positions
blocks_dt[, net_dir_sign := sign(
  vapply(seq_len(.N), function(i)
    positions[blocks_dt$end_idx[i]] - positions[blocks_dt$start_idx[i]],
    numeric(1)
  )
)]
# For exact ties (flat), coerce to expected sign to avoid spurious "mixed"
blocks_dt[net_dir_sign == 0L, net_dir_sign := expected_sign]  # For exact ties (flat), coerce to expected sign to avoid spurious "mixed"
  blocks_dt[net_dir_sign == 0L, net_dir_sign := expected_sign]                  # treat flat as aligned
  
  # ------------------------------------------------------
  # CASE A: All blocks follow the global direction (no inv)
  # ------------------------------------------------------
  all_aligned <- all(blocks_dt$net_dir_sign == expected_sign)  # TRUE if every block aligns
  
  # Initialize output inversion table (may be cleaned/unchanged)
  cleaned_inv <- if (is.null(inv_table)) NULL else data.table::copy(inv_table)  # safe copy
  # Initialize a table for *undefined translocations* inferred at jumps
  trans_table <- data.table::data.table(                                     # empty trans table
    translocation_id = integer(),
    gene_name       = character(),
    role            = character(),
    type            = character()
  )
  
  if (all_aligned) {
    # No true inversions: the “big jumps” between *adjacent blocks* are interpreted
    # as undefined translocations (per your rule #2).
    # We drop any existing inversions inside the window because they are artefacts.
    if (!is.null(cleaned_inv) && nrow(cleaned_inv) > 0L) {
      # Remove any inversion whose gene lives inside this dt_segment window
      in_window <- cleaned_inv$gene_name %in% genes
      cleaned_inv <- cleaned_inv[!in_window]                                  # drop those rows
    }
    # Create one undefined translocation at each *block boundary* (i.e., start of block k>1)
    if (nrow(blocks_dt) >= 2L) {
      # Translocation id counter
      trans_id <- 1L
      for (k in 2:nrow(blocks_dt)) {
        # The "jump" occurs at the start of block k
        start_gene <- genes[blocks_dt$start_idx[k]]                           # gene where block k begins
        # Append a single-row undefined translocation (distinct name!)
        trans_table <- data.table::rbind(
          trans_table,
          data.table::data.table(
            translocation_id = trans_id,                 # new translocation id
            gene_name       = start_gene,                # the boundary gene
            role            = "undefined translocation", # **exact label required**
            type            = "simple translocation"     # same type vocabulary as elsewhere
          ),
          use.names = TRUE
        )
        trans_id <- trans_id + 1L
      }
    }
    
    # Return early: no inversions after reconciliation, only undefined translocations added
    return(list(inv_table = cleaned_inv,
                trans_table = trans_table,
                blocks = blocks_dt))
  }
  
  # ------------------------------------------------------
  # CASE B: Mixed directions across blocks
  # ------------------------------------------------------
  # We keep blocks that align with the global direction unchanged.
  # For blocks *opposite* to the global direction:
  #   - remove any inversions that sit *inside* the block (not at its edges),
  #   - then assign an "undefined" inversion at the block’s END gene (your rule #3).
  # Additionally: inversions that were inferred solely because "type" flipped
  # (simple inversion vs inversion+translocation) inside a block and are not at edges
  # are removed as well (same interior removal catches those).
  
  # If there is no inversion table yet, create an empty one with the expected schema
  if (is.null(cleaned_inv)) {
    cleaned_inv <- data.table::data.table(
      inversion_id = integer(),
      gene_name    = character(),
      role         = character(),
      type         = character()
    )
  } else {
    # Normalize schema (ensure required columns exist)
    need_cols <- c("inversion_id","gene_name","role","type")
    miss_cols <- setdiff(need_cols, names(cleaned_inv))
    if (length(miss_cols) > 0L) stop("inv_table is missing required column(s): ", paste(miss_cols, collapse = ", "))
  }
  
  # Map inversion genes to indices in the current window for fast interior checks
  gene_to_idx <- match(cleaned_inv$gene_name, genes)                          # NA if outside window
  
  # Remove interior inversions for each *opposite* block
  for (k in seq_len(nrow(blocks_dt))) {
    # Skip blocks aligned with the global direction
    if (blocks_dt$net_dir_sign[k] == expected_sign) next
    
    # Current block bounds (inclusive)
    b_start <- blocks_dt$start_idx[k]
    b_end   <- blocks_dt$end_idx[k]
    
    # Identify inversion rows that map inside this block
    inside_block <- !is.na(gene_to_idx) & gene_to_idx >= b_start & gene_to_idx <= b_end
    if (any(inside_block)) {
      # Keep only inversions sitting *at* the block edges; drop interior ones
      at_edge <- gene_to_idx[inside_block] %in% c(b_start, b_end)             # TRUE if at block boundary
      # Build logical mask to DROP interior (inside_block & !at_edge)
      drop_mask <- logical(nrow(cleaned_inv))
      drop_mask[which(inside_block)] <- !at_edge
      if (any(drop_mask)) {
        cleaned_inv <- cleaned_inv[!drop_mask]                                 # drop interior artefacts
        gene_to_idx <- match(cleaned_inv$gene_name, genes)                     # refresh mapping after drop
      }
    }
    
    # Ensure there is an "undefined" inversion at the block *end* (per your rule)
    end_gene <- genes[b_end]                                                  # gene name at block end
    # Check if an undefined at this end already exists
    exists_undef <- nrow(cleaned_inv[gene_name == end_gene & role == "undefined"]) > 0L
    if (!exists_undef) {
      # Compute next inversion_id (compact reindexing)
      next_id <- if (nrow(cleaned_inv) == 0L) 1L else (max(cleaned_inv$inversion_id, na.rm = TRUE) + 1L)
      # Append undefined inversion at block end
      cleaned_inv <- data.table::rbind(
        cleaned_inv,
        data.table::data.table(
          inversion_id = next_id,            # new ID
          gene_name    = end_gene,           # breakpoint gene at block end
          role         = "undefined",        # inversion undefined (NOT translocation)
          type         = "simple inversion"  # keep vocabulary consistent
        ),
        use.names = TRUE
      )
      # Refresh mapping vector (keep it coherent if further blocks need processing)
      gene_to_idx <- match(cleaned_inv$gene_name, genes)
    }
  }
  
  # ---------------------------
  # RETURN reconciled structures
  # ---------------------------
  return(list(inv_table  = cleaned_inv,   # cleaned/adjusted inversions
              trans_table = trans_table,  # undefined translocations added in CASE A (empty here otherwise)
              blocks      = blocks_dt))   # diagnostic blocks (start_idx, end_idx, net_dir_sign)
}
