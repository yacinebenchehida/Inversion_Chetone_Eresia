identify_blocks <- function(dt, inv_table,
                            pos_col_index = 9, gene_col_index = 1,
                            direction = c("positive","negative"), max_gap = 2,
                            resume_confirm = 2) {
  # Ensure data.table objects (inv_table kept for API compatibility)
  dt        <- data.table::as.data.table(dt)          # main window
  inv_table <- data.table::as.data.table(inv_table)   # unused here
  
  # Validate direction argument (kept for interface consistency)
  direction <- match.arg(direction)                   # "positive" | "negative"
  
  # Extract working vectors
  genes     <- as.character(dt[[gene_col_index]])     # gene ids
  positions <- as.numeric(dt[[pos_col_index]])        # coordinates
  
  # Handle trivial sizes
  n <- length(positions)                              # number of genes
  if (n == 0) {
    return(data.table::data.table(
      block_id    = integer(0), start_index = integer(0), end_index = integer(0),
      start_gene  = character(0), end_gene  = character(0), direction = character(0),
      n_genes     = integer(0), mean_step = numeric(0), median_pos = numeric(0),
      genes       = I(list())
    ))
  }
  if (n == 1) {
    return(data.table::data.table(
      block_id    = 1, start_index = 1, end_index = 1,
      start_gene  = genes[1], end_gene = genes[1], direction = "flat",
      n_genes     = 1, mean_step = 0, median_pos = positions[1],
      genes       = list(genes[1])
    ))
  }
  
  # Map delta -> sign in {-1,0,1}
  step_sign <- function(dx) {
    if (is.na(dx)) return(0)                          # treat NA as flat
    if (dx > 0)  return(1)                            # up
    if (dx < 0)  return(-1)                           # down
    return(0)                                         # flat
  }
  
  # Precompute step signs
  diffs <- diff(positions)                            # successive deltas
  signs <- vapply(diffs, step_sign, numeric(1))       # vector of -1/0/1
  
  # All-flat: single flat block
  if (all(signs == 0)) {
    return(data.table::data.table(
      block_id=1, start_index=1, end_index=n,
      start_gene=genes[1], end_gene=genes[n], direction="flat",
      n_genes=n, mean_step=0, median_pos=stats::median(positions),
      genes=list(genes[1:n])
    ))
  }
  
  # Initialize first block direction as first non-zero sign
  first_nonzero <- which(signs != 0)[1]               # index in signs
  cur_sign      <- signs[first_nonzero]               # -1 or +1
  
  # Initialize current block start and anchor (FIRST value rule)
  block_start <- 1                                    # start index of current block
  block_first <- positions[block_start]               # first coordinate of block
  
  # Output collector
  blocks   <- list()                                  # list of data.tables
  block_id <- 1                                       # running block id
  
  # Scan through genes
  i <- 2                                              # current gene index
  while (i <= n) {
    s <- step_sign(positions[i] - positions[i - 1])   # sign of current step
    
    # Flat step → continue
    if (s == 0) { i <- i + 1; next }
    
    # Same direction as current block → continue
    if (s == cur_sign) { i <- i + 1; next }
    
    # Opposite sign → begin a potential bump
    bump_start <- i                                   # first gene of bump
    bump_len   <- 1                                   # bump length (steps)
    lookahead  <- i + 1                               # probe pointer
    
    # Grow bump up to max_gap while steps keep opposing cur_sign
    while (lookahead <= n && bump_len <= max_gap) {
      s_next <- step_sign(positions[lookahead] - positions[lookahead - 1])
      if (s_next == cur_sign) break                   # direction resumed
      bump_len  <- bump_len + 1                       # extend bump
      lookahead <- lookahead + 1                      # advance probe
    }
    
    # If direction resumes within max_gap, require confirmation + anchor rule to MERGE
    if (bump_len <= max_gap && lookahead <= n &&
        step_sign(positions[lookahead] - positions[lookahead - 1]) == cur_sign) {
      
      # Confirm resumption: need >= resume_confirm consecutive steps with cur_sign
      stable_forward <- 1                             # we have at least lookahead step
      j <- lookahead + 1                              # continue confirming
      while (j <= n) {
        s_j <- step_sign(positions[j] - positions[j - 1])
        if (s_j != cur_sign) break                    # confirmation ended
        stable_forward <- stable_forward + 1          # extend confirmation
        if (stable_forward >= resume_confirm) break   # enough confirmation
        j <- j + 1
      }
      
      # Anchor check: resumption point must be beyond FIRST value on correct side
      anchor_ok <- if (cur_sign > 0) {
        positions[lookahead] > block_first            # positive block → above FIRST
      } else {
        positions[lookahead] < block_first            # negative block → below FIRST
      }
      
      # If both pass → MERGE bump, jump scan forward, keep block
      if (stable_forward >= resume_confirm && anchor_ok) {
        i <- lookahead + stable_forward               # skip confirmed run
        next
      }
      # Else → fall through to SPLIT
    }
    
    # SPLIT: close current block right BEFORE the bump
    end_idx <- bump_start - 1                         # last pre-bump gene
    if (end_idx < block_start) end_idx <- block_start # guard
    blk_idx <- block_start:end_idx                    # indices in this block
    
    # Emit finished block
    blocks[[length(blocks) + 1]] <- data.table::data.table(
      block_id    = block_id,
      start_index = block_start,
      end_index   = end_idx,
      start_gene  = genes[block_start],
      end_gene    = genes[end_idx],
      direction   = if (cur_sign > 0) "positive" else "negative",
      n_genes     = length(blk_idx),
      mean_step   = if (length(blk_idx) >= 2) mean(diff(positions[blk_idx])) else 0,
      median_pos  = stats::median(positions[blk_idx]),
      genes       = list(genes[blk_idx])
    )
    
    # Decide where the NEW block should start (singleton-spike suppression)
    # Default: start new block AT bump_start with new direction = s
    new_start <- bump_start                           # candidate start of new block
    new_sign  <- s                                    # candidate direction
    
    # If the very next step immediately flips AGAIN (singleton spike),
    # then skip the spike and start at bump_start+1 with that next sign.
    if (bump_start + 1 <= n) {
      s_after <- step_sign(positions[bump_start + 1] - positions[bump_start])
      if (s_after != 0 && s_after != new_sign) {
        new_start <- bump_start                    # drop singleton spike
        new_sign  <- s_after                          # take the sustained direction
      }
    }
    
    # Start new block at decided position
    block_id   <- block_id + 1                        # next id
    block_start <- new_start                          # new block start index
    block_first <- positions[block_start]             # reset FIRST anchor
    cur_sign    <- new_sign                           # new block direction
    
    # Continue scanning from the next index after new_start
    i <- max(new_start + 1, i + 1)                    # advance safely
  }
  
  # Close the final (open) block to the end
  blk_idx <- block_start:n                            # indices for last block
  blocks[[length(blocks) + 1]] <- data.table::data.table(
    block_id    = block_id,
    start_index = block_start,
    end_index   = n,
    start_gene  = genes[block_start],
    end_gene    = genes[n],
    direction   = if (cur_sign > 0) "positive" else "negative",
    n_genes     = length(blk_idx),
    mean_step   = if (length(blk_idx) >= 2) mean(diff(positions[blk_idx])) else 0,
    median_pos  = stats::median(positions[blk_idx]),
    genes       = list(genes[blk_idx])
  )
  
  # Bind and return
  out <- data.table::rbindlist(blocks, use.names = TRUE)
  data.table::setcolorder(out, c("block_id","start_index","end_index",
                                 "start_gene","end_gene","direction",
                                 "n_genes","mean_step","median_pos","genes"))
  return(out[])
}
