phylo_correction <- function(tree, inv_dt, threshold_ace = 0.9, outgroup = NULL, plot_tree=FALSE, rename_acc_2_sp = NULL) {
  # -------------------------------------------------------------
  # 1. Prepare trait vector for this inversion
  # -------------------------------------------------------------
  inv <- inv_dt$present                         # 0/1 vector
  names(inv) <- inv_dt$GCA                      # names must be GC IDs
  
  tree <- ape::read.tree(tree)                  # tree is a Newick file path
  
  # -------------------------------------------------------------------------
  # 2. Extract species name instead of GC names and root the tree if provided
  # -------------------------------------------------------------------------
  if (!is.null(rename_acc_2_sp)) {

  rename_dt <- read.table(
    rename_acc_2_sp,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  colnames(rename_dt) <- c("accession", "species")

  ## add _1, _2, ... to duplicated species names
  rename_dt$species_unique <- paste0(
    rename_dt$species, "_",
    ave(rename_dt$species, rename_dt$species, FUN = seq_along)
  )

  ## build renaming vector USING species_unique
  rename_vec <- rename_dt$species_unique
  names(rename_vec) <- rename_dt$accession

  ## rename tree tips
  new_labels <- rename_vec[tree$tip.label]
  keep <- !is.na(new_labels)
  tree <- ape::drop.tip(tree, tree$tip.label[!keep])
  tree$tip.label <- new_labels[keep]

  ## rename inversion vector
  names(inv) <- rename_vec[names(inv)]
  inv <- inv[!is.na(names(inv))]
}
  
  
  if (!is.null(outgroup)) {                     # drop outgroup tip(s) if provided
    outgroup_present <- intersect(outgroup, tree$tip.label)
    if (length(outgroup_present) > 0L) {
      tree <- ape::drop.tip(tree, outgroup_present)
    }
  }
  
  # match tree to data
  tree <- ape::keep.tip(tree, intersect(tree$tip.label, names(inv)))
  inv  <- inv[tree$tip.label]                   # reorder to match tree
  inv  <- ifelse(inv == 1, "TRUE", "FALSE")     # convert to "TRUE"/"FALSE"
  
  # ---------------------------
  # 3. Fit Mk model & run ASR
  # ---------------------------
  inv_ard <- fitMk(tree, inv, model = "ARD", pi = "fitzjohn")
  inv_ancr <- ancr(inv_ard)
  
  tree <- attr(inv_ancr, "tree")
  ace  <- inv_ancr$ace
  ace[ace < 0] <- 0

  
  Ntip <- length(tree$tip.label)
  Nnode <- tree$Nnode
  
  # ---------------------------
  # 4. Helper functions
  # ---------------------------
  get_children <- function(node, tree) tree$edge[tree$edge[,1] == node, 2]
  
  get_descendant_tips <- function(node, tree, Ntip) {
    kids <- get_children(node, tree)
    tips <- kids[kids <= Ntip]
    internals <- kids[kids > Ntip]
    for (k in internals) tips <- c(tips, get_descendant_tips(k, tree, Ntip))
    unique(tips)
  }
  
  get_parent <- function(node, tree) {
    idx <- which(tree$edge[,2] == node)
    if (!length(idx)) return(NA_integer_)
    tree$edge[idx,1]
  }
  
  # Terminal monophyly check: opaque tips must be terminal (no branching)
  is_terminal_monophyletic <- function(tips, tree, Ntip) {
    if (length(tips) == 0) return(TRUE)
    if (length(tips) == 1) return(TRUE)
    mrca <- ape::getMRCA(tree, tips)
    if (is.null(mrca) || is.na(mrca)) return(FALSE)
    desc <- get_descendant_tips(mrca, tree, Ntip)
    identical(sort(desc), sort(tips))
  }
  
  # ---------------------------
  # 5. Tip states & ACE
  # ---------------------------
  tip_state <- inv[tree$tip.label]
  inv_col <- which(colnames(ace) == "TRUE")
  ace_inv <- ace[, inv_col]
  
  # ---------------------------
  # 6. Sanity check
  # ---------------------------
  inv_tips_observed <- which(tip_state == "TRUE")
  internal_rows_idx <- which(ace_inv >= threshold_ace)
  inv_nodes_internal <- internal_rows_idx + Ntip
  inv_nodes_full <- c(inv_nodes_internal, inv_tips_observed)
  
  cat("\n================ SANITY CHECK ================\n")
  cat("Observed inverted tips: ", length(inv_tips_observed), "\n")
  if (length(inv_tips_observed) > 0) {
    cat("Labels: ", paste(tree$tip.label[inv_tips_observed], collapse=", "), "\n")
  }
  cat("Internal nodes with ACE≥", threshold_ace, ": ", length(inv_nodes_internal), "\n")
  cat("Combined candidate nodes: ", length(inv_nodes_full), "\n")
  cat("===============================================\n\n")
  
  # ---------------------------
  # 7. Candidate internal nodes (ACE >= threshold + terminal-monophyly)
  # ---------------------------
  candidate_nodes <- integer(0)
  candidate_nodes_info <- list()
  
  for (row_idx in internal_rows_idx) {
    nd <- row_idx + Ntip
    desc_tips <- get_descendant_tips(nd, tree, Ntip)
    opaque_tips <- desc_tips[tip_state[desc_tips] == "FALSE"]
    # Keep node if opaque tips are all terminal
    if (is_terminal_monophyletic(opaque_tips, tree, Ntip)) {
      candidate_nodes <- c(candidate_nodes, nd)
      candidate_nodes_info[[as.character(nd)]] <- list(
        node = nd,
        n_tips = length(desc_tips),
        n_trans = sum(tip_state[desc_tips]=="TRUE"),
        tip_labels = tree$tip.label[desc_tips]
      )
    }
  }
  
  cat("Candidate internal nodes after terminal-monophyly filter: ", length(candidate_nodes), "\n")
  
  # ---------------------------
  # 8. Basal internal nodes (most basal per cluster)
  # ---------------------------
  basal_internal <- integer(0)
  
  for (nd in sort(unique(candidate_nodes))) {
    # check ancestors
    p <- get_parent(nd, tree)
    ancestors <- integer(0)
    while (!is.na(p)) { ancestors <- c(ancestors, p); p <- get_parent(p, tree) }
    if (!any(ancestors %in% candidate_nodes)) basal_internal <- c(basal_internal, nd)
  }
  
  cat("Basal internal nodes: ", length(basal_internal), "\n")
  
  # ---------------------------
  # 9. Tip-only transparent units (not covered by basal nodes)
  # ---------------------------
  is_tip_covered <- function(tip_idx, basal_nodes, tree, Ntip) {
    for (bn in basal_nodes) if (tip_idx %in% get_descendant_tips(bn, tree, Ntip)) return(TRUE)
    FALSE
  }
  
  basal_tip_only <- integer(0)
  for (tip_idx in inv_tips_observed) {
    if (!is_tip_covered(tip_idx, basal_internal, tree, Ntip)) basal_tip_only <- c(basal_tip_only, tip_idx)
  }
  
  # ---------------------------
  # 10. Final counts & summary
  # ---------------------------
  n_internal <- length(basal_internal)
  n_tiponly  <- length(basal_tip_only)
  total_origins <- n_internal + n_tiponly
  
  cat("\n========== INVERSION ORIGIN SUMMARY ==========\n")
  cat("Basal internal origins: ", n_internal, "\n")
  cat("Tip-only single origins: ", n_tiponly, "\n")
  cat("TOTAL independent origins: ", total_origins, "\n")
  
  if (n_internal > 0) {
    cat("Basal internal node IDs: ", paste(basal_internal, collapse=", "), "\n")
  }
  
  if (n_tiponly > 0) {
    cat("Tip-only labels: ", paste(tree$tip.label[basal_tip_only], collapse=", "), "\n")
  }
  cat("===============================================\n\n")
  
  # ---------------------------
  # 11. Plot ASR tree with Transparent highlights
  # ---------------------------
  if(plot_tree) {

    cols <- c("TRUE"="grey", "FALSE"="goldenrod")
    pie_data <- ace[, c("TRUE","FALSE"), drop=FALSE]

    node.cex <- rep(0.2, nrow(pie_data))
    node.cex[pie_data[, "TRUE"] >= threshold_ace] <- 0.5

    tip_labels <- paste0(tree$tip.label, " (", tip_state, ")")
    inv_ancr$ace[inv_ancr$ace < 0] <- 0
    if (!is.null(inv_ancr$lik.anc)) {
      inv_ancr$lik.anc[inv_ancr$lik.anc < 0] <- 0
    }

    inv_ancr$ace <- inv_ancr$ace / rowSums(inv_ancr$ace)

    inv_ancr$tree$tip.label <- tip_labels

    plot(
      inv_ancr,
      args.plotTree=list(size=0.3),
      args.nodelabels=list(pie=pie_data, piecol=cols, cex=0.1),
      args.tiplabels=list(cex=0.1),
      legend=FALSE
    )
    
    mtext(
  paste0("After phylogenetic correction: ", total_origins, " independent origins"),
  side = 1,        # bottom margin
  line = 3,        # distance from plot
  cex  = 1.2
	)
    
    nodelabels(cex=0.4, frame="none", adj=c(1.2, -0.4))
  }
    
    return(total_origins)
}
