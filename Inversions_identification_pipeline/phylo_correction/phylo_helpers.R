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
