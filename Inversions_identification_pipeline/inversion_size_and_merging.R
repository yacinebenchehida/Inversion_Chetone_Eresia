library(data.table)                                                         # load data.table
library(ggplot2)                                                            # load ggplot2

gene_order <- fread("genes_order_busco.txt", header = FALSE)               # load gene order file
setnames(gene_order, "V1", "gene_id")                                      # set column name to gene_id
gene_order[, idx := .I]                                                    # assign index to each gene

dt_clean <- dt[                                                            # clean dt from invalid rows
  start_breakpoint != "Undefined" &
  end_breakpoint   != "Undefined" &
  start_breakpoint != "Not_applicable" &
  end_breakpoint   != "Not_applicable"
]

unique_breakpoints <- unique(dt_clean[, .(start_breakpoint, end_breakpoint)]) # keep unique breakpoint pairs

unique_breakpoints[, c("g1","g2") := tstrsplit(start_breakpoint, "_")]     # left breakpoint: g1 g2
unique_breakpoints[, c("g3","g4") := tstrsplit(end_breakpoint, "_")]       # right breakpoint: g3 g4

unique_breakpoints[, g1_idx := gene_order$idx[match(g1, gene_order$gene_id)]] # index of g1
unique_breakpoints[, g2_idx := gene_order$idx[match(g2, gene_order$gene_id)]] # index of g2
unique_breakpoints[, g3_idx := gene_order$idx[match(g3, gene_order$gene_id)]] # index of g3
unique_breakpoints[, g4_idx := gene_order$idx[match(g4, gene_order$gene_id)]] # index of g4

unique_breakpoints[, left_start  := pmin(g1_idx, g2_idx)]                  # left interval start
unique_breakpoints[, left_end    := pmax(g1_idx, g2_idx)]                  # left interval end
unique_breakpoints[, right_start := pmin(g3_idx, g4_idx)]                  # right interval start
unique_breakpoints[, right_end   := pmax(g3_idx, g4_idx)]                  # right interval end

n <- nrow(unique_breakpoints)                                              # number of inversions
adj <- matrix(FALSE, n, n)                                                 # adjacency matrix

for (i in 1:(n - 1)) {                                                     # loop over first index
  for (j in (i + 1):n) {                                                   # loop over second index
    left_overlap  <- max(unique_breakpoints$left_start[i],
                         unique_breakpoints$left_start[j]) <
                    min(unique_breakpoints$left_end[i],
                         unique_breakpoints$left_end[j])                   # strict overlap on left
    right_overlap <- max(unique_breakpoints$right_start[i],
                         unique_breakpoints$right_start[j]) <
                    min(unique_breakpoints$right_end[i],
                         unique_breakpoints$right_end[j])                  # strict overlap on right
    if (left_overlap && right_overlap) {                                   # require overlap both sides
      adj[i, j] <- TRUE                                                    # mark adjacency i→j
      adj[j, i] <- TRUE                                                    # mark adjacency j→i
    }
  }
}

cluster_id <- rep(NA_integer_, n)                                          # cluster id for each inversion
current_cluster <- 0                                                       # cluster counter

for (i in 1:n) {                                                           # loop over inversions
  if (is.na(cluster_id[i])) {                                              # if not yet assigned
    current_cluster <- current_cluster + 1                                 # start new cluster
    queue <- c(i)                                                          # initialize queue with i
    cluster_id[i] <- current_cluster                                       # assign cluster id
    while (length(queue) > 0) {                                            # BFS over adjacency
      k <- queue[1]                                                        # take first element
      queue <- queue[-1]                                                   # remove from queue
      neighbors <- which(adj[k, ] & is.na(cluster_id))                     # unassigned neighbors
      if (length(neighbors) > 0) {                                         # if neighbors exist
        cluster_id[neighbors] <- current_cluster                           # assign same cluster
        queue <- c(queue, neighbors)                                       # add neighbors to queue
      }
    }
  }
}

unique_breakpoints[, cluster_id := cluster_id]                             # store cluster id in table

cluster_rep <- unique_breakpoints[                                         # cluster-level merged inversion
  ,
  {
    left_start_cluster  <- min(left_start)                                 # min left start
    left_end_cluster    <- max(left_end)                                   # max left end
    right_start_cluster <- min(right_start)                                # min right start
    right_end_cluster   <- max(right_end)                                  # max right end
    g1_merged <- gene_order$gene_id[left_start_cluster]                    # merged g1 gene
    g2_merged <- gene_order$gene_id[left_end_cluster]                      # merged g2 gene
    g3_merged <- gene_order$gene_id[right_start_cluster]                   # merged g3 gene
    g4_merged <- gene_order$gene_id[right_end_cluster]                     # merged g4 gene
    start_idx <- left_start_cluster                                        # start index for ordering
    span_idx  <- right_end_cluster - left_start_cluster                    # total span for ordering
    list(g1 = g1_merged,
         g2 = g2_merged,
         g3 = g3_merged,
         g4 = g4_merged,
         start_idx = start_idx,
         span_idx  = span_idx)
  },
  by = cluster_id
]

setorder(cluster_rep, start_idx, span_idx)                                 # order clusters by start then span

all_inv_list <- list()                                                     # list to collect all inversions
list_i <- 0                                                                 # index for list elements

for (i in seq_len(nrow(cluster_rep))) {                                    # loop over clusters in order
  cid <- cluster_rep$cluster_id[i]                                         # current cluster id
  
  merged_row <- data.table(                                                # merged inversion row
    cluster_id = cid,
    type = "merged",
    g1 = cluster_rep$g1[i],
    g2 = cluster_rep$g2[i],
    g3 = cluster_rep$g3[i],
    g4 = cluster_rep$g4[i]
  )
  
  list_i <- list_i + 1                                                     # increment list index
  all_inv_list[[list_i]] <- merged_row                                     # store merged row
  
  orig_rows <- unique_breakpoints[cluster_id == cid]                       # all original inversions
  orig_rows[, left_for_order := left_start]                                # use left_start for order
  orig_rows[, span_for_order  := right_end - left_start]                   # use span for order
  setorder(orig_rows, left_for_order, span_for_order)                      # order originals inside cluster
  
  orig_dt <- orig_rows[, .(cluster_id = cluster_id,
                           type = "original",
                           g1 = g1,
                           g2 = g2,
                           g3 = g3,
                           g4 = g4)]                                      # keep required columns
  
  list_i <- list_i + 1                                                     # increment list index
  all_inv_list[[list_i]] <- orig_dt                                       # store original rows
}

all_inv <- rbindlist(all_inv_list, use.names = TRUE, fill = TRUE)         # bind all rows into one table
all_inv[, line_id := .I]                                                  # one line per inversion

plot_dt <- melt(                                                          # long format for plot
  all_inv,
  id.vars = c("cluster_id","type","line_id"),                             # identifiers
  measure.vars = c("g1","g2","g3","g4"),                                  # four breakpoint genes
  variable.name = "pos",                                                  # position label
  value.name = "gene_id"                                                  # gene id column
)

plot_dt[, gene_id := factor(gene_id, levels = gene_order$gene_id)]        # enforce x-axis gene order

pdf("overlap_merged_vs_original.pdf", 40, 80)                             # open PDF for plotting

ggplot(plot_dt,
       aes(x = gene_id,
           y = factor(line_id),
           group = line_id,
           color = type)) +                                               # color by type
  geom_line(linewidth = 0.7) +                                            # draw line between 4 genes
  geom_point(size = 3) +                                                  # draw gene points
  scale_color_manual(values = c("original" = "red", "merged" = "blue")) + # set colors
  theme_classic() +                                                       # clean theme
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),       # rotate x labels
    axis.text.y = element_blank(),                                        # hide y text
    axis.ticks.y = element_blank(),                                       # hide y ticks
    axis.title.y = element_blank()                                        # remove y title
  )

dev.off()                                                                 # close PDF device
