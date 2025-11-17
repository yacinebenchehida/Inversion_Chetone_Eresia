library(ggplot2)
library(data.table)
library(dplyr)


##########################
# Test a specific genome #
##########################
setwd("/Users/yacinebenchehida/Desktop/Convergent_evolution/Chetone_histrio/R_package")
# Read BLAST output
dt <- fread("busco/Eresia_pelonia_busco.tsv",header = FALSE)
filtered_dt <- filter_contigs_by_proportion(dt,2, threshold =0.2,outlier = TRUE)
my_dir <- assess_direction(filtered_dt)
blocks <- identify_blocks(dt= filtered_dt,  pos_col_index = 9, gene_col_index = 1, direction = my_dir, max_gap = 2,resume_confirm = 1)

ann_block <- detect_inversions_translo_using_blocks(dt= filtered_dt, direction = my_dir, pos_col_index = 9, gene_col_index = 1, min_consecutive = 2, block_data = blocks)
ann_block <- fix_small_edge_blocks(ann_block)
ann_block <- resolve_smallblock_bridge(block_annotation=ann_block, filtered_dt=filtered_dt)
ann_block <- detect_and_fix_backward_boundary(
  dt = filtered_dt,
  pos_col_index = 10,
  gene_col_index = 1,
  block_data = ann_block)

plot_genes(dt = filtered_dt, block_annotation = ann_block, pos_col_index = 9, gene_col_index = 1, plot_title = "Gene order plot")
check_gene_in_block(block_data=ann_block, gene_order="genes_order.txt", gene_name="HMEL000025-RB")


###########################################
# Plotting busco genes no event detection #
###########################################
setwd("/Users/yacinebenchehida/Desktop/Convergent_evolution/Chetone_histrio/R_package/busco")
busco_files <- list.files(pattern = "*busco.tsv", full.names = TRUE)
get_prefix <- function(x) sub("^(GCA_[0-9.]+).*", "\\1", x)

busco_df <- data.frame(busco = busco_files, key = get_prefix(basename(busco_files)))
noms <- busco_df[2]
busco_df <- busco_df[1]

lapply(seq_len(nrow(busco_df)), function(i) {
  
  busco <- busco_df$busco[i]
  gca_name <- noms[i, 1]
  print(gca_name)
  message(sprintf("[%d/%d] Processing %s...", i, nrow(busco_df), gca_name))
  dt_busco <- fread(busco, header = FALSE)
  
  tryCatch({
    suppressMessages(suppressWarnings({
      filtered_dt_busco <- filter_contigs_by_proportion(
        dt_busco,
        contig_col_index = 2,
        threshold = 0.7,
        outlier = TRUE
      )
      my_dir <- assess_direction(filtered_dt_busco)
      blocks <- identify_blocks(
        dt = filtered_dt_busco,
        pos_col_index = 10,
        gene_col_index = 1,
        direction = my_dir,
        max_gap = 2,
        resume_confirm = 1
      )
      ann_block <- detect_inversions_translo_using_blocks(
        dt = filtered_dt_busco,
        direction = my_dir,
        pos_col_index = 10,
        gene_col_index = 1,
        min_consecutive = 2,
        block_data = blocks
      )
      
      window_start   <- "Hmel215003o:833935-840617"
      cortex         <- "Hmel215003o:1413941-1423696"
      betafructo     <- "Hmel215003o:1335401-1342966"
      extra_gene     <- "Hmel215003o:1785377-1787042"
      window_end     <- "Hmel215003o:2237273-2243166"
      
      gene_names <- as.character(filtered_dt_busco[[1]])
      
      start_idx      <- which(gene_names == window_start)
      end_idx        <- which(gene_names == window_end)
      cortex_idx     <- which(gene_names == cortex)
      betafructo_idx <- which(gene_names == betafructo)
      extra_idx      <- which(gene_names == extra_gene)
      
      rect_df <- data.frame(
        xmin = start_idx - 0.5,
        xmax = end_idx   + 0.5,
        ymin = -Inf,
        ymax = Inf
      )
      
      highlight_cortex     <- filtered_dt_busco[gene_names == cortex, ]
      highlight_betafructo <- filtered_dt_busco[gene_names == betafructo, ]
      highlight_extra      <- filtered_dt_busco[gene_names == extra_gene, ]
      
      pdf(paste0(gca_name, ".pdf"))
      
      p <- plot_genes(
        dt = filtered_dt_busco,
        pos_col_index = 10,
        gene_col_index = 1,
        plot_title = gca_name
      )
      
      p <- p +
        geom_rect(
          data        = rect_df,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          inherit.aes = FALSE,
          fill        = "blue",
          color       = "black",
          alpha       = 0.15
        ) +
        
        # cortex (red)
        geom_point(
          data        = highlight_cortex,
          aes_string(x = names(filtered_dt_busco)[1],
                     y = names(filtered_dt_busco)[10]),
          inherit.aes = FALSE,
          color       = "red",
          size        = 2
        ) +
        
        # betafructo (dark green)
        geom_point(
          data        = highlight_betafructo,
          aes_string(x = names(filtered_dt_busco)[1],
                     y = names(filtered_dt_busco)[10]),
          inherit.aes = FALSE,
          color       = "gold",
          size        = 2
        ) +
        
        # extra gene (purple)
        geom_point(
          data        = highlight_extra,
          aes_string(x = names(filtered_dt_busco)[1],
                     y = names(filtered_dt_busco)[10]),
          inherit.aes = FALSE,
          color       = "purple",
          size        = 2
        )
      
      print(p)
      dev.off()
    }))
  }, error = function(e) {
    message(sprintf("Error processing %s: %s", gca_name, e$message))
    return(NULL)
  }, warning = function(w) {
    message(sprintf("Warning in %s: %s", gca_name, w$message))
    invokeRestart("muffleWarning")
  }, finally = {
    message(sprintf("Finished %s", gca_name))
  })
  
})


###########################
# Pipeline on busco genes #
###########################
setwd("/Users/yacinebenchehida/Desktop/Convergent_evolution/Chetone_histrio/R_package/busco")
busco_files <- list.files(pattern = "*busco.tsv", full.names = TRUE)
get_prefix <- function(x) sub("^(GC[AF]_[0-9.]+).*", "\\1", x)

busco_df <- data.frame(busco = busco_files, key = get_prefix(basename(busco_files)))
noms <- busco_df[2]
busco_df <- busco_df[1]

good_window_counter <- 0 
dropped_counter <- 0 
dropped_counter_error <- 0
colinear_counter <- 0
ivory_in_inversion <- 0
ivory_in_translocation <- 0
ivory_in_inv_trans <- 0

inversion_results_busco <- lapply(seq_len(nrow(busco_df)), function(i) {
  
  busco <- busco_df$busco[i]
  gca_name <- noms[i, 1]
  print(gca_name)
  message(sprintf("[%d/%d] Processing %s...", i, nrow(busco_df), gca_name))
  dt_busco <- fread(busco, header = FALSE)
  
  tryCatch({
    suppressMessages(suppressWarnings({
      
      filtered_dt_busco <- tryCatch(
        filter_contigs_by_proportion(
          dt_busco,
          contig_col_index = 2,
          threshold = 0.7,
          outlier = TRUE,
          min_bit_score = 80
        ),
        error = function(e) {
          dropped_counter_error <<- dropped_counter_error + 1
          message(sprintf("filter_contigs_by_proportion error for %s: %s", gca_name, e$message))
          return(NULL)
        }
      )
      
      if(is.null(filtered_dt_busco) || nrow(filtered_dt_busco) == 0){
        dropped_counter <<- dropped_counter + 1
        print(paste(dropped_counter," dropped windows",sep=""))
      }
      
      if (!is.null(filtered_dt_busco) && nrow(filtered_dt_busco) > 0) {
        
      good_window_counter <<- good_window_counter + 1
      my_dir <- assess_direction(filtered_dt_busco)
      blocks <- identify_blocks(
        dt = filtered_dt_busco,
        pos_col_index = 10,
        gene_col_index = 1,
        direction = my_dir,
        max_gap = 2,
        resume_confirm = 1
      )
      ann_block <- detect_inversions_translo_using_blocks(
        dt = filtered_dt_busco,
        direction = my_dir,
        pos_col_index = 10,
        gene_col_index = 1,
        min_consecutive = 2,
        block_data = blocks
      )
      
      ann_block <- fix_small_edge_blocks(ann_block)
      
      ann_block <- resolve_smallblock_bridge(block_annotation=ann_block, 
        filtered_dt=filtered_dt_busco)
      
      ann_block <- detect_and_fix_backward_boundary(
       dt = filtered_dt_busco,
       pos_col_index = 10,
       gene_col_index = 1,
       block_data = ann_block)
      
      the_event <- events_counter(ann_block)
      
      if(all(the_event=="colinear")){
        colinear_counter <<- colinear_counter + 1
      }
      
      # Multiple event outputs allowed
      if (any(the_event == "inversion")) {
          ann_inv <- ann_block[role == "inversion"]
          is_ivory_in_inversion <- check_gene_in_block(block_data=ann_inv, gene_order="../genes_order_busco.txt", gene_name="Hmel215003o:1413941-1423696")
          
          if(is_ivory_in_inversion==TRUE){
            ivory_in_inversion <<- ivory_in_inversion + 1
          }
        }
      
      if (any(the_event == "translocation")) {
        ann_trans <- ann_block[role == "translocation"]
        is_ivory_in_translocation <- check_gene_in_block(block_data=ann_trans, gene_order="../genes_order_busco.txt", gene_name="Hmel215003o:1413941-1423696")
        
        if(is_ivory_in_translocation==TRUE){
          ivory_in_translocation <<- ivory_in_translocation + 1
        }
      }
      
      if (any(the_event == "inversion+translocation")) {
          ann_inv_trans <- ann_block[role == "inversion+translocation"]
          is_ivory_in_inv_trans <- check_gene_in_block(block_data=ann_inv_trans, gene_order="../genes_order_busco.txt", gene_name="Hmel215003o:1413941-1423696")
          
          if(is_ivory_in_inv_trans==TRUE){
            ivory_in_inv_trans <<- ivory_in_inv_trans + 1
          }
       }
      # RETURN VALUE FOR SUCCESSFUL CASE
      if (!is.null(ann_block) && nrow(ann_block) > 0) {
        pdf(file.path("Pipeline_plots/", paste0(gca_name, ".pdf")),10,10)
        p <- plot_genes(dt = filtered_dt_busco, block_annotation = ann_block, pos_col_index = 10, gene_col_index = 1, plot_title = gca_name)
        print(p)
        dev.off()
        return(list(name = gca_name, data = ann_block))
        
       }
      }
   }
   ))
    
  }, error = function(e) {
    message(sprintf("Error processing %s: %s", gca_name, e$message))
    return(NULL)
  }, warning = function(w) {
    message(sprintf("Warning in %s: %s", gca_name, w$message))
    invokeRestart("muffleWarning")
  }, finally = {
    message(sprintf("Finished %s", gca_name))
  })
  
})

percent_kept <- round((good_window_counter / nrow(busco_df)) * 100, 1)
percent_colinear <- round((colinear_counter / good_window_counter) * 100, 1)
percent_no_colinear <- 100-percent_colinear
total_window_with_ivory <- ivory_in_inversion + ivory_in_translocation + ivory_in_inv_trans
percent_with_ivory <-  round((total_window_with_ivory/(good_window_counter-colinear_counter))*100, 1)


message(sprintf("A total of %d out of %d windows were kept after filtering (%s%%).",
                good_window_counter, nrow(busco_df), percent_kept))

message(sprintf("Of these, %d out of %d are colinear windows (%s%%).",
                colinear_counter, good_window_counter, percent_colinear))

message(sprintf("And %d out of %d contains a structural variant (%s%%).",
                good_window_counter-colinear_counter, good_window_counter, percent_no_colinear))

message(sprintf("Of these, %d out of %d contain ivory (%s%%).",
                total_window_with_ivory, good_window_counter-colinear_counter, percent_with_ivory))

message(sprintf("Of these, %d out of %d are simple inversion (%s%%).",
                ivory_in_inversion, total_window_with_ivory, round((ivory_in_inversion/total_window_with_ivory)*100)))

message(sprintf("And %d out of %d are translocated inversion (%s%%).",
                ivory_in_inv_trans, total_window_with_ivory, round((ivory_in_inv_trans/total_window_with_ivory)*100)))

message(sprintf("And %d out of %d are translocation (%s%%).",
                ivory_in_translocation, total_window_with_ivory, round((ivory_in_translocation/total_window_with_ivory)*100)))

message(sprintf("In total, %d out of %d analysed windows contain ivory (%s%%).",
                total_window_with_ivory, good_window_counter, round((total_window_with_ivory/good_window_counter)*100)))



# Remove NULL results
inversion_results_busco <- Filter(Negate(is.null), inversion_results_busco)

# Combine all inversions
all_inversions <- rbindlist(
  lapply(inversion_results_busco, function(x) {
    x$data[, GCA := x$name]
  })
)

# Build single DT
inversion_results_df_busco <- rbindlist(
  lapply(inversion_results_busco, function(x) {
    df <- x$data
    df[, GCA := x$name]
    setcolorder(df, c("GCA", setdiff(names(df), "GCA")))
    df
  }),
  fill = TRUE
)

dt <- inversion_results_df_busco

# Step 1: gather breakpoints into long format
bp <- melt(
  dt,
  id.vars = c("GCA", "role"),
  measure.vars = c("start_breakpoint", "end_breakpoint"),
  value.name = "breakpoint"
)

# Step 2: filter to only meaningful breakpoints
bp <- bp[
  breakpoint != "Undefined" &
    breakpoint != "Not_applicable" &
    role %in% c("inversion", "inversion+translocation", "translocation")
]

# Step 3: count occurrences of each breakpoint by role
breakpoint_summary <- bp[, .N, by = .(breakpoint, role)]

# ---- 1) READ GENE ORDER ----
genes_order <- fread("../genes_order_busco.txt", header = FALSE)
setnames(genes_order, "V1", "gene")
genes_order[, pos := .I]

# ---- 2) ENSURE data.table TYPE ----
breakpoint_summary <- as.data.table(breakpoint_summary)

# ---- 3) KEEP ONLY EVENT ROLES (NO 'colinear') ----
event_roles <- c("inversion", "inversion+translocation", "translocation")
breakpoint_summary <- breakpoint_summary[role %in% event_roles]

# ---- 4) SPLIT BREAKPOINT INTO GENE1/GENE2 ----
breakpoint_summary[, c("g1","g2") := tstrsplit(breakpoint, "_", fixed = TRUE)]

# ---- 5) MAP GENE POSITIONS FROM genes_order ----
breakpoint_summary <- merge(breakpoint_summary, genes_order, by.x = "g1", by.y = "gene", all.x = TRUE)
setnames(breakpoint_summary, "pos", "pos1")
breakpoint_summary <- merge(breakpoint_summary, genes_order, by.x = "g2", by.y = "gene", all.x = TRUE)
setnames(breakpoint_summary, "pos", "pos2")

# ---- 6) DETERMINE ADJACENCY (NEIGHBOR GENES) ----
breakpoint_summary[, is_adjacent := abs(pos1 - pos2) == 1]

# ---- 7) DEFINE ORDER ALONG GENOME FOR X-AXIS ----
breakpoint_summary[, order_pos := pmin(pos1, pos2)]
setorder(breakpoint_summary, order_pos, breakpoint, role)

# ---- 8) BUILD UNIQUE ORDERED LEVELS (FIXES DUPLICATE-LEVEL ERRORS) ----
ordered_levels <- unique(breakpoint_summary$breakpoint)

# ---- 9) FACTOR FOR X-AXIS USING UNIQUE LEVELS ONLY ----
breakpoint_summary[, breakpoint_f := factor(breakpoint, levels = ordered_levels)]

# ---- MULTI-GENE HIGHLIGHTING (3 GENES, 3 COLORS) ----
highlight_genes <- c(
  betafructo = "Hmel215003o:1335401-1342966",
  cortex = "Hmel215003o:1413941-1423696",
  HMEL00049     = "Hmel215003o:1785377-1787042"
)

breakpoint_summary[, highlight_gene := NA_character_]

breakpoint_summary[
  grepl(highlight_genes["betafructo"], breakpoint, fixed = TRUE),
  highlight_gene := "betafructo"
]

breakpoint_summary[
  grepl(highlight_genes["cortex"], breakpoint, fixed = TRUE),
  highlight_gene := "cortex"
]

breakpoint_summary[
  grepl(highlight_genes["HMEL00049"], breakpoint, fixed = TRUE),
  highlight_gene := "HMEL00049"
]

highlight_cols <- c(
  betafructo = "blue",
  cortex = "red",
  HMEL00049     = "gold"
)

# ---- 10) PER-LABEL COLOR FOR PLOT 2 (BLACK=ADJACENT, RED=NON-ADJACENT) ----
lab_key <- unique(breakpoint_summary[, .(breakpoint, is_adjacent)])[
  match(ordered_levels, unique(breakpoint))
]
label_colors <- ifelse(lab_key$is_adjacent, "black", "red")

# ---- 11) CONSISTENT ROLE COLORS (STACKED BARS) ----
role_cols <- c(
  inversion = "#1b9e77",
  "inversion+translocation" = "#d95f02",
  translocation = "#7570b3"
)

# =====================================================
#   PLOT 1: ONLY ADJACENT BREAKPOINTS (x ordered by genome)
# =====================================================
bp_adj <- breakpoint_summary[is_adjacent == TRUE]
adj_levels <- unique(bp_adj$breakpoint)
bp_adj[, breakpoint_f := factor(breakpoint, levels = adj_levels)]

p1 <- ggplot(bp_adj, aes(x = breakpoint_f, y = N/2, fill = role)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = role_cols) +
  geom_point(
    data = bp_adj[!is.na(highlight_gene)],
    aes(x = breakpoint_f, y = N + 20, colour = highlight_gene),
    size = 1
  ) +
  scale_colour_manual(values = highlight_cols, drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)
  ) +
  labs(x = "Breakpoint (adjacent only)", y = "Count", fill = "Role", colour = "Gene")

ggsave("breakpoints_only_with_ajacent_genes.pdf", p1, width = 21, height = 7, units = "in")

# =====================================================
#   PLOT 2: ALL BREAKPOINTS, NON-ADJACENT LABELS IN RED
# =====================================================
p2 <- ggplot(breakpoint_summary, aes(x = breakpoint_f, y = N, fill = role)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = role_cols) +
  geom_point(
    data = breakpoint_summary[!is.na(highlight_gene)],
    aes(x = breakpoint_f, y = N + 20, colour = highlight_gene),
    size = 1
  ) +
  scale_colour_manual(values = highlight_cols, drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 0.5,
      colour = label_colors
    )
  ) +
  labs(x = "Breakpoint (all; non-adjacent labels in red)", y = "Count", fill = "Role", colour = "Gene")

ggsave("breakpoints_all_with_nonadjacent.pdf", p2, width = 30, height = 8, units = "in")



############################################
# Keep only genomes with ivory in an event #
############################################
setwd("/Users/yacinebenchehida/Desktop/Convergent_evolution/Chetone_histrio/R_package/busco")
busco_files <- list.files(pattern = "*busco.tsv", full.names = TRUE)
get_prefix <- function(x) sub("^(GC[AF]_[0-9.]+).*", "\\1", x)
busco_df <- data.frame(busco = busco_files, key = get_prefix(basename(busco_files)))
noms <- busco_df[2]
busco_df <- busco_df[1]

target_gene <- "Hmel215003o:1413941-1423696"

good_window_counter <- 0 
dropped_counter <- 0 
dropped_counter_error <- 0
colinear_counter <- 0
ivory_in_inversion <- 0
ivory_in_translocation <- 0
ivory_in_inv_trans <- 0

inversion_results_busco <- lapply(seq_len(nrow(busco_df)), function(i) {
  busco <- busco_df$busco[i]
  gca_name <- noms[i, 1]
  dt_busco <- fread(busco, header = FALSE)
  
  tryCatch({
    suppressMessages(suppressWarnings({
      
      filtered_dt_busco <- tryCatch(
        filter_contigs_by_proportion(
          dt_busco,
          contig_col_index = 2,
          threshold = 0.7,
          outlier = TRUE,
          min_bit_score = 80
        ),
        error = function(e) {
          dropped_counter_error <<- dropped_counter_error + 1
          return(NULL)
        }
      )
      
      if(is.null(filtered_dt_busco) || nrow(filtered_dt_busco) == 0){
        dropped_counter <<- dropped_counter + 1
      }
      
      if (!is.null(filtered_dt_busco) && nrow(filtered_dt_busco) > 0) {
        
        good_window_counter <<- good_window_counter + 1
        my_dir <- assess_direction(filtered_dt_busco)
        
        blocks <- identify_blocks(
          dt = filtered_dt_busco,
          pos_col_index = 10,
          gene_col_index = 1,
          direction = my_dir,
          max_gap = 2,
          resume_confirm = 1
        )
        
        ann_block <- detect_inversions_translo_using_blocks(
          dt = filtered_dt_busco,
          direction = my_dir,
          pos_col_index = 10,
          gene_col_index = 1,
          min_consecutive = 2,
          block_data = blocks
        )
        
        ann_block <- fix_small_edge_blocks(ann_block)
        ann_block <- resolve_smallblock_bridge(ann_block, filtered_dt_busco)
        ann_block <- detect_and_fix_backward_boundary(filtered_dt_busco, 10, 1, ann_block)
        
        the_event <- events_counter(ann_block)
        
        keep_this_window <- FALSE
        
        if (any(the_event %in% c("inversion"))) {
          ann_inv <- ann_block[role == "inversion"]
          if (check_gene_in_block(ann_inv, "../genes_order_busco.txt", target_gene)) keep_this_window <- TRUE
        }
        
        if (any(the_event %in% c("translocation"))) {
          ann_trans <- ann_block[role == "translocation"]
          if (check_gene_in_block(ann_trans, "../genes_order_busco.txt", target_gene)) keep_this_window <- TRUE
        }
        
        if (any(the_event %in% c("inversion+translocation"))) {
          ann_inv_trans <- ann_block[role == "inversion+translocation"]
          if (check_gene_in_block(ann_inv_trans, "../genes_order_busco.txt", target_gene)) keep_this_window <- TRUE
        }
        
        if (keep_this_window) {
          return(list(name = gca_name, data = ann_block))
        }
      }
    }))
  }, error = function(e) {
    return(NULL)
  })
})

inversion_results_busco <- Filter(Negate(is.null), inversion_results_busco)

inversion_results_df_busco <- rbindlist(
  lapply(inversion_results_busco, function(x) {
    df <- x$data
    df[, GCA := x$name]
    setcolorder(df, c("GCA", setdiff(names(df), "GCA")))
    df
  }),
  fill = TRUE
)

dt <- inversion_results_df_busco

bp <- melt(
  dt,
  id.vars = c("GCA", "role"),
  measure.vars = c("start_breakpoint", "end_breakpoint"),
  value.name = "breakpoint"
)

bp <- bp[
  breakpoint != "Undefined" &
    breakpoint != "Not_applicable" &
    role %in% c("inversion", "inversion+translocation", "translocation")
]

breakpoint_summary <- bp[, .N, by = .(breakpoint, role)]

genes_order <- fread("../genes_order_busco.txt", header = FALSE)
setnames(genes_order, "V1", "gene")
genes_order[, pos := .I]

breakpoint_summary <- as.data.table(breakpoint_summary)

event_roles <- c("inversion", "inversion+translocation", "translocation")
breakpoint_summary <- breakpoint_summary[role %in% event_roles]

breakpoint_summary[, c("g1","g2") := tstrsplit(breakpoint, "_", fixed = TRUE)]

breakpoint_summary <- merge(breakpoint_summary, genes_order, by.x = "g1", by.y = "gene", all.x = TRUE)
setnames(breakpoint_summary, "pos", "pos1")
breakpoint_summary <- merge(breakpoint_summary, genes_order, by.x = "g2", by.y = "gene", all.x = TRUE)
setnames(breakpoint_summary, "pos", "pos2")

breakpoint_summary[, is_adjacent := abs(pos1 - pos2) == 1]

breakpoint_summary[, order_pos := pmin(pos1, pos2)]
setorder(breakpoint_summary, order_pos, breakpoint, role)

ordered_levels <- unique(breakpoint_summary$breakpoint)
breakpoint_summary[, breakpoint_f := factor(breakpoint, levels = ordered_levels)]

highlight_genes <- c(
  betafructo = "Hmel215003o:1335401-1342966",
  cortex = "Hmel215003o:1413941-1423696",
  HMEL00049 = "Hmel215003o:1785377-1787042"
)

breakpoint_summary[, highlight_gene := NA_character_]

breakpoint_summary[
  grepl(highlight_genes["betafructo"], breakpoint, fixed = TRUE),
  highlight_gene := "betafructo"
]

breakpoint_summary[
  grepl(highlight_genes["cortex"], breakpoint, fixed = TRUE),
  highlight_gene := "cortex"
]

breakpoint_summary[
  grepl(highlight_genes["HMEL00049"], breakpoint, fixed = TRUE),
  highlight_gene := "HMEL00049"
]

highlight_cols <- c(
  betafructo = "blue",
  cortex = "red",
  HMEL00049 = "gold"
)

lab_key <- unique(breakpoint_summary[, .(breakpoint, is_adjacent)])[
  match(ordered_levels, unique(breakpoint))
]
label_colors <- ifelse(lab_key$is_adjacent, "black", "red")

role_cols <- c(
  inversion = "#1b9e77",
  "inversion+translocation" = "#d95f02",
  translocation = "#7570b3"
)

bp_adj <- breakpoint_summary[is_adjacent == TRUE]
adj_levels <- unique(bp_adj$breakpoint)
bp_adj[, breakpoint_f := factor(breakpoint, levels = adj_levels)]

p1 <- ggplot(bp_adj, aes(x = breakpoint_f, y = N/2, fill = role)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = role_cols) +
  geom_point(
    data = bp_adj[!is.na(highlight_gene)],
    aes(x = breakpoint_f, y = N + 20, colour = highlight_gene),
    size = 1
  ) +
  scale_colour_manual(values = highlight_cols, drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)
  ) +
  labs(x = "breakpoint (adjacent only)", y = "count", fill = "role", colour = "gene")

ggsave("breakpoints_adj_only_filtered.pdf", p1, width = 12, height = 6, units = "in")

p2 <- ggplot(breakpoint_summary, aes(x = breakpoint_f, y = N, fill = role)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = role_cols) +
  geom_point(
    data = breakpoint_summary[!is.na(highlight_gene)],
    aes(x = breakpoint_f, y = N + 20, colour = highlight_gene),
    size = 1
  ) +
  scale_colour_manual(values = highlight_cols, drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 0.5,
      colour = label_colors
    )
  ) +
  labs(x = "breakpoint (all; nonadjacent labels red)", y = "count", fill = "role", colour = "gene")

ggsave("breakpoints_all_filtered.pdf", p2, width = 20, height = 6, units = "in")
