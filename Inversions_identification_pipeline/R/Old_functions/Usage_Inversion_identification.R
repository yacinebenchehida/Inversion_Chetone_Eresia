library(ggplot2)
library(data.table)
library(future.apply)

# args <- commandArgs(trailingOnly=TRUE)
folder_path <- "../Inputs"
# List all files starting with GCA
gca_files <- list.files(folder_path, pattern = "^GCA.*", full.names = TRUE)
gca_files <- head(gca_files, 1000)

# Load functions
source("filter_contigs_by_proportion.R")
source("assess_direction.R")
source("detect_inversions.R")
source("wrapper_detect_inversions.R")

# cpus
#plan(multisession, workers = 1)

# Read BLAST output
#dt <- fread(args[1],header = FALSE)
#filtered_dt <- filter_contigs_by_proportion(dt,2, threshold =0.4)
#direction <- assess_direction(filtered_dt)
#whatever <- wrapper_detect_inversions(filtered_dt, min_consecutive=3, plot_window=FALSE, boundary_genes = c("HMEL000058-RA","HMEL000026-RA"))

inversion_results <- lapply(seq_along(gca_files), function(i) { # use future_lapply for parallelisation but speed gain only huge data otherwise even slower due to overhead.
  file <- gca_files[i]
  gca_name <- tools::file_path_sans_ext(basename(file))
  message(sprintf("[%d/%d] Processing %s...", i, length(gca_files), gca_name))
  dt <- fread(file, header = FALSE)
  result <- tryCatch({
    suppressMessages(suppressWarnings({
    filtered_dt <- filter_contigs_by_proportion(dt, contig_col_index = 2, threshold = 0.4)
    
    if (!is.null(filtered_dt) && nrow(filtered_dt) > 0) {
      inv <- wrapper_detect_inversions(
        filtered_dt,
        min_consecutive = 3,
        plot_window = FALSE,
        boundary_genes = c("HMEL000058-RA","HMEL000026-RA")
      )
      
      if (!is.null(inv) && nrow(inv) > 0) {
        gca_name <- tools::file_path_sans_ext(basename(file))
        return(list(name = gca_name, data = inv))
      }
    }
    
    NULL
    }))
  }, error = function(e) {
    message(sprintf("Skipping %s: %s", basename(file), e$message))
    NULL  # skip this file on error
  })
  
  return(result)
})

# Remove NULL results (files with no inversions)
inversion_results <- Filter(Negate(is.null), inversion_results)

# Combine all inversions into one data.table with GCA column
all_inversions <- rbindlist(
  lapply(inversion_results, function(x) {
    x$data[, GCA := x$name]
  })
)

# Check how many species have inversions
# inversion_results


# Convert to a single data.frame with GCA column first
inversion_results_df <- rbindlist(
  lapply(inversion_results, function(x) {
    df <- x$data
    df[, GCA := x$name]   # add GCA column
    setcolorder(df, c("GCA", setdiff(names(df), "GCA"))) # GCA as first col
    return(df)
  }),
  fill = TRUE
)

inversion_results_df$unique_breakpoints <- paste(inversion_results_df$gene_name,inversion_results_df$role,sep="_")

# Convert table to data frame
df <- as.data.frame(table(inversion_results_df$unique_breakpoints))
colnames(df) <- c("Entry", "Frequency")
df$Entry <- factor(df$Entry, levels = df$Entry)
df[df$Frequency > 0,]

# Dot plot with labels on the left
pdf("brekapoint_frequencies.pdf")
ggplot(df, aes(x = Frequency, y = Entry)) +
  geom_point(size = 1.5, color = "red") +
  geom_segment(aes(x = 0, xend = Frequency, y = Entry, yend = Entry),
               color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Frequencies of entries",
       x = "Count", y = "Entry") +
  theme(axis.text.y = element_text(size =8))
dev.off()







library(ggplot2)
library(data.table)

setwd("/Users/yacinebenchehida/Desktop/Convergent_evolution/Chetone_histrio/R_package")
# Read BLAST output
dt <- fread("mmseq2/GCA_910589475.2_best_hits.reordered.tsv",header = FALSE)
dt <- dt[-c((nrow(dt)-2):nrow(dt)),]
dispersion_score(dt, pos_col_index = 9, chaos_excess_score = 2)
filtered_dt <- filter_contigs_by_proportion(dt,2, threshold =0.7,outlier = TRUE)
my_dir <- assess_direction(filtered_dt)

detect_inversions(filtered_dt, pos_col_index = 9, gene_col_index = 1, min_consecutive = 3, direction=my_dir)
print(wrapper_detect_inversions(filtered_dt, min_consecutive=3, plot_window=TRUE))
invers_table <- wrapper_detect_inversions(filtered_dt, min_consecutive=3, plot_window=TRUE)
invers_table <- merge_oversplit_inversions(invers_table, filtered_dt, pos_col_index = 9, gene_col_index = 1, my_dir)
invers_table <- Fix_fake_interruptions(inv_table = invers_table, dt = filtered_dt, pos_col_index = 9, gene_col_index = 1, my_dir, max_gap = 2)

detect_translocations(dt_segment = filtered_dt, inv_table = invers_table, pos_col_index = 9, gene_col_index = 1, min_consecutive = 3, direction = my_dir)











library(ggplot2)
library(data.table)

setwd("/Users/yacinebenchehida/Desktop/Convergent_evolution/Chetone_histrio/R_package")
# Read BLAST output
dt <- fread("Inputs/GCA_963082685.1_ilCycAlbi1.1.txt",header = FALSE)
dispersion_score(dt, pos_col_index = 9, chaos_excess_score = 2)
filtered_dt <- filter_contigs_by_proportion(dt,2, threshold =0.7,outlier = TRUE)
my_dir <- assess_direction(filtered_dt)
invers_table <- wrapper_detect_inversions(filtered_dt, min_consecutive=3, plot_window=TRUE)
invers_table <- merge_oversplit_inversions(invers_table, filtered_dt, pos_col_index = 9, gene_col_index = 1, my_dir)
invers_table <- Fix_fake_interruptions(inv_table = invers_table, dt = filtered_dt, pos_col_index = 9, gene_col_index = 1, my_dir, max_gap = 2)
blocks <- identify_blocks(dt= filtered_dt, inv_table = invers_table, pos_col_index = 9, gene_col_index = 1, direction = my_dir, max_gap = 2,resume_confirm = 2)
identify_blocks(dt= filtered_dt, inv_table = invers_table, pos_col_index = 9, gene_col_index = 1, direction = my_dir, max_gap = 2,resume_confirm = 2)

ann_block <- detect_inversions_translo_using_blocks(dt= filtered_dt, direction = my_dir, pos_col_index = 9, gene_col_index = 1, min_consecutive = 3, block_data = blocks)
plot_genes(dt = filtered_dt, block_annotation = ann_block, pos_col_index = 9, gene_col_index = 1, plot_title = "Gene order plot")





library(ggplot2)
library(data.table)

setwd("/Users/yacinebenchehida/Desktop/Convergent_evolution/Chetone_histrio/R_package")
# Read BLAST output
dt <- fread("busco/GCA_965615795.1_ilMomAlpi1.hap1.1_genomic.fna_busco.tsv",header = FALSE)
dt <- dt[-c((nrow(dt)-10):nrow(dt)),]
filtered_dt <- filter_contigs_by_proportion(dt,2, threshold =0.7,outlier = TRUE)
my_dir <- assess_direction(filtered_dt)
dt <- fread("mmseq2/GCA_965615795.1_best_hits.reordered.tsv",header = FALSE)
dt <- dt[-c((nrow(dt)-10):nrow(dt)),]
filtered_dt <- filter_contigs_by_proportion(dt,2, threshold =0.7,outlier = TRUE)
invers_table <- wrapper_detect_inversions(filtered_dt, min_consecutive=3, plot_window=TRUE)
invers_table <- merge_oversplit_inversions(invers_table, filtered_dt, pos_col_index = 9, gene_col_index = 1, my_dir)
invers_table <- Fix_fake_interruptions(inv_table = invers_table, dt = filtered_dt, pos_col_index = 9, gene_col_index = 1, my_dir, max_gap = 2)
blocks <- identify_blocks(dt= filtered_dt,  pos_col_index = 9, gene_col_index = 1, direction = my_dir, max_gap = 2,resume_confirm = 1)
#identify_blocks(dt= filtered_dt, inv_table = invers_table, pos_col_index = 9, gene_col_index = 1, direction = my_dir, max_gap = 2,resume_confirm = 2)

ann_block <- detect_inversions_translo_using_blocks(dt= filtered_dt, direction = my_dir, pos_col_index = 9, gene_col_index = 1, min_consecutive = 2, block_data = blocks)
ann_block
ann_block <- fix_small_edge_blocks(ann_block)
ann_block <- resolve_smallblock_bridge(block_annotation=ann_block, filtered_dt=filtered_dt)

plot_genes(dt = filtered_dt, block_annotation = ann_block, pos_col_index = 9, gene_col_index = 1, plot_title = "Gene order plot")

events_counter(ann_block)
is_ambiguous(ann_block)


setwd("/Users/yacinebenchehida/Desktop/Convergent_evolution/Chetone_histrio/R_package/")
folder_path <- "/Users/yacinebenchehida/Desktop/Convergent_evolution/Chetone_histrio/R_package/mmseq2/"
get_prefix <- function(x) sub("^(GCA_[0-9.]+).*", "\\1", x)

busco_files <- list.files("busco", full.names = TRUE)
mmseq_files <- list.files("mmseq2", full.names = TRUE)

busco_df <- data.frame(busco = busco_files, key = get_prefix(basename(busco_files)))
mmseq_df <- data.frame(mmseq = mmseq_files, key = get_prefix(basename(mmseq_files)))

pairs <- merge(busco_df, mmseq_df, by = "key")
busco_df <- pairs[2]
mmseq_df <- pairs[3]

good_window_counter <- 0 
dropped_counter <- 0 
complex_window <- 0
colinear_counter <- 0


inversion_results <- lapply(seq_len(nrow(mmseq_df)), function(i) {
  
  file <- mmseq_df[i,]
  busco <- busco_df[i,]
  gca_name <- pairs[i,1]
  message(sprintf("[%d/%d] Processing %s...", i, nrow(mmseq_df), gca_name))
  
  dt <- fread(file, header = FALSE)
  dt <- dt[-c((nrow(dt)-10):nrow(dt)),]
  
  
  dt_busco <- fread(busco, header = FALSE)
  
  result <- tryCatch({
    suppressMessages(suppressWarnings({
      
      filtered_dt <- filter_contigs_by_proportion(dt, contig_col_index = 2, threshold = 0.7, outlier = TRUE)
      filtered_dt_busco <- filter_contigs_by_proportion(dt_busco, contig_col_index = 2, threshold = 0.7, outlier = TRUE)
      
      
      if(is.null(filtered_dt)){
        dropped_counter <<- dropped_counter + 1
        print(paste(dropped_counter," dropped windows",sep=""))
      }
      
      my_dir <- assess_direction(filtered_dt)
      
      if (!is.null(filtered_dt) && nrow(filtered_dt) > 0) {
        
        good_window_counter <<- good_window_counter + 1
        blocks <- identify_blocks(dt = filtered_dt, pos_col_index = 10, gene_col_index = 1, direction = my_dir, max_gap = 2, resume_confirm = 1)
        ann_block <- detect_inversions_translo_using_blocks(dt = filtered_dt, direction = my_dir, pos_col_index = 10, gene_col_index = 1, min_consecutive = 2, block_data = blocks)
        ann_block <-fix_small_edge_blocks(ann_block)
        ann_block <- resolve_smallblock_bridge(block_annotation=ann_block, filtered_dt=filtered_dt)
        the_event <- events_counter(ann_block)
        
        if(all(the_event=="colinear")){
          colinear_counter <<- colinear_counter + 1
        }
        
        # Multiple event outputs allowed
        if (any(the_event == "inversion")) {
          if(is_ambiguous(ann_block)){
            my_dir <- assess_direction(filtered_dt_busco) 
            blocks <- identify_blocks(dt = filtered_dt, pos_col_index = 10, gene_col_index = 1, direction = my_dir, max_gap = 2, resume_confirm = 1)
            ann_block <- detect_inversions_translo_using_blocks(dt = filtered_dt, direction = my_dir, pos_col_index = 10, gene_col_index = 1, min_consecutive = 2, block_data = blocks)
            ann_block <- fix_small_edge_blocks(ann_block)
            ann_block <- resolve_smallblock_bridge(block_annotation=ann_block, filtered_dt=filtered_dt)
            
            pdf(file.path("pdf_outputs/inversions", paste0(gca_name, ".pdf")))
            p  <- plot_genes(dt = filtered_dt, block_annotation = ann_block, pos_col_index = 10, gene_col_index = 1, plot_title = gca_name)
            print(p)
            dev.off()
          }else{
            pdf(file.path("pdf_outputs/inversions", paste0(gca_name, ".pdf")))
            p  <- plot_genes(dt = filtered_dt, block_annotation = ann_block, pos_col_index = 10, gene_col_index = 1, plot_title = gca_name)
            print(p)
            dev.off()
          }
        }
        
        if (any(the_event == "translocation")) {
          pdf(file.path("pdf_outputs/translocations", paste0(gca_name, ".pdf")))
          p <- plot_genes(dt = filtered_dt, block_annotation = ann_block, pos_col_index = 10, gene_col_index = 1, plot_title = gca_name)
          print(p)
          dev.off()
        }
        
        if (any(the_event == "inversion+translocation")) {
          if(is_ambiguous(ann_block)){
            my_dir <- assess_direction(filtered_dt_busco) 
            blocks <- identify_blocks(dt = filtered_dt, pos_col_index = 10, gene_col_index = 1, direction = my_dir, max_gap = 2, resume_confirm = 1)
            ann_block <- detect_inversions_translo_using_blocks(dt = filtered_dt, direction = my_dir, pos_col_index = 10, gene_col_index = 1, min_consecutive = 2, block_data = blocks)
            ann_block <- fix_small_edge_blocks(ann_block)
            ann_block <- resolve_smallblock_bridge(block_annotation=ann_block, filtered_dt=filtered_dt)
            
            pdf(file.path("pdf_outputs/inversions_translocations", paste0(gca_name, ".pdf")))
            p <- plot_genes(dt = filtered_dt, block_annotation = ann_block, pos_col_index = 10, gene_col_index = 1, plot_title = gca_name)
            print(p)
            dev.off()
          }else{
            pdf(file.path("pdf_outputs/inversions_translocations", paste0(gca_name, ".pdf")))
            p <- plot_genes(dt = filtered_dt, block_annotation = ann_block, pos_col_index = 10, gene_col_index = 1, plot_title = gca_name)
            print(p)
            dev.off()
          }
        }
        
        if(nrow(ann_block) > 5){
          complex_window <<- complex_window + 1
          print(paste(complex_window," complex windows",sep=""))
          
          pdf(file.path("pdf_outputs/Complex", paste0(gca_name, ".pdf")))
          p <- plot_genes(dt = filtered_dt, block_annotation = ann_block, pos_col_index = 10, gene_col_index = 1, plot_title = gca_name)
          print(p)
          dev.off()
        }
        
        if (!is.null(ann_block) && nrow(ann_block) > 0) {
          return(list(name = gca_name, data = ann_block))
        }
      }
      
      NULL
    }))
  },
  error = function(e) {
    message(sprintf("Skipping %s: %s", basename(file), e$message))
    NULL
  })
  
  return(result)
})

percent_kept <- round((good_window_counter / nrow(mmseq_df)) * 100, 1)
percent_colinear <- round((colinear_counter / good_window_counter) * 100, 1)

message(sprintf("A total of %d out of %d windows were kept after filtering (%s%%).",
                good_window_counter, nrow(mmseq_df), percent_kept))

message(sprintf("Of these, %d out of %d are colinear windows (%s%%).",
                colinear_counter, good_window_counter, percent_colinear))

# Remove NULL results (files with no inversions)
inversion_results <- Filter(Negate(is.null), inversion_results)

# Combine all inversions into one data.table with GCA column
all_inversions <- rbindlist(
  lapply(inversion_results, function(x) {
    x$data[, GCA := x$name]
  })
)


# Convert to a single data.frame with GCA column first
inversion_results_df <- rbindlist(
  lapply(inversion_results, function(x) {
    df <- x$data
    df[, GCA := x$name]   # add GCA column
    setcolorder(df, c("GCA", setdiff(names(df), "GCA"))) # GCA as first col
    return(df)
  }),
  fill = TRUE
)


dt <- inversion_results_df   # ensure it is a data.table

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
genes_order <- fread("genes_order.txt", header = FALSE)        # read one gene per line
setnames(genes_order, "V1", "gene")                            # rename column to 'gene'
genes_order[, pos := .I]                                       # add 1-based index 'pos' along chromosome

# ---- 2) ENSURE data.table TYPE ----
breakpoint_summary <- as.data.table(breakpoint_summary)        # coerce to data.table

# ---- 3) KEEP ONLY EVENT ROLES (NO 'colinear') ----
event_roles <- c("inversion", "inversion+translocation", "translocation")  # roles to plot
breakpoint_summary <- breakpoint_summary[role %in% event_roles]            # filter to events only

# ---- 4) SPLIT BREAKPOINT INTO GENE1/GENE2 ----
breakpoint_summary[, c("g1","g2") := tstrsplit(breakpoint, "_", fixed = TRUE)]  # split "geneA_geneB"

# ---- 5) MAP GENE POSITIONS FROM genes_order ----
breakpoint_summary <- merge(breakpoint_summary, genes_order, by.x = "g1", by.y = "gene", all.x = TRUE)  # add pos1
setnames(breakpoint_summary, "pos", "pos1")                                 # rename to pos1
breakpoint_summary <- merge(breakpoint_summary, genes_order, by.x = "g2", by.y = "gene", all.x = TRUE)  # add pos2
setnames(breakpoint_summary, "pos", "pos2")                                 # rename to pos2

# ---- 6) DETERMINE ADJACENCY (NEIGHBOR GENES) ----
breakpoint_summary[, is_adjacent := abs(pos1 - pos2) == 1]                  # TRUE if consecutive in genes_order

# ---- 7) DEFINE ORDER ALONG GENOME FOR X-AXIS ----
breakpoint_summary[, order_pos := pmin(pos1, pos2)]                         # use lower index as ordering anchor
setorder(breakpoint_summary, order_pos, breakpoint, role)                   # stable order by genomic position

# ---- 8) BUILD UNIQUE ORDERED LEVELS (FIXES DUPLICATE-LEVEL ERRORS) ----
ordered_levels <- unique(breakpoint_summary$breakpoint)                     # unique breakpoints in desired order

# ---- 9) FACTOR FOR X-AXIS USING UNIQUE LEVELS ONLY ----
breakpoint_summary[, breakpoint_f := factor(breakpoint, levels = ordered_levels)]  # ordered factor for x-axis

# ---- 10) PER-LABEL COLOR FOR PLOT 2 (BLACK=ADJACENT, RED=NON-ADJACENT) ----
# Build a lookup of one row per breakpoint with its adjacency flag in the same order as 'ordered_levels'
lab_key <- unique(breakpoint_summary[, .(breakpoint, is_adjacent)])[match(ordered_levels, unique(breakpoint_summary$breakpoint)), ]
label_colors <- ifelse(lab_key$is_adjacent, "black", "red")                 # vector parallel to x-axis levels

# ---- 11) CONSISTENT ROLE COLORS (STACKED BARS) ----
role_cols <- c("inversion" = "#1b9e77", "inversion+translocation" = "#d95f02", "translocation" = "#7570b3")

# =====================================================
#   PLOT 1: ONLY ADJACENT BREAKPOINTS (x ordered by genome)
# =====================================================
bp_adj <- breakpoint_summary[is_adjacent == TRUE]                            # keep adjacent only
adj_levels <- unique(bp_adj$breakpoint)                                      # unique adjacent labels in genomic order
bp_adj[, breakpoint_f := factor(breakpoint, levels = adj_levels)]            # refactor for this subset only

p1 <- ggplot(bp_adj, aes(x = breakpoint_f, y = N, fill = role)) +            # set aesthetics
  geom_col(width = 0.9) +                                                    # stacked bars
  scale_fill_manual(values = role_cols) +                                    # fixed role colors
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)                        # vertical labels
  ) +
  labs(x = "Breakpoint (adjacent only)", y = "Count", fill = "Role")

ggsave("breakpoints_adjacent_only.pdf", p1, width = 12, height = 6, units = "in")  # save plot 1

# =====================================================
#   PLOT 2: ALL BREAKPOINTS, NON-ADJACENT LABELS IN RED
#   (x-axis ordered by genome; labels colored by adjacency)
# =====================================================
# Use the global 'breakpoint_f' with 'ordered_levels' and 'label_colors'
p2 <- ggplot(breakpoint_summary, aes(x = breakpoint_f, y = N, fill = role)) + # all breakpoints
  geom_col(width = 0.9) +                                                    # stacked bars
  scale_fill_manual(values = role_cols) +                                    # fixed role colors
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust = 0.5, colour = label_colors) # color per label (black/red)
  ) +
  labs(x = "Breakpoint (all; non-adjacent labels in red)", y = "Count", fill = "Role")

ggsave("breakpoints_all_with_nonadjacent_in_red.pdf", p2, width = 12, height = 6, units = "in")  # save plot 2


