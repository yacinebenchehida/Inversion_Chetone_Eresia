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

