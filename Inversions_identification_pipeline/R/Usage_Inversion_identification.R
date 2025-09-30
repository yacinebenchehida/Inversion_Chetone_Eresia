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

inversion_results <- lapply(gca_files, function(file) { # use future_lapply for parallelisation but speed gain only huge data otherwise even slower due to overhead.
  gca_name <- tools::file_path_sans_ext(basename(file))
  message(sprintf("Processing %s...", gca_name))
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
length(inversion_results)
