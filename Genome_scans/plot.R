# Load libraries
library(ggplot2)
library(ggh4x)
library(dplyr)

# Load external arguments 
Input <- commandArgs(trailingOnly=TRUE)

# Check minimum 5 arguments
if(length(Input) < 5) {
  stop("Please provide at least one summary statistic")
}

data_folder <- Input[1]
plot_prefix <- Input[2]
inv_start <- as.numeric(Input[3])
inv_end <- as.numeric(Input[4])
stats_str <- Input[5]

if(nchar(stats_str) == 0) {
  stop("Please provide at least one summary statistic")
}

# Split comma-separated summary statistics into vector
stats <- unlist(strsplit(stats_str, split = ","))

if(length(stats) == 0) {
  stop("Please provide at least one summary statistic")
}

# Function to remove outliers
remove_outliers <- function(df) {
  sd_val <- sd(df[[3]], na.rm = TRUE)
  mean_val <- mean(df[[3]], na.rm = TRUE)
  threshold <- 5 * sd_val
  df_clean <- df[abs(df[[3]] - mean_val) <= threshold, ]
  return(df_clean)
}

# Function for scientific notation on x-axis
scientific_10 <- function(x) {
  ifelse(x > 1000, parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x))), x)
}

# Read and process each requested statistic dynamically
list_df <- list()
for (stat in stats) {
  if (stat == "pi") {
    pi_files <- list.files(data_folder, pattern = "^pi_.*", full.names = TRUE)
    for (file in pi_files) {
      df <- tryCatch({
        read.table(file, fill = TRUE, row.names = NULL)
      }, error = function(e) {
        warning(paste("Failed to read pi file:", file, ":", e$message))
        return(NULL)
      })
      
      if (!is.null(df) && ncol(df) >= 3) {
        colnames(df)[1:3] <- c("Start", "End", "Values")
        df$Statistics <- paste0("pi_", gsub(".*pi_", "", file))
        df <- remove_outliers(df)
        list_df[[length(list_df) + 1]] <- df
      }
    }
  } else {
    filepath <- file.path(data_folder, stat)
    if (!file.exists(filepath)) {
      warning(paste("File for statistic", stat, "not found at", filepath, "- skipping"))
      next
    }
    df <- tryCatch({
      read.table(filepath, fill = TRUE, row.names = NULL)
    }, error = function(e) {
      warning(paste("Failed to read file for statistic", stat, ":", e$message))
      return(NULL)
    })
    
    if (is.null(df) || ncol(df) < 3) {
      warning(paste("File for statistic", stat, "does not have at least 3 columns - skipping"))
      next
    }
    
    colnames(df)[1:3] <- c("Start", "End", "Values")
    df$Statistics <- stat
    df <- remove_outliers(df)
    list_df[[length(list_df) + 1]] <- df
  }
}

# Check if any data loaded
if(length(list_df) == 0) {
  stop("No valid summary statistics files loaded. Please check your input.")
}

# Combine data frames
Comb <- do.call(rbind, list_df)
Comb <- na.omit(Comb)
Comb$Start <- as.numeric(Comb$Start)
Comb$End <- as.numeric(Comb$End)

# Midpoint of windows
Comb$Positions <- (Comb$End + Comb$Start) / 2

# Cut the plot into before, within, and after the inversion
Comb$Block <- cut(Comb$Positions,
                  breaks = c(-Inf, inv_start, inv_end, Inf),
                  labels = c("Left", "Middle", "Right"))

# Calculate median per statistic and block
medians_df <- Comb %>%
  group_by(Statistics, Block) %>%
  summarise(median_val = median(Values, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    x_start = case_when(
      Block == "Left" ~ min(Comb$Positions),
      Block == "Middle" ~ inv_start,
      Block == "Right" ~ inv_end
    ),
    x_end = case_when(
      Block == "Left" ~ inv_start,
      Block == "Middle" ~ inv_end,
      Block == "Right" ~ max(Comb$Positions)
    )
  )

# Adjust PDF height dynamically based on number of statistics (1.5 inch per stat, min 7)
pdf_height <- max(7, length(unique(Comb$Statistics)) * 1.5)

# Plot
pdf(file.path(data_folder, paste0("Results_", plot_prefix, ".pdf")), width = 10, height = pdf_height)
ggplot(Comb, aes(x = Positions, y = Values)) +
  geom_point(size = 2, color = "purple") +
  theme_bw() +
  scale_x_continuous(labels = scientific_10) +
  facet_grid(Statistics ~ ., scales = "free_y") +
  geom_vline(xintercept = c(inv_start, inv_end), linetype = "dashed", size = 1, color = "blue", alpha = 0.5) +
  geom_vline(xintercept = 10069882, linetype = "dashed", size = 1, color = "red", alpha = 0.5) + # Works only for Chetone remove or adjust for other species: this indicates the position of the ivory promotor
  geom_segment(data = medians_df,
               aes(x = x_start, xend = x_end, y = median_val, yend = median_val),
               inherit.aes = FALSE,
               size = 1, color = "black") +
  theme(strip.text = element_text(size = 12, face = "bold"))
dev.off()


# Filter only pi_* data (overlayed so we can compare those)
pi_data <- Comb[grepl("^pi_", Comb$Statistics), ]

# Only proceed if there is pi data
if (nrow(pi_data) > 0) {
  # Create a new PDF
  pdf(file.path(data_folder, paste0("Results_", plot_prefix, "_pi_comparison.pdf")), width = 10, height = 6)

  p<- ggplot(pi_data, aes(x = Positions, y = Values, color = Statistics)) +
    geom_point(size = 2, alpha = 0.6) +
    theme_bw() +
    scale_x_continuous(labels = scientific_10) +
    labs(title = "Pi values across populations",
         x = "Genomic Position",
         y = "Pi",
         color = "Population") +
    geom_vline(xintercept = c(inv_start, inv_end), linetype = "dashed", size = 1, color = "blue", alpha = 0.5) +
    geom_vline(xintercept = 10069882, linetype = "dashed", size = 1, color = "red", alpha = 0.5) # Works only for Chetone remove or adjust for other species: this indicates the position of the ivory promotor
  plot(p)
  dev.off()
}