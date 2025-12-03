#!/usr/bin/env Rscript

# -------------------------------------------------------------
# load required packages
# -------------------------------------------------------------
library(data.table)
library(ggplot2)

# -------------------------------------------------------------
# read gene order
# -------------------------------------------------------------
gene_order <- fread("genes_order_busco.txt", header = FALSE)[[1]]

# -------------------------------------------------------------
# read all result files
# -------------------------------------------------------------
files <- list.files(
    pattern = "^results_.*\\.txt$",
    full.names = TRUE
)

dt_list <- lapply(files, fread)
dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# -------------------------------------------------------------
# enforce correct gene order
# -------------------------------------------------------------
dt$gene_name <- factor(dt$gene_name, levels = gene_order)
setorder(dt, gene_name)

# -------------------------------------------------------------
# plot histogram of frequencies
# -------------------------------------------------------------
p <- ggplot(dt, aes(x = gene_name, y = frequency)) +
    geom_col() +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major.x = element_blank()
    ) +
    xlab("Busco gene") +
    ylab("Prop. of species where the gene falls inside an inversion") 

ggsave("SV_frequency_histogram.pdf", p, width = 20, height = 6)


