#!/usr/bin/env Rscript

# -------------------------------------------------------------
# load command line arguments: each argument is a BUSCO gene name
# -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)             # batch of gene names
if (length(args) == 0) stop("No genes provided")     # require at least one gene

# -------------------------------------------------------------
# list all species BUSCO files to process
# -------------------------------------------------------------
busco_files <- list.files(
    "/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Results/busco",
    pattern = "busco\\.tsv$",
    full.names = TRUE,
    recursive = TRUE
)
# -------------------------------------------------------------
# load required R packages
# -------------------------------------------------------------
library(ggplot2)
library(data.table)
library(dplyr)

# -------------------------------------------------------------
# Load functions
# -------------------------------------------------------------
source("/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Scripts/R_functions/filter_contigs_by_proportion.R")             
source("/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Scripts/R_functions/assess_direction.R")
source("/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Scripts/R_functions/identify_blocks.R")
source("/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Scripts/R_functions/detect_inversions_translo_using_blocks.R")
source("/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Scripts/R_functions/fix_small_edge_blocks.R")
source("/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Scripts/R_functions/resolve_smallblock_bridge.R")
source("/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Scripts/R_functions/detect_and_fix_backward_boundary.R")
source("/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Scripts/R_functions/remove_too_fuzzy_breakpoints.R")
source("/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Scripts/R_functions/check_gene_in_block.R")
source("/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Scripts/R_functions/pct_genes_in_scaffold.R")
source("/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Scripts/R_functions/is_big_enough.R")
source("/mnt/scratch/projects/biol-specgen-2018/yacine/Inversions/Inversions_identification_pipeline/Scripts/R_functions/dispersion_score.R")

# -------------------------------------------------------------
# master loop over genes in this batch
# -------------------------------------------------------------
for (gene_name in args) {                            # process each gene in batch

    # initialize counters for this gene
    total_present <- 0                               # how many species have this gene
    total_with_sv <- 0                               # how many species have gene in SV

    # ---------------------------------------------------------
    # process each species one by one
    # ---------------------------------------------------------
    results_list <- lapply(seq_along(busco_files), function(i) {
        
        busco_path <- busco_files[i]
        print(busco_path)
        message(sprintf("[%d/%d] Processing %s...", i, length(busco_files), busco_path))

        dt <- fread(busco_path, header = FALSE)      # BUSCO data for one species

        if (!(gene_name %in% dt$V1)) return(NULL)

        filtered_dt <- tryCatch(
            filter_contigs_by_proportion(
                dt,
                contig_col_index = 2,
                threshold = 0.7,
                outlier = TRUE,
                min_bit_score = 80,
                gene_order = "genes_order_busco.txt"
            ),
            error = function(e) NULL
        )
        if (is.null(filtered_dt) || nrow(filtered_dt) == 0) return(NULL)

        direction <- assess_direction(filtered_dt)   # forward or reverse orientation

        # Identify blocks
        blocks <- identify_blocks(
            dt = filtered_dt,
            pos_col_index = 9,
            gene_col_index = 1,
            direction = direction,
            max_gap = 2,
            resume_confirm = 1
        )

        # Initial block annotation
        ann_block <- detect_inversions_translo_using_blocks(
            dt = filtered_dt,
            direction = direction,
            pos_col_index = 9,
            gene_col_index = 1,
            min_consecutive = 2,
            block_data = blocks
        )

        # Further block refinements
        ann_block <- fix_small_edge_blocks(ann_block)
        ann_block <- resolve_smallblock_bridge(ann_block, filtered_dt)
        ann_block <- detect_and_fix_backward_boundary(
            dt = filtered_dt,
            pos_col_index = 9,
            gene_col_index = 1,
            block_data = ann_block
        )

        ann_block <- remove_too_fuzzy_breakpoints(
            block_annotation = ann_block,
            gene_order = "genes_order_busco.txt",
            max_genes_dist = 3
        )

         # -----------------------------------------------------
        # check whether this gene is inside an SV
        # -----------------------------------------------------
        in_sv <- check_gene_in_block(
            block_data = ann_block,
            gene_order = "genes_order_busco.txt",
            gene_name = gene_name
        )

        return(list(present = TRUE, in_sv = isTRUE(in_sv)))
    })

    present_vec <- sapply(results_list, function(x) if (is.null(x)) FALSE else x$present)
    sv_vec      <- sapply(results_list, function(x) if (is.null(x)) FALSE else x$in_sv)

    total_present <- sum(present_vec)
    total_with_sv <- sum(sv_vec)

    # -------------------------------------------------------------
    # compute frequency for this gene
    # -------------------------------------------------------------
    if (total_present == 0) {
        frequency <- 0
    } else {
        frequency <- total_with_sv / total_present
    }

    # -------------------------------------------------------------
    # output is one file per gene
    # -------------------------------------------------------------
    outfile <- paste0("results_", gene_name, ".txt")      # filename for this gene

    fwrite(
        data.table(
            gene_name = gene_name,
            total_with_SV = total_with_sv,
            frequency = frequency
        ),
        file = outfile,
        sep = "\t",
        col.names = TRUE
    )
}
