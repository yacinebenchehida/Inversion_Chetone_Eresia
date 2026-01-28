library(ggplot2)
library(ggpattern)

setwd("/Users/ybc502/Desktop/Inversion_project/Sliding_windows/")

# Load datasets
Chetone <- read.table("Da_chetone.txt")
Chetone$species <- "Chetone"

Eresia <- read.table("Da_eresia.txt")
Eresia$species <- "Eresia"

Numata <- read.table("Da_numata.txt")
Numata$species <- "Numata"

# Set gene coordinates for each species (example)
gene_coords <- list(
  Chetone = 51039870,   
  Eresia  = 3346772, 
  Numata  = 1595580    
)

# Apply filters for the full inversion 
Eresia <- Eresia[Eresia$V1 >= 2808389 & Eresia$V1 <= 3352577,]
Numata <- Numata[Numata$V1 >= 1346000 & Numata$V1 <= 1776500,]

# Create full inversion datasets
Chetone_full <- Chetone
Eresia_full  <- Eresia
Numata_full  <- Numata

# Create truncated datasets up to gene
Chetone_gene <- Chetone[Chetone$V1 <= gene_coords$Chetone, ]
Eresia_gene  <- Eresia[Eresia$V1 <= gene_coords$Eresia, ]
Numata_gene  <- Numata[Numata$V1 <= gene_coords$Numata, ]

# Add region_type column
Chetone_full$region <- "full_inversion"
Eresia_full$region  <- "full_inversion"
Numata_full$region  <- "full_inversion"

Chetone_gene$region <- "common_genes"
Eresia_gene$region  <- "common_genes"
Numata_gene$region  <- "common_genes"

# Combine all datasets
All_sp <- rbind(Chetone_full, Eresia_full, Numata_full,
                Chetone_gene, Eresia_gene, Numata_gene)

# Rename columns
colnames(All_sp) <- c("start", "end", "value", "species", "region")
All_sp$region <- factor(All_sp$region, levels = c("full_inversion", "common_genes"))  # enforce region order



ggplot(All_sp, aes(
  x = species,
  y = value,
  group = interaction(species, region, sep = "-"),
  color = interaction(species, region, sep = "-"),
  fill  = interaction(species, region, sep = "-"))) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    alpha = 0.6,
    width = 0.6
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width  = 0.8
    ),
    size = 1.5,
    alpha = 0.6
  ) +
  scale_color_manual(
    name = "Species_Region",
    values = c(
      "Chetone-full_inversion" = "#009E73",
      "Chetone-common_genes"   = "#009E73",
      "Numata-full_inversion"  = "#56B4E9",
      "Numata-common_genes"    = "#56B4E9",
      "Eresia-full_inversion"  = "#E69F00",
      "Eresia-common_genes"    = "#E69F00"
    )
  ) +
  scale_fill_manual(
    name = "Species_Region",
    values = c(
      "Chetone-full_inversion" = "#009E73",
      "Numata-full_inversion"  = "#56B4E9",
      "Eresia-full_inversion"  = "#E69F00",
      "Chetone-common_genes"   = NA,
      "Numata-common_genes"    = NA,
      "Eresia-common_genes"    = NA
    ),
    na.translate = FALSE
  ) +
  theme_classic(base_size = 14) +
  labs(x = "Species", y = "Da") +
  scale_x_discrete(labels = c(
    Chetone = expression(italic("Chetone histrio")),
    Numata  = expression(italic("Heliconius numata")),
    Eresia  = expression(italic("Eresia pelonia"))
  ))
