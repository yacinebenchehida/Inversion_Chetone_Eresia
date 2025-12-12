library(ape)
library(tidyverse)
library(ggtree)

tree <- read.tree("../busco_phylogeny/busco_nj_phylo.nwk")

meta <- read_tsv(
  "../busco_phylogeny/meta_info.txt",
  col_names = c("accession", "species", "assembly_type", "ploidy", "n_contigs", "N50", "assembly_size", "other")
)

rename_vec <- meta$species
names(rename_vec) <- meta$accession

tree$tip.label <- rename_vec[tree$tip.label]

missing <- setdiff(tree$tip.label, meta$accession)
missing

pdf("../busco_phylogeny/complete_tree.pdf",20,200)
ggtree(tree) +
  geom_tiplab()
dev.off()
