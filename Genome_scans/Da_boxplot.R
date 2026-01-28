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
