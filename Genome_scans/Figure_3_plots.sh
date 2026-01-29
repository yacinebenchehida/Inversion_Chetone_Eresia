############################################################
## Genome scan plotting script (annotated, with medians)   ##
## Nested facets: 4 rows (stats) x 3 columns (species)     ##
############################################################

## Set working directory (adjust if needed)
setwd("/Users/ybc502/Desktop/Inversion_project/Sliding_windows/genome_scans")  # set base directory

## Load required libraries
library(ggplot2)  # plotting
library(dplyr)    # data manipulation
library(scales)   # scientific_format for axis labels
library(ggh4x)    # facet_grid2 + facetted_pos_scales

############################################################
## User-defined configuration                              ##
############################################################

folders <- c("chetone", "numata", "eresia")  # species folders to iterate over

stats_files <- c("Fst", "Dxy", "r2", "pi_hom_anc", "pi_hom_inv")  # expected files per folder

inv_coords <- list(  # inversion coordinates per species
  chetone = c(start = 50673046, end = 51693436),  # chetone inversion start and end
  numata  = c(start = 1346000,  end = 1776500),   # numata inversion start and end
  eresia  = c(start = 2808389,  end = 3352577)    # eresia inversion start and end
)

stat_colors <- c(  # point color for non-pi statistics (one per species)
  eresia  = "#E69F00",  # non-pi point color for eresia
  numata  = "#56B4E9",  # non-pi point color for numata
  chetone = "#009E73"   # non-pi point color for chetone
)

pi_point_colors <- c(  # pi point colors (same across species)
  pi_hom_anc = "#F8766D",  # pi population 1 point color
  pi_hom_inv = "#00BFC4"   # pi population 2 point color
)

y_scales <- list(  # y-axis ranges per Stat row (must match Stat factor levels below)
  scale_y_continuous(limits = c(0, 0.45)),   # Fst y limits
  scale_y_continuous(limits = c(0, 0.035)),  # Dxy y limits
  scale_y_continuous(limits = c(0, 0.022)),  # pi y limits
  scale_y_continuous(limits = c(0, 0.08))   # r2 y limits
)

stat_labeller <- as_labeller(  # facet strip label formatting for Stat
  c(  # mapping from Stat values to plotmath strings
    Fst = "italic(F)[ST]",  # italic F with ST subscript
    Dxy = "italic(D)[xy]",  # italic D with xy subscript
    r2  = "italic(r)^2",    # italic r squared
    pi  = "italic(pi)"      # italic pi (Greek letter via plotmath)
  ),  # end mapping
  label_parsed  # parse plotmath expressions
)  # end labeller

species_labeller <- as_labeller(c(
  chetone = "italic(Chetone~histrio)",
  numata  = "italic(Heliconius~numata)",
  eresia  = "italic(Eresia~pelonia)"
), label_parsed)


############################################################
## Helper functions                                        ##
############################################################

scientific_10 <- function(x) {  # function for scientific notation on x-axis
  ifelse(  # conditional formatting
    x > 1000,  # apply scientific formatting only for large coordinates
    parse(text = gsub("e\\+*", " %*% 10^", scales::scientific_format()(x))),  # convert 1e+07 to 1 * 10^7 expression
    x  # otherwise keep raw values
  )  # end ifelse
}  # end function

############################################################
## Build inversion lookup table (vectorised, robust)       ##
############################################################

inv_df <- data.frame(  # build a lookup table for joins
  Species   = names(inv_coords),  # species names from inv_coords
  inv_start = as.numeric(sapply(inv_coords, `[[`, "start")),  # extract start for each species
  inv_end   = as.numeric(sapply(inv_coords, `[[`, "end")),    # extract end for each species
  stringsAsFactors = FALSE  # keep strings as character
)  # end inv_df

############################################################
## Load data for all species                               ##
############################################################

all_dfs <- list()  # initialise list to store per-species combined data frames

for (sp in folders) {  # loop over species folders

  dfs <- list()  # initialise list to store per-file data frames for this species

  for (f in stats_files) {  # loop over expected statistic files

    filepath <- file.path(sp, f)  # build full path to file
    if (!file.exists(filepath)) next  # skip if file does not exist

    df <- read.table(filepath)  # read table from file
    colnames(df)[1:3] <- c("Start", "End", "Value")  # enforce standard column names

    df$Start <- as.numeric(df$Start)  # ensure numeric Start
    df$End <- as.numeric(df$End)      # ensure numeric End
    df$Value <- as.numeric(df$Value)  # ensure numeric Value

    df$Position <- (df$Start + df$End) / 2  # compute midpoint position
    df$Species <- sp  # store species name for column facets

    if (grepl("^pi_", f)) {  # check if file is a pi file
      df$Stat <- "pi"  # collapse both pi files into a single Stat facet
      df$Group <- f    # store which pi group this is (pi_hom_anc or pi_hom_inv)
    } else {  # otherwise this is a non-pi statistic
      df$Stat <- f     # Stat is the filename (Fst, Dxy, r2)
      df$Group <- NA   # Group is not used for non-pi
    }  # end if

    dfs[[length(dfs) + 1]] <- df  # append df to species list
  }  # end inner loop

  Comb_sp <- bind_rows(dfs)  # combine all stats for this species
  Comb_sp <- Comb_sp[!is.na(Comb_sp$Position) & !is.na(Comb_sp$Value), ]  # drop rows missing Position or Value

  all_dfs[[length(all_dfs) + 1]] <- Comb_sp  # append species combined df to global list
}  # end outer loop

Comb <- bind_rows(all_dfs)  # combine all species into one data frame

############################################################
## Add inversion coords, facet ordering, and blocks         ##
############################################################

Comb <- left_join(Comb, inv_df, by = "Species")  # join inv_start and inv_end onto each row

Comb$Stat <- factor(Comb$Stat, levels = c("Fst", "Dxy", "pi", "r2"))  # enforce row facet order

Comb$Block <- dplyr::case_when(  # define Left/Middle/Right blocks per row
  Comb$Position < Comb$inv_start ~ "Left",  # left of inversion start
  Comb$Position <= Comb$inv_end ~ "Middle", # within inversion region (inclusive)
  TRUE ~ "Right"  # right of inversion end
)  # end case_when
Comb$Block <- factor(Comb$Block, levels = c("Left", "Middle", "Right"))  # set block order

############################################################
## Build a single color key so one color scale works        ##
############################################################

Comb$ColorKey <- dplyr::if_else(  # create a single key for color mapping
  Comb$Stat == "pi",  # if pi panel
  as.character(Comb$Group),  # use pi group as color key
  as.character(Comb$Species) # otherwise use species as color key
)  # end if_else

color_values <- c(  # define colors for both species and pi groups
  chetone = stat_colors[["chetone"]],  # species color for chetone non-pi
  numata  = stat_colors[["numata"]],   # species color for numata non-pi
  eresia  = stat_colors[["eresia"]],   # species color for eresia non-pi
  pi_hom_anc = pi_point_colors[["pi_hom_anc"]],  # pi group color for pi_hom_anc
  pi_hom_inv = pi_point_colors[["pi_hom_inv"]]   # pi group color for pi_hom_inv
)  # end color_values

############################################################
## Species-specific x ranges for median segments            ##
############################################################

species_ranges <- Comb %>%  # start from combined data
  group_by(Species) %>%  # group by species
  summarise(  # compute per-species plot ranges and inversion coords
    min_pos = min(Position, na.rm = TRUE),  # minimum midpoint position for species
    max_pos = max(Position, na.rm = TRUE),  # maximum midpoint position for species
    inv_start = first(inv_start),  # inversion start for species
    inv_end   = first(inv_end),    # inversion end for species
    .groups = "drop"  # drop grouping
  )  # end summarise

############################################################
## Compute medians per block (species-aware)               ##
############################################################

med_nonpi <- Comb %>%  # start from combined data
  filter(Stat != "pi") %>%  # keep non-pi panels
  group_by(Species, Stat, Block) %>%  # group by species, stat, and block
  summarise(median_val = median(Value, na.rm = TRUE), .groups = "drop") %>%  # compute median
  left_join(species_ranges, by = "Species") %>%  # add min_pos, max_pos, inv_start, inv_end for segments
  mutate(  # compute segment start and end per block
    x_start = case_when(  # define segment x start
      Block == "Left"   ~ min_pos,     # left block starts at min position
      Block == "Middle" ~ inv_start,   # middle block starts at inversion start
      Block == "Right"  ~ inv_end      # right block starts at inversion end
    ),  # end case_when
    x_end = case_when(  # define segment x end
      Block == "Left"   ~ inv_start,   # left block ends at inversion start
      Block == "Middle" ~ inv_end,     # middle block ends at inversion end
      Block == "Right"  ~ max_pos      # right block ends at max position
    )  # end case_when
  )  # end mutate

med_pi <- Comb %>%  # start from combined data
  filter(Stat == "pi") %>%  # keep only pi panel data
  group_by(Species, Stat, Group, Block) %>%  # group by species, pi group, and block
  summarise(median_val = median(Value, na.rm = TRUE), .groups = "drop") %>%  # compute median
  left_join(species_ranges, by = "Species") %>%  # add min_pos, max_pos, inv_start, inv_end for segments
  mutate(  # compute segment start and end per block
    x_start = case_when(  # define segment x start
      Block == "Left"   ~ min_pos,     # left block starts at min position
      Block == "Middle" ~ inv_start,   # middle block starts at inversion start
      Block == "Right"  ~ inv_end      # right block starts at inversion end
    ),  # end case_when
    x_end = case_when(  # define segment x end
      Block == "Left"   ~ inv_start,   # left block ends at inversion start
      Block == "Middle" ~ inv_end,     # middle block ends at inversion end
      Block == "Right"  ~ max_pos      # right block ends at max position
    )  # end case_when
  )  # end mutate

############################################################
## Build vline data for species-specific inversion lines   ##
############################################################

vlines_df <- inv_df %>%  # start from inversion lookup
  tidyr::pivot_longer(  # reshape start/end into long format
    cols = c(inv_start, inv_end),  # columns to pivot
    names_to = "which",  # name of the new identifier column
    values_to = "xint"   # name of the x-intercept column
  ) %>%  # end pivot
  select(Species, xint)  # keep only needed columns

############################################################
## Plot: 4 rows (Stat) x 3 cols (Species)                  ##
############################################################
pdf("Results_all_species_nested.pdf", width = 15, height = 8)  # open one PDF for all species

p <- ggplot() +  # start plot

  geom_point(  # plot points for all statistics
    data = Comb,  # use combined data
    aes(x = Position, y = Value, color = ColorKey),  # map x, y, and unified color key
    size = 2,  # point size
    alpha = dplyr::if_else(Comb$Stat == "pi", 0.4, 0.7)  # alpha differs for pi vs non-pi
  ) +  # end geom_point

  geom_segment(  # draw median segments for non-pi stats
    data = med_nonpi,  # use non-pi median table
    aes(x = x_start, xend = x_end, y = median_val, yend = median_val),  # segment aesthetics
    inherit.aes = FALSE,  # do not inherit aesthetics
    color = "black",  # median segment color for non-pi
    size = 1  # segment thickness
  ) +  # end geom_segment

  geom_segment(  # draw median segments for pi (colored by group)
    data = med_pi,  # use pi median table
    aes(x = x_start, xend = x_end, y = median_val, yend = median_val, color = Group),  # segment aesthetics
    inherit.aes = FALSE,  # do not inherit aesthetics
    size = 1  # segment thickness
  ) +  # end geom_segment

  geom_vline(  # draw inversion boundaries per species column
    data = vlines_df,  # use per-species vline data
    aes(xintercept = xint),  # map x-intercept
    inherit.aes = FALSE,  # do not inherit aesthetics
    linetype = "dashed",  # dashed line
    color = "blue",  # line color
    alpha = 0.5  # line transparency
  ) +  # end geom_vline

  ggh4x::facet_grid2(  # nested grid facet (rows and columns)
    Stat ~ Species,  # 4 rows by Stat and 3 columns by Species
    scales = "free",  # free y across facets (limits forced by facetted_pos_scales below)
    labeller = labeller(Stat = stat_labeller, Species = species_labeller)  # formatted Stat and species strip labels
  ) +  # end facet_grid2

  ggh4x::facetted_pos_scales(  # apply facet-specific y scales by row order
    y = y_scales  # scales must match Stat factor order: Fst, Dxy, pi, r2
  ) +  # end facetted_pos_scales

  scale_x_continuous(labels = scientific_10) +  # scientific labels on x-axis

  scale_color_manual(
  values = color_values,
  breaks = c("pi_hom_anc", "pi_hom_inv"),
  labels = c(
    "pi_hom_anc" = expression(pi~"(ancestral)"),
    "pi_hom_inv" = expression(pi~"(inversion)")
  ),
  na.translate = FALSE
) +  # end scale_color_manual

  theme_bw() +  # use bw theme
  theme(strip.text = element_text(size = 12, face = "bold")) +  # facet strip styling
  labs(x = "Genomic position", y = NULL, color = NULL)  # axis labels and legend title

print(p)  # render plot into the PDF device
dev.off()  # close PDF device
