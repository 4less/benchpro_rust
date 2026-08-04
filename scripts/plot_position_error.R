#!/usr/bin/env Rscript
# Plot per-sample position error metrics from benchpro strain position_error.tsv.
#
# Usage:
#   Rscript scripts/plot_position_error.R <position_error.tsv> <output_dir>

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_position_error.R <position_error.tsv> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_tsv(input_file, col_types = cols(
  ID         = col_character(),
  Species    = col_character(),
  Sample     = col_character(),
  Gene       = col_character(),
  ErrorNum   = col_double(),
  ErrorRate  = col_double(),
  GeneLength = col_double()
))

# Summarise per (ID, Species, Sample): sum genes, recompute rate.
summary_df <- df |>
  group_by(ID, Species, Sample) |>
  summarise(
    GenomeLength = sum(GeneLength, na.rm = TRUE),
    ErrorNum     = sum(ErrorNum,   na.rm = TRUE),
    ErrorRate    = ErrorNum / GenomeLength,
    .groups = "drop"
  )

# One plot per (ID, Species) pair.
pairs <- summary_df |> distinct(ID, Species)

for (i in seq_len(nrow(pairs))) {
  id_val      <- pairs$ID[i]
  species_val <- pairs$Species[i]

  sub <- summary_df |>
    filter(ID == id_val, Species == species_val) |>
    pivot_longer(
      cols      = c(ErrorNum, ErrorRate, GenomeLength),
      names_to  = "Metric",
      values_to = "Value"
    ) |>
    mutate(Metric = factor(Metric, levels = c("ErrorNum", "ErrorRate", "GenomeLength")))

  p <- ggplot(sub, aes(x = Metric, y = Value)) +
    geom_boxplot(outlier.shape = NA, width = 0.4, colour = "grey40") +
    geom_jitter(aes(colour = Metric), width = 0.15, size = 1.8, alpha = 0.7,
                show.legend = FALSE) +
    facet_wrap(~Metric, scales = "free_y", nrow = 1) +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    labs(
      title    = id_val,
      subtitle = species_val,
      x        = NULL,
      y        = "Value"
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.text       = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    )

  safe_id      <- gsub("[^A-Za-z0-9_-]", "_", id_val)
  safe_species <- gsub("[^A-Za-z0-9_-]", "_", species_val)
  out_file     <- file.path(output_dir, paste0(safe_id, "__", safe_species, ".pdf"))
  ggsave(out_file, p, width = 9, height = 5)
  message("Saved: ", out_file)
}
