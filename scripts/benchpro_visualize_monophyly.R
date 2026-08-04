#!/usr/bin/env Rscript
# Visualize sliding-window mean monophyly vs. closest-neighbor similarity.
# One line per Tool, summarised across all species.
#
# Usage:
#   Rscript visualize_monophyly.R <genome_proximity.tsv> [output.pdf]
#
# Arguments:
#   genome_proximity.tsv  *.genome_proximity.tsv produced by 'benchpro strain'
#   output.pdf            Optional output path (default: same base name, .pdf)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
usage <- "
Usage: benchpro_visualize_monophyly <genome_proximity.tsv> [output.pdf]

  genome_proximity.tsv  *.genome_proximity.tsv produced by 'benchpro strain'
  output.pdf            Output path (default: same name with .pdf extension)

Description:
  Plots sliding-window mean monophyly score against the similarity of each
  genome to its closest neighbour in the gold-standard tree.  One line is
  drawn per Tool (derived from the ID column by stripping the species suffix,
  or from an explicit 'Tool' column if present in the file).  Genomes from
  all species are pooled per tool and the window mean is recomputed across
  all of them before plotting.

Options:
  -h, --help  Show this help message and exit
"

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || args[1] %in% c("-h", "--help")) {
  cat(usage)
  quit(status = if (length(args) == 0) 1L else 0L)
}

input_file  <- args[1]
output_file <- if (length(args) >= 2) args[2] else sub("\\.tsv$", ".pdf", input_file)

# ── Load ──────────────────────────────────────────────────────────────────────
df <- read.table(input_file, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, check.names = FALSE)

# Use 'Tool' column if present; otherwise strip the species suffix from ID
# (everything from '_s__' onward) to get the tool name.
if (!"Tool" %in% colnames(df)) {
  df$Tool <- sub("_s__.*$", "", df$ID)
}

# ── Re-compute sliding window across all species per Tool ─────────────────────
# The Rust output contains per-job (Tool × Species) window means.  Here we
# pool all genomes belonging to the same Tool, sort by similarity, and apply
# a centred window of size 30 so the curve summarises across species.
window_mean <- function(x, k = 30) {
  n    <- length(x)
  half <- k %/% 2
  vapply(seq_along(x), function(i) {
    mean(x[max(1L, i - half):min(n, i + half)])
  }, numeric(1))
}

df <- df %>%
  group_by(Tool) %>%
  arrange(ClosestNeighborSimilarity, .by_group = TRUE) %>%
  mutate(WindowMean = window_mean(MonophylyScore, k = 30)) %>%
  ungroup()

# ── X-axis transformation ─────────────────────────────────────────────────────
# Map similarity s to -log10(1 - s) for log-spaced visual separation;
# tick labels show the original negative-similarity value (−s).
df$x_trans <- -log10(1 - df$ClosestNeighborSimilarity)

break_sims   <- c(0.98, 0.99, 0.995, 0.998, 0.9999, 0.99999, 0.999999)
break_x      <- -log10(1 - break_sims)
break_labels <- paste0("\u2212",
                       c("0.98", "0.99", "0.995",
                         "0.998", "0.9999", "0.99999", "0.999999"))

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(df, aes(x = x_trans, y = WindowMean, colour = Tool, group = Tool)) +
  geom_line(linewidth = 0.6, alpha = 0.35, na.rm = TRUE) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              linewidth = 1.1, na.rm = TRUE) +
  scale_x_continuous(
    name   = "Similarity to closest neighboring genome in tree",
    breaks = break_x,
    labels = break_labels
  ) +
  scale_y_continuous(
    name   = "Monophyly per genome",
    limits = c(NA, 1)
  ) +
  labs(
    title  = "Mean monophyly per genome",
    colour = "Tool"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.tag         = element_text(face = "bold", size = 14),
    plot.title       = element_text(face = "bold"),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 30, hjust = 1)
  )

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave(output_file, p, width = 6, height = 4.5)
message("Written: ", output_file)
