#!/usr/bin/env Rscript

# Supplementary Figure S1: patient-level expression support for the frozen
# GSE206785 Fibroblast-to-immune CellChat interactions. This is an expression
# availability audit, not a patient-specific CellChat probability analysis.

options(stringsAsFactors = FALSE, warn = 1)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
figure_lib <- file.path(project, "runtime", "R_figure")
.libPaths(c(figure_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(ragg)
  library(svglite)
  library(ggtext)
})

source_dir <- file.path(project, "figures", "data", "Supplementary_Figure_S1_donor_support")
supp_dir <- file.path(project, "figures", "supplementary")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

support_file <- file.path(source_dir, "panel_a_donor_receiver_support.tsv")
summary_file <- file.path(source_dir, "panel_b_donor_receiver_support_summary.tsv")
if (!file.exists(support_file) || !file.exists(summary_file)) {
  stop("Missing Supplementary Figure S1 source data")
}

donor_support <- fread(support_file)
support_summary <- fread(summary_file)
stopifnot(nrow(donor_support) == 72L, nrow(support_summary) == 7L)

# Publication axes, ticks, labels and summary markers use true black for
# print contrast; only semantic lineage accents remain coloured.
ink <- "#000000"
grey <- "#6F7880"
grid_grey <- "#E9EDF0"
accent_dark <- "#98495E"
lineage_cols <- c(
  "B cells" = "#89A8CC", "CD4+ T" = "#A9BBD7", "CD8+ T" = "#6F8FB8",
  "NK/ILC" = "#91B88F", "Mast" = "#BFA6D4", "Myeloid" = "#74A69B",
  "Plasma" = "#9DA4CE"
)

theme_pub <- function(base_size = 6.2) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = ink),
      axis.ticks = element_line(linewidth = 0.35, colour = ink),
      axis.title = element_text(size = base_size, colour = ink),
      axis.text = element_text(size = base_size - 0.4, colour = ink),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(3, 3, 3, 3, unit = "pt")
    )
}

myeloid_markdown <- function(labels) {
  ifelse(
    grepl("^Myeloid", labels),
    paste0("<b><span style='color:", accent_dark, "'>", labels, "</span></b>"),
    labels
  )
}

support_order <- support_summary[order(median_weighted_support), target]
support_labels <- setNames(
  paste0(support_summary$target, " (", support_summary$eligible_donors, ")"),
  support_summary$target
)
donor_support[, target_plot := factor(target, levels = support_order)]
donor_support[, highlighted := target == "Myeloid"]
donor_medians <- donor_support[, .(
  median_support = median(pooled_weight_supported_fraction)
), by = target]
donor_medians[, target_plot := factor(target, levels = support_order)]

p_a <- ggplot(
  donor_support,
  aes(pooled_weight_supported_fraction * 100, target_plot, fill = target)
) +
  geom_jitter(
    aes(colour = highlighted), shape = 21, size = 2.2, stroke = 0.50,
    height = 0.12, width = 0, alpha = 0.80
  ) +
  geom_point(
    data = donor_medians,
    aes(median_support * 100, target_plot),
    inherit.aes = FALSE, shape = 23, size = 2.9,
    fill = ink, colour = "white", stroke = 0.45
  ) +
  scale_fill_manual(values = lineage_cols, guide = "none") +
  scale_colour_manual(values = c(`FALSE` = "white", `TRUE` = accent_dark), guide = "none") +
  scale_x_continuous(
    limits = c(60, 100.8), breaks = c(60, 80, 100),
    labels = function(x) paste0(x, "%")
  ) +
  scale_y_discrete(labels = function(x) myeloid_markdown(support_labels[x])) +
  labs(x = "Donor-supported pooled weight", y = NULL) +
  theme_pub(5.8) +
  theme(
    axis.text.y = element_markdown(size = 5.2),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(colour = grid_grey, linewidth = 0.28)
  )

eligibility <- copy(support_summary)
eligibility[, target_plot := factor(target, levels = support_order)]
eligibility[, highlighted := target == "Myeloid"]

p_b <- ggplot(eligibility, aes(eligible_donors, target_plot, fill = target)) +
  geom_col(aes(colour = highlighted), width = 0.64, linewidth = 0.55) +
  geom_text(
    aes(label = eligible_donors), hjust = -0.18,
    family = "sans", size = 1.85, colour = ink
  ) +
  scale_fill_manual(values = lineage_cols, guide = "none") +
  scale_colour_manual(values = c(`FALSE` = "white", `TRUE` = accent_dark), guide = "none") +
  scale_x_continuous(limits = c(0, 15), breaks = c(0, 5, 10, 15), expand = c(0, 0)) +
  scale_y_discrete(labels = function(x) myeloid_markdown(x)) +
  labs(x = "Eligible tumour patients", y = NULL) +
  theme_pub(5.8) +
  theme(
    axis.text.y = element_markdown(size = 5.2),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(colour = grid_grey, linewidth = 0.28)
  )

figure_s1 <- plot_grid(
  p_a, p_b, nrow = 1, rel_widths = c(1.18, 0.82),
  labels = c("a", "b"), label_fontfamily = "sans",
  label_fontface = "bold", label_size = 9,
  label_x = c(0.004, 0.004), label_y = c(0.996, 0.996),
  hjust = 0, vjust = 1, align = "h", axis = "tb"
)

width_mm <- 150
height_mm <- 68
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
stem <- file.path(supp_dir, "Supplementary_Figure_S1_donor_support")

agg_png(paste0(stem, ".png"), width = width_in, height = height_in,
        units = "in", res = 400, background = "white")
print(figure_s1)
dev.off()

agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in,
         units = "in", res = 600, compression = "lzw", background = "white")
print(figure_s1)
dev.off()

cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in,
          family = "sans", onefile = TRUE)
print(figure_s1)
dev.off()

svglite(paste0(stem, ".svg"), width = width_in, height = height_in,
        bg = "white", system_fonts = list(sans = "Arial"))
print(figure_s1)
dev.off()

cat("Supplementary Figure S1 exported successfully\n")
cat("Canvas: ", width_mm, " x ", height_mm, " mm\n", sep = "")
cat("Myeloid eligible patients: ", support_summary[target == "Myeloid", eligible_donors],
    "; median weighted support: ",
    sprintf("%.6f", support_summary[target == "Myeloid", median_weighted_support]),
    "\n", sep = "")
