#!/usr/bin/env Rscript

# Figure 1: unbiased ecosystem-level discovery of a recurrent, pathway-distinct
# Fibroblast -> Myeloid communication branch in gastric cancer.
#
# Drawing, assembly, preview and export are R-only. PLAT, F13 and fine myeloid
# states are intentionally absent. Figure 1 is discovery-cohort only; donor
# robustness is exported separately as Supplementary Figure S1 and independent
# cohort validation is deferred until that cohort is formally annotated.

options(stringsAsFactors = FALSE, warn = 1)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
figure_lib <- file.path(project, "runtime", "R_figure")
cellchat_lib <- "I:/ECM-new/scRNA/20.rebuild_GSE163558/R_libs_4.4"
.libPaths(c(figure_lib, cellchat_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(ragg)
  library(svglite)
  library(ggrepel)
  library(ggtext)
  library(scales)
})

source_dir <- file.path(project, "figures", "data", "Figure1_GSE206785_unbiased")
panel_dir <- file.path(project, "figures", "panels", "Figure1_GSE206785_restructured")
main_dir <- file.path(project, "figures", "main")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  "panel_d_global_umap.tsv.gz",
  "panel_e_cellchat_group_coverage.tsv",
  "panel_e_cellchat_broad_network.tsv",
  "panel_e_cellchat_fibroblast_receiver_summary.tsv",
  "panel_d_cellchat_fibroblast_pathway_receiver.tsv",
  "panel_c_cellchat_signaling_roles.tsv"
)
missing_files <- required_files[!file.exists(file.path(source_dir, required_files))]
if (length(missing_files)) {
  stop("Missing restructured Figure 1 Source Data: ", paste(missing_files, collapse = ", "))
}

read_source <- function(name) fread(file.path(source_dir, name))

umap <- read_source("panel_d_global_umap.tsv.gz")
coverage <- read_source("panel_e_cellchat_group_coverage.tsv")
network <- read_source("panel_e_cellchat_broad_network.tsv")
fibroblast_receiver <- read_source("panel_e_cellchat_fibroblast_receiver_summary.tsv")
pathways <- read_source("panel_d_cellchat_fibroblast_pathway_receiver.tsv")
roles <- read_source("panel_c_cellchat_signaling_roles.tsv")

# ------------------------------ assertions ---------------------------------

stopifnot(nrow(umap) == 8337L)
stopifnot(nrow(network) == 64L)
stopifnot(nrow(roles) == 8L)
stopifnot(nrow(pathways) == 140L)
stopifnot(uniqueN(pathways$pathway_name) == 20L)
stopifnot(uniqueN(pathways$target) == 7L)

fib_role <- roles[cell_type == "Fibroblast"]
if (nrow(fib_role) != 1L || fib_role$outgoing_weight_rank != 1L ||
    abs(fib_role$outgoing_weight - 24.9233803671) > 1e-6) {
  stop("Fibroblast global outgoing role changed")
}
immune_targets <- c("B cells", "CD4+ T", "CD8+ T", "NK/ILC", "Mast", "Myeloid", "Plasma")
fibroblast_receiver[, `:=`(immune_count_rank = NA_integer_, immune_weight_rank = NA_integer_)]
fibroblast_receiver[target %chin% immune_targets, `:=`(
  immune_count_rank = frank(-interaction_count, ties.method = "min"),
  immune_weight_rank = frank(-interaction_weight, ties.method = "min")
)]
myeloid_discovery <- fibroblast_receiver[target == "Myeloid"]
if (nrow(myeloid_discovery) != 1L || myeloid_discovery$immune_count_rank != 2L ||
    myeloid_discovery$immune_weight_rank != 4L) {
  stop("GSE206785 Fibroblast -> Myeloid ranks changed")
}
# ------------------------------- palette -----------------------------------

ink <- "#000000"
grey <- "#6F7880"
light_grey <- "#E6EAED"
grid_grey <- "#E9EDF0"
accent <- "#B95F72"
accent_dark <- "#98495E"
accent_light <- "#F3E2E6"
teal_dark <- "#2F7566"

lineage_cols <- c(
  "B cells" = "#89A8CC",
  "B_cell" = "#89A8CC",
  "CD4+ T" = "#A9BBD7",
  "CD8+ T" = "#6F8FB8",
  "T_cell" = "#8EA7C8",
  "Endothelial" = "#82BEB7",
  "Epithelial" = "#D8A98F",
  "Fibroblast" = "#C98393",
  "Mast" = "#BFA6D4",
  "Mural" = "#CDBB83",
  "Myeloid" = "#74A69B",
  "NK/ILC" = "#91B88F",
  "NK" = "#91B88F",
  "Plasma" = "#9DA4CE"
)

theme_pub <- function(base_size = 6.4) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = ink),
      axis.ticks = element_line(linewidth = 0.35, colour = ink),
      axis.title = element_text(size = base_size, colour = ink),
      axis.text = element_text(size = base_size - 0.4, colour = ink),
      legend.title = element_text(size = base_size - 0.2, colour = ink),
      legend.text = element_text(size = base_size - 0.6, colour = ink),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(2, 2, 2, 2, unit = "pt")
    )
}
theme_set(theme_pub())

panel_label <- function(plot, label, x = 0.004, y = 0.996) {
  ggdraw(plot) +
    draw_label(
      label, x = x, y = y, hjust = 0, vjust = 1,
      fontfamily = "sans", fontface = "bold", size = 9, colour = ink
    )
}

myeloid_markdown <- function(labels) {
  ifelse(
    labels == "Myeloid",
    paste0("<b><span style='color:", accent_dark, "'>Myeloid</span></b>"),
    labels
  )
}

# ----------------------------- panel a: UMAP -------------------------------

umap[, cell_type := factor(cell_type, levels = names(lineage_cols)[
  !duplicated(names(lineage_cols)) & names(lineage_cols) %chin% unique(umap$cell_type)
])]
centroids <- umap[, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by = cell_type]
centroids[cell_type == "CD8+ T", `:=`(UMAP1 = UMAP1 - 0.7, UMAP2 = UMAP2 - 1.0)]
centroids[cell_type == "NK/ILC", `:=`(UMAP1 = UMAP1 + 0.7, UMAP2 = UMAP2 + 0.75)]
centroids[cell_type == "CD4+ T", UMAP2 := UMAP2 - 0.45]
centroids[cell_type == "B cells", `:=`(UMAP1 = UMAP1 - 0.55, UMAP2 = UMAP2 + 0.25)]

p_a <- ggplot(umap, aes(UMAP1, UMAP2, colour = cell_type)) +
  geom_point(size = 0.24, alpha = 0.72, stroke = 0) +
  geom_label(
    data = centroids, aes(label = cell_type),
    family = "sans", fontface = "bold", size = 1.8,
    linewidth = 0, label.padding = unit(0.32, "mm"),
    fill = alpha("white", 0.74), colour = ink, show.legend = FALSE
  ) +
  scale_colour_manual(values = lineage_cols, drop = FALSE) +
  coord_equal(expand = FALSE) +
  theme_void(base_family = "sans") +
  theme(legend.position = "none", plot.margin = margin(2, 2, 1, 1, unit = "pt"))

# ---------------------- panel b: classic CellChat circles ------------------

groups <- c(
  "Fibroblast", "B cells", "CD4+ T", "CD8+ T",
  "NK/ILC", "Mast", "Myeloid", "Plasma"
)
theta <- seq(pi / 2, pi / 2 - 2 * pi, length.out = length(groups) + 1L)[-length(groups) - 1L]
nodes <- data.table(
  cell_type = groups,
  x = cos(theta), y = sin(theta), node_order = seq_along(groups)
)
nodes <- merge(nodes, coverage, by = "cell_type", all.x = TRUE, sort = FALSE)
nodes[, `:=`(
  label_x = 1.22 * x,
  label_y = 1.22 * y,
  hjust = fifelse(x > 0.15, 0, fifelse(x < -0.15, 1, 0.5))
)]

edges <- copy(network[source != target])
edges <- merge(
  edges,
  nodes[, .(source = cell_type, x_source = x, y_source = y,
            source_node_order = node_order)],
  by = "source", all.x = TRUE
)
edges <- merge(
  edges,
  nodes[, .(target = cell_type, x_target = x, y_target = y,
            target_node_order = node_order)],
  by = "target", all.x = TRUE
)

cellchat_circle <- function(value_col, short_label) {
  z <- copy(edges[get(value_col) > 0])
  z[, scaled_metric := rescale(get(value_col), to = c(0.10, 1.00))]
  z[, edge_alpha := rescale(get(value_col), to = c(0.18, 0.68))]
  forward <- z[source_node_order < target_node_order]
  reverse <- z[source_node_order > target_node_order]

  ggplot() +
    geom_curve(
      data = forward,
      aes(x = x_source, y = y_source, xend = x_target, yend = y_target,
          linewidth = scaled_metric, colour = source, alpha = edge_alpha),
      curvature = 0.18,
      arrow = arrow(length = unit(0.58, "mm"), type = "closed")
    ) +
    geom_curve(
      data = reverse,
      aes(x = x_source, y = y_source, xend = x_target, yend = y_target,
          linewidth = scaled_metric, colour = source, alpha = edge_alpha),
      curvature = -0.18,
      arrow = arrow(length = unit(0.58, "mm"), type = "closed")
    ) +
    geom_point(
      data = nodes, aes(x, y, fill = cell_type),
      shape = 21, size = 5.5, colour = "white", stroke = 0.62
    ) +
    geom_text(
      data = nodes,
      aes(label_x, label_y, label = cell_type, hjust = hjust),
      family = "sans", size = 1.85, colour = ink
    ) +
    annotate(
      "text", x = 0, y = 1.48, label = short_label,
      family = "sans", fontface = "bold", size = 2.0, colour = ink
    ) +
    scale_fill_manual(values = lineage_cols, guide = "none") +
    scale_colour_manual(values = lineage_cols, guide = "none") +
    scale_alpha_identity() +
    scale_linewidth_identity() +
    coord_equal(xlim = c(-1.48, 1.48), ylim = c(-1.40, 1.54), clip = "off") +
    theme_void(base_family = "sans") +
    theme(plot.margin = margin(1, 8, 1, 8, unit = "pt"))
}

p_b_count <- cellchat_circle("interaction_count", "Count")
p_b_weight <- cellchat_circle("interaction_weight", "Weight")
p_b <- plot_grid(p_b_count, p_b_weight, nrow = 1, rel_widths = c(1, 1), align = "h")

# ---------------------- panel c: global signaling roles --------------------

p_c <- ggplot(
  roles,
  aes(incoming_weight, outgoing_weight, fill = cell_type, size = outgoing_count)
) +
  geom_point(shape = 21, colour = "white", stroke = 0.65, alpha = 0.95) +
  geom_point(
    data = roles[cell_type == "Fibroblast"],
    shape = 21, fill = NA, colour = accent_dark, stroke = 0.9
  ) +
  geom_text_repel(
    aes(label = cell_type), family = "sans", size = 1.75, colour = ink,
    min.segment.length = 0, segment.size = 0.25, segment.color = "#AAB2B8",
    box.padding = 0.25, point.padding = 0.16, max.overlaps = Inf, seed = 20260804
  ) +
  scale_fill_manual(values = lineage_cols, guide = "none") +
  scale_size_continuous(range = c(2.4, 6.8), breaks = c(50, 150, 250),
                        name = "Outgoing\ninteractions") +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.18))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.10))) +
  labs(x = "Incoming weight", y = "Outgoing weight") +
  theme_pub(5.8) +
  theme(
    panel.grid.major = element_line(colour = grid_grey, linewidth = 0.30),
    legend.position = c(0.76, 0.72),
    legend.key.height = unit(2.8, "mm"),
    plot.margin = margin(4, 5, 2, 3, unit = "pt")
  )

# ---------------- panel d: Fibroblast -> immune receivers ------------------

receiver_bar <- function(value_col, axis_label, digits = 0L) {
  z <- copy(fibroblast_receiver[target %chin% immune_targets])
  setorderv(z, value_col, order = -1L)
  z[, target_plot := factor(target, levels = rev(target))]
  z[, highlighted := target == "Myeloid"]
  z[, value_label := if (digits == 0L) {
    sprintf("%.0f", get(value_col))
  } else {
    sprintf(paste0("%.", digits, "f"), get(value_col))
  }]

  ggplot(z, aes(target_plot, .data[[value_col]], fill = target)) +
    geom_col(aes(colour = highlighted), width = 0.67, linewidth = 0.55) +
    geom_text(
      aes(label = value_label), hjust = -0.15,
      family = "sans", size = 1.70, colour = ink
    ) +
    scale_fill_manual(values = lineage_cols, guide = "none") +
    scale_colour_manual(values = c(`FALSE` = "white", `TRUE` = accent_dark), guide = "none") +
    scale_x_discrete(labels = myeloid_markdown) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    coord_flip(clip = "off") +
    labs(x = NULL, y = axis_label) +
    theme_pub(5.7) +
    theme(
      axis.text.y = element_markdown(size = 5.05),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(colour = grid_grey, linewidth = 0.30),
      plot.margin = margin(3, 8, 2, 1, unit = "pt")
    )
}

p_d_count <- receiver_bar("interaction_count", "Interaction count", 0L)
p_d_weight <- receiver_bar("interaction_weight", "Interaction weight", 2L)
p_d <- plot_grid(p_d_count, p_d_weight, nrow = 1, rel_widths = c(1, 1), align = "h")

# -------- panel e: CellChat-style overall pathway signaling patterns -------

pathway_order <- unique(pathways[order(-pathway_total_weight, pathway_name), pathway_name])
pathway_plot <- copy(pathways)
pathway_plot[, pathway_name := factor(pathway_name, levels = rev(pathway_order))]
pathway_plot[, target := factor(target, levels = immune_targets)]

myeloid_only_pathways <- pathways[
  pathway_receivers == 1L & target == "Myeloid" & interaction_weight > 0,
  unique(pathway_name)
]
pathway_labels <- setNames(
  ifelse(
    pathway_order %chin% myeloid_only_pathways,
    paste0("<b><span style='color:", accent_dark, "'>", pathway_order, "</span></b>"),
    pathway_order
  ),
  pathway_order
)

receiver_totals <- pathways[, .(
  receiver_total_weight = sum(interaction_weight)
), by = target]
receiver_totals[, target := factor(target, levels = immune_targets)]

p_e_top <- ggplot(receiver_totals, aes(target, receiver_total_weight, fill = target)) +
  geom_col(width = 0.66, colour = "white", linewidth = 0.28) +
  scale_fill_manual(values = lineage_cols, guide = "none") +
  scale_y_continuous(
    limits = c(0, max(receiver_totals$receiver_total_weight) * 1.10),
    breaks = c(0, 2, 4), expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "Total weight") +
  theme_pub(5.2) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.grid.major.y = element_line(colour = grid_grey, linewidth = 0.27),
    plot.margin = margin(1, 1, 1, 2, unit = "pt")
  )

p_e_heat_base <- ggplot(pathway_plot, aes(target, pathway_name, fill = relative_strength)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradientn(
    colours = c("#F5F8F7", "#D5E9E1", "#97C9B4", "#4D9B78", "#1E684B"),
    limits = c(0, 1), breaks = c(0, 0.5, 1), name = "Relative strength"
  ) +
  scale_x_discrete() +
  scale_y_discrete(labels = pathway_labels) +
  labs(x = NULL, y = "CellChat pathway") +
  theme_pub(5.5) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_markdown(size = 4.65),
    axis.line = element_blank(), axis.ticks = element_blank(),
    legend.position = "right",
    legend.direction = "horizontal",
    legend.key.height = unit(1.35, "mm"),
    legend.key.width = unit(11, "mm"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    plot.margin = margin(0, 1, 0, 2, unit = "pt")
  ) +
  guides(fill = guide_colourbar(
    direction = "horizontal", title.position = "top", title.hjust = 0,
    barwidth = unit(18, "mm"), barheight = unit(1.35, "mm")
  ))

relative_legend <- get_legend(
  p_e_heat_base +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 5.1),
      legend.text = element_text(size = 4.8),
      legend.box.margin = margin(0, 0, 0, 0, unit = "pt")
    )
)
p_e_heat <- p_e_heat_base + theme(legend.position = "none")

receiver_strip <- data.table(
  target = factor(immune_targets, levels = immune_targets),
  strip_y = 1
)
p_e_strip <- ggplot(receiver_strip, aes(target, strip_y, fill = target)) +
  geom_tile(width = 0.97, height = 0.22, colour = "white", linewidth = 0.30) +
  scale_fill_manual(values = lineage_cols, guide = "none") +
  scale_x_discrete(labels = myeloid_markdown) +
  scale_y_continuous(limits = c(0.86, 1.14), expand = c(0, 0)) +
  labs(x = "Immune receiver", y = NULL) +
  theme_void(base_family = "sans") +
  theme(
    axis.text.x = element_markdown(
      angle = 38, hjust = 1, vjust = 1, size = 5.0,
      margin = margin(t = 1.5, unit = "pt")
    ),
    axis.title.x = element_text(size = 5.5, colour = ink, margin = margin(t = 2, unit = "pt")),
    plot.margin = margin(0, 1, 3, 2, unit = "pt")
  )

pathway_totals <- unique(pathways[, .(pathway_name, pathway_total_weight)])
pathway_totals[, pathway_name := factor(pathway_name, levels = rev(pathway_order))]

p_e_abs <- ggplot(pathway_totals, aes(pathway_total_weight, pathway_name)) +
  geom_col(width = 0.64, fill = "#AFC7BE", colour = "white", linewidth = 0.25) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Pathway weight", y = NULL) +
  theme_pub(5.3) +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.x = element_line(colour = grid_grey, linewidth = 0.28),
    plot.margin = margin(0, 3, 0, 0, unit = "pt")
  )

p_e_left <- plot_grid(
  p_e_top, p_e_heat, p_e_strip,
  ncol = 1, rel_heights = c(0.25, 1.00, 0.18),
  align = "v", axis = "lr"
)
p_e_right <- plot_grid(
  ggdraw(), p_e_abs, relative_legend,
  ncol = 1, rel_heights = c(0.25, 1.00, 0.18)
)
p_e <- plot_grid(
  p_e_left, p_e_right,
  nrow = 1, rel_widths = c(1.00, 0.24),
  align = "h", axis = "tb"
)

# ----------------------------- composition ---------------------------------

row1 <- plot_grid(
  panel_label(p_a, "a"), panel_label(p_b, "b"),
  nrow = 1, rel_widths = c(0.78, 1.32), align = "h"
)
row2 <- plot_grid(
  panel_label(p_c, "c"), panel_label(p_d, "d"),
  nrow = 1, rel_widths = c(0.78, 1.32), align = "h"
)
figure1 <- plot_grid(
  row1, row2, panel_label(p_e, "e"),
  ncol = 1,
  rel_heights = c(0.96, 0.90, 1.48),
  align = "v"
)

# ------------------------------- export ------------------------------------

width_mm <- 183
height_mm <- 182
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
stem <- file.path(main_dir, "Figure1_GSE206785_unbiased")
alias_stem <- file.path(main_dir, "Figure1_GSE206785_restructured")

agg_png(
  paste0(stem, ".png"), width = width_in, height = height_in,
  units = "in", res = 400, background = "white"
)
print(figure1)
dev.off()

agg_tiff(
  paste0(stem, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw", background = "white"
)
print(figure1)
dev.off()

cairo_pdf(
  paste0(stem, ".pdf"), width = width_in, height = height_in,
  family = "sans", onefile = TRUE
)
print(figure1)
dev.off()

svglite(
  paste0(stem, ".svg"), width = width_in, height = height_in,
  bg = "white", system_fonts = list(sans = "Arial")
)
print(figure1)
dev.off()

for (ext in c(".png", ".tiff", ".pdf", ".svg")) {
  if (!file.copy(paste0(stem, ext), paste0(alias_stem, ext), overwrite = TRUE)) {
    stop("Failed to update restructured alias for ", ext)
  }
}

cat("Restructured Figure 1 exported successfully\n")
cat("Canvas: ", width_mm, " x ", height_mm, " mm\n", sep = "")
cat("Fibroblast outgoing weight: ", sprintf("%.6f", fib_role$outgoing_weight),
    "; rank ", fib_role$outgoing_weight_rank, "/", nrow(roles), "\n", sep = "")
cat("GSE206785 Fibroblast -> Myeloid: count rank ",
    myeloid_discovery$immune_count_rank, "/7; weight rank ",
    myeloid_discovery$immune_weight_rank, "/7\n", sep = "")
cat("Main panels: a-e; donor support moved to Supplementary Figure S1; ",
    "GSE183904 deferred until formal annotation\n", sep = "")
