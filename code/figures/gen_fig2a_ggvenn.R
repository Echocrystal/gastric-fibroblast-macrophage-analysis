# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) four-set Venn -> project-native ggvenn production panel -> native run
# RULE: this script renders the complete 15-region panel used by Figure 2a.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggvenn)
  library(ragg)
  library(grid)
})

# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
categorical <- c("#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666")
categorical_extended <- c(
  "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
  "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999"
)
diverging  <- c("#2166AC", "#F7F7F7", "#B2182B")
sequential <- c("#F7FBFF", "#6BAED6", "#08306B")
accent_red <- "#B2182B"
grey       <- "#999999"
black      <- "#222222"

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_file <- file.path(
  project, "figures", "data", "Figure1_GSE206785_unbiased",
  "panel_a_four_set_activity_membership.tsv.gz"
)
panel_dir <- file.path(project, "figures", "panels", "Figure2_PLAT_nomination")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

membership <- fread(source_file)
set_cols <- c(
  "Matrisome", "LR prior", "Tumour fibroblast expression", "NicheNet activity"
)
stopifnot(all(c("gene_symbol", set_cols) %in% names(membership)))
for (nm in set_cols) set(membership, j = nm, value = as.logical(membership[[nm]]))

sets_plot <- lapply(set_cols, function(nm) membership[get(nm), gene_symbol])
names(sets_plot) <- c(
  "Matrisome (513)",
  "LR prior (1,226)",
  "Fibroblast expr. (5,872)",
  "NicheNet (87)"
)

stopifnot(
  length(sets_plot[[1]]) == 513L,
  length(sets_plot[[2]]) == 1226L,
  length(sets_plot[[3]]) == 5872L,
  length(sets_plot[[4]]) == 87L,
  length(Reduce(intersect, sets_plot)) == 32L
)

venn_outline <- "#59637C"
ink <- "#000000"

panel_a <- ggvenn(
  sets_plot,
  show_elements = FALSE,
  show_percentage = FALSE,
  fill_color = c("#D7AAB4", "#E6BEA4", "#B9D0AF", "#AFC8E5"),
  fill_alpha = 0.42,
  stroke_color = venn_outline,
  stroke_alpha = 0.82,
  stroke_size = 0.45,
  set_name_color = ink,
  set_name_size = 2.50,
  text_color = ink,
  text_size = 3.00,
  show_outside = "none",
  auto_scale = FALSE
) +
  coord_fixed(xlim = c(-2.72, 2.72), ylim = c(-1.78, 1.62), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(3, 3, 3, 3, unit = "pt"))

# ggvenn 0.1.10 has no comma-separator option.  Modify only the already
# computed label strings in the built ggplot object; geometry and membership
# calculations remain entirely native to ggvenn.
panel_a_built <- ggplot_build(panel_a)
panel_a_built$data[[4]]$label <- format(
  as.integer(panel_a_built$data[[4]]$label), big.mark = ",",
  scientific = FALSE, trim = TRUE
)
panel_a_grob <- ggplot_gtable(panel_a_built)

draw_panel_a <- function() {
  grid.newpage()
  grid.draw(panel_a_grob)
}

width_mm <- 84
height_mm <- 58
png_file <- file.path(panel_dir, "Figure2A_ggvenn.png")
pdf_file <- file.path(panel_dir, "Figure2A_ggvenn.pdf")
svg_file <- file.path(panel_dir, "Figure2A_ggvenn.svg")

agg_png(png_file, width = width_mm, height = height_mm, units = "mm",
        res = 800, background = "white")
draw_panel_a()
dev.off()

cairo_pdf(pdf_file, width = width_mm / 25.4, height = height_mm / 25.4,
          family = "Arial", onefile = TRUE)
draw_panel_a()
dev.off()

svg(svg_file, width = width_mm / 25.4, height = height_mm / 25.4,
    family = "Arial", bg = "white", onefile = TRUE)
draw_panel_a()
dev.off()

message("Rendered native ggvenn Figure 2a: ", png_file)
