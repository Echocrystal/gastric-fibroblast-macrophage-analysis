#!/usr/bin/env Rscript

# Deterministic aggregation of the frozen broad-lineage CellChat interaction
# table for the ecosystem-first Figure 1. No probabilities are recomputed.

options(stringsAsFactors = FALSE, warn = 1)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
figure_lib <- file.path(project, "runtime", "R_figure")
.libPaths(c(figure_lib, .libPaths()))

suppressPackageStartupMessages(library(data.table))

fig_dir <- file.path(project, "figures", "data", "Figure1_GSE206785_unbiased")
input <- file.path(fig_dir, "panel_e_cellchat_interactions.tsv.gz")
output <- file.path(fig_dir, "panel_d_cellchat_fibroblast_pathway_receiver.tsv")
stopifnot(file.exists(input))

comm <- fread(input)
required <- c("source", "target", "prob", "pathway_name", "annotation")
if (!all(required %in% names(comm))) stop("Unexpected CellChat interaction schema")

immune_targets <- c("B cells", "CD4+ T", "CD8+ T", "NK/ILC", "Mast", "Myeloid", "Plasma")
observed <- comm[
  source == "Fibroblast" & target %chin% immune_targets,
  .(
    interaction_count = .N,
    interaction_weight = sum(prob),
    annotation = annotation[1L]
  ),
  by = .(pathway_name, target)
]

pathway_annotation <- unique(observed[, .(pathway_name, annotation)])
out <- CJ(
  pathway_name = sort(unique(observed$pathway_name)),
  target = immune_targets,
  unique = TRUE
)
out <- merge(out, observed, by = c("pathway_name", "target"), all.x = TRUE, sort = FALSE)
out <- merge(
  out, pathway_annotation, by = "pathway_name", all.x = TRUE,
  suffixes = c("", "_pathway"), sort = FALSE
)
out[is.na(interaction_count), interaction_count := 0L]
out[is.na(interaction_weight), interaction_weight := 0]
out[is.na(annotation), annotation := annotation_pathway]
out[, annotation_pathway := NULL]
out[, `:=`(
  pathway_total_weight = sum(interaction_weight),
  pathway_receivers = sum(interaction_weight > 0),
  relative_strength = if (max(interaction_weight) > 0) {
    interaction_weight / max(interaction_weight)
  } else {
    rep(0, .N)
  }
), by = pathway_name]
out[, target_order := match(target, immune_targets)]
setorder(out, -pathway_total_weight, pathway_name, target_order)
out[, target_order := NULL]

if (nrow(out) != 20L * 7L) stop("Unexpected pathway-by-receiver dimensions")
myeloid_only <- out[
  pathway_receivers == 1L & target == "Myeloid" & interaction_weight > 0,
  pathway_name
]
expected_myeloid_only <- c("ANGPTL", "ANNEXIN", "COMPLEMENT", "CSF", "THY1")
if (!all(expected_myeloid_only %chin% myeloid_only)) {
  stop("Expected Myeloid-selective pathways changed")
}

fwrite(out, output, sep = "\t")
cat("CellChat pathway summary exported: ", nrow(out), " rows\n", sep = "")
cat("Myeloid-only pathways: ", paste(sort(myeloid_only), collapse = ", "), "\n", sep = "")
