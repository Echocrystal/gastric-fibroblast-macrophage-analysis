#!/usr/bin/env Rscript

# CellChat broad-lineage communication analysis for Figure 1. The input is the
# frozen donor-balanced tumour-cell h5ad object. This script estimates CellChat
# mass-action communication probabilities among Fibroblast and seven broad
# immune lineages and exports count/weight matrices for an R-drawn classic
# circle network. It does not use PLAT status to select cells, groups or edges.

options(stringsAsFactors = FALSE, warn = 1)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
figure_lib <- file.path(project, "runtime", "R_figure")
cellchat_lib <- "I:/ECM-new/scRNA/20.rebuild_GSE163558/R_libs_4.4"
.libPaths(c(figure_lib, cellchat_lib, .libPaths()))

suppressPackageStartupMessages({
  library(hdf5r)
  library(Matrix)
  library(CellChat)
  library(data.table)
})

input_h5ad <- file.path(
  project, "analysis", "results", "figure1_gse206785_unbiased",
  "GSE206785_donor_balanced_all_lineage.h5ad"
)
out_dir <- file.path(project, "analysis", "results", "figure1_gse206785_unbiased")
fig_dir <- file.path(project, "figures", "data", "Figure1_GSE206785_unbiased")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_h5ad))
set.seed(20260803L)

read_h5ad_csr <- function(path) {
  f <- H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)

  if (!all(c("data", "indices", "indptr") %in% names(f[["X"]]))) {
    stop("The h5ad X slot is not encoded as CSR")
  }
  data <- f[["X/data"]][]
  indices <- as.integer(f[["X/indices"]][])
  indptr <- as.integer(f[["X/indptr"]][])
  shape <- as.integer(h5attr(f[["X"]], "shape"))
  if (length(shape) != 2L || length(indptr) != shape[1] + 1L) {
    stop("Invalid h5ad CSR dimensions")
  }
  x_cell_gene <- new(
    "dgRMatrix", p = indptr, j = indices, x = as.numeric(data), Dim = shape
  )

  obs_index_name <- as.character(h5attr(f[["obs"]], "_index"))
  var_index_name <- as.character(h5attr(f[["var"]], "_index"))
  cell_ids <- as.character(f[[paste0("obs/", obs_index_name)]][])
  genes <- as.character(f[[paste0("var/", var_index_name)]][])
  if (length(cell_ids) != shape[1] || length(genes) != shape[2]) {
    stop("h5ad observation/variable names do not match X")
  }

  read_obs_column <- function(name) {
    node <- f[[paste0("obs/", name)]]
    if (inherits(node, "H5Group")) {
      categories <- as.character(node[["categories"]][])
      codes <- as.integer(node[["codes"]][])
      out <- rep(NA_character_, length(codes))
      valid <- codes >= 0L
      out[valid] <- categories[codes[valid] + 1L]
      out
    } else {
      as.character(node[])
    }
  }

  meta <- data.frame(
    cell_id = cell_ids,
    cell_type = read_obs_column("cell_type"),
    Patient = read_obs_column("Patient"),
    Annotation = read_obs_column("Annotation"),
    row.names = cell_ids,
    check.names = FALSE
  )
  expr <- as(t(x_cell_gene), "CsparseMatrix")
  rownames(expr) <- genes
  colnames(expr) <- cell_ids
  list(expr = expr, meta = meta)
}

message("Reading frozen donor-balanced h5ad")
obj <- read_h5ad_csr(input_h5ad)
groups <- c(
  "Fibroblast", "B cells", "CD4+ T", "CD8+ T",
  "NK/ILC", "Mast", "Myeloid", "Plasma"
)
keep <- obj$meta$cell_type %in% groups
expr <- obj$expr[, keep, drop = FALSE]
meta <- obj$meta[keep, , drop = FALSE]
meta$cell_type <- factor(meta$cell_type, levels = groups)
meta$samples <- meta$Patient
if (any(table(meta$cell_type) < 10L)) stop("A CellChat group has fewer than 10 cells")

coverage <- as.data.table(meta)[, .(
  cells = .N,
  patients = uniqueN(Patient)
), by = cell_type]
coverage[, group_order := match(cell_type, groups)]
setorder(coverage, group_order)
coverage[, group_order := NULL]
fwrite(coverage, file.path(fig_dir, "panel_e_cellchat_group_coverage.tsv"), sep = "\t")

message("Running CellChat v", as.character(packageVersion("CellChat")))
cellchat <- createCellChat(object = expr, meta = meta, group.by = "cell_type")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(
  cellchat,
  type = "triMean",
  trim = 0.1,
  raw.use = TRUE,
  population.size = FALSE,
  nboot = 100
)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file.path(out_dir, "GSE206785_CellChat_broad_fibroblast_immune.rds"))

matrix_to_long <- function(mat, value_name) {
  z <- as.data.table(as.table(mat))
  setnames(z, c("source", "target", value_name))
  z[, source := as.character(source)]
  z[, target := as.character(target)]
  z[]
}

count_long <- matrix_to_long(cellchat@net$count, "interaction_count")
weight_long <- matrix_to_long(cellchat@net$weight, "interaction_weight")
network <- merge(count_long, weight_long, by = c("source", "target"), all = TRUE)
network[is.na(interaction_count), interaction_count := 0]
network[is.na(interaction_weight), interaction_weight := 0]
network[, source_order := match(source, groups)]
network[, target_order := match(target, groups)]
setorder(network, source_order, target_order)
fwrite(network, file.path(fig_dir, "panel_e_cellchat_broad_network.tsv"), sep = "\t")

comm <- as.data.table(subsetCommunication(cellchat))
fwrite(comm, file.path(fig_dir, "panel_e_cellchat_interactions.tsv.gz"), sep = "\t")

# Receiver-resolved pathway summary for a transparent comparison of absolute
# fibroblast outgoing strength and pathway selectivity. Relative strength is
# normalized within each pathway across the seven immune receivers; the raw
# CellChat probability sum remains available in interaction_weight.
immune_targets <- setdiff(groups, "Fibroblast")
pathway_observed <- comm[
  source == "Fibroblast" & target %chin% immune_targets,
  .(
    interaction_count = .N,
    interaction_weight = sum(prob),
    annotation = annotation[1L]
  ),
  by = .(pathway_name, target)
]
pathway_annotation <- unique(pathway_observed[, .(pathway_name, annotation)])
pathway_grid <- CJ(
  pathway_name = sort(unique(pathway_observed$pathway_name)),
  target = immune_targets,
  unique = TRUE
)
pathway_grid <- merge(
  pathway_grid, pathway_observed,
  by = c("pathway_name", "target"), all.x = TRUE, sort = FALSE
)
pathway_grid <- merge(
  pathway_grid, pathway_annotation,
  by = "pathway_name", all.x = TRUE, suffixes = c("", "_pathway"), sort = FALSE
)
pathway_grid[is.na(interaction_count), interaction_count := 0L]
pathway_grid[is.na(interaction_weight), interaction_weight := 0]
pathway_grid[is.na(annotation), annotation := annotation_pathway]
pathway_grid[, annotation_pathway := NULL]
pathway_grid[, `:=`(
  pathway_total_weight = sum(interaction_weight),
  pathway_receivers = sum(interaction_weight > 0),
  relative_strength = if (max(interaction_weight) > 0) {
    interaction_weight / max(interaction_weight)
  } else {
    rep(0, .N)
  }
), by = pathway_name]
pathway_grid[, target_order := match(target, immune_targets)]
setorder(pathway_grid, -pathway_total_weight, pathway_name, target_order)
pathway_grid[, target_order := NULL]
fwrite(
  pathway_grid,
  file.path(fig_dir, "panel_d_cellchat_fibroblast_pathway_receiver.tsv"),
  sep = "\t"
)

sender_summary <- network[source == "Fibroblast", .(
  source, target, interaction_count, interaction_weight
)]
sender_summary[, count_rank := frank(-interaction_count, ties.method = "min")]
sender_summary[, weight_rank := frank(-interaction_weight, ties.method = "min")]
fwrite(sender_summary, file.path(fig_dir, "panel_e_cellchat_fibroblast_receiver_summary.tsv"), sep = "\t")

manifest <- data.table(
  field = c(
    "dataset", "input", "groups", "database", "expression_summary",
    "population_size", "bootstrap", "minimum_cells", "claim_boundary"
  ),
  value = c(
    "GSE206785", "donor-balanced tumour cells", paste(groups, collapse = "; "),
    "CellChatDB.human", "trimean", "FALSE", "100", "10",
    "mass-action communication probability; not binding, spatial contact or causality"
  )
)
fwrite(manifest, file.path(audit_dir, "Figure1_GSE206785_CellChat_manifest.tsv"), sep = "\t")

print(coverage)
print(sender_summary[order(-interaction_weight)])
cat("Nonzero broad CellChat edges:", network[interaction_weight > 0, .N], "\n")
