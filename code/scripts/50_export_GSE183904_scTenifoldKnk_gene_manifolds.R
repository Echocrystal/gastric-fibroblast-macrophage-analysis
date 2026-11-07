# Export paired WT/KO gene-manifold coordinates from the frozen
# scTenifoldKnk objects. This script performs data conversion only; all
# Figure 7 rendering is performed by the Python/matplotlib figure script.

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
result_dir <- file.path(
  project,
  "analysis/results/virtual_knockout/GSE183904/C1QC_APOE_macrophage_v1"
)
out_file <- file.path(
  project,
  "figures/data/Figure7C_paired_gene_manifold_source.tsv"
)

targets <- c("ANXA2", "ITGAM", "ITGB2", "LRP1")
seeds <- c(2026072739L, 2026072740L, 2026072741L)

extract_one <- function(target, seed) {
  rds_file <- file.path(
    result_dir,
    sprintf("scTenifoldKnk_seed%d__%s.rds", seed, target)
  )
  if (!file.exists(rds_file)) {
    stop("Missing frozen scTenifoldKnk object: ", rds_file)
  }

  fit <- readRDS(rds_file)
  alignment <- fit[["manifoldAlignment"]]
  if (!is.matrix(alignment) || ncol(alignment) != 2L) {
    stop("Unexpected manifoldAlignment structure in: ", rds_file)
  }

  wt <- alignment[grepl("^X_", rownames(alignment)), , drop = FALSE]
  ko <- alignment[grepl("^Y_", rownames(alignment)), , drop = FALSE]
  rownames(wt) <- sub("^X_", "", rownames(wt))
  rownames(ko) <- sub("^Y_", "", rownames(ko))
  if (!setequal(rownames(wt), rownames(ko))) {
    stop("WT/KO gene sets do not match in: ", rds_file)
  }
  ko <- ko[rownames(wt), , drop = FALSE]

  displacement <- sqrt(rowSums((ko - wt)^2))
  data.frame(
    target = target,
    seed = seed,
    gene = rownames(wt),
    wt_nlma1 = wt[, 1],
    wt_nlma2 = wt[, 2],
    ko_nlma1 = ko[, 1],
    ko_nlma2 = ko[, 2],
    delta_nlma1 = ko[, 1] - wt[, 1],
    delta_nlma2 = ko[, 2] - wt[, 2],
    euclidean_displacement = displacement,
    displacement_rank = rank(-displacement, ties.method = "min"),
    stringsAsFactors = FALSE
  )
}

pieces <- lapply(targets, function(target) {
  lapply(seeds, function(seed) extract_one(target, seed))
})
manifold <- do.call(rbind, unlist(pieces, recursive = FALSE))

expected_rows <- length(targets) * length(seeds) * 774L
if (nrow(manifold) != expected_rows) {
  stop("Unexpected exported row count: ", nrow(manifold), " != ", expected_rows)
}
if (anyDuplicated(manifold[c("target", "seed", "gene")])) {
  stop("Duplicate target/seed/gene rows in exported manifold table")
}

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.table(
  manifold,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Wrote", nrow(manifold), "paired gene-manifold rows to\n", out_file, "\n")
