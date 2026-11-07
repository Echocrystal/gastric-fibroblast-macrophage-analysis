#!/usr/bin/env Rscript

# Export the real paired WT/ITGB2-KO gene manifold from receiver-screen v2.
# Each point is a gene. No post-knockout cells are simulated.

options(stringsAsFactors = FALSE, warn = 1)
project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
result_dir <- file.path(project,
  "analysis/results/virtual_knockout/GSE183904/C1QC_APOE_macrophage_receptor_screen_v2")
out_dir <- file.path(project,
  "figures/source_data/Figure7_receiver_screen_v2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

target <- "ITGB2"
seeds <- c(2026072739L, 2026072740L, 2026072741L)

extract_one <- function(seed) {
  rds_file <- file.path(result_dir,
    sprintf("scTenifoldKnk_seed%d__%s.rds", seed, target))
  if (!file.exists(rds_file)) stop("Missing v2 fit: ", rds_file)
  fit <- readRDS(rds_file)
  alignment <- fit[["manifoldAlignment"]]
  if (!is.matrix(alignment) || ncol(alignment) != 2L) {
    stop("Unexpected manifoldAlignment structure: ", rds_file)
  }
  wt <- alignment[grepl("^X_", rownames(alignment)), , drop = FALSE]
  ko <- alignment[grepl("^Y_", rownames(alignment)), , drop = FALSE]
  rownames(wt) <- sub("^X_", "", rownames(wt))
  rownames(ko) <- sub("^Y_", "", rownames(ko))
  if (!setequal(rownames(wt), rownames(ko))) stop("WT/KO gene mismatch")
  ko <- ko[rownames(wt), , drop = FALSE]
  displacement <- sqrt(rowSums((ko - wt)^2))
  dr <- fit$diffRegulation
  dr <- dr[match(rownames(wt), dr$gene), , drop = FALSE]
  if (anyNA(dr$gene) || !identical(dr$gene, rownames(wt))) {
    stop("Differential-regulation/manifold mismatch")
  }
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
    diffreg_distance = dr$distance,
    diffreg_p_adj = dr$p.adj,
    stringsAsFactors = FALSE
  )
}

manifold <- do.call(rbind, lapply(seeds, extract_one))
if (nrow(manifold) != 3L * 775L ||
    anyDuplicated(manifold[c("target", "seed", "gene")])) {
  stop("Unexpected v2 manifold export dimensions")
}
out_file <- file.path(out_dir, "Figure7D_ITGB2_v2_gene_manifold.tsv")
write.table(manifold, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote", nrow(manifold), "paired gene-manifold rows to", out_file, "\n")

