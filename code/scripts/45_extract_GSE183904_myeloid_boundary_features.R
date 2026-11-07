#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

processed_dir <- file.path(project, "analysis", "processed", "GSE183904_fine_states")
sce_dir <- file.path(
  source_project, "data", "intermediate", "scrna", "GSE183904", "sample_sce"
)
meta_file <- file.path(
  processed_dir, "GSE183904_Myeloid_fine_state_cell_metadata.tsv.gz"
)
out_file <- file.path(
  processed_dir, "GSE183904_Myeloid_boundary_marker_cell_features.tsv.gz"
)

genes <- c(
  "FCN1", "AREG", "S100A8", "S100A9",
  "SPP1", "APOC1", "MARCO",
  "C1QA", "C1QB", "C1QC", "APOE", "CTSD", "LGMN", "ACP5",
  "IL1B", "CXCL8", "NFKBIA",
  "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "FCER1A", "CD207"
)

target <- fread(
  meta_file,
  select = c("cell_id", "sample_id", "Patient_ID", "cluster", "state_label", "UMAP1", "UMAP2")
)
if (nrow(target) != 9002L) stop("Expected 9,002 myeloid cells")
if (anyDuplicated(target$cell_id)) stop("Duplicated cell IDs")

for (gene in genes) {
  target[[paste0("expr_", gene)]] <- NA_real_
  target[[paste0("detect_", gene)]] <- NA
}

samples <- unique(target$sample_id)
for (i in seq_along(samples)) {
  sample_id <- samples[[i]]
  message(sprintf("[%d/%d] extracting boundary markers from %s", i, length(samples), sample_id))
  sce_file <- file.path(sce_dir, paste0(sample_id, "_secondary_qc_sce.rds"))
  if (!file.exists(sce_file)) stop("Missing SCE: ", sce_file)
  sce <- readRDS(sce_file)
  counts_mat <- SummarizedExperiment::assay(sce, "counts")

  missing_genes <- setdiff(genes, rownames(counts_mat))
  if (length(missing_genes)) {
    stop("Missing genes in ", sample_id, ": ", paste(missing_genes, collapse = ", "))
  }

  idx_target <- which(target$sample_id == sample_id)
  idx_sce <- match(target$cell_id[idx_target], colnames(counts_mat))
  if (anyNA(idx_sce)) stop("Cell alignment failure in ", sample_id)

  lib <- Matrix::colSums(counts_mat[, idx_sce, drop = FALSE])
  if (any(!is.finite(lib) | lib <= 0)) stop("Invalid library size in ", sample_id)

  raw <- as.matrix(counts_mat[genes, idx_sce, drop = FALSE])
  expr <- log1p(t(raw) * (1e4 / lib))
  detect <- t(raw > 0)
  for (j in seq_along(genes)) {
    gene <- genes[[j]]
    target[[paste0("expr_", gene)]][idx_target] <- expr[, j]
    target[[paste0("detect_", gene)]][idx_target] <- detect[, j]
  }

  rm(sce, counts_mat, raw, expr, detect)
  gc(verbose = FALSE)
}

expr_cols <- paste0("expr_", genes)
detect_cols <- paste0("detect_", genes)
if (anyNA(target[, ..expr_cols]) || anyNA(target[, ..detect_cols])) {
  stop("Incomplete marker extraction")
}

setcolorder(target, c(
  "cell_id", "sample_id", "Patient_ID", "cluster", "state_label", "UMAP1", "UMAP2",
  as.vector(rbind(expr_cols, detect_cols))
))
fwrite(target, out_file, sep = "\t", compress = "gzip")

cat("Wrote:", out_file, "\n")
cat("Cells:", nrow(target), "\n")
cat("Patients:", uniqueN(target$Patient_ID), "\n")
cat("Clusters:", paste(sort(unique(target$cluster)), collapse = ","), "\n")

