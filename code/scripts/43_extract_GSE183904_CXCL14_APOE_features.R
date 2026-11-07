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
out_file <- file.path(processed_dir, "GSE183904_CXCL14_APOE_cell_features.tsv.gz")

read_target_cells <- function(lineage) {
  f <- file.path(
    processed_dir,
    paste0("GSE183904_", lineage, "_fine_state_cell_metadata.tsv.gz")
  )
  x <- fread(f, select = c("cell_id", "sample_id", "lineage"))
  if (anyDuplicated(x$cell_id)) stop("Duplicated cell IDs in ", lineage)
  x
}

target <- rbindlist(list(
  read_target_cells("Fibroblast"),
  read_target_cells("Myeloid")
))
target[, `:=`(
  expr_CXCL14 = NA_real_, detect_CXCL14 = NA,
  expr_APOE = NA_real_, detect_APOE = NA
)]

genes <- c("CXCL14", "APOE")
samples <- unique(target$sample_id)

for (i in seq_along(samples)) {
  sample_id <- samples[[i]]
  message(sprintf("[%d/%d] extracting CXCL14/APOE from %s", i, length(samples), sample_id))
  sce_file <- file.path(sce_dir, paste0(sample_id, "_secondary_qc_sce.rds"))
  if (!file.exists(sce_file)) stop("Missing SCE: ", sce_file)
  sce <- readRDS(sce_file)
  counts_mat <- SummarizedExperiment::assay(sce, "counts")

  idx_target <- which(target$sample_id == sample_id)
  idx_sce <- match(target$cell_id[idx_target], colnames(counts_mat))
  if (anyNA(idx_sce)) stop("Cell alignment failure in ", sample_id)

  lib <- Matrix::colSums(counts_mat[, idx_sce, drop = FALSE])
  if (any(!is.finite(lib) | lib <= 0)) stop("Invalid library size in ", sample_id)

  for (gene in genes) {
    if (!(gene %in% rownames(counts_mat))) stop("Missing gene ", gene, " in ", sample_id)
    raw <- as.numeric(counts_mat[gene, idx_sce, drop = TRUE])
    expr <- log1p(raw * 1e4 / lib)
    target[[paste0("expr_", gene)]][idx_target] <- expr
    target[[paste0("detect_", gene)]][idx_target] <- raw > 0
  }
  rm(sce, counts_mat)
  gc(verbose = FALSE)
}

if (anyNA(target$expr_CXCL14) || anyNA(target$expr_APOE)) {
  stop("Incomplete feature extraction")
}

setcolorder(target, c(
  "cell_id", "sample_id", "lineage",
  "expr_CXCL14", "detect_CXCL14", "expr_APOE", "detect_APOE"
))
fwrite(target, out_file, sep = "\t", compress = "gzip")

cat("Wrote:", out_file, "\n")
cat("Cells:", nrow(target), "\n")
cat("Fibroblast CXCL14 detection:",
    mean(target[lineage == "Fibroblast"]$detect_CXCL14), "\n")
cat("Myeloid APOE detection:",
    mean(target[lineage == "Myeloid"]$detect_APOE), "\n")
