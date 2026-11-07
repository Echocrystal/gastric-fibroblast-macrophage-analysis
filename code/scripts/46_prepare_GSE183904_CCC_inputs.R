#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260801L)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
})

manifest_file <- file.path(
  source_project, "data", "processed", "scrna", "GSE183904",
  "GSE183904_sample_patient_clinical_manifest.tsv"
)
coarse_file <- file.path(
  source_project, "data", "processed", "scrna", "GSE183904",
  "GSE183904_secondary_qc_coarse_annotation_metadata.tsv.gz"
)
sce_dir <- file.path(
  source_project, "data", "intermediate", "scrna", "GSE183904", "sample_sce"
)
fine_dir <- file.path(project, "analysis", "processed", "GSE183904_fine_states")
out_dir <- file.path(project, "analysis", "processed", "communication", "GSE183904")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(manifest_file)[
  include_primary_replication == TRUE & specimen_class == "Primary_Tumor"
]
stopifnot(nrow(manifest) == 26L, uniqueN(manifest$Patient_ID) == 26L)

coarse <- fread(coarse_file)[
  sample_id %in% manifest$sample_id & passes_secondary_qc_singlet == TRUE
]
coarse_keep <- c(
  "Epithelial", "Endothelial", "Fibroblast", "Mural", "T_cell",
  "NK", "B_cell", "Plasma", "Myeloid", "Mast"
)
coarse <- coarse[coarse_label %in% coarse_keep]
setorder(coarse, sample_id, coarse_label, cell_id)

fib <- fread(file.path(
  fine_dir, "GSE183904_Fibroblast_fine_state_cell_metadata.tsv.gz"
))
mye <- fread(file.path(
  fine_dir, "GSE183904_Myeloid_fine_state_cell_metadata.tsv.gz"
))

fib <- fib[cluster %in% c(1L, 2L, 3L)]
fib[, ccc_group := fcase(
  cluster == 1L, "F1_PLATrich",
  cluster == 2L, "F2_MuralControl",
  cluster == 3L, "F3_AntigenPresentation"
)]
mye_all5 <- mye[cluster %in% c(1L, 2L, 3L, 4L, 5L)]
mye_all5[, ccc_group := fcase(
  cluster == 1L, "M1_Unresolved",
  cluster == 2L, "M2_SPP1_APOE_sensitivity",
  cluster == 3L, "M3_C1QC_APOE_target",
  cluster == 4L, "M4_NeutrophilLike",
  cluster == 5L, "M5_C1QC_boundary"
)]
mye <- mye_all5[cluster != 2L]
fine <- rbindlist(list(
  fib[, .(cell_id, sample_id, Patient_ID, lineage, cluster, state_label, ccc_group)],
  mye[, .(cell_id, sample_id, Patient_ID, lineage, cluster, state_label, ccc_group)]
))
if (anyDuplicated(fine$cell_id)) stop("Duplicated fine-state cell identifiers")
fine_all5 <- rbindlist(list(
  fib[, .(cell_id, sample_id, Patient_ID, lineage, cluster, state_label, ccc_group)],
  mye_all5[, .(cell_id, sample_id, Patient_ID, lineage, cluster, state_label, ccc_group)]
))
if (anyDuplicated(fine_all5$cell_id)) stop("Duplicated all-five-state cell identifiers")

deterministic_cap <- function(dt, strata, cap = 100L) {
  z <- copy(dt)
  z[, random_key := runif(.N)]
  setorderv(z, c(strata, "random_key"))
  z <- z[, head(.SD, cap), by = strata]
  z[, random_key := NULL]
  z
}

global_balanced <- deterministic_cap(
  coarse[, .(cell_id, sample_id, Patient_ID, coarse_label)],
  c("Patient_ID", "coarse_label"), 100L
)
fine_balanced <- deterministic_cap(
  fine, c("Patient_ID", "ccc_group"), 100L
)
fine_all5_balanced <- deterministic_cap(
  fine_all5, c("Patient_ID", "ccc_group"), 100L
)

write_sparse_bundle <- function(selection, prefix, group_col) {
  expected <- file.path(out_dir, c(
    paste0(prefix, "_raw_counts_genes_by_cells.mtx"),
    paste0(prefix, "_genes.tsv"),
    paste0(prefix, "_cells.tsv"),
    paste0(prefix, "_group_summary.tsv")
  ))
  if (all(file.exists(expected)) && all(file.info(expected)$size > 0)) {
    message(prefix, ": complete bundle already exists; skipping")
    return(invisible(NULL))
  }
  mats <- list()
  metas <- list()
  reference_genes <- NULL
  for (sid in manifest$sample_id) {
    wanted <- selection[sample_id == sid]
    if (!nrow(wanted)) next
    message(prefix, ": load ", sid, " (", nrow(wanted), " selected cells)")
    sce <- readRDS(file.path(sce_dir, paste0(sid, "_secondary_qc_sce.rds")))
    idx <- match(wanted$cell_id, colnames(sce))
    if (anyNA(idx)) stop("Cell alignment failure in ", sid)
    mat <- counts(sce)[, idx, drop = FALSE]
    colnames(mat) <- wanted$cell_id
    if (is.null(reference_genes)) {
      reference_genes <- rownames(mat)
    } else if (!identical(reference_genes, rownames(mat))) {
      stop("Gene-order mismatch in ", sid)
    }
    mats[[sid]] <- mat
    metas[[sid]] <- wanted
    rm(sce, mat)
    gc(verbose = FALSE)
  }
  mat <- do.call(cbind, mats)
  meta <- rbindlist(metas)
  idx <- match(colnames(mat), meta$cell_id)
  if (anyNA(idx)) stop("Metadata alignment failure for ", prefix)
  meta <- meta[idx]
  if (!identical(colnames(mat), meta$cell_id)) stop("Order mismatch for ", prefix)
  Matrix::writeMM(mat, file.path(out_dir, paste0(prefix, "_raw_counts_genes_by_cells.mtx")))
  fwrite(data.table(gene = rownames(mat)), file.path(out_dir, paste0(prefix, "_genes.tsv")), sep = "\t")
  fwrite(meta, file.path(out_dir, paste0(prefix, "_cells.tsv")), sep = "\t")
  summary <- meta[, .(
    n_cells = .N,
    n_patients = uniqueN(Patient_ID)
  ), by = group_col]
  setnames(summary, group_col, "group")
  fwrite(summary, file.path(out_dir, paste0(prefix, "_group_summary.tsv")), sep = "\t")
  message(prefix, ": wrote ", nrow(mat), " genes x ", ncol(mat), " cells")
}

write_sparse_bundle(global_balanced, "GSE183904_global_donor_balanced", "coarse_label")
write_sparse_bundle(fine_balanced, "GSE183904_fine_donor_balanced", "ccc_group")
write_sparse_bundle(fine, "GSE183904_fine_full", "ccc_group")
write_sparse_bundle(fine_all5_balanced, "GSE183904_fine_all5_donor_balanced", "ccc_group")
write_sparse_bundle(fine_all5, "GSE183904_fine_all5_full", "ccc_group")

fwrite(
  coarse[, .(n_cells = .N, n_patients = uniqueN(Patient_ID)), by = coarse_label],
  file.path(out_dir, "GSE183904_global_full_group_summary.tsv"), sep = "\t"
)

cat("Communication input preparation complete.\n")
