#!/usr/bin/env Rscript

# Lossless data-conversion step for the Python-only Figure 2 atlas workflow.
# This script does not normalize, select features, cluster, embed, or draw.

options(stringsAsFactors = FALSE, warn = 1)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(rhdf5)
  library(jsonlite)
  library(digest)
})

manifest_file <- file.path(
  source_project, "data", "processed", "scrna", "GSE183904",
  "GSE183904_sample_patient_clinical_manifest.tsv"
)
cell_meta_file <- file.path(
  source_project, "data", "processed", "scrna", "GSE183904",
  "GSE183904_secondary_qc_coarse_annotation_metadata.tsv.gz"
)
sce_dir <- file.path(
  source_project, "data", "intermediate", "scrna", "GSE183904", "sample_sce"
)
out_dir <- file.path(project, "analysis", "processed", "GSE183904_global_atlas")
raw_dir <- file.path(out_dir, "raw_h5")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(manifest_file)[
  include_primary_replication == TRUE & specimen_class == "Primary_Tumor"
]
cell_meta <- fread(cell_meta_file)

if (nrow(manifest) != 26L || uniqueN(manifest$Patient_ID) != 26L) {
  stop("Expected 26 patient-independent primary gastric tumours")
}

metadata_parts <- vector("list", nrow(manifest))
reference_genes <- NULL
sample_audit <- vector("list", nrow(manifest))

for (k in seq_len(nrow(manifest))) {
  info <- manifest[k]
  message(sprintf("[%02d/26] exporting %s", k, info$sample_id))

  sce_file <- file.path(sce_dir, paste0(info$sample_id, "_secondary_qc_sce.rds"))
  sce <- readRDS(sce_file)
  cm <- cell_meta[sample_id == info$sample_id]
  idx <- match(colnames(sce), cm$cell_id)
  if (anyNA(idx)) stop("Cell alignment failure for ", info$sample_id)
  cm <- cm[idx]

  keep <- which(
    cm$passes_secondary_qc_singlet &
      cm$specimen_class == "Primary_Tumor" &
      cm$Patient_ID == info$Patient_ID
  )
  if (!length(keep)) stop("No eligible primary-tumour singlets for ", info$sample_id)

  mat <- counts(sce)[, keep, drop = FALSE]
  mat <- as(mat, "dgCMatrix")
  genes <- rownames(mat)
  if (is.null(reference_genes)) {
    reference_genes <- genes
    fwrite(data.table(gene_symbol = reference_genes),
           file.path(out_dir, "GSE183904_global_atlas_genes.tsv"), sep = "\t")
  } else if (!identical(reference_genes, genes)) {
    stop("Gene-order mismatch in ", info$sample_id)
  }

  h5_file <- file.path(raw_dir, paste0(info$sample_id, "_counts_csc.h5"))
  if (file.exists(h5_file)) file.remove(h5_file)
  h5createFile(h5_file)
  h5write(as.integer(mat@p), h5_file, "indptr")
  h5write(as.integer(mat@i), h5_file, "indices")
  h5write(as.numeric(mat@x), h5_file, "data")
  h5write(as.integer(dim(mat)), h5_file, "shape_genes_by_cells")
  h5write(as.character(colnames(mat)), h5_file, "cell_id")
  H5close()

  part <- cm[keep, .(
    cell_id, sample_id, Patient_ID, specimen_class, coarse_label,
    sum_counts, detected_genes, mito_percent,
    secondary_scDblFinder_score, annotation_margin
  )]
  part[, local_matrix_column := seq_len(.N)]
  part[, h5_file := basename(h5_file)]
  metadata_parts[[k]] <- part

  sample_audit[[k]] <- data.table(
    sample_id = info$sample_id,
    Patient_ID = info$Patient_ID,
    n_cells = ncol(mat),
    n_genes = nrow(mat),
    nonzero_entries = length(mat@x),
    source_sce_sha256 = digest(sce_file, algo = "sha256", file = TRUE),
    exported_h5_sha256 = digest(h5_file, algo = "sha256", file = TRUE)
  )

  rm(sce, cm, mat, part)
  gc(verbose = FALSE)
}

metadata <- rbindlist(metadata_parts)
sample_summary <- rbindlist(sample_audit)
fwrite(metadata, file.path(out_dir, "GSE183904_global_atlas_cells.tsv.gz"), sep = "\t")
fwrite(sample_summary, file.path(out_dir, "GSE183904_global_atlas_export_summary.tsv"), sep = "\t")

if (nrow(metadata) != 109788L) {
  stop("Unexpected primary-tumour singlet count: ", nrow(metadata))
}
if (uniqueN(metadata$cell_id) != nrow(metadata)) stop("Duplicate cell IDs after export")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  purpose = "lossless raw-count conversion for Python-only global atlas analysis",
  dataset = "GSE183904",
  primary_tumours = nrow(manifest),
  patients = uniqueN(metadata$Patient_ID),
  cells = nrow(metadata),
  genes = length(reference_genes),
  nonzero_entries = sum(sample_summary$nonzero_entries),
  transformations = list(
    normalization = FALSE,
    feature_selection = FALSE,
    dimensionality_reduction = FALSE,
    clustering = FALSE,
    visualization = FALSE
  ),
  input_sha256 = list(
    manifest = digest(manifest_file, algo = "sha256", file = TRUE),
    cell_metadata = digest(cell_meta_file, algo = "sha256", file = TRUE)
  )
)
write_json(
  audit,
  file.path(audit_dir, "GSE183904_global_atlas_raw_export.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA
)

cat(sprintf(
  "Export complete: %d cells, %d genes, %d samples, %s nonzero entries\n",
  nrow(metadata), length(reference_genes), nrow(sample_summary),
  format(sum(sample_summary$nonzero_entries), big.mark = ",", scientific = FALSE)
))
