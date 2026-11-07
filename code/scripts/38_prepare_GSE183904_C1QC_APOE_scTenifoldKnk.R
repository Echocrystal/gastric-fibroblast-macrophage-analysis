#!/usr/bin/env Rscript

# Freeze a donor-capped C1QC+APOE+ macrophage matrix and compact GRN universe
# for an exploratory scTenifoldKnk receiver-node screen.

options(stringsAsFactors = FALSE, warn = 1)
seed <- 2026072738L
set.seed(seed)

args <- commandArgs(trailingOnly = TRUE)
screen_version <- if (length(args)) args[[1L]] else "v1"
if (!screen_version %in% c("v1", "receptor_screen_v2")) {
  stop("Unknown screen version: ", screen_version)
}

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(jsonlite)
  library(digest)
})

manifest_file <- file.path(source_project, "data", "processed", "scrna", "GSE183904",
                           "GSE183904_sample_patient_clinical_manifest.tsv")
sce_dir <- file.path(source_project, "data", "intermediate", "scrna", "GSE183904", "sample_sce")
myeloid_meta_file <- file.path(project, "analysis", "processed", "GSE183904_fine_states",
                               "GSE183904_Myeloid_fine_state_cell_metadata.tsv.gz")
program_file <- file.path(project, "resources", "gene_sets", "myeloid_receiver_programs_frozen.tsv")
analysis_label <- if (screen_version == "v1") {
  "C1QC_APOE_macrophage_v1"
} else {
  "C1QC_APOE_macrophage_receptor_screen_v2"
}
out_dir <- file.path(project, "analysis", "processed", "virtual_knockout", "GSE183904",
                     analysis_label)
resource_dir <- file.path(project, "resources", "gene_sets")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(manifest_file)[include_primary_replication == TRUE &
                                 specimen_class == "Primary_Tumor"]
meta <- fread(myeloid_meta_file)[cluster == 3L]
if (nrow(meta) != 1977L || uniqueN(meta$Patient_ID) != 26L) {
  stop("Corrected C1QC+APOE+ state size or patient coverage changed")
}

donor_cap <- 30L
selected_rows <- list()
for (i in seq_len(nrow(manifest))) {
  info <- manifest[i]
  cm <- meta[sample_id == info$sample_id]
  if (!nrow(cm)) stop("Missing corrected receiver state: ", info$sample_id)
  if (nrow(cm) > donor_cap) {
    set.seed(seed + i)
    cm <- cm[sort(sample.int(nrow(cm), donor_cap, replace = FALSE))]
  }
  selected_rows[[info$sample_id]] <- cm
}
selected_meta <- rbindlist(selected_rows)
if (any(selected_meta[, .N, by = Patient_ID]$N > donor_cap)) stop("Donor cap failure")

matrix_parts <- list()
reference_genes <- NULL
for (i in seq_len(nrow(manifest))) {
  info <- manifest[i]
  message(sprintf("[%d/26] virtual-KO extraction: %s", i, info$sample_id))
  sce <- readRDS(file.path(sce_dir, paste0(info$sample_id, "_secondary_qc_sce.rds")))
  x <- counts(sce)
  if (is.null(reference_genes)) {
    reference_genes <- rownames(x)
  } else if (!identical(reference_genes, rownames(x))) {
    stop("Gene-order mismatch: ", info$sample_id)
  }
  cm <- selected_meta[sample_id == info$sample_id]
  idx <- match(cm$cell_id, colnames(x))
  if (anyNA(idx)) stop("Cell/count alignment failure: ", info$sample_id)
  matrix_parts[[info$sample_id]] <- x[, idx, drop = FALSE]
  rm(sce, x)
  gc(verbose = FALSE)
}
x <- do.call(cbind, matrix_parts)
colnames(x) <- selected_meta$cell_id
rownames(x) <- reference_genes
if (ncol(x) != nrow(selected_meta) || anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
  stop("Extracted virtual-KO matrix integrity failure")
}

detected <- x > 0
cell_fraction <- as.numeric(Matrix::rowMeans(detected))
patient_count <- integer(nrow(x))
for (patient in unique(selected_meta$Patient_ID)) {
  idx <- which(selected_meta$Patient_ID == patient)
  patient_count <- patient_count + as.integer(Matrix::rowSums(detected[, idx, drop = FALSE]) > 0)
}
lib <- Matrix::colSums(x)
if (any(lib <= 0)) stop("Zero-library cell in receiver state")
norm <- x %*% Matrix::Diagonal(x = 1e4 / lib)
norm@x <- log1p(norm@x)
lognorm_mean <- as.numeric(Matrix::rowMeans(norm))
lognorm_second <- as.numeric(Matrix::rowMeans(norm^2))
lognorm_variance <- pmax(0, (lognorm_second - lognorm_mean^2) * ncol(norm) / (ncol(norm) - 1))

metrics <- data.table(
  gene_symbol = rownames(x),
  cell_fraction = cell_fraction,
  patient_count = patient_count,
  lognorm_mean = lognorm_mean,
  lognorm_variance = lognorm_variance
)
metrics[, passes_detection_filter := cell_fraction >= 0.05 & patient_count >= 3L]
metrics[, technical := grepl("^(MT-|RPL[0-9]|RPS[0-9]|HB[ABDEGQZ])", gene_symbol,
                             ignore.case = TRUE)]
eligible <- which(metrics$passes_detection_filter & !metrics$technical &
                  is.finite(metrics$lognorm_variance))
rank_mean <- rank(metrics$lognorm_mean[eligible], ties.method = "first")
bin <- pmin(20L, ceiling(rank_mean / (length(eligible) / 20)))
hvg_score <- rep(NA_real_, nrow(metrics))
for (b in sort(unique(bin))) {
  idx_local <- which(bin == b)
  values <- log1p(metrics$lognorm_variance[eligible[idx_local]])
  s <- sd(values)
  z <- if (!is.finite(s) || s <= 0) rep(0, length(values)) else (values - mean(values)) / s
  hvg_score[eligible[idx_local]] <- z
}
metrics[, hvg_score := hvg_score]
hvg_genes <- metrics[!is.na(hvg_score)][order(-hvg_score, gene_symbol)][1:300]$gene_symbol

focus_programs <- c(
  "GOBP_PHAGOCYTOSIS",
  "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
  "ID_INFLAMMATORY_MONO",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "GOBP_CELL_MATRIX_ADHESION",
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX"
)
targets <- if (screen_version == "v1") {
  c("ANXA2", "ITGAM", "ITGB2", "LRP1")
} else {
  c("CD74", "CCRL2", "CXCR4", "ITGB2", "LRP1", "ITGAM", "CMKLR1", "ANXA2", "TREM2")
}
mechanism_guard <- c("ILK", "PLAUR", "PTK2", "SRC", "TLN1", "VCL", "FERMT3",
                     "RHOA", "RAC1", "RELA", "JUN", "FOS", "LGMN", "CTSD", "ACP5")
forced <- unique(c(targets, mechanism_guard))
if (!all(forced %chin% rownames(x))) stop("Forced virtual-KO genes missing")

programs <- fread(program_file)[program %chin% focus_programs]
eligible_program_genes <- intersect(
  programs$gene_symbol,
  metrics[passes_detection_filter == TRUE]$gene_symbol
)
universe <- sort(unique(c(hvg_genes, eligible_program_genes, forced)))
universe_metrics <- metrics[match(universe, gene_symbol)]
universe_metrics[, `:=`(
  receiver_hvg = gene_symbol %chin% hvg_genes,
  frozen_program = gene_symbol %chin% eligible_program_genes,
  forced_target_or_guard = gene_symbol %chin% forced
)]

filtered_programs <- programs[
  gene_symbol %chin% universe & !gene_symbol %chin% c(targets, "PLAT", "ILK")
]
if (!setequal(unique(filtered_programs$program), focus_programs)) {
  stop("Program family lost during universe filtering")
}

resource_prefix <- if (screen_version == "v1") {
  "GSE183904_C1QC_APOE_scTenifoldKnk"
} else {
  "GSE183904_C1QC_APOE_receptor_screen_v2_scTenifoldKnk"
}
universe_file <- file.path(resource_dir, paste0(resource_prefix, "_universe.tsv"))
filtered_program_file <- file.path(resource_dir,
                                   paste0(resource_prefix, "_programs_axis_excluded.tsv"))
fwrite(universe_metrics, universe_file, sep = "\t")
fwrite(filtered_programs, filtered_program_file, sep = "\t")

selected_x <- x[match(universe, rownames(x)), , drop = FALSE]
saveRDS(selected_x, file.path(out_dir, "GSE183904_C1QC_APOE_sanity_counts.rds"),
        compress = "xz")
fwrite(selected_meta,
       file.path(out_dir, "GSE183904_C1QC_APOE_sanity_metadata.tsv"), sep = "\t")

patient_sizes <- selected_meta[, .N, by = Patient_ID]
audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE183904",
  state = "C1QC_APOE_macrophage",
  screen_version = screen_version,
  analysis_label = analysis_label,
  source_cluster = 3L,
  cells_before_cap = nrow(meta),
  cells_after_cap = ncol(selected_x),
  patients = uniqueN(selected_meta$Patient_ID),
  donor_cap = donor_cap,
  minimum_cells_per_patient_after_cap = min(patient_sizes$N),
  maximum_cells_per_patient_after_cap = max(patient_sizes$N),
  seed = seed,
  selection_is_outcome_blind = TRUE,
  hvg_genes = length(hvg_genes),
  hvg_method = "log1p(CP10K) variance z-scored within 20 mean-expression rank bins",
  detection_filter = list(minimum_cell_fraction = 0.05, minimum_patients = 3L),
  frozen_programs = focus_programs,
  targets = targets,
  target_tiers = if (screen_version == "receptor_screen_v2") list(
    communication_derived = c("CD74", "CCRL2", "CXCR4", "ITGB2"),
    direct_PLAT_prior_sensitivity = c("LRP1", "ITGAM"),
    chemerin_mechanism_comparator = "CMKLR1",
    surface_scaffold_comparator = "ANXA2",
    macrophage_receptor_calibrator = "TREM2"
  ) else NULL,
  mechanism_guard = mechanism_guard,
  universe_genes = length(universe),
  program_rows_after_axis_exclusion = nrow(filtered_programs),
  differential_regulation_semantics = "unsigned manifold distance",
  files = list(
    count_matrix = file.path(out_dir, "GSE183904_C1QC_APOE_sanity_counts.rds"),
    metadata = file.path(out_dir, "GSE183904_C1QC_APOE_sanity_metadata.tsv"),
    universe = universe_file,
    programs = filtered_program_file
  ),
  sha256 = list(
    count_matrix = digest(file.path(out_dir, "GSE183904_C1QC_APOE_sanity_counts.rds"),
                          algo = "sha256", file = TRUE),
    metadata = digest(file.path(out_dir, "GSE183904_C1QC_APOE_sanity_metadata.tsv"),
                      algo = "sha256", file = TRUE),
    universe = digest(universe_file, algo = "sha256", file = TRUE),
    programs = digest(filtered_program_file, algo = "sha256", file = TRUE)
  ),
  warning = paste(
    "This donor-capped compact-universe layer is an exploratory receiver-node screen.",
    "Virtual knockout output does not establish extracellular binding, signaling direction or biological causality."
  )
)
audit_file <- if (screen_version == "v1") {
  "GSE183904_C1QC_APOE_scTenifoldKnk_universe.json"
} else {
  "GSE183904_C1QC_APOE_receptor_screen_v2_scTenifoldKnk_universe.json"
}
write_json(audit,
           file.path(audit_dir, audit_file),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("Prepared matrix:", nrow(selected_x), "genes x", ncol(selected_x), "cells\n")
cat("Patients:", uniqueN(selected_meta$Patient_ID),
    "cell range after cap:", min(patient_sizes$N), "-", max(patient_sizes$N), "\n")
cat("Program rows after axis exclusion:", nrow(filtered_programs), "\n")
