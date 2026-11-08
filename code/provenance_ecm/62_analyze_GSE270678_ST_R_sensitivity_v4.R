#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else
  "I:/ECM-new/ECM-public"
legacy_lib <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026/runtime/R_library"
project_lib <- file.path(project_dir, "environment", "R_library")
.libPaths(c(project_lib, legacy_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(rhdf5)
  library(digest)
})

set.seed(270679L)

model_dir <- file.path(project_dir, "models", "spatial", "v4")
lock_file <- file.path(model_dir, "GSE270678_project_compartment_analysis_lock_v4.rds")
lock_hash_file <- file.path(model_dir, "GSE270678_project_compartment_lock_hashes_v4.tsv")
sample_manifest_file <- file.path(project_dir, "manifests", "v4", "GSE270678_spatial_sample_manifest.tsv")
patient_manifest_file <- file.path(project_dir, "manifests", "v4", "GSE270678_spatial_patient_manifest.tsv")
verification_file <- file.path(project_dir, "manifests", "v4", "GSE270678_ST_R_file_verification.tsv")
main_patient_file <- file.path(project_dir, "results", "spatial", "v4",
                               "GSE270678_project_compartments",
                               "GSE270678_primary_endpoints_patient_level_v4.tsv")
input_dir <- file.path(project_dir, "data", "raw", "spatial", "GSE270678", "ST_R_public_files")
required <- c(lock_file, lock_hash_file, sample_manifest_file, patient_manifest_file,
              verification_file, main_patient_file)
if (!all(file.exists(required))) stop("Missing ST_R inputs/lock/main result")

result_dir <- file.path(project_dir, "results", "spatial", "v4", "GSE270678_ST_R_sensitivity")
processed_dir <- file.path(project_dir, "data", "processed", "spatial", "GSE270678", "v4")
audit_dir <- file.path(project_dir, "audit", "v4")
for (d in c(result_dir, processed_dir, audit_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

lock_hashes <- fread(lock_hash_file)
lock_inputs <- lock_hashes[grepl("GSE270678_project_compartment_(gene_sets|endpoints|analysis_lock)_v4|tcga_ecm_centroid|GSE183904_coarse|scrna_spatial_programs|GSE270678_spatial_sample_manifest", artifact)]
for (i in seq_len(nrow(lock_inputs))) {
  observed <- digest(lock_inputs$artifact[[i]], algo = "sha256", file = TRUE)
  if (!identical(tolower(observed), tolower(lock_inputs$sha256[[i]]))) {
    stop("Pre-expression lock hash mismatch: ", lock_inputs$artifact[[i]])
  }
}

lock <- readRDS(lock_file)
gene_sets_dt <- as.data.table(lock$gene_sets)
gene_sets <- split(gene_sets_dt$gene_symbol, gene_sets_dt$score)
endpoints <- as.data.table(lock$endpoints)[tier == "primary"]
samples_all <- fread(sample_manifest_file)
samples <- samples_all[platform_role == "ST_R"]
patients <- fread(patient_manifest_file)[main_cohort_donor == TRUE]
eligible_donors <- patients$donor_id
verification <- fread(verification_file)
if (nrow(samples) != 39L || uniqueN(samples$donor_id) != 17L) stop("Unexpected ST_R cohort")
if (nrow(verification) != 156L || !all(verification$size_verified == TRUE)) stop("ST_R files not verified")
if (!all(samples$donor_id[samples$donor_id != "D30"] %chin% eligible_donors)) stop("Unexpected extension donor")

robust_qc <- function(sum_counts, detected, mito) {
  ls <- log10(sum_counts + 1); ld <- log10(detected + 1)
  ms <- median(ls); mads <- mad(ls, constant = 1.4826)
  md <- median(ld); madd <- mad(ld, constant = 1.4826)
  mm <- median(mito); madm <- mad(mito, constant = 1.4826)
  if (!is.finite(mads) || mads == 0) mads <- IQR(ls) / 1.349
  if (!is.finite(madd) || madd == 0) madd <- IQR(ld) / 1.349
  if (!is.finite(madm) || madm == 0) madm <- IQR(mito) / 1.349
  if (!is.finite(mads)) mads <- 0
  if (!is.finite(madd)) madd <- 0
  if (!is.finite(madm)) madm <- 0
  list(min_sum = max(0, 10^(ms - 4 * mads) - 1),
       min_detected = max(0, 10^(md - 4 * madd) - 1),
       max_mito = min(100, mm + 5 * madm))
}

read_10x_h5 <- function(path) {
  data <- h5read(path, "/matrix/data")
  indices <- h5read(path, "/matrix/indices")
  indptr <- h5read(path, "/matrix/indptr")
  shape <- as.integer(h5read(path, "/matrix/shape"))
  barcodes <- as.character(h5read(path, "/matrix/barcodes"))
  genes <- as.character(h5read(path, "/matrix/features/name"))
  feature_type <- tryCatch(as.character(h5read(path, "/matrix/features/feature_type")),
                           error = function(e) rep("Gene Expression", length(genes)))
  x <- sparseMatrix(i = as.integer(indices) + 1L, p = as.integer(indptr),
                    x = as.numeric(data), dims = shape, giveCsparse = TRUE)
  keep <- feature_type == "Gene Expression"
  list(counts = x[keep, , drop = FALSE], genes = genes[keep], barcodes = barcodes)
}

read_positions <- function(path) {
  z <- fread(path, header = FALSE)
  if (tolower(as.character(z[[1L]][1L])) == "barcode") {
    z <- fread(path, header = TRUE)
  }
  if (ncol(z) < 6L) stop("Unexpected positions file: ", path)
  setnames(z, names(z)[1:6], c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  z[, .(barcode = as.character(barcode), in_tissue = as.integer(in_tissue),
        array_row = as.integer(array_row), array_col = as.integer(array_col),
        pxl_row = as.numeric(pxl_row), pxl_col = as.numeric(pxl_col))]
}

spatial_lag <- function(values, array_row, array_col) {
  offsets <- rbind(c(-1, -1), c(-1, 1), c(0, -2), c(0, 2), c(1, -1), c(1, 1))
  keys <- paste(array_row, array_col, sep = ":")
  idx <- do.call(cbind, lapply(seq_len(nrow(offsets)), function(i) {
    match(paste(array_row + offsets[i, 1], array_col + offsets[i, 2], sep = ":"), keys)
  }))
  out <- rep(NA_real_, length(values))
  for (i in seq_along(values)) {
    j <- idx[i, ]; j <- j[is.finite(j)]
    if (length(j)) out[[i]] <- mean(values[j], na.rm = TRUE)
  }
  out
}

spot_outputs <- list(); qc_outputs <- list(); slide_outputs <- list(); coverage_outputs <- list()

for (ii in seq_len(nrow(samples))) {
  info <- samples[ii]
  message(sprintf("[%02d/%02d] %s (%s)", ii, nrow(samples), info$sample_id, info$donor_id))
  h5_file <- file.path(input_dir, info$filtered_h5)
  pos_file <- file.path(input_dir, info$positions)
  verify_rows <- verification[Name %chin% c(info$filtered_h5, info$positions, info$scalefactors)]
  if (nrow(verify_rows) != 3L || !all(file.exists(verify_rows$destination)) ||
      !all(as.numeric(file.info(verify_rows$destination)$size) == verify_rows$expected_size_bytes)) {
    stop("ST_R size verification failed for ", info$sample_id)
  }
  h <- read_10x_h5(h5_file)
  pos <- read_positions(pos_file)[in_tissue == 1L]
  common <- intersect(h$barcodes, pos$barcode)
  if (length(common) < 100L) stop("Fewer than 100 in-tissue barcodes")
  counts <- h$counts[, match(common, h$barcodes), drop = FALSE]
  pos <- pos[match(common, barcode)]
  sum_counts <- as.numeric(colSums(counts))
  detected <- as.numeric(colSums(counts > 0))
  mito_rows <- grepl("^MT-", h$genes)
  mito <- 100 * if (any(mito_rows)) as.numeric(colSums(counts[mito_rows, , drop = FALSE])) else 0
  mito <- mito / pmax(sum_counts, 1)
  thresholds <- robust_qc(sum_counts, detected, mito)
  keep <- sum_counts >= thresholds$min_sum & detected >= thresholds$min_detected & mito <= thresholds$max_mito
  if (sum(keep) < 100L) stop("Fewer than 100 QC spots")
  counts <- counts[, keep, drop = FALSE]
  pos <- pos[keep]
  sum_counts <- sum_counts[keep]

  requested <- unique(gene_sets_dt$gene_symbol)
  rel_counts <- vapply(requested, function(g) {
    rows <- which(h$genes == g)
    if (!length(rows)) return(rep(0, ncol(counts)))
    as.numeric(colSums(counts[rows, , drop = FALSE]))
  }, numeric(ncol(counts)))
  rel_counts <- t(rel_counts); rownames(rel_counts) <- requested
  norm <- log1p(sweep(rel_counts, 2L, 1e4 / pmax(sum_counts, 1), "*"))
  z <- t(scale(t(norm), center = TRUE, scale = TRUE)); z[!is.finite(z)] <- 0
  scores <- as.matrix(vapply(names(gene_sets), function(nm) {
    colMeans(z[gene_sets[[nm]], , drop = FALSE])
  }, numeric(ncol(z))))
  spot <- cbind(pos, data.table(sample_id = info$sample_id, main_slice_id = info$main_slice_id,
                                donor_id = info$donor_id, stage = info$stage_geo,
                                lauren = info$lauren_geo,
                                sampling_location = info$sampling_location_geo),
                as.data.table(scores))
  spot[, frozen_ECM_direction := FROZEN_ECM_C2_HIGH50 - FROZEN_ECM_C1_HIGH50]
  spot[, epithelial_dominance := MARKER_Epithelial - MARKER_Fibroblast]
  spot[, stromal_dominance := MARKER_Fibroblast - MARKER_Epithelial]
  spot[, epithelial_rich := epithelial_dominance >= quantile(epithelial_dominance, .75, na.rm = TRUE)]
  spot[, stromal_rich := stromal_dominance >= quantile(stromal_dominance, .75, na.rm = TRUE)]
  q_ecm <- quantile(spot$frozen_ECM_direction, c(.25, .5, .75), na.rm = TRUE)
  spot[, `:=`(
    ECM_bottom_quartile = frozen_ECM_direction <= q_ecm[[1]],
    low_ECM_context = frozen_ECM_direction < q_ecm[[2]],
    high_ECM_context = frozen_ECM_direction >= q_ecm[[2]],
    ECM_top_quartile = frozen_ECM_direction >= q_ecm[[3]],
    lag_MARKER_T_cell = spatial_lag(MARKER_T_cell, array_row, array_col)
  )]

  min_group <- as.integer(lock$definitions$minimum_contrast_group_spots)
  add_simple <- function(id, g1, g2, variable) {
    ok <- nrow(g1) >= min_group && nrow(g2) >= min_group
    data.table(sample_id = info$sample_id, donor_id = info$donor_id, endpoint_id = id,
               n_group1 = nrow(g1), n_group2 = nrow(g2),
               slide_effect = if (ok) median(g1[[variable]], na.rm = TRUE) - median(g2[[variable]], na.rm = TRUE) else NA_real_,
               minimum_group_rule_met = ok)
  }
  epi_top <- spot[epithelial_rich & ECM_top_quartile]
  epi_bottom <- spot[epithelial_rich & ECM_bottom_quartile]
  slide_outputs[[length(slide_outputs) + 1L]] <- add_simple(
    "P1_epi_high_vs_low_ECM_same_T", epi_top, epi_bottom, "MARKER_T_cell")
  slide_outputs[[length(slide_outputs) + 1L]] <- add_simple(
    "P2_epi_high_vs_low_ECM_neighbor_T", epi_top, epi_bottom, "lag_MARKER_T_cell")

  stromal_high <- spot[stromal_rich & high_ECM_context]
  epi_high <- spot[epithelial_rich & high_ECM_context]
  stromal_low <- spot[stromal_rich & low_ECM_context]
  epi_low <- spot[epithelial_rich & low_ECM_context]
  add_did <- function(id, variable) {
    ns <- c(nrow(stromal_high), nrow(epi_high), nrow(stromal_low), nrow(epi_low))
    ok <- all(ns >= min_group)
    high_gap <- median(stromal_high[[variable]], na.rm = TRUE) - median(epi_high[[variable]], na.rm = TRUE)
    low_gap <- median(stromal_low[[variable]], na.rm = TRUE) - median(epi_low[[variable]], na.rm = TRUE)
    data.table(sample_id = info$sample_id, donor_id = info$donor_id, endpoint_id = id,
               n_stromal_high = ns[[1]], n_epithelial_high = ns[[2]],
               n_stromal_low = ns[[3]], n_epithelial_low = ns[[4]],
               slide_effect = if (ok) high_gap - low_gap else NA_real_,
               minimum_group_rule_met = ok)
  }
  slide_outputs[[length(slide_outputs) + 1L]] <- add_did(
    "P3_ECM_specific_T_confinement_DiD", "MARKER_T_cell")
  slide_outputs[[length(slide_outputs) + 1L]] <- add_did(
    "P4_ECM_specific_cytotoxic_confinement_DiD", "PROGRAM_T_CELL_CYTOTOXIC")

  qc_outputs[[length(qc_outputs) + 1L]] <- data.table(
    sample_id = info$sample_id, donor_id = info$donor_id,
    primary_cohort_donor = info$donor_id %chin% eligible_donors,
    in_tissue_spots = length(common), QC_passing_spots = sum(keep), QC_fraction = mean(keep),
    min_counts = thresholds$min_sum, min_detected = thresholds$min_detected,
    max_mito = thresholds$max_mito)
  coverage_one <- gene_sets_dt[, .(
    requested_genes = .N,
    genes_detected_after_QC = sum(vapply(gene_symbol, function(g) any(rel_counts[g, ] > 0), logical(1)))
  ), by = score]
  coverage_one[, `:=`(sample_id = info$sample_id, donor_id = info$donor_id)]
  coverage_outputs[[length(coverage_outputs) + 1L]] <- coverage_one
  spot_outputs[[length(spot_outputs) + 1L]] <- spot
  rm(h, counts, rel_counts, norm, z, scores, spot); gc(verbose = FALSE)
}

spots <- rbindlist(spot_outputs, fill = TRUE)
qc <- rbindlist(qc_outputs)
coverage <- rbindlist(coverage_outputs)
slide <- rbindlist(slide_outputs, fill = TRUE)
slide_primary_donors <- slide[donor_id %chin% eligible_donors]
patient_st_r <- slide_primary_donors[is.finite(slide_effect), .(
  n_evaluable_ST_R_slides = .N,
  ST_R_patient_effect = median(slide_effect),
  min_ST_R_slide_effect = min(slide_effect), max_ST_R_slide_effect = max(slide_effect)
), by = .(donor_id, endpoint_id)]

bootstrap_median <- function(x, B = 10000L, seed = 1L) {
  if (length(x) < 2L) return(c(NA_real_, NA_real_))
  set.seed(seed); b <- replicate(B, median(sample(x, length(x), replace = TRUE)))
  quantile(b, c(.025, .975), names = FALSE)
}
meta_st_r <- patient_st_r[, {
  x <- ST_R_patient_effect
  ci <- bootstrap_median(x, 10000L, 270679L + .GRP)
  wt <- suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE, correct = FALSE))
  .(n_donors = .N, median_patient_effect = median(x), mean_patient_effect = mean(x),
    ci_low = ci[[1]], ci_high = ci[[2]], p_value = wt$p.value,
    positive_donors = sum(x > 0), negative_donors = sum(x < 0), zero_donors = sum(x == 0))
}, by = endpoint_id]
meta_st_r <- merge(meta_st_r, endpoints[, .(endpoint_id, supportive_direction)], by = "endpoint_id", all.x = TRUE)
meta_st_r[, p_fdr_BH_four_primary := p.adjust(p_value, method = "BH")]

main_patient <- fread(main_patient_file)[, .(donor_id, endpoint_id, main_patient_effect = patient_effect)]
concordance_patient <- merge(main_patient, patient_st_r, by = c("donor_id", "endpoint_id"), all = FALSE)
concordance_meta <- concordance_patient[, {
  ok <- is.finite(main_patient_effect) & is.finite(ST_R_patient_effect)
  x <- main_patient_effect[ok]; y <- ST_R_patient_effect[ok]
  nonzero_both <- x != 0 & y != 0
  .(n_overlap_donors = length(x),
    spearman_rho = if (length(x) >= 5L) suppressWarnings(cor(x, y, method = "spearman")) else NA_real_,
    same_strict_direction = sum(sign(x[nonzero_both]) == sign(y[nonzero_both])),
    nonzero_in_both = sum(nonzero_both),
    either_zero = sum(!nonzero_both),
    median_main = median(x), median_ST_R = median(y),
    median_ST_R_minus_main = median(y - x))
}, by = endpoint_id]

fwrite(spots, file.path(processed_dir, "GSE270678_ST_R_project_compartment_spot_scores_v4.tsv.gz"), sep = "\t", quote = FALSE)
fwrite(qc, file.path(result_dir, "GSE270678_ST_R_slide_QC_v4.tsv"), sep = "\t", quote = FALSE)
fwrite(coverage, file.path(result_dir, "GSE270678_ST_R_gene_set_coverage_v4.tsv"), sep = "\t", quote = FALSE)
fwrite(slide, file.path(result_dir, "GSE270678_ST_R_primary_endpoints_slide_level_v4.tsv"), sep = "\t", quote = FALSE)
fwrite(patient_st_r, file.path(result_dir, "GSE270678_ST_R_primary_endpoints_patient_level_v4.tsv"), sep = "\t", quote = FALSE)
fwrite(meta_st_r, file.path(result_dir, "GSE270678_ST_R_primary_endpoints_patient_equal_meta_v4.tsv"), sep = "\t", quote = FALSE)
fwrite(concordance_patient, file.path(result_dir, "GSE270678_main_vs_ST_R_patient_concordance_v4.tsv"), sep = "\t", quote = FALSE)
fwrite(concordance_meta, file.path(result_dir, "GSE270678_main_vs_ST_R_concordance_summary_v4.tsv"), sep = "\t", quote = FALSE)

audit_lines <- c(
  "# GSE270678 ST_R replicate/extension sensitivity audit",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("- ST_R slides processed: ", uniqueN(qc$sample_id), "."),
  paste0("- Donors represented in ST_R: ", uniqueN(qc$donor_id), "."),
  paste0("- Main-cohort donors represented in ST_R: ", uniqueN(qc[primary_cohort_donor == TRUE, donor_id]), "."),
  "- Extension-only donor D30 is retained in QC outputs but excluded from patient-level sensitivity inference.",
  "- The same pre-expression lock, QC formula, gene sets, scoring, compartment definitions, six-neighbour geometry, 20-spot rule, four endpoints, and BH family were reused without change.",
  "- ST_R slides are not additional independent patients. Slide effects were aggregated within donor before inference.",
  "- Main-versus-ST_R concordance is descriptive because ST_R contains both technical/section replicates and extension locations.",
  "- No new spatial axis was searched."
)
writeLines(audit_lines, file.path(audit_dir, "GSE270678_ST_R_sensitivity_audit_v4.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "GSE270678_ST_R_sensitivity_sessionInfo_v4.txt"), useBytes = TRUE)

message("Completed GSE270678 ST_R sensitivity.")
print(meta_st_r)
print(concordance_meta)
