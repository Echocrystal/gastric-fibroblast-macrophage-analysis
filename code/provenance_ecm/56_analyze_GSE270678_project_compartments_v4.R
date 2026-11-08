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

set.seed(270678L)

model_dir <- file.path(project_dir, "models", "spatial", "v4")
lock_file <- file.path(model_dir, "GSE270678_project_compartment_analysis_lock_v4.rds")
lock_hash_file <- file.path(model_dir, "GSE270678_project_compartment_lock_hashes_v4.tsv")
sample_manifest_file <- file.path(project_dir, "manifests", "v4", "GSE270678_spatial_sample_manifest.tsv")
patient_manifest_file <- file.path(project_dir, "manifests", "v4", "GSE270678_spatial_patient_manifest.tsv")
verification_file <- file.path(project_dir, "manifests", "v4", "GSE270678_main_analysis_file_verification.tsv")
input_dir <- file.path(project_dir, "data", "raw", "spatial", "GSE270678", "main_public_files")
required <- c(lock_file, lock_hash_file, sample_manifest_file, patient_manifest_file, verification_file)
if (!all(file.exists(required))) {
  stop("Missing locked model/manifests or completed direct-download verification: ",
       paste(required[!file.exists(required)], collapse = ", "))
}

result_dir <- file.path(project_dir, "results", "spatial", "v4", "GSE270678_project_compartments")
processed_dir <- file.path(project_dir, "data", "processed", "spatial", "GSE270678", "v4")
audit_dir <- file.path(project_dir, "audit", "v4")
for (d in c(result_dir, processed_dir, audit_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Verify every pre-expression-lock artefact before loading expression data.
lock_hashes <- fread(lock_hash_file)
lock_inputs <- lock_hashes[grepl("GSE270678_project_compartment_(gene_sets|endpoints|analysis_lock)_v4|tcga_ecm_centroid|GSE183904_coarse|scrna_spatial_programs|GSE270678_spatial_sample_manifest", artifact)]
for (i in seq_len(nrow(lock_inputs))) {
  observed <- digest(lock_inputs$artifact[[i]], algo = "sha256", file = TRUE)
  if (!identical(tolower(observed), tolower(lock_inputs$sha256[[i]]))) {
    stop("Pre-expression lock hash mismatch: ", lock_inputs$artifact[[i]])
  }
}

lock <- readRDS(lock_file)
if (!identical(lock$outcome_or_expression_loaded, FALSE)) stop("Lock is not pre-expression")
if (!identical(lock$version, "V4.1_pre_expression_lock")) stop("Unexpected lock version")
gene_sets_dt <- as.data.table(lock$gene_sets)
endpoints <- as.data.table(lock$endpoints)
gene_sets <- split(gene_sets_dt$gene_symbol, gene_sets_dt$score)

samples <- fread(sample_manifest_file)[platform_role == "ST"]
patients <- fread(patient_manifest_file)[main_cohort_donor == TRUE]
verification <- fread(verification_file)
if (nrow(samples) != 33L || uniqueN(samples$donor_id) != 19L) stop("Unexpected main spatial cohort")
if (nrow(verification) != 99L || !all(verification$size_verified == TRUE)) {
  stop("The 99 main-analysis files have not all passed exact-size verification")
}

robust_qc <- function(sum_counts, detected, mito) {
  ls <- log10(sum_counts + 1)
  ld <- log10(detected + 1)
  ms <- median(ls); mads <- mad(ls, constant = 1.4826)
  md <- median(ld); madd <- mad(ld, constant = 1.4826)
  mm <- median(mito); madm <- mad(mito, constant = 1.4826)
  if (!is.finite(mads) || mads == 0) mads <- IQR(ls) / 1.349
  if (!is.finite(madd) || madd == 0) madd <- IQR(ld) / 1.349
  if (!is.finite(madm) || madm == 0) madm <- IQR(mito) / 1.349
  if (!is.finite(mads)) mads <- 0
  if (!is.finite(madd)) madd <- 0
  if (!is.finite(madm)) madm <- 0
  list(
    min_sum = max(0, 10^(ms - 4 * mads) - 1),
    min_detected = max(0, 10^(md - 4 * madd) - 1),
    max_mito = min(100, mm + 5 * madm)
  )
}

read_10x_h5 <- function(path) {
  data <- h5read(path, "/matrix/data")
  indices <- h5read(path, "/matrix/indices")
  indptr <- h5read(path, "/matrix/indptr")
  shape <- as.integer(h5read(path, "/matrix/shape"))
  barcodes <- as.character(h5read(path, "/matrix/barcodes"))
  genes <- as.character(h5read(path, "/matrix/features/name"))
  feature_type <- tryCatch(
    as.character(h5read(path, "/matrix/features/feature_type")),
    error = function(e) rep("Gene Expression", length(genes))
  )
  x <- sparseMatrix(
    i = as.integer(indices) + 1L,
    p = as.integer(indptr),
    x = as.numeric(data),
    dims = shape,
    giveCsparse = TRUE
  )
  keep <- feature_type == "Gene Expression"
  x <- x[keep, , drop = FALSE]
  list(counts = x, genes = genes[keep], barcodes = barcodes)
}

read_positions <- function(path) {
  z <- fread(path, header = FALSE)
  if (tolower(as.character(z[[1L]][1L])) == "barcode") {
    z <- fread(path, header = TRUE)
    setnames(z, names(z)[1:6], c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  } else {
    if (ncol(z) < 6L) stop("Unexpected Visium positions file: ", path)
    setnames(z, names(z)[1:6], c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  }
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
  n_neighbors <- integer(length(values))
  for (i in seq_along(values)) {
    j <- idx[i, ]
    j <- j[is.finite(j)]
    n_neighbors[[i]] <- length(j)
    if (length(j)) out[[i]] <- mean(values[j], na.rm = TRUE)
  }
  list(value = out, n_neighbors = n_neighbors)
}

rank_partial <- function(x, y, covariates = NULL) {
  ok <- is.finite(x) & is.finite(y)
  if (!is.null(covariates)) {
    covariates <- as.matrix(covariates)
    if (nrow(covariates) != length(x)) {
      stop("Covariate row count does not match x/y length: ", nrow(covariates), " versus ", length(x))
    }
    ok <- ok & complete.cases(covariates)
  }
  if (sum(ok) < 50L) return(c(rho = NA_real_, n = sum(ok)))
  xr <- rank(x[ok], ties.method = "average")
  yr <- rank(y[ok], ties.method = "average")
  if (!is.null(covariates)) {
    cv <- covariates[ok, , drop = FALSE]
    cv <- apply(cv, 2L, rank, ties.method = "average")
    if (is.null(dim(cv))) cv <- matrix(cv, ncol = 1L)
    design <- cbind(1, cv)
    xr <- lm.fit(design, xr)$residuals
    yr <- lm.fit(design, yr)$residuals
  }
  c(rho = suppressWarnings(cor(xr, yr, method = "pearson")), n = sum(ok))
}

score_outputs <- list()
qc_outputs <- list()
coverage_outputs <- list()
primary_slide_outputs <- list()
secondary_slide_outputs <- list()

for (ii in seq_len(nrow(samples))) {
  info <- samples[ii]
  message(sprintf("[%02d/%02d] %s (%s)", ii, nrow(samples), info$sample_id, info$donor_id))
  h5_file <- file.path(input_dir, info$filtered_h5)
  pos_file <- file.path(input_dir, info$positions)
  verify_names <- c(info$filtered_h5, info$positions, info$scalefactors)
  verify_rows <- verification[Name %chin% verify_names]
  if (nrow(verify_rows) != 3L || !all(file.exists(verify_rows$destination)) ||
      !all(as.numeric(file.info(verify_rows$destination)$size) == verify_rows$expected_size_bytes)) {
    stop("Input-size verification failed for ", info$sample_id)
  }

  h <- read_10x_h5(h5_file)
  pos <- read_positions(pos_file)[in_tissue == 1L]
  common <- intersect(h$barcodes, pos$barcode)
  if (length(common) < 100L) stop("Fewer than 100 in-tissue barcodes for ", info$sample_id)
  j <- match(common, h$barcodes)
  counts <- h$counts[, j, drop = FALSE]
  pos <- pos[match(common, barcode)]
  if (!all(pos$barcode == common)) stop("Barcode/position order mismatch")

  sum_counts <- as.numeric(colSums(counts))
  detected <- as.numeric(colSums(counts > 0))
  mito_rows <- grepl("^MT-", h$genes, ignore.case = FALSE)
  mito_counts <- if (any(mito_rows)) as.numeric(colSums(counts[mito_rows, , drop = FALSE])) else rep(0, ncol(counts))
  mito <- 100 * mito_counts / pmax(sum_counts, 1)
  thresholds <- robust_qc(sum_counts, detected, mito)
  keep <- sum_counts >= thresholds$min_sum & detected >= thresholds$min_detected & mito <= thresholds$max_mito
  if (sum(keep) < 100L) stop("Fewer than 100 QC-passing spots for ", info$sample_id)
  counts <- counts[, keep, drop = FALSE]
  pos <- pos[keep]
  sum_counts <- sum_counts[keep]
  detected <- detected[keep]
  mito <- mito[keep]

  requested <- unique(gene_sets_dt$gene_symbol)
  rel_counts <- vapply(requested, function(g) {
    rows <- which(h$genes == g)
    if (!length(rows)) return(rep(0, ncol(counts)))
    as.numeric(colSums(counts[rows, , drop = FALSE]))
  }, numeric(ncol(counts)))
  rel_counts <- t(rel_counts)
  rownames(rel_counts) <- requested
  norm <- log1p(sweep(rel_counts, 2L, 1e4 / pmax(sum_counts, 1), "*"))
  z <- t(scale(t(norm), center = TRUE, scale = TRUE))
  gene_variable <- apply(norm, 1L, sd) > 0
  z[!is.finite(z)] <- 0

  score_matrix <- vapply(names(gene_sets), function(nm) {
    colMeans(z[gene_sets[[nm]], , drop = FALSE])
  }, numeric(ncol(z)))
  score_matrix <- as.matrix(score_matrix)
  spot <- cbind(
    pos,
    data.table(
      sample_id = info$sample_id,
      donor_id = info$donor_id,
      stage = info$stage_geo,
      lauren = info$lauren_geo,
      sampling_location = info$sampling_location_geo,
      total_counts = sum_counts,
      detected_features = detected,
      percent_mito = mito
    ),
    as.data.table(score_matrix)
  )
  spot[, frozen_ECM_direction := FROZEN_ECM_C2_HIGH50 - FROZEN_ECM_C1_HIGH50]
  spot[, epithelial_dominance := MARKER_Epithelial - MARKER_Fibroblast]
  spot[, stromal_dominance := MARKER_Fibroblast - MARKER_Epithelial]
  spot[, epithelial_rich := epithelial_dominance >= quantile(epithelial_dominance, 0.75, na.rm = TRUE)]
  spot[, stromal_rich := stromal_dominance >= quantile(stromal_dominance, 0.75, na.rm = TRUE)]
  q_ecm <- quantile(spot$frozen_ECM_direction, c(0.25, 0.5, 0.75), na.rm = TRUE)
  spot[, ECM_bottom_quartile := frozen_ECM_direction <= q_ecm[[1L]]]
  spot[, high_ECM_context := frozen_ECM_direction >= q_ecm[[2L]]]
  spot[, low_ECM_context := frozen_ECM_direction < q_ecm[[2L]]]
  spot[, ECM_top_quartile := frozen_ECM_direction >= q_ecm[[3L]]]
  lag_t <- spatial_lag(spot$MARKER_T_cell, spot$array_row, spot$array_col)
  lag_c <- spatial_lag(spot$PROGRAM_T_CELL_CYTOTOXIC, spot$array_row, spot$array_col)
  lag_x <- spatial_lag(spot$PROGRAM_T_CELL_EXHAUSTION, spot$array_row, spot$array_col)
  spot[, `:=`(
    lag_MARKER_T_cell = lag_t$value,
    lag_PROGRAM_T_CELL_CYTOTOXIC = lag_c$value,
    lag_PROGRAM_T_CELL_EXHAUSTION = lag_x$value,
    n_first_ring_neighbors = lag_t$n_neighbors
  )]

  qc_outputs[[length(qc_outputs) + 1L]] <- data.table(
    sample_id = info$sample_id, donor_id = info$donor_id,
    in_tissue_intersection_spots = length(common), QC_passing_spots = sum(keep),
    QC_fraction = mean(keep), min_counts = thresholds$min_sum,
    min_detected_features = thresholds$min_detected, max_percent_mito = thresholds$max_mito,
    median_counts_after_QC = median(sum_counts), median_features_after_QC = median(detected),
    median_percent_mito_after_QC = median(mito)
  )
  coverage_one <- gene_sets_dt[, .(
    requested_genes = .N,
    genes_present_in_H5 = sum(gene_symbol %chin% h$genes),
    genes_detected_after_QC = sum(vapply(gene_symbol, function(g) any(rel_counts[g, ] > 0), logical(1))),
    genes_variable_after_normalisation = sum(gene_variable[gene_symbol])
  ), by = score]
  coverage_one[, `:=`(sample_id = info$sample_id, donor_id = info$donor_id)]
  setcolorder(coverage_one, c("sample_id", "donor_id", "score"))
  coverage_outputs[[length(coverage_outputs) + 1L]] <- coverage_one

  min_group <- as.integer(lock$definitions$minimum_contrast_group_spots)
  get_primary <- function(endpoint_id, group1, group2, variable) {
    ok <- nrow(group1) >= min_group && nrow(group2) >= min_group
    data.table(
      sample_id = info$sample_id, donor_id = info$donor_id,
      endpoint_id = endpoint_id, n_group1 = nrow(group1), n_group2 = nrow(group2),
      median_group1 = if (nrow(group1)) median(group1[[variable]], na.rm = TRUE) else NA_real_,
      median_group2 = if (nrow(group2)) median(group2[[variable]], na.rm = TRUE) else NA_real_,
      slide_effect = if (ok) median(group1[[variable]], na.rm = TRUE) - median(group2[[variable]], na.rm = TRUE) else NA_real_,
      minimum_group_rule_met = ok
    )
  }

  epi_top <- spot[epithelial_rich == TRUE & ECM_top_quartile == TRUE]
  epi_bottom <- spot[epithelial_rich == TRUE & ECM_bottom_quartile == TRUE]
  primary_slide_outputs[[length(primary_slide_outputs) + 1L]] <- get_primary(
    "P1_epi_high_vs_low_ECM_same_T", epi_top, epi_bottom, "MARKER_T_cell")
  primary_slide_outputs[[length(primary_slide_outputs) + 1L]] <- get_primary(
    "P2_epi_high_vs_low_ECM_neighbor_T", epi_top, epi_bottom, "lag_MARKER_T_cell")

  stromal_high <- spot[stromal_rich == TRUE & high_ECM_context == TRUE]
  epi_high <- spot[epithelial_rich == TRUE & high_ECM_context == TRUE]
  stromal_low <- spot[stromal_rich == TRUE & low_ECM_context == TRUE]
  epi_low <- spot[epithelial_rich == TRUE & low_ECM_context == TRUE]
  get_did <- function(endpoint_id, variable) {
    ns <- c(nrow(stromal_high), nrow(epi_high), nrow(stromal_low), nrow(epi_low))
    ok <- all(ns >= min_group)
    high_gap <- if (all(ns[1:2] > 0)) median(stromal_high[[variable]], na.rm = TRUE) - median(epi_high[[variable]], na.rm = TRUE) else NA_real_
    low_gap <- if (all(ns[3:4] > 0)) median(stromal_low[[variable]], na.rm = TRUE) - median(epi_low[[variable]], na.rm = TRUE) else NA_real_
    data.table(
      sample_id = info$sample_id, donor_id = info$donor_id, endpoint_id = endpoint_id,
      n_stromal_high = ns[[1]], n_epithelial_high = ns[[2]],
      n_stromal_low = ns[[3]], n_epithelial_low = ns[[4]],
      high_ECM_stromal_minus_epithelial_gap = high_gap,
      low_ECM_stromal_minus_epithelial_gap = low_gap,
      slide_effect = if (ok) high_gap - low_gap else NA_real_,
      minimum_group_rule_met = ok
    )
  }
  primary_slide_outputs[[length(primary_slide_outputs) + 1L]] <- get_did(
    "P3_ECM_specific_T_confinement_DiD", "MARKER_T_cell")
  primary_slide_outputs[[length(primary_slide_outputs) + 1L]] <- get_did(
    "P4_ECM_specific_cytotoxic_confinement_DiD", "PROGRAM_T_CELL_CYTOTOXIC")

  specs <- list(
    list(id = "S1_epi_ECM_same_T_rho", x = "frozen_ECM_direction", y = "MARKER_T_cell",
         subset = spot$epithelial_rich, cov = NULL),
    list(id = "S2_epi_ECM_neighbor_T_rho", x = "frozen_ECM_direction", y = "lag_MARKER_T_cell",
         subset = spot$epithelial_rich, cov = NULL),
    list(id = "S3_allspot_ECM_T_partial_rho", x = "frozen_ECM_direction", y = "MARKER_T_cell",
         subset = rep(TRUE, nrow(spot)), cov = c("MARKER_Epithelial", "MARKER_Fibroblast", "MARKER_Endothelial")),
    list(id = "S4_allspot_ECM_neighbor_T_partial_rho", x = "frozen_ECM_direction", y = "lag_MARKER_T_cell",
         subset = rep(TRUE, nrow(spot)), cov = c("MARKER_Epithelial", "MARKER_Fibroblast", "MARKER_Endothelial")),
    list(id = "S5_fibroblast_T_partial_rho", x = "MARKER_Fibroblast", y = "MARKER_T_cell",
         subset = rep(TRUE, nrow(spot)), cov = c("frozen_ECM_direction", "MARKER_Epithelial")),
    list(id = "S6_fibroblast_neighbor_T_partial_rho", x = "MARKER_Fibroblast", y = "lag_MARKER_T_cell",
         subset = rep(TRUE, nrow(spot)), cov = c("frozen_ECM_direction", "MARKER_Epithelial")),
    list(id = "S7_endothelial_T_partial_rho", x = "MARKER_Endothelial", y = "MARKER_T_cell",
         subset = rep(TRUE, nrow(spot)), cov = c("frozen_ECM_direction", "MARKER_Epithelial", "MARKER_Fibroblast")),
    list(id = "S8_endothelial_neighbor_T_partial_rho", x = "MARKER_Endothelial", y = "lag_MARKER_T_cell",
         subset = rep(TRUE, nrow(spot)), cov = c("frozen_ECM_direction", "MARKER_Epithelial", "MARKER_Fibroblast"))
  )
  for (s in specs) {
    q <- spot[s$subset == TRUE]
    cv <- if (is.null(s$cov)) NULL else as.matrix(q[, s$cov, with = FALSE])
    tt <- rank_partial(q[[s$x]], q[[s$y]], cv)
    secondary_slide_outputs[[length(secondary_slide_outputs) + 1L]] <- data.table(
      sample_id = info$sample_id, donor_id = info$donor_id,
      endpoint_id = s$id, n_spots = as.integer(tt[["n"]]), slide_rho = tt[["rho"]]
    )
  }

  score_outputs[[length(score_outputs) + 1L]] <- spot
  rm(h, counts, rel_counts, norm, z, score_matrix, spot)
  gc(verbose = FALSE)
}

spots <- rbindlist(score_outputs, fill = TRUE)
qc <- rbindlist(qc_outputs, fill = TRUE)
coverage <- rbindlist(coverage_outputs, fill = TRUE)
primary_slide <- rbindlist(primary_slide_outputs, fill = TRUE)
secondary_slide <- rbindlist(secondary_slide_outputs, fill = TRUE)

primary_patient <- primary_slide[is.finite(slide_effect), .(
  n_evaluable_slides = .N,
  patient_effect = median(slide_effect),
  min_slide_effect = min(slide_effect),
  max_slide_effect = max(slide_effect)
), by = .(donor_id, endpoint_id)]
primary_patient <- merge(primary_patient, patients[, .(
  donor_id, stage_geo, lauren_geo, n_main_ST, main_locations
)], by = "donor_id", all.x = TRUE)

bootstrap_median <- function(x, B = 10000L, seed = 1L) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(low = NA_real_, high = NA_real_))
  set.seed(seed)
  b <- replicate(B, median(sample(x, length(x), replace = TRUE)))
  quantile(b, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

primary_meta <- primary_patient[, {
  x <- patient_effect[is.finite(patient_effect)]
  ci <- bootstrap_median(x, 10000L, 270678L + .GRP)
  wt <- if (length(x) >= 2L) suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE, correct = FALSE)) else NULL
  .(n_donors = length(x), median_patient_effect = median(x), mean_patient_effect = mean(x),
    ci_low = ci[[1]], ci_high = ci[[2]], p_value = if (is.null(wt)) NA_real_ else wt$p.value,
    positive_donors = sum(x > 0), negative_donors = sum(x < 0), zero_donors = sum(x == 0))
}, by = endpoint_id]
primary_meta <- merge(primary_meta, endpoints[tier == "primary", .(
  endpoint_id, estimand, supportive_direction
)], by = "endpoint_id", all.x = TRUE)
primary_meta[, p_fdr_BH_four_primary := p.adjust(p_value, method = "BH")]
primary_meta[, supportive_direction_met := fifelse(
  supportive_direction == "negative", median_patient_effect < 0,
  fifelse(supportive_direction == "positive", median_patient_effect > 0, NA)
)]

secondary_patient <- secondary_slide[is.finite(slide_rho), .(
  n_evaluable_slides = .N,
  patient_rho = tanh(mean(atanh(pmax(pmin(slide_rho, 0.999999), -0.999999)))),
  median_slide_rho = median(slide_rho)
), by = .(donor_id, endpoint_id)]

bootstrap_fisher_rho <- function(x, B = 10000L, seed = 1L) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(low = NA_real_, high = NA_real_))
  z <- atanh(pmax(pmin(x, 0.999999), -0.999999))
  set.seed(seed)
  b <- replicate(B, tanh(mean(sample(z, length(z), replace = TRUE))))
  quantile(b, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

secondary_meta <- secondary_patient[, {
  x <- patient_rho[is.finite(patient_rho)]
  ci <- bootstrap_fisher_rho(x, 10000L, 270778L + .GRP)
  wt <- if (length(x) >= 2L) suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE, correct = FALSE)) else NULL
  .(n_donors = length(x), pooled_patient_equal_rho = tanh(mean(atanh(pmax(pmin(x, 0.999999), -0.999999)))),
    median_patient_rho = median(x), ci_low = ci[[1]], ci_high = ci[[2]],
    p_value = if (is.null(wt)) NA_real_ else wt$p.value,
    positive_donors = sum(x > 0), negative_donors = sum(x < 0), zero_donors = sum(x == 0))
}, by = endpoint_id]
secondary_meta <- merge(secondary_meta, endpoints[tier == "secondary", .(
  endpoint_id, estimand, supportive_direction
)], by = "endpoint_id", all.x = TRUE)
secondary_meta[, p_fdr_BH_eight_secondary := p.adjust(p_value, method = "BH")]

primary_patient[, stage_group := fifelse(stage_geo %chin% c("Stage I", "Stage II"), "Stage I-II", "Stage III")]
strata_long <- rbindlist(list(
  primary_patient[, .(stratum_family = "Lauren", stratum = lauren_geo,
                      endpoint_id, donor_id, patient_effect)],
  primary_patient[, .(stratum_family = "Stage", stratum = stage_group,
                      endpoint_id, donor_id, patient_effect)]
))
primary_strata <- strata_long[, .(
  n_donors = uniqueN(donor_id),
  median_patient_effect = median(patient_effect),
  positive_donors = sum(patient_effect > 0),
  negative_donors = sum(patient_effect < 0),
  minimum_five_donors_rule_met = uniqueN(donor_id) >= 5L
), by = .(stratum_family, stratum, endpoint_id)]
primary_strata[n_donors < 5L, `:=`(
  median_patient_effect = NA_real_, positive_donors = NA_integer_, negative_donors = NA_integer_
)]

fwrite(spots, file.path(processed_dir, "GSE270678_main_ST_project_compartment_spot_scores_v4.tsv.gz"),
       sep = "\t", quote = FALSE)
fwrite(qc, file.path(result_dir, "GSE270678_main_ST_slide_QC_v4.tsv"), sep = "\t", quote = FALSE)
fwrite(coverage, file.path(result_dir, "GSE270678_main_ST_gene_set_coverage_v4.tsv"), sep = "\t", quote = FALSE)
fwrite(primary_slide, file.path(result_dir, "GSE270678_primary_endpoints_slide_level_v4.tsv"), sep = "\t", quote = FALSE)
fwrite(primary_patient, file.path(result_dir, "GSE270678_primary_endpoints_patient_level_v4.tsv"), sep = "\t", quote = FALSE)
fwrite(primary_meta, file.path(result_dir, "GSE270678_primary_endpoints_patient_equal_meta_v4.tsv"), sep = "\t", quote = FALSE)
fwrite(primary_strata, file.path(result_dir, "GSE270678_primary_endpoints_stage_Lauren_sensitivity_v4.tsv"), sep = "\t", quote = FALSE)
fwrite(secondary_slide, file.path(result_dir, "GSE270678_secondary_axes_slide_level_v4.tsv"), sep = "\t", quote = FALSE)
fwrite(secondary_patient, file.path(result_dir, "GSE270678_secondary_axes_patient_level_v4.tsv"), sep = "\t", quote = FALSE)
fwrite(secondary_meta, file.path(result_dir, "GSE270678_secondary_axes_patient_equal_meta_v4.tsv"), sep = "\t", quote = FALSE)

saveRDS(list(
  lock_sha256 = digest(lock_file, algo = "sha256", file = TRUE),
  primary_meta = primary_meta,
  secondary_meta = secondary_meta,
  qc = qc,
  gene_set_coverage = coverage,
  author_region_status = "NOT EVALUABLE",
  labels = "project-reconstructed expression compartments"
), file.path(model_dir, "GSE270678_project_compartment_results_v4.rds"), compress = "xz")

audit_lines <- c(
  "# GSE270678 project-reconstructed compartment analysis run audit",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("- Main slides analysed: ", uniqueN(qc$sample_id), "."),
  paste0("- Main donors: ", uniqueN(samples$donor_id), "."),
  paste0("- QC-passing spots: ", sum(qc$QC_passing_spots), "."),
  "- All pre-expression-lock SHA-256 values were verified before expression matrices were read.",
  "- The source-study LAR/SMSR/ICR/NCR barcode labels remain NOT EVALUABLE and were not imputed.",
  "- Source-study region DE genes were not loaded or used for classification.",
  "- Compartments are project-reconstructed expression compartments based on previously frozen markers.",
  "- Slide-level results were aggregated within donor before inference; spots and slides were not counted as independent patients.",
  "- Four primary endpoints use one BH family; eight secondary correlation axes use a separate BH family.",
  "- Stage/Lauren sensitivity summaries are descriptive and suppressed for strata with fewer than five donors."
)
writeLines(audit_lines, file.path(audit_dir, "GSE270678_project_compartment_analysis_run_v4.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "GSE270678_project_compartment_analysis_sessionInfo_v4.txt"), useBytes = TRUE)

message("Completed GSE270678 project-compartment analysis.")
print(primary_meta)
