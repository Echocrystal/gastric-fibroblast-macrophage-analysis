#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else
  "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- if (length(args) >= 2L) normalizePath(args[[2]], winslash = "/", mustWork = TRUE) else
  "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(digest)
  library(jsonlite)
})

set.seed(25195005L)

raw_dir <- file.path(source_project, "data", "raw", "spatial", "GSE251950")
sample_root <- file.path(raw_dir, "samples")
out_dir <- file.path(project_dir, "analysis", "results", "spatial", "GSE251950_candidate_axes")
processed_dir <- file.path(project_dir, "analysis", "processed", "spatial", "GSE251950_candidate_axes")
audit_dir <- file.path(project_dir, "analysis", "audit")
figure_dir <- file.path(project_dir, "figures", "spatial_candidate_axes")
for (d in c(out_dir, processed_dir, audit_dir, figure_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

permutation_B <- 199L
primary_detection_fraction <- 0.05
exploratory_detection_fraction <- 0.01

axis_manifest <- data.table(
  axis_id = c(
    "PLAU_PLAUR", "ANGPTL2_ITGA5_ITGB1", "TNC_ITGA5_ITGB1", "HLAF_LILRB2",
    "CXCL12_CXCR4", "SPP1_CD44", "TGFB1_TGFBR2", "CSF1_CSF1R", "CCL2_CCR2"
  ),
  display_axis = c(
    "PLAU-PLAUR", "ANGPTL2-ITGA5/ITGB1", "TNC-ITGA5/ITGB1", "HLA-F-LILRB2",
    "CXCL12-CXCR4", "SPP1-CD44", "TGFB1-TGFBR2", "CSF1-CSF1R", "CCL2-CCR2"
  ),
  ligand = c("PLAU", "ANGPTL2", "TNC", "HLA-F", "CXCL12", "SPP1", "TGFB1", "CSF1", "CCL2"),
  receptor_components = c("PLAUR", "ITGA5|ITGB1", "ITGA5|ITGB1", "LILRB2",
                          "CXCR4", "CD44", "TGFBR2", "CSF1R", "CCR2"),
  axis_class = c(
    "direct_ECM_proteolysis", "ECM_adhesion_mechanosensing", "core_ECM_adhesion",
    "ECM_associated_immunoregulation", rep("published_benchmark", 5)
  ),
  candidate_type = c(rep("candidate", 4), rep("benchmark", 5)),
  proposed_receiver = c(rep("Myeloid", 4), "mixed_immune", "Myeloid", "mixed", "Myeloid", "Myeloid"),
  single_cell_expected_ecm_sign = c(rep(-1L, 4), rep(NA_integer_, 5))
)

manifest <- fread(file.path(source_project, "data", "processed", "spatial", "GSE251950",
                            "GSE251950_spatial_sample_manifest.tsv"))
if (nrow(manifest) != 10L || sum(manifest$primary_analysis) != 9L || uniqueN(manifest$patient_id) != 9L) {
  stop("Unexpected frozen GSE251950 manifest")
}

frozen <- fread(file.path(source_project, "models", "ecotype",
                          "tcga_ecm_centroid_gene_contrasts.tsv"))
contrast <- frozen[order(-centroid_contrast_C2_minus_C1)]
ecm_high50 <- head(contrast$gene_symbol, 50L)
ecm_low50 <- tail(contrast$gene_symbol, 50L)

program_manifest <- fread(file.path(source_project, "data", "raw", "gene_sets",
                                    "scrna_spatial_programs_locked.tsv"))
marker_manifest <- fread(file.path(source_project, "data", "raw", "gene_sets",
                                   "GSE183904_coarse_annotation_markers_locked.tsv"))
programs <- split(program_manifest$gene_symbol, program_manifest$program)
markers <- split(marker_manifest$gene_symbol, marker_manifest$coarse_label)
required_programs <- c("CAF_ACTIVATED_MATRIX", "MYELOID_INFLAMMATORY",
                       "MYELOID_SPP1_TREM2", "TLS_B_SUPPORT")
required_markers <- c("Fibroblast", "Myeloid", "Epithelial", "Endothelial",
                      "T_cell", "B_cell", "Plasma")
if (!all(required_programs %chin% names(programs)) || !all(required_markers %chin% names(markers))) {
  stop("Frozen spatial program/marker manifest is missing required sets")
}

receptor_genes <- unique(unlist(strsplit(axis_manifest$receptor_components, "|", fixed = TRUE)))
candidate_genes <- unique(c(axis_manifest$ligand, receptor_genes))
score_sets <- c(
  list(FROZEN_C2_HIGH50 = ecm_high50, FROZEN_C1_HIGH50 = ecm_low50),
  setNames(programs[required_programs], paste0("PROGRAM_", required_programs)),
  setNames(markers[required_markers], paste0("MARKER_", required_markers))
)
needed_genes <- unique(c(candidate_genes, unlist(score_sets), programs$CAF_ACTIVATED_MATRIX))

robust_qc <- function(sum_counts, detected, mito) {
  ls <- log10(sum_counts + 1)
  ld <- log10(detected + 1)
  ms <- median(ls); mads <- mad(ls, constant = 1.4826)
  md <- median(ld); madd <- mad(ld, constant = 1.4826)
  mm <- median(mito); madm <- mad(mito, constant = 1.4826)
  if (!is.finite(mads) || mads == 0) mads <- IQR(ls) / 1.349
  if (!is.finite(madd) || madd == 0) madd <- IQR(ld) / 1.349
  if (!is.finite(madm) || madm == 0) madm <- IQR(mito) / 1.349
  list(
    min_sum = max(0, 10^(ms - 4 * mads) - 1),
    min_detected = max(0, 10^(md - 4 * madd) - 1),
    max_mito = min(100, mm + 5 * madm)
  )
}

spatial_lag <- function(values, array_row, array_col) {
  offsets <- rbind(c(-1, -1), c(-1, 1), c(0, -2), c(0, 2), c(1, -1), c(1, 1))
  keys <- paste(array_row, array_col, sep = ":")
  neighbors <- lapply(seq_len(nrow(offsets)), function(i) {
    match(paste(array_row + offsets[i, 1], array_col + offsets[i, 2], sep = ":"), keys)
  })
  idx <- do.call(cbind, neighbors)
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

zscore_rows <- function(mat) {
  z <- t(scale(t(as.matrix(mat)), center = TRUE, scale = TRUE))
  z[!is.finite(z)] <- 0
  z
}

score_gene_sets <- function(norm, sets) {
  z <- zscore_rows(norm[intersect(unique(unlist(sets)), rownames(norm)), , drop = FALSE])
  scores <- vapply(sets, function(gs) {
    present <- intersect(gs, rownames(z))
    if (length(present) < 2L) return(rep(NA_real_, ncol(z)))
    colMeans(z[present, , drop = FALSE])
  }, numeric(ncol(z)))
  if (is.null(dim(scores))) scores <- matrix(scores, ncol = 1L)
  colnames(scores) <- names(sets)
  list(
    scores = scores,
    coverage = data.table(
      score = names(sets),
      available_genes = vapply(sets, function(gs) length(intersect(gs, rownames(norm))), integer(1)),
      requested_genes = lengths(sets)
    )
  )
}

safe_quantile_bin <- function(x, n_bins = 4L) {
  r <- frank(x, ties.method = "average", na.last = "keep")
  out <- ceiling(n_bins * r / max(r, na.rm = TRUE))
  out[!is.finite(out)] <- NA_integer_
  pmax(1L, pmin(n_bins, as.integer(out)))
}

rank_partial_residuals <- function(x, y, covariates = NULL) {
  ok <- is.finite(x) & is.finite(y)
  if (!is.null(covariates)) {
    ok <- ok & apply(covariates, 1L, function(v) all(is.finite(v)))
  }
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 20L || uniqueN(x) < 3L || uniqueN(y) < 3L) {
    return(list(ok = ok, x_resid = numeric(), y_resid = numeric(), rho = NA_real_, p = NA_real_, n = length(x)))
  }
  xr <- rank(x, ties.method = "average")
  yr <- rank(y, ties.method = "average")
  if (!is.null(covariates)) {
    cmat <- covariates[ok, , drop = FALSE]
    design <- cbind(1, apply(cmat, 2L, rank, ties.method = "average"))
    xr <- lm.fit(design, xr)$residuals
    yr <- lm.fit(design, yr)$residuals
  }
  test <- suppressWarnings(cor.test(xr, yr, method = "pearson"))
  list(ok = ok, x_resid = xr, y_resid = yr, rho = unname(test$estimate),
       p = test$p.value, n = length(xr))
}

stratified_permutation_p <- function(fit, strata, B = permutation_B) {
  if (!is.finite(fit$rho) || length(fit$x_resid) < 20L) return(NA_real_)
  strata <- strata[fit$ok]
  strata[is.na(strata)] <- "missing"
  groups <- split(seq_along(strata), strata)
  null <- numeric(B)
  for (b in seq_len(B)) {
    yp <- fit$y_resid
    for (idx in groups) {
      if (length(idx) > 1L) yp[idx] <- sample(yp[idx], length(idx), replace = FALSE)
    }
    null[[b]] <- suppressWarnings(cor(fit$x_resid, yp, method = "pearson"))
  }
  (1 + sum(abs(null) >= abs(fit$rho), na.rm = TRUE)) / (B + 1)
}

complete_receptor_expression <- function(norm, components) {
  genes <- strsplit(components, "|", fixed = TRUE)[[1L]]
  if (!all(genes %chin% rownames(norm))) {
    return(list(expression = rep(NA_real_, ncol(norm)), detected = rep(FALSE, ncol(norm)),
                available = FALSE, available_components = intersect(genes, rownames(norm))))
  }
  x <- as.matrix(norm[genes, , drop = FALSE])
  if (length(genes) == 1L) {
    expression <- as.numeric(x[1L, ])
  } else {
    expression <- apply(x, 2L, min)
  }
  list(expression = expression, detected = apply(x > 0, 2L, all),
       available = TRUE, available_components = genes)
}

patient_equal_weight_summary <- function(rho, n, B = 5000L) {
  ok <- is.finite(rho) & n > 3L
  r <- pmax(-0.999999, pmin(0.999999, rho[ok]))
  k <- length(r)
  if (!k) {
    return(data.table(k = 0L, pooled_rho = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                      p_value = NA_real_, median_slide_rho = NA_real_, positive_slides = 0L,
                      negative_slides = 0L, min_loo_rho = NA_real_, max_loo_rho = NA_real_))
  }
  pooled <- tanh(mean(atanh(r)))
  if (k > 1L) {
    boot <- replicate(B, tanh(mean(atanh(sample(r, k, replace = TRUE)))))
    ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
    w <- suppressWarnings(wilcox.test(r, mu = 0, exact = FALSE, correct = FALSE))
    loo <- vapply(seq_len(k), function(i) tanh(mean(atanh(r[-i]))), numeric(1))
  } else {
    ci <- c(NA_real_, NA_real_)
    w <- list(p.value = NA_real_)
    loo <- NA_real_
  }
  data.table(
    k = k, pooled_rho = pooled, ci_low = as.numeric(ci[[1L]]), ci_high = as.numeric(ci[[2L]]),
    p_value = w$p.value, median_slide_rho = median(r), positive_slides = sum(r > 0),
    negative_slides = sum(r < 0), min_loo_rho = min(loo, na.rm = TRUE),
    max_loo_rho = max(loo, na.rm = TRUE)
  )
}

endpoint_specs <- data.table(
  endpoint = c(
    "ligand_receptor_same_raw", "ligand_receptor_same_adjusted",
    "ligand_receptor_neighbor_raw", "ligand_receptor_neighbor_adjusted",
    "matrix_axis_same_adjusted", "matrix_axis_neighbor_adjusted",
    "ecm_direction_axis_same_adjusted", "ecm_direction_axis_neighbor_adjusted",
    "fibroblast_source_attribution", "activated_CAF_source_attribution",
    "myeloid_receiver_attribution"
  ),
  expected_sign = c(1L, 1L, 1L, 1L, 1L, 1L, rep(NA_integer_, 2), 1L, 1L, 1L)
)

spot_axis_outputs <- list()
slide_stats_outputs <- list()
coverage_outputs <- list()
qc_outputs <- list()

for (i in seq_len(nrow(manifest))) {
  info <- manifest[i]
  message(sprintf("[%d/%d] Candidate-axis spatial analysis: %s", i, nrow(manifest), info$sample_title))
  folder <- file.path(sample_root, info$folder_id)
  prefix <- file.path(folder, info$folder_id)
  features <- fread(paste0(prefix, "_features.tsv.gz"), header = FALSE)
  barcodes <- fread(paste0(prefix, "_barcodes.tsv.gz"), header = FALSE)$V1
  raw_counts <- as(readMM(gzfile(paste0(prefix, "_matrix.mtx.gz"))), "dgCMatrix")
  if (nrow(raw_counts) != nrow(features) || ncol(raw_counts) != length(barcodes)) {
    stop("Matrix dimension mismatch for ", info$folder_id)
  }
  symbols <- features$V2
  valid <- !is.na(symbols) & nzchar(symbols)
  symbols_valid <- symbols[valid]
  unique_symbols <- sort(unique(symbols_valid))
  mapper <- sparseMatrix(i = match(symbols_valid, unique_symbols), j = which(valid), x = 1,
                         dims = c(length(unique_symbols), nrow(raw_counts)))
  counts <- mapper %*% raw_counts
  rownames(counts) <- unique_symbols
  colnames(counts) <- barcodes
  rm(raw_counts, mapper); gc(verbose = FALSE)

  positions <- fread(paste0(prefix, "_tissue_positions_list.csv"), header = FALSE)
  setnames(positions, c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  positions <- positions[match(barcodes, barcode)]
  if (anyNA(positions$barcode)) stop("Missing spatial coordinates for ", info$folder_id)

  total <- Matrix::colSums(counts)
  detected <- Matrix::colSums(counts > 0)
  mt <- grepl("^MT-", rownames(counts))
  mito <- 100 * Matrix::colSums(counts[mt, , drop = FALSE]) / pmax(total, 1)
  thr <- robust_qc(total, detected, mito)
  pass <- positions$in_tissue == 1L & total >= thr$min_sum &
    detected >= thr$min_detected & mito <= thr$max_mito
  if (sum(pass) < 200L) stop("Too few QC spots in ", info$folder_id)
  counts <- counts[, pass, drop = FALSE]
  positions <- positions[pass]
  total <- total[pass]
  detected <- detected[pass]
  mito <- mito[pass]

  present <- intersect(needed_genes, rownames(counts))
  norm <- counts[present, , drop = FALSE] %*% Diagonal(x = 1e4 / pmax(total, 1))
  rownames(norm) <- present
  colnames(norm) <- colnames(counts)
  norm@x <- log1p(norm@x)
  scored <- score_gene_sets(norm, score_sets)
  spot_base <- cbind(
    data.table(
      GSM = info$GSM, sample_title = info$sample_title, patient_id = info$patient_id,
      specimen_type = info$specimen_type, primary_analysis = info$primary_analysis,
      barcode = colnames(counts), total_UMI = total, detected_genes = detected,
      mito_percent = mito
    ),
    positions[, .(array_row, array_col, pxl_row, pxl_col)],
    as.data.table(scored$scores)
  )
  spot_base[, frozen_ECM_direction := FROZEN_C2_HIGH50 - FROZEN_C1_HIGH50]
  lag_myeloid <- spatial_lag(spot_base$MARKER_Myeloid, spot_base$array_row, spot_base$array_col)
  spot_base[, lag_MARKER_Myeloid := lag_myeloid$value]
  spot_base[, n_first_ring_neighbors := lag_myeloid$n_neighbors]

  for (a in seq_len(nrow(axis_manifest))) {
    spec <- axis_manifest[a]
    ligand_available <- spec$ligand %chin% rownames(norm)
    ligand_expr <- if (ligand_available) as.numeric(norm[spec$ligand, ]) else rep(NA_real_, ncol(norm))
    ligand_detected <- if (ligand_available) ligand_expr > 0 else rep(FALSE, ncol(norm))
    receptor <- complete_receptor_expression(norm, spec$receptor_components)
    receptor_lag <- spatial_lag(receptor$expression, spot_base$array_row, spot_base$array_col)$value
    same_axis_score <- sqrt(pmax(ligand_expr, 0) * pmax(receptor$expression, 0))
    neighbor_axis_score <- sqrt(pmax(ligand_expr, 0) * pmax(receptor_lag, 0))

    caf_genes_no_ligand <- setdiff(programs$CAF_ACTIVATED_MATRIX,
                                   c(spec$ligand, strsplit(spec$receptor_components, "|", fixed = TRUE)[[1L]]))
    caf_present <- intersect(caf_genes_no_ligand, rownames(norm))
    caf_z <- zscore_rows(norm[caf_present, , drop = FALSE])
    caf_score_no_ligand <- if (length(caf_present) >= 2L) colMeans(caf_z) else rep(NA_real_, ncol(norm))

    axis_spot <- spot_base[, .(
      GSM, sample_title, patient_id, specimen_type, primary_analysis, barcode,
      array_row, array_col, pxl_row, pxl_col, total_UMI, detected_genes,
      n_first_ring_neighbors, frozen_ECM_direction, MARKER_Fibroblast,
      MARKER_Myeloid, lag_MARKER_Myeloid, PROGRAM_CAF_ACTIVATED_MATRIX,
      PROGRAM_MYELOID_INFLAMMATORY, PROGRAM_MYELOID_SPP1_TREM2, PROGRAM_TLS_B_SUPPORT
    )]
    axis_spot[, `:=`(
      axis_id = spec$axis_id,
      display_axis = spec$display_axis,
      ligand = spec$ligand,
      receptor_components = spec$receptor_components,
      ligand_expression = ligand_expr,
      receptor_complete_expression = receptor$expression,
      neighbor_receptor_complete_expression = receptor_lag,
      same_axis_score = same_axis_score,
      neighbor_axis_score = neighbor_axis_score,
      CAF_activated_matrix_no_axis_genes = caf_score_no_ligand
    )]
    spot_axis_outputs[[length(spot_axis_outputs) + 1L]] <- axis_spot

    ligand_fraction <- mean(ligand_detected, na.rm = TRUE)
    receptor_fraction <- mean(receptor$detected, na.rm = TRUE)
    informative_primary <- ligand_fraction >= primary_detection_fraction &
      receptor_fraction >= primary_detection_fraction & receptor$available
    informative_exploratory <- ligand_fraction >= exploratory_detection_fraction &
      receptor_fraction >= exploratory_detection_fraction & receptor$available

    coverage_outputs[[length(coverage_outputs) + 1L]] <- data.table(
      GSM = info$GSM, sample_title = info$sample_title, patient_id = info$patient_id,
      specimen_type = info$specimen_type, primary_analysis = info$primary_analysis,
      axis_id = spec$axis_id, display_axis = spec$display_axis,
      ligand = spec$ligand, ligand_available = ligand_available,
      ligand_detected_fraction = ligand_fraction,
      receptor_components = spec$receptor_components,
      receptor_components_available = paste(receptor$available_components, collapse = "|"),
      complete_receptor_available = receptor$available,
      complete_receptor_detected_fraction = receptor_fraction,
      informative_primary = informative_primary,
      informative_exploratory = informative_exploratory,
      CAF_program_genes_after_axis_exclusion = length(caf_present)
    )

    log_umi <- log1p(spot_base$total_UMI)
    common_cov <- cbind(log_umi, spot_base$MARKER_Fibroblast, spot_base$MARKER_Myeloid,
                        spot_base$frozen_ECM_direction)
    neighbor_cov <- cbind(log_umi, spot_base$MARKER_Fibroblast,
                          spot_base$lag_MARKER_Myeloid, spot_base$frozen_ECM_direction)
    abundance_cov <- cbind(log_umi, spot_base$MARKER_Fibroblast, spot_base$MARKER_Myeloid)
    neighbor_abundance_cov <- cbind(log_umi, spot_base$MARKER_Fibroblast,
                                    spot_base$lag_MARKER_Myeloid)
    source_cov <- cbind(log_umi, spot_base$MARKER_Myeloid,
                        spot_base$MARKER_Epithelial, spot_base$MARKER_Endothelial)
    receiver_cov <- cbind(log_umi, spot_base$MARKER_Fibroblast,
                          spot_base$MARKER_Epithelial, spot_base$MARKER_Endothelial)
    same_strata <- paste(safe_quantile_bin(spot_base$MARKER_Fibroblast),
                         safe_quantile_bin(spot_base$MARKER_Myeloid), sep = ":")
    neighbor_strata <- paste(safe_quantile_bin(spot_base$MARKER_Fibroblast),
                             safe_quantile_bin(spot_base$lag_MARKER_Myeloid), sep = ":")

    endpoint_inputs <- list(
      ligand_receptor_same_raw = list(x = ligand_expr, y = receptor$expression, cov = NULL,
                                      strata = same_strata, permute = FALSE),
      ligand_receptor_same_adjusted = list(x = ligand_expr, y = receptor$expression,
                                           cov = common_cov, strata = same_strata, permute = TRUE),
      ligand_receptor_neighbor_raw = list(x = ligand_expr, y = receptor_lag, cov = NULL,
                                          strata = neighbor_strata, permute = FALSE),
      ligand_receptor_neighbor_adjusted = list(x = ligand_expr, y = receptor_lag,
                                               cov = neighbor_cov, strata = neighbor_strata, permute = TRUE),
      matrix_axis_same_adjusted = list(x = caf_score_no_ligand, y = same_axis_score,
                                       cov = abundance_cov, strata = same_strata, permute = FALSE),
      matrix_axis_neighbor_adjusted = list(x = caf_score_no_ligand, y = neighbor_axis_score,
                                           cov = neighbor_abundance_cov, strata = neighbor_strata, permute = FALSE),
      ecm_direction_axis_same_adjusted = list(x = spot_base$frozen_ECM_direction, y = same_axis_score,
                                              cov = abundance_cov, strata = same_strata, permute = FALSE),
      ecm_direction_axis_neighbor_adjusted = list(x = spot_base$frozen_ECM_direction,
                                                  y = neighbor_axis_score,
                                                  cov = neighbor_abundance_cov,
                                                  strata = neighbor_strata, permute = FALSE),
      fibroblast_source_attribution = list(x = ligand_expr, y = spot_base$MARKER_Fibroblast,
                                           cov = source_cov, strata = same_strata, permute = FALSE),
      activated_CAF_source_attribution = list(x = ligand_expr, y = caf_score_no_ligand,
                                              cov = source_cov, strata = same_strata, permute = FALSE),
      myeloid_receiver_attribution = list(x = receptor$expression, y = spot_base$MARKER_Myeloid,
                                          cov = receiver_cov, strata = same_strata, permute = FALSE)
    )

    for (endpoint in names(endpoint_inputs)) {
      ep <- endpoint_inputs[[endpoint]]
      fit <- rank_partial_residuals(ep$x, ep$y, ep$cov)
      empirical_p <- if (ep$permute && informative_exploratory) {
        stratified_permutation_p(fit, ep$strata, permutation_B)
      } else {
        NA_real_
      }
      slide_stats_outputs[[length(slide_stats_outputs) + 1L]] <- data.table(
        GSM = info$GSM, sample_title = info$sample_title, patient_id = info$patient_id,
        specimen_type = info$specimen_type, primary_analysis = info$primary_analysis,
        axis_id = spec$axis_id, display_axis = spec$display_axis,
        axis_class = spec$axis_class, candidate_type = spec$candidate_type,
        endpoint = endpoint, n_spots = fit$n, rho = fit$rho,
        spot_level_p_descriptive = fit$p,
        composition_stratified_permutation_p = empirical_p,
        informative_primary = informative_primary,
        informative_exploratory = informative_exploratory,
        ligand_detected_fraction = ligand_fraction,
        complete_receptor_detected_fraction = receptor_fraction
      )
    }
  }

  scored$coverage[, `:=`(GSM = info$GSM, sample_title = info$sample_title)]
  qc_outputs[[i]] <- data.table(
    GSM = info$GSM, sample_title = info$sample_title,
    input_filtered_barcodes = length(barcodes), retained_spots = nrow(spot_base),
    min_UMI = thr$min_sum, min_detected = thr$min_detected, max_mito_percent = thr$max_mito
  )
  rm(counts, norm, spot_base, scored); gc(verbose = FALSE)
}

spot_axes <- rbindlist(spot_axis_outputs, fill = TRUE)
slide_stats <- rbindlist(slide_stats_outputs, fill = TRUE)
coverage <- rbindlist(coverage_outputs, fill = TRUE)
qc_summary <- rbindlist(qc_outputs, fill = TRUE)

spot_path <- file.path(processed_dir, "GSE251950_candidate_axis_spot_scores.tsv.gz")
slide_path <- file.path(out_dir, "GSE251950_candidate_axis_slide_statistics.tsv")
coverage_path <- file.path(out_dir, "GSE251950_candidate_axis_gene_coverage.tsv")
qc_path <- file.path(out_dir, "GSE251950_candidate_axis_spatial_QC.tsv")
fwrite(spot_axes, spot_path, sep = "\t", na = "")
fwrite(slide_stats, slide_path, sep = "\t", na = "")
fwrite(coverage, coverage_path, sep = "\t", na = "")
fwrite(qc_summary, qc_path, sep = "\t", na = "")

primary_stats <- slide_stats[
  primary_analysis == TRUE & informative_exploratory == TRUE & is.finite(rho)
]
meta <- primary_stats[, patient_equal_weight_summary(rho, n_spots),
                      by = .(axis_id, display_axis, axis_class, candidate_type, endpoint)]
meta[, fdr := p.adjust(p_value, method = "BH"), by = .(endpoint, candidate_type)]
meta[, direction_fraction := pmax(positive_slides, negative_slides) / pmax(k, 1L)]
meta[, majority_direction := fifelse(positive_slides >= negative_slides, "positive", "negative")]

meta_path <- file.path(out_dir, "GSE251950_candidate_axis_patient_level_meta.tsv")
fwrite(meta, meta_path, sep = "\t", na = "")

wide <- dcast(
  meta[candidate_type == "candidate"],
  axis_id + display_axis + axis_class ~ endpoint,
  value.var = c("k", "pooled_rho", "ci_low", "ci_high", "p_value", "fdr",
                "positive_slides", "negative_slides", "min_loo_rho", "max_loo_rho")
)

get_col <- function(prefix, endpoint) {
  x <- paste(prefix, endpoint, sep = "_")
  if (x %chin% names(wide)) wide[[x]] else rep(NA, nrow(wide))
}
wide[, spatial_same_rho := get_col("pooled_rho", "ligand_receptor_same_adjusted")]
wide[, spatial_neighbor_rho := get_col("pooled_rho", "ligand_receptor_neighbor_adjusted")]
wide[, matrix_same_rho := get_col("pooled_rho", "matrix_axis_same_adjusted")]
wide[, matrix_neighbor_rho := get_col("pooled_rho", "matrix_axis_neighbor_adjusted")]
wide[, ecm_same_rho := get_col("pooled_rho", "ecm_direction_axis_same_adjusted")]
wide[, source_CAF_rho := get_col("pooled_rho", "activated_CAF_source_attribution")]
wide[, receiver_myeloid_rho := get_col("pooled_rho", "myeloid_receiver_attribution")]
wide[, same_positive_fraction := get_col("positive_slides", "ligand_receptor_same_adjusted") /
       pmax(get_col("k", "ligand_receptor_same_adjusted"), 1)]
wide[, neighbor_positive_fraction := get_col("positive_slides", "ligand_receptor_neighbor_adjusted") /
       pmax(get_col("k", "ligand_receptor_neighbor_adjusted"), 1)]
wide[, matrix_positive_fraction := get_col("positive_slides", "matrix_axis_same_adjusted") /
       pmax(get_col("k", "matrix_axis_same_adjusted"), 1)]
wide[, min_informative_slides := pmin(
  get_col("k", "ligand_receptor_same_adjusted"),
  get_col("k", "ligand_receptor_neighbor_adjusted"),
  get_col("k", "matrix_axis_same_adjusted"), na.rm = TRUE
)]
wide[, same_endpoint_supported :=
       same_positive_fraction >= 0.60 & spatial_same_rho > 0 &
       (get_col("ci_low", "ligand_receptor_same_adjusted") > 0 |
          get_col("fdr", "ligand_receptor_same_adjusted") <= 0.10)]
wide[, neighbor_endpoint_supported :=
       neighbor_positive_fraction >= 0.60 & spatial_neighbor_rho > 0 &
       (get_col("ci_low", "ligand_receptor_neighbor_adjusted") > 0 |
          get_col("fdr", "ligand_receptor_neighbor_adjusted") <= 0.10)]
wide[, spatial_endpoint_supported := same_endpoint_supported | neighbor_endpoint_supported]
wide[, matrix_endpoint_supported :=
       matrix_positive_fraction >= 0.60 & matrix_same_rho > 0 &
       (get_col("ci_low", "matrix_axis_same_adjusted") > 0 |
          get_col("fdr", "matrix_axis_same_adjusted") <= 0.10)]
wide[, source_attribution_supported :=
       source_CAF_rho > 0 &
       (get_col("ci_low", "activated_CAF_source_attribution") > 0 |
          get_col("fdr", "activated_CAF_source_attribution") <= 0.10)]
wide[, receiver_attribution_supported :=
       receiver_myeloid_rho > 0 &
       (get_col("ci_low", "myeloid_receiver_attribution") > 0 |
          get_col("fdr", "myeloid_receiver_attribution") <= 0.10)]
wide[, evidence_tier := fcase(
  min_informative_slides >= 5L &
    spatial_endpoint_supported & matrix_endpoint_supported &
    source_attribution_supported & receiver_attribution_supported,
  "confirmatory_spatial_candidate",
  min_informative_slides >= 2L &
    pmax(same_positive_fraction, neighbor_positive_fraction, na.rm = TRUE) >= 0.50 &
    matrix_positive_fraction >= 0.50 &
    pmax(spatial_same_rho, spatial_neighbor_rho, na.rm = TRUE) > 0,
  "exploratory_spatial_candidate",
  min_informative_slides >= 2L,
  "spatially_tested_not_supported",
  default = "technically_uninformative"
)]

integrated_path <- file.path(out_dir, "GSE251950_candidate_axis_integrated_gate.tsv")
wide[, tier_order := match(evidence_tier, c(
  "confirmatory_spatial_candidate", "exploratory_spatial_candidate",
  "spatially_tested_not_supported", "technically_uninformative"
))]
setorder(wide, tier_order, -matrix_same_rho, -spatial_neighbor_rho)
fwrite(wide, integrated_path, sep = "\t", na = "")

plot_endpoints <- c("ligand_receptor_same_adjusted", "ligand_receptor_neighbor_adjusted",
                    "matrix_axis_same_adjusted", "ecm_direction_axis_same_adjusted")
plot_data <- meta[candidate_type == "candidate" & endpoint %chin% plot_endpoints]
plot_data[, endpoint_label := factor(endpoint, levels = plot_endpoints,
  labels = c("Same-spot LR", "Neighbor LR", "Matrix-program coupling", "Frozen ECM-direction coupling"))]
plot_data[, display_axis := factor(display_axis, levels = rev(axis_manifest[candidate_type == "candidate"]$display_axis))]
p <- ggplot(plot_data, aes(pooled_rho, display_axis, color = endpoint_label)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0,
                position = position_dodge(width = 0.58), linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.58), size = 2.2) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Equal-weight pooled within-slide partial Spearman rho", y = NULL,
       color = "Spatial endpoint",
       title = "Frozen candidate-axis spatial and ECM-remodelling gate",
       subtitle = "Nine independent primary gastric tumors; slide/patient is the inference unit") +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
ggsave(file.path(figure_dir, "GSE251950_candidate_axis_spatial_gate.pdf"), p,
       width = 8.4, height = 4.8, device = cairo_pdf)
ggsave(file.path(figure_dir, "GSE251950_candidate_axis_spatial_gate.png"), p,
       width = 8.4, height = 4.8, dpi = 300, bg = "white")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE251950",
  primary_tumor_slides = sum(manifest$primary_analysis),
  sensitivity_metastasis_slides = sum(!manifest$primary_analysis),
  candidates = axis_manifest[candidate_type == "candidate"]$display_axis,
  benchmarks = axis_manifest[candidate_type == "benchmark"]$display_axis,
  fixed_neighborhood = "exact first-ring Visium hex-neighbor mean (six offsets)",
  primary_inference_unit = "patient/slide; spot-level P values are descriptive",
  receptor_complex_rule = "for multi-subunit receptors, spot expression is the minimum log-normalized expression across required components",
  matrix_program_rule = "frozen CAF_ACTIVATED_MATRIX score with each axis ligand and receptor components removed",
  composition_adjustment = "log UMI, frozen fibroblast marker, frozen myeloid marker; neighbor analyses use lagged myeloid marker",
  permutation = paste0(permutation_B, " composition-stratified within-slide residual permutations for adjusted LR endpoints"),
  permutation_limitation = "composition-stratified permutation does not fully preserve spatial autocorrelation; patient/slide synthesis is primary",
  confirmatory_detection_fraction = primary_detection_fraction,
  exploratory_detection_fraction = exploratory_detection_fraction,
  interpretation_warning = paste(
    "Spatial association is not receptor activation or causality.",
    "ECM-program coupling is not proof that ECM is upstream.",
    "Exploratory candidates may be retained without global FDR significance if effect direction and biological support are disclosed."
  ),
  spot_output_sha256 = digest(file = spot_path, algo = "sha256"),
  slide_output_sha256 = digest(file = slide_path, algo = "sha256"),
  meta_output_sha256 = digest(file = meta_path, algo = "sha256"),
  integrated_output_sha256 = digest(file = integrated_path, algo = "sha256")
)
write_json(audit, file.path(audit_dir, "GSE251950_candidate_axis_spatial_manifest.json"),
           pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(audit, pretty = TRUE, auto_unbox = TRUE), "\n")
