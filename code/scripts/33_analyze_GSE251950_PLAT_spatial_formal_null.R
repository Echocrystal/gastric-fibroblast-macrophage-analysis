#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
set.seed(2026072633L)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(jsonlite)
  library(digest)
})

sample_root <- file.path(source_project, "data", "raw", "spatial", "GSE251950", "samples")
frozen_spot_file <- file.path(
  project, "analysis", "processed", "spatial", "GSE251950_candidate_axes",
  "GSE251950_candidate_axis_spot_scores.tsv.gz"
)
program_file <- file.path(project, "resources", "gene_sets", "myeloid_receiver_programs_frozen.tsv")
matched_file <- file.path(
  project, "analysis", "results", "axis_scan", "GSE206785",
  "matched_validation_finalists20_ecm", "GSE206785_top_ECM_ligand_pair_ablation.tsv"
)
out_dir <- file.path(project, "analysis", "results", "PLAT_axis", "GSE251950", "formal_null")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# These endpoints and the positive direction were frozen before this formal-null run.
primary_programs <- c(
  "GOBP_CELL_MATRIX_ADHESION",
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX"
)
ligand <- "PLAT"
receptors <- c("ANXA2", "ITGAM", "ITGB2", "LRP1", "ILK")
competitors <- c("CCL2", "MIF", "THBS1", "SPON2", "CTHRC1")
block_targets <- c(75L, 150L, 300L)
primary_block_target <- 150L
n_permutations <- 9999L

matched_table <- fread(matched_file)
matched_controls <- unique(matched_table[
  variant == "all_20_patients" & candidate == ligand,
  matched_ligand
])
if (length(matched_controls) != 20L) {
  stop("Expected 20 frozen ECM-matched controls for PLAT; observed ", length(matched_controls))
}
expected_controls <- c(
  "SERPINF1", "MMP2", "NID1", "A2M", "MFGE8", "BGN", "TIMP2", "COL5A1",
  "SPON2", "COL12A1", "LAMA4", "FN1", "MMP14", "COL5A2", "COL6A1", "FBLN1",
  "COL4A1", "SERPING1", "LAMB1", "TIMP1"
)
if (!setequal(matched_controls, expected_controls)) {
  stop("Frozen matched-control set differs from the preregistered 20-ligand set")
}
matched_controls <- expected_controls
matched_panel <- c(ligand, matched_controls)

manifest <- fread(file.path(
  source_project, "data", "processed", "spatial", "GSE251950",
  "GSE251950_spatial_sample_manifest.tsv"
))[primary_analysis == TRUE]
if (nrow(manifest) != 9L || uniqueN(manifest$patient_id) != 9L) {
  stop("Expected nine patient-independent primary spatial sections")
}

base_cols <- c(
  "GSM", "sample_title", "patient_id", "barcode", "array_row", "array_col",
  "total_UMI", "frozen_ECM_direction", "MARKER_Fibroblast", "MARKER_Myeloid",
  "lag_MARKER_Myeloid"
)
frozen <- fread(frozen_spot_file)
base_spots <- unique(frozen[primary_analysis == TRUE, ..base_cols])

# To prevent direct part-whole correlation, both primary program scores exclude
# PLAT, every frozen matched ligand and all receptor/competitor transcripts.
score_exclusions <- unique(c(matched_panel, competitors, receptors))
program_manifest <- fread(program_file)[
  program %chin% primary_programs & !gene_symbol %chin% score_exclusions
]
programs <- split(program_manifest$gene_symbol, program_manifest$program)
if (!setequal(names(programs), primary_programs)) stop("Primary program definitions are incomplete")

spatial_lag <- function(values, row, col) {
  offsets <- rbind(c(-1, -1), c(-1, 1), c(0, -2), c(0, 2), c(1, -1), c(1, 1))
  keys <- paste(row, col, sep = ":")
  idx <- do.call(cbind, lapply(seq_len(nrow(offsets)), function(i) {
    match(paste(row + offsets[i, 1], col + offsets[i, 2], sep = ":"), keys)
  }))
  vapply(seq_along(values), function(i) {
    j <- idx[i, ]
    j <- j[is.finite(j)]
    if (!length(j)) NA_real_ else mean(values[j], na.rm = TRUE)
  }, numeric(1))
}

zscore_rows <- function(mat) {
  z <- t(scale(t(as.matrix(mat))))
  z[!is.finite(z)] <- 0
  z
}

rank_residual_pair <- function(x, y, covariates) {
  ok <- is.finite(x) & is.finite(y) & apply(covariates, 1L, function(v) all(is.finite(v)))
  if (sum(ok) < 50L || uniqueN(x[ok]) < 3L || uniqueN(y[ok]) < 3L) return(NULL)
  xr <- rank(x[ok], ties.method = "average")
  yr <- rank(y[ok], ties.method = "average")
  cr <- apply(covariates[ok, , drop = FALSE], 2L, rank, ties.method = "average")
  design <- cbind(1, cr)
  rx <- lm.fit(design, xr)$residuals
  ry <- lm.fit(design, yr)$residuals
  list(ok = ok, rx = rx, ry = ry, rho = unname(cor(rx, ry)), n = sum(ok))
}

patient_meta <- function(rho, B = 9999L) {
  r <- pmax(-0.999999, pmin(0.999999, rho[is.finite(rho)]))
  k <- length(r)
  if (!k) return(data.table(
    k = 0L, pooled_rho = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    wilcoxon_p = NA_real_, sign_test_one_sided_p = NA_real_,
    positive_patients = 0L, negative_patients = 0L, median_patient_rho = NA_real_
  ))
  pooled <- tanh(mean(atanh(r)))
  boot <- replicate(B, tanh(mean(atanh(sample(r, k, replace = TRUE)))))
  ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
  wp <- if (k > 1L) suppressWarnings(wilcox.test(r, mu = 0, exact = FALSE)$p.value) else NA_real_
  pos <- sum(r > 0)
  sign_p <- pbinom(pos - 1L, size = k, prob = 0.5, lower.tail = FALSE)
  data.table(
    k = k, pooled_rho = pooled, ci_low = ci[[1]], ci_high = ci[[2]], wilcoxon_p = wp,
    sign_test_one_sided_p = sign_p, positive_patients = pos,
    negative_patients = sum(r < 0), median_patient_rho = median(r)
  )
}

make_restricted_groups <- function(base, ok, target_size, seed) {
  d <- copy(base[ok])
  coords <- scale(cbind(d$array_row, d$array_col))
  k <- max(2L, as.integer(round(nrow(d) / target_size)))
  k <- min(k, max(2L, floor(nrow(d) / 20L)))
  set.seed(seed)
  km <- kmeans(coords, centers = k, nstart = 20L, iter.max = 100L)
  d[, spatial_block := km$cluster]
  d[, fibro_state := as.integer(MARKER_Fibroblast > median(MARKER_Fibroblast, na.rm = TRUE))]
  d[, ecm_state := as.integer(frozen_ECM_direction > median(frozen_ECM_direction, na.rm = TRUE))]
  d[, stratum := paste(spatial_block, fibro_state, ecm_state, sep = ":")]
  sizes <- d[, .N, by = stratum]
  small <- sizes[N < 4L, stratum]
  d[stratum %chin% small, stratum := paste0(spatial_block, ":small")]
  # A serpentine coordinate order supports circular shifts that preserve the
  # within-stratum residual sequence and therefore more spatial autocorrelation
  # than unrestricted random shuffling.
  d[, row_bin := frank(array_row, ties.method = "dense"), by = stratum]
  d[, path_col := fifelse(row_bin %% 2L == 0L, -array_col, array_col)]
  d[, path_order := frankv(.SD, cols = c("array_row", "path_col"),
                            ties.method = "first"),
    by = stratum, .SDcols = c("array_row", "path_col")]
  groups <- split(seq_len(nrow(d)), d$stratum)
  groups <- lapply(groups, function(idx) idx[order(d$path_order[idx])])
  list(
    groups = groups,
    qc = data.table(
      target_spots_per_block = target_size, n_spots = nrow(d), n_blocks = uniqueN(d$spatial_block),
      n_strata = length(groups), median_stratum_n = median(lengths(groups)),
      min_stratum_n = min(lengths(groups)), max_stratum_n = max(lengths(groups)),
      permutable_fraction = mean(lengths(groups)[match(d$stratum, names(groups))] > 1L)
    )
  )
}

# Exact, fast null for a restricted circular-shift randomization. The numerator
# of Pearson's r is a sum of within-stratum circular cross-products; the global
# mean and variance are invariant to within-stratum shifts.
restricted_shift_null <- function(rx, ry, groups, B, seed) {
  rx <- rx - mean(rx)
  ry <- ry - mean(ry)
  denom <- sqrt(sum(rx^2) * sum(ry^2))
  if (!is.finite(denom) || denom <= 0) return(rep(NA_real_, B))
  set.seed(seed)
  numerators <- numeric(B)
  for (idx in groups) {
    m <- length(idx)
    xg <- rx[idx]
    yg <- ry[idx]
    if (m <= 1L) {
      numerators <- numerators + sum(xg * yg)
      next
    }
    dots <- vapply(0:(m - 1L), function(s) {
      shifted <- if (s == 0L) xg else c(xg[(s + 1L):m], xg[seq_len(s)])
      sum(shifted * yg)
    }, numeric(1))
    # Excluding the zero shift ensures every permutable stratum is perturbed.
    sampled_shift <- sample.int(m - 1L, size = B, replace = TRUE) + 1L
    numerators <- numerators + dots[sampled_shift]
  }
  numerators / denom
}

slide_primary <- list()
slide_matched <- list()
strata_qc <- list()
residual_store <- list()
coverage_rows <- list()

for (i in seq_len(nrow(manifest))) {
  info <- manifest[i]
  message(sprintf("[%d/9] formal spatial null input: %s", i, info$sample_title))
  base <- copy(base_spots[GSM == info$GSM])
  setorder(base, barcode)
  folder <- file.path(sample_root, info$folder_id)
  prefix <- file.path(folder, info$folder_id)
  features <- fread(paste0(prefix, "_features.tsv.gz"), header = FALSE)
  barcodes <- fread(paste0(prefix, "_barcodes.tsv.gz"), header = FALSE)$V1
  raw <- as(readMM(gzfile(paste0(prefix, "_matrix.mtx.gz"))), "dgCMatrix")
  symbols <- features$V2
  valid <- !is.na(symbols) & nzchar(symbols)
  unique_symbols <- sort(unique(symbols[valid]))
  mapper <- sparseMatrix(
    i = match(symbols[valid], unique_symbols), j = which(valid), x = 1,
    dims = c(length(unique_symbols), nrow(raw))
  )
  counts <- mapper %*% raw
  rownames(counts) <- unique_symbols
  colnames(counts) <- barcodes
  idx <- match(base$barcode, colnames(counts))
  if (anyNA(idx)) stop("Frozen spot mismatch: ", info$GSM)
  counts <- counts[, idx, drop = FALSE]

  needed <- intersect(unique(c(matched_panel, competitors, unlist(programs))), rownames(counts))
  norm <- counts[needed, , drop = FALSE] %*% Diagonal(x = 1e4 / pmax(base$total_UMI, 1))
  rownames(norm) <- needed
  colnames(norm) <- base$barcode
  norm@x <- log1p(norm@x)
  missing_panel <- setdiff(unique(c(matched_panel, competitors)), rownames(norm))
  if (length(missing_panel)) stop("Required spatial ligand missing in ", info$GSM, ": ", paste(missing_panel, collapse = ", "))
  z <- zscore_rows(norm)

  base_cov <- cbind(
    log1p(base$total_UMI), base$MARKER_Fibroblast,
    base$lag_MARKER_Myeloid, base$frozen_ECM_direction
  )
  competitor_matrix <- t(as.matrix(norm[competitors, , drop = FALSE]))
  incremental_cov <- cbind(base_cov, competitor_matrix)

  for (program_id in primary_programs) {
    present <- intersect(programs[[program_id]], rownames(z))
    informative <- present[apply(z[present, , drop = FALSE], 1L, sd) > 0]
    coverage_rows[[length(coverage_rows) + 1L]] <- data.table(
      GSM = info$GSM, patient_id = info$patient_id, program = program_id,
      defined_after_common_exclusion = length(programs[[program_id]]),
      available = length(present), informative = length(informative)
    )
    if (length(informative) < 5L) stop("Insufficient program coverage in ", info$GSM)
    program_score <- colMeans(z[informative, , drop = FALSE])
    lag_program <- spatial_lag(program_score, base$array_row, base$array_col)

    plat <- as.numeric(norm[ligand, ])
    inc_fit <- rank_residual_pair(plat, lag_program, incremental_cov)
    base_fit <- rank_residual_pair(plat, lag_program, base_cov)
    if (is.null(inc_fit) || is.null(base_fit)) stop("PLAT model failed in ", info$GSM)
    slide_primary[[length(slide_primary) + 1L]] <- data.table(
      GSM = info$GSM, sample_title = info$sample_title, patient_id = info$patient_id,
      program = program_id,
      model = c("base_adjusted", "all_competitors_adjusted"),
      n_spots = c(base_fit$n, inc_fit$n), rho = c(base_fit$rho, inc_fit$rho),
      PLAT_detected_fraction = mean(plat > 0)
    )
    residual_store[[paste(info$GSM, program_id, sep = "|")]] <- list(
      info = info, base = base, fit = inc_fit
    )

    for (candidate in matched_panel) {
      candidate_x <- as.numeric(norm[candidate, ])
      fit <- rank_residual_pair(candidate_x, lag_program, base_cov)
      slide_matched[[length(slide_matched) + 1L]] <- data.table(
        GSM = info$GSM, sample_title = info$sample_title, patient_id = info$patient_id,
        program = program_id, ligand = candidate, n_spots = if (is.null(fit)) NA_integer_ else fit$n,
        rho = if (is.null(fit)) NA_real_ else fit$rho,
        detected_fraction = mean(candidate_x > 0)
      )
    }
  }
  rm(raw, mapper, counts, norm, z)
  gc(verbose = FALSE)
}

slide_primary_dt <- rbindlist(slide_primary)
slide_matched_dt <- rbindlist(slide_matched)
coverage_dt <- rbindlist(coverage_rows)

primary_meta <- slide_primary_dt[, patient_meta(rho), by = .(program, model)]
primary_meta[, wilcoxon_FDR_BH := p.adjust(wilcoxon_p, "BH"), by = model]

matched_meta <- slide_matched_dt[, patient_meta(rho), by = .(program, ligand)]
matched_meta[, is_PLAT := ligand == "PLAT"]
matched_meta[, rank_descending := frank(-pooled_rho, ties.method = "min"), by = program]
matched_comparison <- matched_meta[, {
  plat_rho <- pooled_rho[ligand == "PLAT"]
  controls <- pooled_rho[ligand != "PLAT"]
  plat_pos <- positive_patients[ligand == "PLAT"]
  control_pos <- positive_patients[ligand != "PLAT"]
  .(
    PLAT_pooled_rho = plat_rho,
    PLAT_rank_of_21 = rank_descending[ligand == "PLAT"],
    controls_exceeding_or_equal_PLAT = sum(controls >= plat_rho),
    matched_ligand_empirical_p = (1 + sum(controls >= plat_rho)) / (1 + length(controls)),
    PLAT_positive_patients = plat_pos,
    controls_with_at_least_as_many_positive_patients = sum(control_pos >= plat_pos)
  )
}, by = program]

# Formal architecture-restricted null for the preregistered incremental PLAT model.
null_meta_rows <- list()
null_draw_rows <- list()
for (target in block_targets) {
  message("Restricted circular-shift null; target block size = ", target)
  null_by_program <- setNames(vector("list", length(primary_programs)), primary_programs)
  observed_by_program <- setNames(vector("list", length(primary_programs)), primary_programs)
  for (program_id in primary_programs) {
    patient_null <- matrix(NA_real_, nrow = nrow(manifest), ncol = n_permutations)
    patient_observed <- numeric(nrow(manifest))
    for (i in seq_len(nrow(manifest))) {
      info <- manifest[i]
      obj <- residual_store[[paste(info$GSM, program_id, sep = "|")]]
      restricted <- make_restricted_groups(
        obj$base, obj$fit$ok, target,
        seed = 2026072633L + target * 100L + i
      )
      q <- copy(restricted$qc)
      q[, `:=`(GSM = info$GSM, patient_id = info$patient_id, program = program_id)]
      strata_qc[[length(strata_qc) + 1L]] <- q
      patient_observed[i] <- obj$fit$rho
      patient_null[i, ] <- restricted_shift_null(
        obj$fit$rx, obj$fit$ry, restricted$groups, n_permutations,
        seed = 2026072633L + target * 1000L + i * 10L + match(program_id, primary_programs)
      )
    }
    obs_pooled <- tanh(mean(atanh(pmax(-0.999999, pmin(0.999999, patient_observed)))))
    null_pooled <- apply(patient_null, 2L, function(v) {
      tanh(mean(atanh(pmax(-0.999999, pmin(0.999999, v)))))
    })
    obs_positive <- sum(patient_observed > 0)
    null_positive <- colSums(patient_null > 0, na.rm = TRUE)
    null_meta_rows[[length(null_meta_rows) + 1L]] <- data.table(
      program = program_id, target_spots_per_block = target,
      n_patients = length(patient_observed), observed_pooled_rho = obs_pooled,
      observed_positive_patients = obs_positive,
      null_pooled_mean = mean(null_pooled), null_pooled_sd = sd(null_pooled),
      null_pooled_q025 = quantile(null_pooled, 0.025),
      null_pooled_q975 = quantile(null_pooled, 0.975),
      restricted_empirical_p_one_sided = (1 + sum(null_pooled >= obs_pooled)) / (n_permutations + 1),
      restricted_empirical_p_two_sided = (1 + sum(abs(null_pooled) >= abs(obs_pooled))) / (n_permutations + 1),
      null_positive_patients_mean = mean(null_positive),
      null_positive_patients_q95 = quantile(null_positive, 0.95),
      positive_concordance_empirical_p = (1 + sum(null_positive >= obs_positive)) / (n_permutations + 1),
      n_permutations = n_permutations
    )
    null_draw_rows[[length(null_draw_rows) + 1L]] <- data.table(
      program = program_id, target_spots_per_block = target,
      permutation = seq_len(n_permutations), pooled_rho = null_pooled,
      positive_patients = null_positive
    )
  }
}

null_meta_dt <- rbindlist(null_meta_rows)
null_meta_dt[, restricted_FDR_BH_one_sided := p.adjust(restricted_empirical_p_one_sided, "BH"),
             by = target_spots_per_block]
null_draw_dt <- rbindlist(null_draw_rows)
strata_qc_dt <- rbindlist(strata_qc)

fwrite(slide_primary_dt, file.path(out_dir, "GSE251950_PLAT_formal_null_slide_statistics.tsv"), sep = "\t")
fwrite(primary_meta, file.path(out_dir, "GSE251950_PLAT_formal_null_patient_meta.tsv"), sep = "\t")
fwrite(coverage_dt, file.path(out_dir, "GSE251950_PLAT_formal_null_program_coverage.tsv"), sep = "\t")
fwrite(strata_qc_dt, file.path(out_dir, "GSE251950_PLAT_formal_null_strata_qc.tsv"), sep = "\t")
fwrite(null_meta_dt, file.path(out_dir, "GSE251950_PLAT_restricted_permutation_summary.tsv"), sep = "\t")
fwrite(null_draw_dt, file.path(out_dir, "GSE251950_PLAT_restricted_permutation_draws.tsv.gz"), sep = "\t")
fwrite(slide_matched_dt, file.path(out_dir, "GSE251950_frozen_ECM_matched_ligand_slide_statistics.tsv"), sep = "\t")
fwrite(matched_meta, file.path(out_dir, "GSE251950_frozen_ECM_matched_ligand_patient_meta.tsv"), sep = "\t")
fwrite(matched_comparison, file.path(out_dir, "GSE251950_PLAT_vs_frozen_ECM_matched_ligands.tsv"), sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE251950 / Lee et al. Gut 2025, DOI 10.1136/gutjnl-2024-332901",
  analysis_role = "secondary spatial reanalysis; not an independent validation cohort",
  primary_endpoints = primary_programs,
  positive_direction_frozen = TRUE,
  primary_model = paste(
    "within-section partial rank correlation of PLAT with first-ring neighboring program, adjusted for",
    "log UMI, fibroblast marker, lagged myeloid marker, frozen ECM direction, CCL2, MIF, THBS1, SPON2 and CTHRC1"
  ),
  program_score_exclusions = score_exclusions,
  matched_ligands = matched_controls,
  matched_ligand_model = paste(
    "common base-adjusted model for PLAT and all 20 controls; program scores exclude the entire ligand panel",
    "to prevent direct part-whole correlation"
  ),
  spatial_null = list(
    method = paste(
      "within-slide restricted circular shift of rank-model PLAT residuals inside coordinate k-means blocks",
      "crossed with median-split fibroblast and frozen-ECM states; strata smaller than four spots pooled",
      "within the same spatial block"
    ),
    primary_target_spots_per_block = primary_block_target,
    sensitivity_target_spots_per_block = setdiff(block_targets, primary_block_target),
    permutations = n_permutations,
    note = paste(
      "No author pathology-region labels were available. These are coordinate- and expression-derived",
      "architecture restrictions and must not be described as pathologist-annotated tumor/stroma strata."
    )
  ),
  inference_unit = "nine patient-independent primary-tumor sections; spot-level tests are not inferential units",
  claim_boundary = paste(
    "A positive result supports a small, patient-consistent, PLAT-linked spatial association beyond the specified",
    "competitors and spatial null. It does not establish fibroblast secretion, extracellular protease activity,",
    "receptor binding, signaling direction, ECM degradation caused by PLAT, or biological causality."
  ),
  input_sha256 = as.list(setNames(
    vapply(c(frozen_spot_file, program_file, matched_file), digest, character(1), algo = "sha256", file = TRUE),
    basename(c(frozen_spot_file, program_file, matched_file))
  ))
)
write_json(
  audit, file.path(audit_dir, "GSE251950_PLAT_spatial_formal_null.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA
)

cat("\nPrimary incremental model:\n")
print(primary_meta[model == "all_competitors_adjusted"])
cat("\nRestricted permutation summary:\n")
print(null_meta_dt[order(target_spots_per_block, restricted_empirical_p_one_sided)])
cat("\nFrozen ECM-matched ligand comparison:\n")
print(matched_comparison)
