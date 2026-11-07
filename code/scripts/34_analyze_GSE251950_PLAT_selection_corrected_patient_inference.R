#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
set.seed(2026072634L)

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
out_dir <- file.path(project, "analysis", "results", "PLAT_axis", "GSE251950", "selection_corrected")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# This is the complete 11-program family searched in script 32. No endpoint is
# called preregistered or independently frozen. Family-wise inference below
# explicitly corrects for having searched this entire family on the same data.
program_family <- c(
  "ID_INFLAMMATORY_MONO", "ID_MONO_FCN1_VCAN",
  "GOBP_CELL_MATRIX_ADHESION", "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_MYELOID_LEUKOCYTE_MIGRATION", "GOBP_PHAGOCYTOSIS",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
  "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS", "REACTOME_INTEGRIN_SIGNALING"
)
ligand <- "PLAT"
receptors <- c("ANXA2", "ITGAM", "ITGB2", "LRP1", "ILK")
competitors <- c("CCL2", "MIF", "THBS1", "SPON2", "CTHRC1")

matched_table <- fread(matched_file)
matched_controls <- unique(matched_table[
  variant == "all_20_patients" & candidate == ligand,
  matched_ligand
])
if (length(matched_controls) != 20L) stop("Expected 20 frozen PLAT matched controls")
score_exclusions <- unique(c(ligand, matched_controls, competitors, receptors))

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
program_manifest <- fread(program_file)[
  program %chin% program_family & !gene_symbol %chin% score_exclusions
]
programs <- split(program_manifest$gene_symbol, program_manifest$program)
if (!setequal(names(programs), program_family)) stop("The complete 11-program family is unavailable")

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

rank_partial <- function(x, y, covariates) {
  ok <- is.finite(x) & is.finite(y) & apply(covariates, 1L, function(v) all(is.finite(v)))
  if (sum(ok) < 50L || uniqueN(x[ok]) < 3L || uniqueN(y[ok]) < 3L) {
    return(c(n = sum(ok), rho = NA_real_))
  }
  xr <- rank(x[ok], ties.method = "average")
  yr <- rank(y[ok], ties.method = "average")
  cr <- apply(covariates[ok, , drop = FALSE], 2L, rank, ties.method = "average")
  design <- cbind(1, cr)
  rx <- lm.fit(design, xr)$residuals
  ry <- lm.fit(design, yr)$residuals
  c(n = sum(ok), rho = unname(cor(rx, ry)))
}

slide_rows <- list()
coverage_rows <- list()

for (i in seq_len(nrow(manifest))) {
  info <- manifest[i]
  message(sprintf("[%d/9] selection-corrected patient inference: %s", i, info$sample_title))
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

  needed <- intersect(unique(c(ligand, competitors, unlist(programs))), rownames(counts))
  norm <- counts[needed, , drop = FALSE] %*% Diagonal(x = 1e4 / pmax(base$total_UMI, 1))
  rownames(norm) <- needed
  colnames(norm) <- base$barcode
  norm@x <- log1p(norm@x)
  if (!all(c(ligand, competitors) %chin% rownames(norm))) stop("Required ligand missing in ", info$GSM)
  z <- zscore_rows(norm)
  plat <- as.numeric(norm[ligand, ])
  competitor_matrix <- t(as.matrix(norm[competitors, , drop = FALSE]))
  covariates <- cbind(
    log1p(base$total_UMI), base$MARKER_Fibroblast,
    base$lag_MARKER_Myeloid, base$frozen_ECM_direction,
    competitor_matrix
  )

  for (program_id in program_family) {
    present <- intersect(programs[[program_id]], rownames(z))
    informative <- present[apply(z[present, , drop = FALSE], 1L, sd) > 0]
    coverage_rows[[length(coverage_rows) + 1L]] <- data.table(
      GSM = info$GSM, patient_id = info$patient_id, program = program_id,
      defined_after_common_exclusion = length(programs[[program_id]]),
      available = length(present), informative = length(informative)
    )
    if (length(informative) < 5L) next
    score <- colMeans(z[informative, , drop = FALSE])
    lag_score <- spatial_lag(score, base$array_row, base$array_col)
    fit <- rank_partial(plat, lag_score, covariates)
    slide_rows[[length(slide_rows) + 1L]] <- data.table(
      GSM = info$GSM, sample_title = info$sample_title, patient_id = info$patient_id,
      program = program_id, n_spots = as.integer(fit[["n"]]), rho = fit[["rho"]],
      PLAT_detected_fraction = mean(plat > 0)
    )
  }
  rm(raw, mapper, counts, norm, z)
  gc(verbose = FALSE)
}

slide <- rbindlist(slide_rows)
coverage <- rbindlist(coverage_rows)
if (slide[, uniqueN(patient_id), by = program][, any(V1 != 9L)]) stop("Incomplete patient coverage")

wide <- dcast(slide, patient_id ~ program, value.var = "rho")
setorder(wide, patient_id)
rho_matrix <- as.matrix(wide[, ..program_family])
rho_matrix[rho_matrix > 0.999999] <- 0.999999
rho_matrix[rho_matrix < -0.999999] <- -0.999999
z_matrix <- atanh(rho_matrix)

# Exhaust all 2^9 patient-level sign flips. A single patient's sign is flipped
# jointly across all 11 endpoints, preserving endpoint dependence. The maximum
# absolute Fisher-z mean supplies selection-aware two-sided FWER control.
sign_grid <- as.matrix(expand.grid(rep(list(c(-1, 1)), nrow(z_matrix))))
observed_stat <- abs(colMeans(z_matrix))
permutation_stats <- t(apply(sign_grid, 1L, function(s) abs(colMeans(z_matrix * s))))
maxT <- apply(permutation_stats, 1L, max)

meta_rows <- lapply(seq_along(program_family), function(j) {
  program_id <- program_family[j]
  r <- rho_matrix[, j]
  pooled <- tanh(mean(atanh(r)))
  set.seed(2026072634L + j)
  boot <- replicate(9999L, tanh(mean(atanh(sample(r, length(r), replace = TRUE)))))
  ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
  pos <- sum(r > 0)
  one_tail <- pbinom(pos - 1L, size = length(r), prob = 0.5, lower.tail = FALSE)
  two_sided_sign <- min(1, 2 * min(one_tail, pbinom(pos, size = length(r), prob = 0.5)))
  data.table(
    program = program_id, k = length(r), pooled_rho = pooled,
    ci_low = ci[[1]], ci_high = ci[[2]], positive_patients = pos,
    negative_patients = sum(r < 0), median_patient_rho = median(r),
    wilcoxon_two_sided_p = suppressWarnings(wilcox.test(r, mu = 0, exact = FALSE)$p.value),
    exact_sign_two_sided_p = two_sided_sign,
    patient_signflip_unadjusted_p = mean(permutation_stats[, j] >= observed_stat[j]),
    patient_signflip_maxT_FWER_p = mean(maxT >= observed_stat[j]),
    observed_abs_mean_fisher_z = observed_stat[j]
  )
})
meta <- rbindlist(meta_rows)
meta[, wilcoxon_Bonferroni_11 := p.adjust(wilcoxon_two_sided_p, "bonferroni")]
meta[, wilcoxon_FDR_BH_11 := p.adjust(wilcoxon_two_sided_p, "BH")]
meta[, exact_sign_Bonferroni_11 := p.adjust(exact_sign_two_sided_p, "bonferroni")]
meta[, exact_sign_FDR_BH_11 := p.adjust(exact_sign_two_sided_p, "BH")]
meta[, maxT_FWER_gate_0_05 := patient_signflip_maxT_FWER_p <= 0.05]

maxT_null <- data.table(
  sign_pattern = seq_len(nrow(sign_grid)),
  max_abs_mean_fisher_z = maxT
)

fwrite(slide, file.path(out_dir, "GSE251950_PLAT_all11_selection_corrected_slide_statistics.tsv"), sep = "\t")
fwrite(meta, file.path(out_dir, "GSE251950_PLAT_all11_selection_corrected_patient_meta.tsv"), sep = "\t")
fwrite(coverage, file.path(out_dir, "GSE251950_PLAT_all11_selection_corrected_program_coverage.tsv"), sep = "\t")
fwrite(maxT_null, file.path(out_dir, "GSE251950_PLAT_all11_patient_signflip_maxT_null.tsv"), sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  reason = paste(
    "Corrects the data-dependent selection in script 32, which examined all 11 programs on the same nine sections",
    "before cell-matrix adhesion and ECM degradation were promoted."
  ),
  endpoint_family = program_family,
  inference = paste(
    "patient-equal two-sided sign-flip test over all 2^9 patterns; patient signs flipped jointly across all 11",
    "endpoints; max absolute mean Fisher-z statistic controls FWER for endpoint selection"
  ),
  robustness = c("two-sided exact sign test", "two-sided Wilcoxon test", "nonparametric patient bootstrap CI"),
  program_score_exclusions = score_exclusions,
  warning = paste(
    "This corrects the 11-endpoint selection but remains a secondary reanalysis of one published nine-section",
    "dataset. It does not make the spatial cohort independent and does not establish PLAT specificity or causality."
  ),
  input_sha256 = as.list(setNames(
    vapply(c(frozen_spot_file, program_file, matched_file), digest, character(1), algo = "sha256", file = TRUE),
    basename(c(frozen_spot_file, program_file, matched_file))
  ))
)
write_json(
  audit, file.path(audit_dir, "GSE251950_PLAT_selection_corrected_patient_inference.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA
)

cat("\nSelection-corrected patient-level inference across all 11 originally searched programs:\n")
print(meta[order(patient_signflip_maxT_FWER_p, -abs(pooled_rho))])
