#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
set.seed(2026072632L)

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
out_dir <- file.path(project, "analysis", "results", "PLAT_axis", "GSE251950")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

focus_programs <- c(
  "ID_INFLAMMATORY_MONO", "ID_MONO_FCN1_VCAN",
  "GOBP_CELL_MATRIX_ADHESION", "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_MYELOID_LEUKOCYTE_MIGRATION", "GOBP_PHAGOCYTOSIS",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
  "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS", "REACTOME_INTEGRIN_SIGNALING"
)
ligand <- "PLAT"
receptors <- c("ANXA2", "ITGAM", "ITGB2")
competitors <- c("CCL2", "MIF", "THBS1", "SPON2", "CTHRC1")
excluded <- c(ligand, receptors, competitors, "ILK", "LRP1")
program_manifest <- fread(program_file)[program %chin% focus_programs & !gene_symbol %chin% excluded]
programs <- split(program_manifest$gene_symbol, program_manifest$program)

spatial_lag <- function(values, row, col) {
  offsets <- rbind(c(-1, -1), c(-1, 1), c(0, -2), c(0, 2), c(1, -1), c(1, 1))
  keys <- paste(row, col, sep = ":")
  idx <- do.call(cbind, lapply(seq_len(nrow(offsets)), function(i) {
    match(paste(row + offsets[i, 1], col + offsets[i, 2], sep = ":"), keys)
  }))
  vapply(seq_along(values), function(i) {
    j <- idx[i, ]; j <- j[is.finite(j)]
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
    return(c(n = sum(ok), rho = NA_real_, p = NA_real_))
  }
  xr <- rank(x[ok], ties.method = "average")
  yr <- rank(y[ok], ties.method = "average")
  cr <- apply(covariates[ok, , drop = FALSE], 2L, rank, ties.method = "average")
  rx <- lm.fit(cbind(1, cr), xr)$residuals
  ry <- lm.fit(cbind(1, cr), yr)$residuals
  tst <- suppressWarnings(cor.test(rx, ry, method = "pearson"))
  c(n = sum(ok), rho = unname(tst$estimate), p = tst$p.value)
}

patient_meta <- function(rho, B = 9999L) {
  r <- pmax(-0.999999, pmin(0.999999, rho[is.finite(rho)]))
  k <- length(r)
  if (!k) return(data.table(k = 0L, pooled_rho = NA_real_, ci_low = NA_real_,
                            ci_high = NA_real_, sign_p = NA_real_, positive_patients = 0L,
                            negative_patients = 0L, median_patient_rho = NA_real_))
  pooled <- tanh(mean(atanh(r)))
  boot <- replicate(B, tanh(mean(atanh(sample(r, k, replace = TRUE)))))
  ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
  sign_p <- if (k > 1L) suppressWarnings(wilcox.test(r, mu = 0, exact = FALSE)$p.value) else NA_real_
  data.table(
    k = k, pooled_rho = pooled, ci_low = ci[[1]], ci_high = ci[[2]], sign_p = sign_p,
    positive_patients = sum(r > 0), negative_patients = sum(r < 0),
    median_patient_rho = median(r)
  )
}

slide_rows <- list()
coverage_rows <- list()

for (i in seq_len(nrow(manifest))) {
  info <- manifest[i]
  message(sprintf("[%d/9] PLAT spatial incremental analysis: %s", i, info$sample_title))
  base <- copy(base_spots[GSM == info$GSM]); setorder(base, barcode)
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
  rownames(counts) <- unique_symbols; colnames(counts) <- barcodes
  idx <- match(base$barcode, colnames(counts))
  if (anyNA(idx)) stop("Frozen spot mismatch: ", info$GSM)
  counts <- counts[, idx, drop = FALSE]

  needed <- intersect(unique(c(ligand, receptors, competitors, unlist(programs))), rownames(counts))
  norm <- counts[needed, , drop = FALSE] %*% Diagonal(x = 1e4 / pmax(base$total_UMI, 1))
  rownames(norm) <- needed; colnames(norm) <- base$barcode; norm@x <- log1p(norm@x)
  if (!all(c(ligand, receptors, competitors) %chin% rownames(norm))) {
    stop("Required PLAT/receptor/competitor gene missing in ", info$GSM)
  }
  z <- zscore_rows(norm)
  plat <- as.numeric(norm[ligand, ])
  competitor_matrix <- t(as.matrix(norm[competitors, , drop = FALSE]))
  receptor_norm <- as.matrix(norm[receptors, , drop = FALSE])
  receptor_z <- colMeans(z[receptors, , drop = FALSE])
  receptor_lag <- spatial_lag(receptor_z, base$array_row, base$array_col)
  plat_z <- as.numeric(scale(plat))
  receptor_lag_z <- as.numeric(scale(receptor_lag))
  plat_z[!is.finite(plat_z)] <- 0
  receptor_lag_z[!is.finite(receptor_lag_z)] <- 0
  axis_neighbor <- rowMeans(cbind(plat_z, receptor_lag_z), na.rm = TRUE)
  axis_neighbor[!is.finite(axis_neighbor)] <- 0

  base_cov <- cbind(
    log1p(base$total_UMI), base$MARKER_Fibroblast,
    base$lag_MARKER_Myeloid, base$frozen_ECM_direction
  )
  all_competitor_cov <- cbind(base_cov, competitor_matrix)

  for (program_id in names(programs)) {
    present <- intersect(programs[[program_id]], rownames(z))
    informative <- present[apply(z[present, , drop = FALSE], 1L, sd) > 0]
    coverage_rows[[length(coverage_rows) + 1L]] <- data.table(
      GSM = info$GSM, patient_id = info$patient_id, program = program_id,
      defined = length(programs[[program_id]]), available = length(present),
      informative = length(informative)
    )
    if (length(informative) < 5L) next
    program_score <- colMeans(z[informative, , drop = FALSE])
    lag_program <- spatial_lag(program_score, base$array_row, base$array_col)
    endpoints <- list(
      PLAT_neighbor_program_base_adjusted = rank_partial(plat, lag_program, base_cov),
      PLAT_neighbor_program_all_competitors_adjusted = rank_partial(
        plat, lag_program, all_competitor_cov
      ),
      receptor_neighbor_program_base_adjusted = rank_partial(
        receptor_lag, lag_program, base_cov
      ),
      PLAT_receptor_axis_neighbor_program_adjusted = rank_partial(
        axis_neighbor, lag_program, base_cov
      )
    )
    for (competitor in competitors) {
      endpoints[[paste0("PLAT_neighbor_program_plus_", competitor)]] <- rank_partial(
        plat, lag_program, cbind(base_cov, competitor_matrix[, competitor])
      )
    }
    for (endpoint in names(endpoints)) {
      ep <- endpoints[[endpoint]]
      slide_rows[[length(slide_rows) + 1L]] <- data.table(
        GSM = info$GSM, sample_title = info$sample_title, patient_id = info$patient_id,
        program = program_id, endpoint = endpoint, n_spots = as.integer(ep[["n"]]),
        rho = ep[["rho"]], spot_p_descriptive = ep[["p"]],
        PLAT_detected_fraction = mean(plat > 0),
        complete_receptor_detected_fraction = mean(apply(receptor_norm > 0, 2L, all))
      )
    }
  }
  rm(raw, mapper, counts, norm, z); gc(verbose = FALSE)
}

slide <- rbindlist(slide_rows)
coverage <- rbindlist(coverage_rows)
meta <- slide[, patient_meta(rho), by = .(program, endpoint)]
meta[, positive_fraction := positive_patients / pmax(k, 1L)]
meta[, endpoint_FDR_BH := p.adjust(sign_p, "BH"), by = endpoint]

fwrite(slide, file.path(out_dir, "GSE251950_PLAT_spatial_incremental_slide_statistics.tsv"), sep = "\t")
fwrite(meta, file.path(out_dir, "GSE251950_PLAT_spatial_incremental_patient_meta.tsv"), sep = "\t")
fwrite(coverage, file.path(out_dir, "GSE251950_PLAT_spatial_program_coverage.tsv"), sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE251950 / Lee et al. Gut 2025, DOI 10.1136/gutjnl-2024-332901",
  samples = 9L,
  ligand = ligand,
  receptor_transcript_proxy = receptors,
  competitors = competitors,
  neighborhood = "Visium first-ring six-neighbor mean",
  covariates = c("log UMI", "fibroblast marker", "lagged myeloid marker", "frozen ECM direction"),
  primary_incremental_endpoint = "PLAT_neighbor_program_all_competitors_adjusted",
  inference_unit = "patient/section; spot-level P values are descriptive only",
  warning = paste(
    "This is secondary analysis of the nine-section discovery dataset used by Lee et al.;",
    "it is not an independent spatial experiment. Association does not establish sender identity,",
    "surface receptor assembly, PLAT secretion/activity, signaling direction, or causality."
  ),
  input_sha256 = as.list(setNames(vapply(c(frozen_spot_file, program_file), digest, character(1),
                                          algo = "sha256", file = TRUE),
                                  basename(c(frozen_spot_file, program_file))))
)
write_json(audit, file.path(audit_dir, "GSE251950_PLAT_spatial_incremental.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("Primary incremental endpoint:\n")
print(meta[endpoint == "PLAT_neighbor_program_all_competitors_adjusted"][order(endpoint_FDR_BH)])
