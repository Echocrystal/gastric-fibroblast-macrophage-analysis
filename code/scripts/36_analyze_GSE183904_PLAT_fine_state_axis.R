#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
set.seed(2026072736L)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(edgeR)
  library(jsonlite)
  library(digest)
})

manifest_file <- file.path(source_project, "data", "processed", "scrna", "GSE183904",
                           "GSE183904_sample_patient_clinical_manifest.tsv")
sce_dir <- file.path(source_project, "data", "intermediate", "scrna", "GSE183904", "sample_sce")
fib_meta_file <- file.path(project, "analysis", "processed", "GSE183904_fine_states",
                           "GSE183904_Fibroblast_fine_state_cell_metadata.tsv.gz")
myeloid_meta_file <- file.path(project, "analysis", "processed", "GSE183904_fine_states",
                               "GSE183904_Myeloid_fine_state_cell_metadata.tsv.gz")
label_file <- file.path(project, "analysis", "audit", "GSE183904_fine_state_label_review.tsv")
prereg_file <- file.path(project, "analysis", "audit",
                         "GSE183904_fine_state_axis_association_preregistration.json")
amendment_file <- file.path(project, "analysis", "audit",
                            "GSE183904_fine_state_axis_association_amendment_01.json")
program_file <- file.path(project, "resources", "gene_sets", "myeloid_receiver_programs_frozen.tsv")
ecm_file <- file.path(source_project, "results", "scrna", "GSE183904",
                      "GSE183904_primary_tumor_frozen_ECM_scores.tsv")
out_dir <- file.path(project, "analysis", "results", "PLAT_axis", "GSE183904", "fine_state_axis")
processed_dir <- file.path(project, "analysis", "processed", "GSE183904_fine_states")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(manifest_file, fib_meta_file, myeloid_meta_file, label_file,
            prereg_file, amendment_file, program_file, ecm_file)
if (!all(file.exists(inputs))) stop("Missing fine-state axis input")

manifest <- fread(manifest_file)[include_primary_replication == TRUE &
                                   specimen_class == "Primary_Tumor"]
fib <- fread(fib_meta_file)
myeloid <- fread(myeloid_meta_file)
labels <- fread(label_file)
program_manifest <- fread(program_file)
ecm <- fread(ecm_file)[, .(sample_id, ECM_direction)]
if (nrow(manifest) != 26L || uniqueN(manifest$Patient_ID) != 26L) {
  stop("Expected 26 patient-independent primary tumors")
}

fib_map <- setNames(labels[lineage == "Fibroblast"]$refined_label,
                    labels[lineage == "Fibroblast"]$cluster)
myeloid_map <- setNames(labels[lineage == "Myeloid"]$refined_label,
                        labels[lineage == "Myeloid"]$cluster)
fib[, refined_label := unname(fib_map[as.character(cluster)])]
myeloid[, refined_label := unname(myeloid_map[as.character(cluster)])]
if (fib[, anyNA(refined_label)] || myeloid[, anyNA(refined_label)]) {
  stop("Incomplete refined-state mapping")
}

sender_state <- "CXCL14_POSTN_PLAT_rich_fibroblast"
receiver_state <- "C1QC_APOE_macrophage"
focus_programs <- c(
  "GOBP_PHAGOCYTOSIS",
  "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
  "ID_INFLAMMATORY_MONO",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "GOBP_CELL_MATRIX_ADHESION",
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX"
)
axis_genes <- c("PLAT", "ANXA2", "ITGAM", "ITGB2", "LRP1", "ILK")
program_manifest <- program_manifest[program %chin% focus_programs &
                                       !gene_symbol %chin% axis_genes]
programs <- split(program_manifest$gene_symbol, program_manifest$program)
if (!setequal(names(programs), focus_programs)) stop("Incomplete frozen program family")
lineage_counts <- rbindlist(list(
  fib[, .(lineage_cells = .N), by = .(sample_id, Patient_ID)][, lineage := "Fibroblast"],
  myeloid[, .(lineage_cells = .N), by = .(sample_id, Patient_ID)][, lineage := "Myeloid"]
))
state_counts <- rbindlist(list(
  fib[, .(state_cells = .N), by = .(sample_id, Patient_ID, refined_label)][
    , lineage := "Fibroblast"],
  myeloid[, .(state_cells = .N), by = .(sample_id, Patient_ID, refined_label)][
    , lineage := "Myeloid"]
))
state_counts <- merge(state_counts, lineage_counts,
                      by = c("sample_id", "Patient_ID", "lineage"), all.x = TRUE)
state_counts[, state_fraction := state_cells / lineage_cells]
fwrite(state_counts, file.path(out_dir, "GSE183904_fine_state_counts_by_patient.tsv"), sep = "\t")

aggregate_rows <- list()
for (i in seq_len(nrow(manifest))) {
  info <- manifest[i]
  message(sprintf("[%d/26] fine-state pseudobulk: %s", i, info$sample_id))
  sce <- readRDS(file.path(sce_dir, paste0(info$sample_id, "_secondary_qc_sce.rds")))
  counts_mat <- counts(sce)
  present <- rownames(counts_mat)
  if (!all(axis_genes %chin% present)) stop("Missing axis genes: ", info$sample_id)

  for (lineage_id in c("Fibroblast", "Myeloid")) {
    cm <- if (lineage_id == "Fibroblast") fib[sample_id == info$sample_id] else
      myeloid[sample_id == info$sample_id]
    idx <- match(cm$cell_id, colnames(counts_mat))
    if (anyNA(idx)) stop("Cell/count alignment failure: ", info$sample_id)
    groups <- split(idx, cm$refined_label)
    groups <- groups[names(groups) %chin% c(sender_state, receiver_state)]
    for (state_id in names(groups)) {
      cells <- groups[[state_id]]
      sums <- Matrix::rowSums(counts_mat[present, cells, drop = FALSE])
      detected <- Matrix::rowSums(counts_mat[present, cells, drop = FALSE] > 0)
      aggregate_rows[[length(aggregate_rows) + 1L]] <- data.table(
        sample_id = info$sample_id, Patient_ID = info$Patient_ID,
        lineage = lineage_id, refined_label = state_id, n_cells = length(cells),
        gene_symbol = present, count_sum = as.numeric(sums),
        detected_fraction = as.numeric(detected) / length(cells)
      )
    }
  }
  rm(sce, counts_mat)
  gc(verbose = FALSE)
}
aggregated <- rbindlist(aggregate_rows)
fwrite(aggregated, file.path(processed_dir,
                             "GSE183904_fine_state_full_transcriptome_pseudobulk.tsv.gz"), sep = "\t")

normalize_state <- function(state_id, min_cells = 10L) {
  dt <- aggregated[refined_label == state_id & n_cells >= min_cells]
  if (!nrow(dt)) stop("No eligible pseudobulk state: ", state_id)
  wide <- dcast(dt, gene_symbol ~ sample_id, value.var = "count_sum", fill = 0)
  genes <- wide$gene_symbol
  mat <- as.matrix(wide[, -1L])
  rownames(mat) <- genes
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  dge <- calcNormFactors(DGEList(counts = mat), method = "TMM")
  lcpm <- cpm(dge, log = TRUE, prior.count = 0.5)
  list(
    expression = lcpm,
    normalization = data.table(
      sample_id = colnames(lcpm), target_panel_library_size = dge$samples$lib.size,
      TMM_factor = dge$samples$norm.factors,
      effective_library_size = dge$samples$lib.size * dge$samples$norm.factors
    ),
    cell_counts = unique(dt[, .(sample_id, Patient_ID, n_cells)])
  )
}

safe_z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

sender_norm <- normalize_state(sender_state, 10L)
receiver_norm <- normalize_state(receiver_state, 10L)
fwrite(sender_norm$normalization, file.path(out_dir, "GSE183904_sender_state_TMM_factors.tsv"),
       sep = "\t")
fwrite(receiver_norm$normalization, file.path(out_dir, "GSE183904_receiver_state_TMM_factors.tsv"),
       sep = "\t")

sender_expr <- data.table(
  sample_id = colnames(sender_norm$expression),
  sender_PLAT_logCPM = as.numeric(sender_norm$expression["PLAT", ])
)
sender_expr[, sender_PLAT_z := safe_z(sender_PLAT_logCPM)]

receiver_z <- t(apply(receiver_norm$expression, 1L, safe_z))
rownames(receiver_z) <- rownames(receiver_norm$expression)
colnames(receiver_z) <- colnames(receiver_norm$expression)
program_rows <- list()
for (program_id in focus_programs) {
  present <- intersect(programs[[program_id]], rownames(receiver_z))
  informative <- present[apply(receiver_z[present, , drop = FALSE], 1L,
                               function(v) any(is.finite(v)))]
  if (length(informative) < 5L) {
    stop("Insufficient receiver-state program coverage: ", program_id)
  }
  program_rows[[program_id]] <- data.table(
    sample_id = colnames(receiver_z), program = program_id,
    program_score = colMeans(receiver_z[informative, , drop = FALSE], na.rm = TRUE),
    n_program_genes = length(informative)
  )
}
program_scores <- rbindlist(program_rows)

sample_features <- merge(data.table(sample_id = manifest$sample_id,
                                    Patient_ID = manifest$Patient_ID),
                         ecm, by = "sample_id", all.x = TRUE)
fib_total <- lineage_counts[lineage == "Fibroblast",
                            .(sample_id, fibroblast_cells = lineage_cells)]
my_total <- lineage_counts[lineage == "Myeloid",
                           .(sample_id, myeloid_cells = lineage_cells)]
sample_features <- merge(sample_features, fib_total, by = "sample_id", all.x = TRUE)
sample_features <- merge(sample_features, my_total, by = "sample_id", all.x = TRUE)

sender_counts <- state_counts[refined_label == sender_state,
                              .(sample_id, sender_state_cells = state_cells,
                                sender_state_fraction = state_fraction)]
receiver_counts <- state_counts[refined_label == receiver_state,
                                .(sample_id, receiver_state_cells = state_cells,
                                  receiver_state_fraction = state_fraction)]
sample_features <- merge(sample_features, sender_counts, by = "sample_id", all.x = TRUE)
sample_features <- merge(sample_features, receiver_counts, by = "sample_id", all.x = TRUE)
sample_features[is.na(sender_state_cells), sender_state_cells := 0L]
sample_features[is.na(sender_state_fraction), sender_state_fraction := 0]
sample_features[is.na(receiver_state_cells), receiver_state_cells := 0L]
sample_features[is.na(receiver_state_fraction), receiver_state_fraction := 0]
sample_features <- merge(sample_features, sender_expr, by = "sample_id", all.x = TRUE)

receiver_cells <- myeloid[refined_label == receiver_state]
codetection <- receiver_cells[, .(
  n_receiver_state_cells = .N,
  receiver_triple_cells = sum(ANXA2_ITGAM_ITGB2_triple_detected),
  receiver_triple_fraction = mean(ANXA2_ITGAM_ITGB2_triple_detected),
  receiver_ANXA2_LRP1_cells = sum(ANXA2_LRP1_double_detected),
  receiver_ANXA2_LRP1_fraction = mean(ANXA2_LRP1_double_detected),
  receiver_ANXA2_mean = mean(expr_ANXA2),
  receiver_ITGAM_detected_fraction = mean(detect_ITGAM),
  receiver_ITGB2_detected_fraction = mean(detect_ITGB2)
), by = .(sample_id, Patient_ID)]
sample_features <- merge(sample_features, codetection,
                         by = c("sample_id", "Patient_ID"), all.x = TRUE)
sample_features[, sender_state_fraction_asin := asin(sqrt(sender_state_fraction))]
sample_features[, receiver_state_fraction_asin := asin(sqrt(receiver_state_fraction))]
sample_features[, receiver_triple_fraction_asin := asin(sqrt(receiver_triple_fraction))]
fwrite(sample_features, file.path(out_dir, "GSE183904_fine_state_axis_patient_features.tsv"),
       sep = "\t")
fwrite(program_scores, file.path(out_dir, "GSE183904_receiver_state_program_scores.tsv"),
       sep = "\t")

partial_test <- function(x, y, covariates, B = 9999L, seed = 1L) {
  ok <- is.finite(x) & is.finite(y) & complete.cases(covariates)
  if (sum(ok) < 10L) {
    return(data.table(n = sum(ok), rho = NA_real_, raw_p = NA_real_,
                      adjusted_rho = NA_real_, partial_df = NA_integer_,
                      adjusted_asymptotic_p = NA_real_,
                      adjusted_permutation_p = NA_real_))
  }
  xx <- x[ok]
  yy <- y[ok]
  cc <- as.matrix(covariates[ok, , drop = FALSE])
  raw <- suppressWarnings(cor.test(xx, yy, method = "spearman", exact = FALSE))
  xr <- lm.fit(cbind(1, cc), rank(xx, ties.method = "average"))$residuals
  yr <- lm.fit(cbind(1, cc), rank(yy, ties.method = "average"))$residuals
  obs <- unname(cor(xr, yr))
  partial_df <- length(xx) - ncol(cc) - 2L
  if (partial_df <= 0L) stop("Insufficient partial-correlation degrees of freedom")
  partial_t <- obs * sqrt(partial_df / (1 - obs^2))
  ap <- 2 * pt(-abs(partial_t), df = partial_df)
  set.seed(seed)
  null <- replicate(B, cor(xr, sample(yr, replace = FALSE)))
  pp <- (1 + sum(abs(null) >= abs(obs))) / (B + 1)
  data.table(n = sum(ok), rho = unname(raw$estimate), raw_p = raw$p.value,
             adjusted_rho = obs, partial_df = partial_df,
             adjusted_asymptotic_p = ap,
             adjusted_permutation_p = pp)
}

partial_loo <- function(x, y, covariates) {
  ok <- is.finite(x) & is.finite(y) & complete.cases(covariates)
  x <- x[ok]
  y <- y[ok]
  covariates <- covariates[ok, , drop = FALSE]
  if (length(x) < 11L) return(c(min = NA, max = NA, positive_fraction = NA))
  vals <- vapply(seq_along(x), function(drop) {
    keep <- seq_along(x) != drop
    xr <- lm.fit(cbind(1, covariates[keep, , drop = FALSE]),
                 rank(x[keep], ties.method = "average"))$residuals
    yr <- lm.fit(cbind(1, covariates[keep, , drop = FALSE]),
                 rank(y[keep], ties.method = "average"))$residuals
    cor(xr, yr)
  }, numeric(1))
  c(min = min(vals), max = max(vals), positive_fraction = mean(vals > 0))
}

program_assoc_rows <- list()
exposures <- c("sender_PLAT_z", "sender_state_fraction_asin")
for (program_id in focus_programs) {
  ps <- program_scores[program == program_id]
  dt <- merge(sample_features, ps, by = "sample_id")
  for (exposure_id in exposures) {
    expression_model <- exposure_id == "sender_PLAT_z"
    cov <- cbind(
      ECM_direction = dt$ECM_direction,
      log_fibroblast_cells = log1p(dt$fibroblast_cells),
      log_myeloid_cells = log1p(dt$myeloid_cells),
      log_receiver_state_cells = log1p(dt$receiver_state_cells)
    )
    if (expression_model) {
      cov <- cbind(cov, log_sender_state_cells = log1p(dt$sender_state_cells))
    }
    tst <- partial_test(dt[[exposure_id]], dt$program_score, cov, 9999L,
                        2026072736L + match(program_id, focus_programs) * 10L +
                          match(exposure_id, exposures))
    loo <- partial_loo(dt[[exposure_id]], dt$program_score, cov)
    program_assoc_rows[[paste(program_id, exposure_id, sep = "::")]] <- cbind(
      data.table(exposure = exposure_id, receiver_state = receiver_state,
                 program = program_id),
      tst,
      data.table(loo_adjusted_rho_min = loo[["min"]],
                 loo_adjusted_rho_max = loo[["max"]],
                 loo_positive_fraction = loo[["positive_fraction"]])
    )
  }
}
program_associations <- rbindlist(program_assoc_rows)
program_associations[, adjusted_permutation_FDR_BH :=
                       p.adjust(adjusted_permutation_p, "BH"), by = exposure]
fwrite(program_associations,
       file.path(out_dir, "GSE183904_fine_state_PLAT_program_associations.tsv"),
       sep = "\t")

state_endpoint_rows <- list()
state_endpoints <- c("receiver_state_fraction_asin", "receiver_triple_fraction_asin")
for (exposure_id in exposures) {
  for (endpoint_id in state_endpoints) {
    dt <- copy(sample_features)
    cov <- cbind(
      ECM_direction = dt$ECM_direction,
      log_fibroblast_cells = log1p(dt$fibroblast_cells),
      log_myeloid_cells = log1p(dt$myeloid_cells)
    )
    if (exposure_id == "sender_PLAT_z") {
      cov <- cbind(cov, log_sender_state_cells = log1p(dt$sender_state_cells))
    }
    if (endpoint_id == "receiver_triple_fraction_asin") {
      cov <- cbind(cov, log_receiver_state_cells = log1p(dt$receiver_state_cells))
    }
    tst <- partial_test(dt[[exposure_id]], dt[[endpoint_id]], cov, 9999L,
                        2026072836L + match(exposure_id, exposures) * 10L +
                          match(endpoint_id, state_endpoints))
    loo <- partial_loo(dt[[exposure_id]], dt[[endpoint_id]], cov)
    state_endpoint_rows[[paste(exposure_id, endpoint_id, sep = "::")]] <- cbind(
      data.table(exposure = exposure_id, endpoint = endpoint_id),
      tst,
      data.table(loo_adjusted_rho_min = loo[["min"]],
                 loo_adjusted_rho_max = loo[["max"]],
                 loo_positive_fraction = loo[["positive_fraction"]])
    )
  }
}
state_endpoints_result <- rbindlist(state_endpoint_rows)
state_endpoints_result[, adjusted_permutation_FDR_BH :=
                         p.adjust(adjusted_permutation_p, "BH"), by = exposure]
fwrite(state_endpoints_result,
       file.path(out_dir, "GSE183904_fine_state_PLAT_receiver_state_associations.tsv"),
       sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE183904", sender_state = sender_state, receiver_state = receiver_state,
  primary_exposure = "TMM logCPM PLAT within cluster-defined sender state; state defined without PLAT",
  receiver_state_definition = "patient-centered cluster defined without ANXA2/ITGAM/ITGB2/LRP1",
  patients = 26L, program_family = focus_programs,
  adjustment = c("ECM direction", "log lineage cell counts",
                 "log sender-state cells for sender-expression models",
                 "log receiver-state cells for receiver-program models"),
  inference_unit = "patient/sample",
  interpretation = paste(
    "Exploratory fine-state follow-up in the same cohort used for coarse replication.",
    "Residual-permutation P values are approximate Monte Carlo sensitivity values.",
    "Association and transcript co-detection do not establish secretion, receptor assembly,",
    "ligand binding, signaling direction or causality."
  ),
  input_sha256 = as.list(setNames(vapply(inputs, digest, character(1),
                                          algo = "sha256", file = TRUE),
                                  basename(inputs)))
)
write_json(audit, file.path(audit_dir, "GSE183904_PLAT_fine_state_axis.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nFine-state PLAT program associations:\n")
print(program_associations[order(exposure, adjusted_permutation_p)])
cat("\nFine-state PLAT receiver-state associations:\n")
print(state_endpoints_result[order(exposure, adjusted_permutation_p)])
cat("\nReceiver-state co-detection summary:\n")
print(codetection[, .(patients = .N, median_state_cells = median(n_receiver_state_cells),
                      median_triple_fraction = median(receiver_triple_fraction),
                      median_ANXA2_LRP1_fraction = median(receiver_ANXA2_LRP1_fraction))])
