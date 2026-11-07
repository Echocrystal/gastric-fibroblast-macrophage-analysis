#!/usr/bin/env Rscript

# Patient-blocked characterization of the corrected recurrent C1QC+APOE+
# macrophage state and specificity analysis for its abundance association.

options(stringsAsFactors = FALSE, warn = 1)
set.seed(2026072737L)

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
label_file <- file.path(project, "analysis", "audit", "GSE183904_fine_state_label_review.tsv")
program_file <- file.path(project, "resources", "gene_sets", "myeloid_receiver_programs_frozen.tsv")
feature_file <- file.path(project, "analysis", "results", "PLAT_axis", "GSE183904",
                          "fine_state_axis", "GSE183904_fine_state_axis_patient_features.tsv")
state_count_file <- file.path(project, "analysis", "results", "PLAT_axis", "GSE183904",
                              "fine_state_axis", "GSE183904_fine_state_counts_by_patient.tsv")
out_dir <- file.path(project, "analysis", "results", "PLAT_axis", "GSE183904",
                     "receiver_state_characterization")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(manifest_file, myeloid_meta_file, label_file, program_file,
            feature_file, state_count_file)
if (!all(file.exists(inputs))) stop("Missing receiver-state characterization input")

manifest <- fread(manifest_file)[include_primary_replication == TRUE &
                                 specimen_class == "Primary_Tumor"]
myeloid <- fread(myeloid_meta_file)
labels <- fread(label_file)
program_manifest <- fread(program_file)
features <- fread(feature_file)
state_counts <- fread(state_count_file)
if (nrow(manifest) != 26L || uniqueN(manifest$Patient_ID) != 26L) {
  stop("Expected 26 patient-independent primary tumors")
}

myeloid_map <- setNames(labels[lineage == "Myeloid"]$refined_label,
                        labels[lineage == "Myeloid"]$cluster)
myeloid[, refined_label := unname(myeloid_map[as.character(cluster)])]
if (myeloid[, anyNA(refined_label)]) stop("Incomplete corrected myeloid-state mapping")

target_state <- "C1QC_APOE_macrophage"
recurrent_states <- labels[lineage == "Myeloid" & recurrent_gate == TRUE]$refined_label
comparator_states <- setdiff(recurrent_states, target_state)
if (!target_state %chin% recurrent_states || length(comparator_states) != 3L) {
  stop("Unexpected recurrent-state family")
}

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
program_genes <- sort(unique(unlist(programs)))

cell_score_rows <- list()
coverage_rows <- list()
for (i in seq_len(nrow(manifest))) {
  info <- manifest[i]
  message(sprintf("[%d/26] receiver-state programs: %s", i, info$sample_id))
  sce <- readRDS(file.path(sce_dir, paste0(info$sample_id, "_secondary_qc_sce.rds")))
  counts_mat <- counts(sce)
  cm <- myeloid[sample_id == info$sample_id]
  idx <- match(cm$cell_id, colnames(counts_mat))
  if (anyNA(idx)) stop("Cell/count alignment failure: ", info$sample_id)
  present <- intersect(program_genes, rownames(counts_mat))
  if (length(present) < 20L) stop("Insufficient program coverage: ", info$sample_id)

  lib <- Matrix::colSums(counts_mat[, idx, drop = FALSE])
  if (any(!is.finite(lib) | lib <= 0)) stop("Invalid myeloid library size: ", info$sample_id)
  norm <- counts_mat[present, idx, drop = FALSE] %*%
    Matrix::Diagonal(x = 1e4 / lib)
  norm@x <- log1p(norm@x)
  dense <- as.matrix(norm)
  mu <- rowMeans(dense)
  sigma <- apply(dense, 1L, sd)
  informative <- is.finite(sigma) & sigma > 0
  z <- matrix(0, nrow = nrow(dense), ncol = ncol(dense),
              dimnames = dimnames(dense))
  z[informative, ] <- (dense[informative, , drop = FALSE] - mu[informative]) /
    sigma[informative]

  score_dt <- cm[, .(cell_id, sample_id, Patient_ID, cluster, refined_label)]
  for (program_id in focus_programs) {
    genes <- intersect(programs[[program_id]], rownames(z)[informative])
    if (length(genes) < 5L) stop("Insufficient informative genes: ", program_id)
    score_dt[, (program_id) := colMeans(z[genes, , drop = FALSE])]
    coverage_rows[[length(coverage_rows) + 1L]] <- data.table(
      sample_id = info$sample_id,
      program = program_id,
      n_frozen_axis_excluded = uniqueN(programs[[program_id]]),
      n_present_informative = length(genes)
    )
  }
  cell_score_rows[[info$sample_id]] <- score_dt
  rm(sce, counts_mat, norm, dense, z)
  gc(verbose = FALSE)
}

cell_scores <- rbindlist(cell_score_rows)
coverage <- rbindlist(coverage_rows)
fwrite(coverage, file.path(out_dir, "GSE183904_receiver_program_coverage_by_patient.tsv"), sep = "\t")

long_scores <- melt(
  cell_scores,
  id.vars = c("cell_id", "sample_id", "Patient_ID", "cluster", "refined_label"),
  measure.vars = focus_programs,
  variable.name = "program",
  value.name = "patient_centered_cell_score"
)
patient_state <- long_scores[, .(
  n_cells = .N,
  mean_program_score = mean(patient_centered_cell_score),
  median_program_score = median(patient_centered_cell_score)
), by = .(sample_id, Patient_ID, refined_label, program)]
fwrite(patient_state,
       file.path(out_dir, "GSE183904_myeloid_program_scores_by_patient_state.tsv"),
       sep = "\t")

# Target versus all other recurrent states, paired within patient.
target_scores <- long_scores[refined_label == target_state, .(
  target_cells = .N,
  target_mean = mean(patient_centered_cell_score)
), by = .(sample_id, Patient_ID, program)]
other_scores <- long_scores[refined_label %chin% comparator_states, .(
  comparator_cells = .N,
  comparator_mean = mean(patient_centered_cell_score)
), by = .(sample_id, Patient_ID, program)]
paired <- merge(target_scores, other_scores,
                by = c("sample_id", "Patient_ID", "program"), all = FALSE)
paired <- paired[target_cells >= 5L & comparator_cells >= 20L]
paired[, delta_target_minus_other_recurrent := target_mean - comparator_mean]
fwrite(paired,
       file.path(out_dir, "GSE183904_C1QC_APOE_vs_other_recurrent_paired_differences.tsv"),
       sep = "\t")

program_contrast <- paired[, {
  d <- delta_target_minus_other_recurrent
  wt <- suppressWarnings(wilcox.test(d, mu = 0, alternative = "two.sided", exact = FALSE))
  positives <- sum(d > 0)
  negatives <- sum(d < 0)
  nonzero <- positives + negatives
  sign_p <- if (nonzero > 0L) binom.test(positives, nonzero, p = 0.5)$p.value else NA_real_
  .(
    n_patients = .N,
    median_delta = median(d),
    mean_delta = mean(d),
    positive_patients = positives,
    negative_patients = negatives,
    paired_wilcoxon_p = wt$p.value,
    exact_sign_p = sign_p
  )
}, by = program]
program_contrast[, paired_wilcoxon_FDR_BH := p.adjust(paired_wilcoxon_p, "BH")]
program_contrast[, exact_sign_FDR_BH := p.adjust(exact_sign_p, "BH")]
fwrite(program_contrast,
       file.path(out_dir, "GSE183904_C1QC_APOE_program_contrasts.tsv"), sep = "\t")

partial_test <- function(x, y, covariates, B = 9999L, seed = 1L) {
  ok <- is.finite(x) & is.finite(y) & complete.cases(covariates)
  xx <- x[ok]
  yy <- y[ok]
  cc <- as.matrix(covariates[ok, , drop = FALSE])
  if (length(xx) < 10L) {
    return(data.table(n = length(xx), raw_rho = NA_real_, raw_p = NA_real_,
                      adjusted_rho = NA_real_, adjusted_permutation_p = NA_real_))
  }
  raw <- suppressWarnings(cor.test(xx, yy, method = "spearman", exact = FALSE))
  xr <- lm.fit(cbind(1, cc), rank(xx, ties.method = "average"))$residuals
  yr <- lm.fit(cbind(1, cc), rank(yy, ties.method = "average"))$residuals
  obs <- unname(cor(xr, yr))
  set.seed(seed)
  null <- replicate(B, cor(xr, sample(yr, replace = FALSE)))
  pp <- (1 + sum(abs(null) >= abs(obs))) / (B + 1)
  data.table(n = length(xx), raw_rho = unname(raw$estimate), raw_p = raw$p.value,
             adjusted_rho = obs, adjusted_permutation_p = pp)
}

partial_loo <- function(x, y, covariates) {
  ok <- is.finite(x) & is.finite(y) & complete.cases(covariates)
  x <- x[ok]
  y <- y[ok]
  covariates <- covariates[ok, , drop = FALSE]
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

# Specificity of the sender-state abundance association across all recurrent
# myeloid states. Fractions use the same myeloid denominator and are treated
# as compositional descriptive endpoints.
myeloid_counts <- state_counts[lineage == "Myeloid" & refined_label %chin% recurrent_states]
abundance_rows <- list()
for (state_id in recurrent_states) {
  yy <- myeloid_counts[refined_label == state_id,
                       .(sample_id, receiver_fraction = state_fraction)]
  dt <- merge(features, yy, by = "sample_id", all.x = TRUE)
  dt[is.na(receiver_fraction), receiver_fraction := 0]
  y <- asin(sqrt(dt$receiver_fraction))
  cov <- cbind(
    ECM_direction = dt$ECM_direction,
    log_fibroblast_cells = log1p(dt$fibroblast_cells),
    log_myeloid_cells = log1p(dt$myeloid_cells)
  )
  tst <- partial_test(dt$sender_state_fraction_asin, y, cov, 9999L,
                      2026072737L + match(state_id, recurrent_states))
  loo <- partial_loo(dt$sender_state_fraction_asin, y, cov)
  abundance_rows[[state_id]] <- cbind(
    data.table(receiver_state = state_id), tst,
    data.table(
      loo_adjusted_rho_min = loo[["min"]],
      loo_adjusted_rho_max = loo[["max"]],
      loo_positive_fraction = loo[["positive_fraction"]]
    )
  )
}
abundance_specificity <- rbindlist(abundance_rows)
abundance_specificity[, adjusted_permutation_FDR_BH :=
                        p.adjust(adjusted_permutation_p, "BH")]
fwrite(abundance_specificity,
       file.path(out_dir, "GSE183904_sender_abundance_vs_recurrent_myeloid_states.tsv"),
       sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE183904",
  target_state = target_state,
  comparator_states = comparator_states,
  program_scoring = paste(
    "log1p(CP10K) using full-library size; gene-wise z score within patient across all myeloid cells;",
    "cell-level program mean; patient-state mean; axis genes excluded"
  ),
  program_contrast = "paired target-state minus pooled other recurrent-state score within patient",
  minimum_cells = list(target_state = 5L, pooled_comparator = 20L),
  abundance_specificity = paste(
    "arcsin-sqrt state fractions among myeloid cells; rank residual correlation;",
    "ECM direction and lineage cell-count adjustment; 9999 residual permutations"
  ),
  inference_unit = "patient/sample",
  interpretation = paste(
    "Program contrasts characterize a cluster discovered from the same transcriptomes and are descriptive.",
    "State-fraction associations are compositional and do not establish direction or causality."
  ),
  input_sha256 = as.list(setNames(vapply(inputs, digest, character(1),
                                        algo = "sha256", file = TRUE), basename(inputs)))
)
write_json(audit,
           file.path(audit_dir, "GSE183904_C1QC_APOE_receiver_characterization.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nC1QC+APOE+ program contrasts:\n")
print(program_contrast[order(paired_wilcoxon_p)])
cat("\nSender-state abundance specificity:\n")
print(abundance_specificity[order(adjusted_permutation_p)])
