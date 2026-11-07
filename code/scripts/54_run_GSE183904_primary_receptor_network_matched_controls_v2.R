#!/usr/bin/env Rscript

# Run one prespecified seed for all unique controls frozen for the four primary
# communication-derived receptors, then calibrate each target against its own
# 20 matched controls.

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(
  file.path(project, "runtime", "R_virtual_knockout"),
  "I:/ECM-new/scRNA/20.rebuild_GSE163558/R_libs_4.4",
  file.path(source_project, "runtime", "R_library"),
  .libPaths()
))

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(scTenifoldKnk)
  library(scTenifoldNet)
  library(jsonlite)
  library(digest)
})

analysis_label <- "C1QC_APOE_macrophage_receptor_screen_v2"
input_dir <- file.path(project, "analysis", "processed", "virtual_knockout", "GSE183904",
                       analysis_label)
result_dir <- file.path(project, "analysis", "results", "virtual_knockout", "GSE183904",
                        analysis_label)
control_dir <- file.path(result_dir, "primary_receptor_network_matched_controls_seed2026072739")
resource_dir <- file.path(project, "resources", "gene_sets")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(control_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(input_dir, "GSE183904_C1QC_APOE_sanity_counts.rds")
program_file <- file.path(resource_dir,
  "GSE183904_C1QC_APOE_receptor_screen_v2_scTenifoldKnk_programs_axis_excluded.tsv")
control_file <- file.path(resource_dir,
  "GSE183904_primary_receptors_v2_scTenifoldKnk_network_matched_controls.tsv")
x <- readRDS(matrix_file)
program_manifest <- fread(program_file)
control_manifest <- fread(control_file)

primary_targets <- c("CD74", "CCRL2", "CXCR4", "ITGB2")
all_screen_targets <- c(primary_targets, "LRP1", "ITGAM", "CMKLR1", "ANXA2", "TREM2")
axis_genes <- unique(c(all_screen_targets, "PLAT", "ILK"))
controls <- unique(control_manifest$control_gene)
controls_to_run <- controls
fit_only_worker <- FALSE
worker_label <- NULL
if (length(args) >= 1L && nzchar(args[[1L]])) {
  bounds <- as.integer(strsplit(args[[1L]], ":", fixed = TRUE)[[1L]])
  if (length(bounds) != 2L || anyNA(bounds) || bounds[[1L]] < 1L ||
      bounds[[2L]] > length(controls) || bounds[[1L]] > bounds[[2L]]) {
    stop("Worker range must be START:END within the frozen unique-control order")
  }
  controls_to_run <- controls[seq.int(bounds[[1L]], bounds[[2L]])]
  fit_only_worker <- TRUE
  worker_label <- paste0(bounds[[1L]], "_", bounds[[2L]])
}
seed <- 2026072739L
programs <- split(program_manifest$gene_symbol, program_manifest$program)
if (!all(primary_targets %chin% control_manifest$target) ||
    control_manifest[, uniqueN(control_gene), by = target][, any(V1 != 20L)] ||
    anyDuplicated(control_manifest, by = c("target", "control_gene")) ||
    !all(controls %in% rownames(x))) {
  stop("Frozen primary-receptor control manifest changed")
}

parameters <- list(
  qc = FALSE, nc_lambda = 0, nc_nNet = 3L, nc_nCells = 350L,
  nc_nComp = 3L, nc_scaleScores = TRUE, nc_symmetric = FALSE,
  nc_q = 0.90, td_K = 3L, td_maxIter = 1000L, td_maxError = 1e-5,
  td_nDecimal = 3L, ma_nDim = 2L, nCores = 4L
)

enrich_one <- function(dr, ko_gene, role) {
  excluded <- unique(c(axis_genes, ko_gene))
  dr <- dr[!gene %in% excluded]
  dr[, distance_percentile := frank(distance, ties.method = "average") / .N]
  rows <- lapply(names(programs), function(program_id) {
    genes <- setdiff(unique(programs[[program_id]]), excluded)
    in_set <- dr$gene %in% genes
    n_in <- sum(in_set)
    n_out <- sum(!in_set)
    if (n_in < 5L || n_out < 5L) stop("Degenerate program during control calibration")
    wt <- suppressWarnings(wilcox.test(dr$distance[in_set], dr$distance[!in_set],
                                       alternative = "greater", exact = FALSE))
    data.table(
      ko_gene = ko_gene, role = role, program = program_id,
      n_in_universe = n_in,
      mean_distance_percentile = mean(dr$distance_percentile[in_set]),
      enrichment_auc = as.numeric(wt$statistic) / (n_in * n_out),
      wilcoxon_p = wt$p.value
    )
  })
  out <- rbindlist(rows)
  out[, wilcoxon_FDR_within_KO := p.adjust(wilcoxon_p, "BH")]
  out
}

reference_file <- file.path(result_dir,
  paste0("scTenifoldKnk_seed", seed, "__", primary_targets[[1L]], ".rds"))
reference_fit <- readRDS(reference_file)
wt_reference <- as.matrix(reference_fit$tensorNetworks$WT)
rm(reference_fit)

target_rows <- list()
for (target in primary_targets) {
  fit_file <- file.path(result_dir, paste0("scTenifoldKnk_seed", seed, "__", target, ".rds"))
  fit <- readRDS(fit_file)
  if (max(abs(as.matrix(fit$tensorNetworks$WT) - wt_reference)) > 1e-12) {
    stop("Primary target WT differs: ", target)
  }
  target_rows[[target]] <- enrich_one(as.data.table(fit$diffRegulation), target, "target")
  rm(fit)
}
target_enrichment <- rbindlist(target_rows)

control_rows <- list()
control_summary_rows <- list()
for (control in controls_to_run) {
  message("Running matched control ", control)
  fit_file <- file.path(control_dir, paste0("scTenifoldKnk_control__", control, ".rds"))
  if (file.exists(fit_file)) {
    fit <- readRDS(fit_file)
    elapsed <- NA_real_
  } else {
    set.seed(seed)
    started <- Sys.time()
    fit <- scTenifoldKnk::scTenifoldKnk(
      countMatrix = x, qc = parameters$qc, gKO = control,
      nc_lambda = parameters$nc_lambda, nc_nNet = parameters$nc_nNet,
      nc_nCells = parameters$nc_nCells, nc_nComp = parameters$nc_nComp,
      nc_scaleScores = parameters$nc_scaleScores,
      nc_symmetric = parameters$nc_symmetric, nc_q = parameters$nc_q,
      td_K = parameters$td_K, td_maxIter = parameters$td_maxIter,
      td_maxError = parameters$td_maxError, td_nDecimal = parameters$td_nDecimal,
      ma_nDim = parameters$ma_nDim, nCores = parameters$nCores
    )
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  }
  wt_difference <- max(abs(as.matrix(fit$tensorNetworks$WT) - wt_reference))
  if (!is.finite(wt_difference) || wt_difference > 1e-12) {
    stop("Matched control WT differs from primary-target WT: ", control)
  }
  dr <- as.data.table(fit$diffRegulation)
  control_rows[[control]] <- enrich_one(dr, control, "matched_control")
  control_row <- dr[gene == control]
  control_summary_rows[[control]] <- data.table(
    control_gene = control, elapsed_seconds = elapsed,
    wt_max_abs_difference_vs_target = wt_difference,
    target_distance_rank = frank(-dr$distance, ties.method = "min")[which(dr$gene == control)],
    target_distance = control_row$distance,
    n_diffreg_FDR_0_05 = sum(dr$p.adj < 0.05, na.rm = TRUE)
  )
  if (!file.exists(fit_file)) saveRDS(fit, fit_file, compress = "xz")
  rm(fit, dr)
  gc(verbose = FALSE)
}
control_enrichment <- rbindlist(control_rows)
control_summary <- rbindlist(control_summary_rows)

if (fit_only_worker) {
  fwrite(control_summary, file.path(control_dir,
    paste0("worker_", worker_label, "_run_summary.tsv")), sep = "\t")
  cat("Worker", worker_label, "completed", length(controls_to_run), "fits\n")
  quit(save = "no", status = 0L)
}

comparison_rows <- list()
for (target in primary_targets) {
  target_id <- target
  target_controls <- control_manifest[target == target_id, control_gene]
  for (program_id in names(programs)) {
    target_row <- target_enrichment[ko_gene == target & program == program_id]
    null <- control_enrichment[ko_gene %chin% target_controls & program == program_id]
    if (nrow(target_row) != 1L || nrow(null) != 20L) stop("Comparison row mismatch")
    key <- paste(target, program_id, sep = "__")
    comparison_rows[[key]] <- data.table(
      target = target, seed = seed, program = program_id,
      target_auc = target_row$enrichment_auc,
      matched_control_median_auc = median(null$enrichment_auc),
      matched_control_minimum_auc = min(null$enrichment_auc),
      matched_control_maximum_auc = max(null$enrichment_auc),
      target_minus_control_median_auc = target_row$enrichment_auc - median(null$enrichment_auc),
      controls_at_least_target = sum(null$enrichment_auc >= target_row$enrichment_auc),
      target_empirical_upper_tail_p =
        (1 + sum(null$enrichment_auc >= target_row$enrichment_auc)) / (1 + nrow(null)),
      target_percentile_among_21 =
        rank(c(null$enrichment_auc, target_row$enrichment_auc),
             ties.method = "average")[21] / 21
    )
  }
}
comparison <- rbindlist(comparison_rows)
comparison[, empirical_FDR_BH_within_target :=
             p.adjust(target_empirical_upper_tail_p, "BH"), by = target]

fwrite(target_enrichment, file.path(control_dir,
  "GSE183904_primary_receptor_target_program_enrichment.tsv"), sep = "\t")
fwrite(control_enrichment, file.path(control_dir,
  "GSE183904_primary_receptor_matched_control_program_enrichment.tsv"), sep = "\t")
fwrite(control_summary, file.path(control_dir,
  "GSE183904_primary_receptor_matched_control_run_summary.tsv"), sep = "\t")
fwrite(comparison, file.path(control_dir,
  "GSE183904_primary_receptors_vs_network_matched_controls.tsv"), sep = "\t")

manifest <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE183904", state = "C1QC_APOE_macrophage",
  targets = primary_targets, seed = seed, unique_controls = controls,
  controls_frozen_without_v2_KO_outcomes = TRUE,
  same_seed_common_WT = all(control_summary$wt_max_abs_difference_vs_target <= 1e-12),
  parameters = parameters,
  multiplicity = "BH across six empirical P values within each target",
  empirical_resolution = "20 matched controls; minimum upper-tail P = 1/21",
  differential_regulation_semantics = "unsigned manifold distance",
  input_sha256 = list(
    matrix = digest(matrix_file, algo = "sha256", file = TRUE),
    programs = digest(program_file, algo = "sha256", file = TRUE),
    controls = digest(control_file, algo = "sha256", file = TRUE)
  ),
  interpretation = paste(
    "Matched controls calibrate relative model perturbation at one prespecified seed.",
    "They do not establish ligand binding, receptor activation, PLAT dependence or causality."
  )
)
write_json(manifest, file.path(audit_dir,
  "GSE183904_primary_receptors_v2_scTenifoldKnk_matched_controls.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nPrimary receptors versus network-matched controls:\n")
print(comparison[order(target, target_empirical_upper_tail_p)])
