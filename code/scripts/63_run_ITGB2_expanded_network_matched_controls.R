#!/usr/bin/env Rscript

# Post-hoc exploratory calibration of ITGB2 against 60 controls frozen by
# 62_freeze_ITGB2_expanded_network_matched_controls.py. The target, programmes,
# scTenifoldKnk parameters and seed are unchanged from the v2 receptor screen.

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
input_dir <- file.path(project, "analysis", "processed", "virtual_knockout",
                       "GSE183904", analysis_label)
result_dir <- file.path(project, "analysis", "results", "virtual_knockout",
                        "GSE183904", analysis_label)
old_control_dir <- file.path(result_dir,
  "primary_receptor_network_matched_controls_seed2026072739")
expanded_dir <- file.path(result_dir,
  "ITGB2_expanded_60_network_matched_controls_seed2026072739_posthoc")
resource_dir <- file.path(project, "resources", "gene_sets")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(expanded_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(input_dir, "GSE183904_C1QC_APOE_sanity_counts.rds")
program_file <- file.path(resource_dir,
  "GSE183904_C1QC_APOE_receptor_screen_v2_scTenifoldKnk_programs_axis_excluded.tsv")
control_file <- file.path(resource_dir,
  "GSE183904_ITGB2_expanded_60_network_matched_controls_posthoc.tsv")
x <- readRDS(matrix_file)
program_manifest <- fread(program_file)
control_manifest <- fread(control_file)

target <- "ITGB2"
seed <- 2026072739L
all_screen_targets <- c("CD74", "CCRL2", "CXCR4", "ITGB2", "LRP1", "ITGAM",
                        "CMKLR1", "ANXA2", "TREM2")
axis_genes <- unique(c(all_screen_targets, "PLAT", "ILK"))
controls <- control_manifest[order(expanded_control_rank), gene_symbol]
programs <- split(program_manifest$gene_symbol, program_manifest$program)

if (length(controls) != 60L || anyDuplicated(controls) ||
    !identical(control_manifest$expanded_control_rank, seq_len(60L)) ||
    !all(controls %in% rownames(x))) {
  stop("Expanded ITGB2 control manifest failed integrity checks")
}

parameters <- list(
  qc = FALSE, nc_lambda = 0, nc_nNet = 3L, nc_nCells = 350L,
  nc_nComp = 3L, nc_scaleScores = TRUE, nc_symmetric = FALSE,
  nc_q = 0.90, td_K = 3L, td_maxIter = 1000L, td_maxError = 1e-5,
  td_nDecimal = 3L, ma_nDim = 2L, nCores = 4L
)

fit_path <- function(gene) {
  old <- file.path(old_control_dir, paste0("scTenifoldKnk_control__", gene, ".rds"))
  new <- file.path(expanded_dir, paste0("scTenifoldKnk_control__", gene, ".rds"))
  if (file.exists(old)) return(old)
  new
}

enrich_one <- function(dr, ko_gene, role) {
  excluded <- unique(c(axis_genes, ko_gene))
  dr <- dr[!gene %chin% excluded]
  dr[, distance_percentile := frank(distance, ties.method = "average") / .N]
  rows <- lapply(names(programs), function(program_id) {
    genes <- setdiff(unique(programs[[program_id]]), excluded)
    in_set <- dr$gene %in% genes
    n_in <- sum(in_set)
    n_out <- sum(!in_set)
    if (n_in < 5L || n_out < 5L) stop("Degenerate programme during calibration")
    wt <- suppressWarnings(wilcox.test(
      dr$distance[in_set], dr$distance[!in_set],
      alternative = "greater", exact = FALSE
    ))
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

target_file <- file.path(result_dir,
  paste0("scTenifoldKnk_seed", seed, "__", target, ".rds"))
target_fit <- readRDS(target_file)
wt_reference <- as.matrix(target_fit$tensorNetworks$WT)
target_enrichment <- enrich_one(as.data.table(target_fit$diffRegulation), target, "target")
rm(target_fit)

controls_to_run <- controls
fit_only_worker <- FALSE
worker_label <- NULL
if (length(args) >= 1L && nzchar(args[[1L]])) {
  bounds <- as.integer(strsplit(args[[1L]], ":", fixed = TRUE)[[1L]])
  if (length(bounds) != 2L || anyNA(bounds) || bounds[[1L]] < 1L ||
      bounds[[2L]] > length(controls) || bounds[[1L]] > bounds[[2L]]) {
    stop("Worker range must be START:END within 1:60")
  }
  controls_to_run <- controls[seq.int(bounds[[1L]], bounds[[2L]])]
  fit_only_worker <- TRUE
  worker_label <- paste0(bounds[[1L]], "_", bounds[[2L]])
}

control_rows <- list()
run_rows <- list()
for (control in controls_to_run) {
  message("Running/reading expanded ITGB2 control ", control)
  existing_path <- fit_path(control)
  if (file.exists(existing_path)) {
    fit <- readRDS(existing_path)
    elapsed <- NA_real_
    newly_run <- FALSE
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
    newly_run <- TRUE
  }
  wt_difference <- max(abs(as.matrix(fit$tensorNetworks$WT) - wt_reference))
  if (!is.finite(wt_difference) || wt_difference > 1e-12) {
    stop("Matched-control WT differs from target WT: ", control)
  }
  dr <- as.data.table(fit$diffRegulation)
  control_rows[[control]] <- enrich_one(dr, control, "matched_control")
  run_rows[[control]] <- data.table(
    control_gene = control,
    elapsed_seconds = elapsed,
    newly_run = newly_run,
    fit_file = existing_path,
    wt_max_abs_difference_vs_target = wt_difference
  )
  if (newly_run) saveRDS(fit, existing_path, compress = "xz")
  rm(fit, dr)
  gc(verbose = FALSE)
}

run_summary <- rbindlist(run_rows)
if (fit_only_worker) {
  fwrite(run_summary, file.path(expanded_dir,
    paste0("worker_", worker_label, "_run_summary.tsv")), sep = "\t")
  cat("Worker", worker_label, "completed", length(controls_to_run), "controls\n")
  quit(save = "no", status = 0L)
}

control_enrichment <- rbindlist(control_rows)
comparison_rows <- lapply(names(programs), function(program_id) {
  target_row <- target_enrichment[program == program_id]
  null <- control_enrichment[program == program_id]
  if (nrow(target_row) != 1L || nrow(null) != 60L) stop("Comparison row mismatch")
  data.table(
    target = target, seed = seed, program = program_id,
    target_auc = target_row$enrichment_auc,
    matched_control_median_auc = median(null$enrichment_auc),
    matched_control_minimum_auc = min(null$enrichment_auc),
    matched_control_maximum_auc = max(null$enrichment_auc),
    target_minus_control_median_auc =
      target_row$enrichment_auc - median(null$enrichment_auc),
    controls_at_least_target = sum(null$enrichment_auc >= target_row$enrichment_auc),
    target_empirical_upper_tail_p =
      (1 + sum(null$enrichment_auc >= target_row$enrichment_auc)) / (1 + nrow(null)),
    target_percentile_among_61 =
      rank(c(null$enrichment_auc, target_row$enrichment_auc),
           ties.method = "average")[61] / 61
  )
})
comparison <- rbindlist(comparison_rows)
comparison[, empirical_FDR_BH_six_programmes :=
             p.adjust(target_empirical_upper_tail_p, "BH")]

fwrite(target_enrichment, file.path(expanded_dir,
  "GSE183904_ITGB2_expanded_target_program_enrichment.tsv"), sep = "\t")
fwrite(control_enrichment, file.path(expanded_dir,
  "GSE183904_ITGB2_expanded_matched_control_program_enrichment.tsv"), sep = "\t")
fwrite(run_summary, file.path(expanded_dir,
  "GSE183904_ITGB2_expanded_matched_control_run_summary.tsv"), sep = "\t")
fwrite(comparison, file.path(expanded_dir,
  "GSE183904_ITGB2_vs_expanded_60_network_matched_controls.tsv"), sep = "\t")

manifest <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  analysis_status = "post_hoc_exploratory",
  target = target,
  seed = seed,
  controls = controls,
  controls_frozen_before_expanded_KO_outcomes = TRUE,
  controls_per_target = length(controls),
  parameters = parameters,
  multiplicity = "BH across six empirical P values",
  empirical_resolution = "60 matched controls; minimum upper-tail P = 1/61",
  differential_regulation_semantics = "unsigned manifold distance",
  input_sha256 = list(
    matrix = digest(matrix_file, algo = "sha256", file = TRUE),
    programmes = digest(program_file, algo = "sha256", file = TRUE),
    controls = digest(control_file, algo = "sha256", file = TRUE)
  ),
  interpretation = paste(
    "Post-hoc expanded controls improve empirical-null resolution only.",
    "They do not establish ligand binding, receptor activation or causality."
  )
)
write_json(manifest, file.path(audit_dir,
  "GSE183904_ITGB2_expanded_60_matched_controls_posthoc_results.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)

print(comparison[order(target_empirical_upper_tail_p)])
