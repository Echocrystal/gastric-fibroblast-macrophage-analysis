#!/usr/bin/env Rscript

# Same-seed network/expression-matched control KOs for the exploratory ITGB2
# scTenifoldKnk result. All program scores exclude the perturbed gene.

options(stringsAsFactors = FALSE, warn = 1)

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

input_dir <- file.path(project, "analysis", "processed", "virtual_knockout", "GSE183904",
                       "C1QC_APOE_macrophage_v1")
result_dir <- file.path(project, "analysis", "results", "virtual_knockout", "GSE183904",
                        "C1QC_APOE_macrophage_v1")
control_dir <- file.path(result_dir, "ITGB2_network_matched_controls_seed2026072739")
resource_dir <- file.path(project, "resources", "gene_sets")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(control_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(input_dir, "GSE183904_C1QC_APOE_sanity_counts.rds")
program_file <- file.path(resource_dir,
                          "GSE183904_C1QC_APOE_scTenifoldKnk_programs_axis_excluded.tsv")
control_file <- file.path(resource_dir,
                          "GSE183904_ITGB2_scTenifoldKnk_network_matched_controls.tsv")
x <- readRDS(matrix_file)
program_manifest <- fread(program_file)
control_manifest <- fread(control_file)

target <- "ITGB2"
controls <- control_manifest$control_gene
seed <- 2026072739L
axis_genes <- c("ANXA2", "ITGAM", "ITGB2", "LRP1", "PLAT", "ILK")
programs <- split(program_manifest$gene_symbol, program_manifest$program)
if (length(controls) != 20L || anyDuplicated(controls) || !all(controls %in% rownames(x))) {
  stop("Frozen ITGB2 control set changed")
}

parameters <- list(
  qc = FALSE,
  nc_lambda = 0,
  nc_nNet = 3L,
  nc_nCells = 350L,
  nc_nComp = 3L,
  nc_scaleScores = TRUE,
  nc_symmetric = FALSE,
  nc_q = 0.90,
  td_K = 3L,
  td_maxIter = 1000L,
  td_maxError = 1e-5,
  td_nDecimal = 3L,
  ma_nDim = 2L,
  nCores = 4L
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
    wt <- suppressWarnings(wilcox.test(dr$distance[in_set], dr$distance[!in_set],
                                       alternative = "greater", exact = FALSE))
    data.table(
      ko_gene = ko_gene,
      role = role,
      program = program_id,
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

target_fit_file <- file.path(result_dir,
                             paste0("scTenifoldKnk_seed", seed, "__", target, ".rds"))
target_fit <- readRDS(target_fit_file)
wt_reference <- as.matrix(target_fit$tensorNetworks$WT)
target_dr <- as.data.table(target_fit$diffRegulation)
target_enrichment <- enrich_one(target_dr, target, "target")
rm(target_fit)

control_enrichment_rows <- list()
control_summary_rows <- list()
wt_differences <- numeric()
for (control in controls) {
  message("Running matched control ", control)
  fit_file <- file.path(control_dir, paste0("scTenifoldKnk_control__", control, ".rds"))
  if (file.exists(fit_file)) {
    message("Resuming completed fit: ", basename(fit_file))
    fit <- readRDS(fit_file)
    elapsed <- NA_real_
  } else {
    set.seed(seed)
    started <- Sys.time()
    fit <- scTenifoldKnk::scTenifoldKnk(
      countMatrix = x,
      qc = parameters$qc,
      gKO = control,
      nc_lambda = parameters$nc_lambda,
      nc_nNet = parameters$nc_nNet,
      nc_nCells = parameters$nc_nCells,
      nc_nComp = parameters$nc_nComp,
      nc_scaleScores = parameters$nc_scaleScores,
      nc_symmetric = parameters$nc_symmetric,
      nc_q = parameters$nc_q,
      td_K = parameters$td_K,
      td_maxIter = parameters$td_maxIter,
      td_maxError = parameters$td_maxError,
      td_nDecimal = parameters$td_nDecimal,
      ma_nDim = parameters$ma_nDim,
      nCores = parameters$nCores
    )
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  }
  wt <- as.matrix(fit$tensorNetworks$WT)
  wt_difference <- max(abs(wt - wt_reference))
  if (!is.finite(wt_difference) || wt_difference > 1e-12) {
    stop("Matched control WT differs from target WT: ", control)
  }
  dr <- as.data.table(fit$diffRegulation)
  control_enrichment_rows[[control]] <- enrich_one(dr, control, "matched_control")
  control_row <- dr[gene == control]
  control_summary_rows[[control]] <- data.table(
    control_gene = control,
    elapsed_seconds = elapsed,
    wt_max_abs_difference_vs_target = wt_difference,
    target_distance_rank = frank(-dr$distance, ties.method = "min")[which(dr$gene == control)],
    target_distance = control_row$distance,
    n_diffreg_FDR_0_05 = sum(dr$p.adj < 0.05, na.rm = TRUE)
  )
  wt_differences <- c(wt_differences, wt_difference)
  if (!file.exists(fit_file)) saveRDS(fit, fit_file, compress = "xz")
  rm(fit, wt, dr)
  gc(verbose = FALSE)
}

control_enrichment <- rbindlist(control_enrichment_rows)
control_summary <- rbindlist(control_summary_rows)
comparison_rows <- lapply(names(programs), function(program_id) {
  target_row <- target_enrichment[program == program_id]
  null <- control_enrichment[program == program_id]
  if (nrow(target_row) != 1L || nrow(null) != 20L) stop("Comparison row mismatch")
  data.table(
    target = target,
    seed = seed,
    program = program_id,
    target_auc = target_row$enrichment_auc,
    matched_control_median_auc = median(null$enrichment_auc),
    matched_control_minimum_auc = min(null$enrichment_auc),
    matched_control_maximum_auc = max(null$enrichment_auc),
    target_minus_control_median_auc = target_row$enrichment_auc - median(null$enrichment_auc),
    controls_at_least_target = sum(null$enrichment_auc >= target_row$enrichment_auc),
    target_empirical_upper_tail_p =
      (1 + sum(null$enrichment_auc >= target_row$enrichment_auc)) / (1 + nrow(null)),
    target_percentile_among_21 =
      rank(c(null$enrichment_auc, target_row$enrichment_auc), ties.method = "average")[21] / 21
  )
})
comparison <- rbindlist(comparison_rows)
comparison[, empirical_FDR_BH := p.adjust(target_empirical_upper_tail_p, "BH")]

fwrite(target_enrichment,
       file.path(control_dir, "GSE183904_ITGB2_target_program_enrichment.tsv"), sep = "\t")
fwrite(control_enrichment,
       file.path(control_dir, "GSE183904_ITGB2_matched_control_program_enrichment.tsv"), sep = "\t")
fwrite(control_summary,
       file.path(control_dir, "GSE183904_ITGB2_matched_control_run_summary.tsv"), sep = "\t")
fwrite(comparison,
       file.path(control_dir, "GSE183904_ITGB2_vs_network_matched_controls.tsv"), sep = "\t")

manifest <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE183904",
  state = "C1QC_APOE_macrophage",
  target = target,
  seed = seed,
  controls = controls,
  controls_frozen_without_KO_outcomes = TRUE,
  same_seed_common_WT = all(wt_differences <= 1e-12),
  parameters = parameters,
  comparison = paste(
    "Upper-tail empirical P uses the target plus 20 frozen network/expression/program-membership-matched controls;",
    "the perturbed gene and all axis genes are excluded from every program score."
  ),
  multiplicity = "BH across six target-versus-control empirical P values",
  differential_regulation_semantics = "unsigned manifold distance",
  interpretation = paste(
    "This is a one-seed matched-control specificity test following a four-target, three-seed exploratory screen.",
    "It quantifies relative model perturbation and does not establish extracellular binding, signaling direction or biological causality."
  ),
  input_sha256 = list(
    matrix = digest(matrix_file, algo = "sha256", file = TRUE),
    programs = digest(program_file, algo = "sha256", file = TRUE),
    controls = digest(control_file, algo = "sha256", file = TRUE),
    target_fit = digest(target_fit_file, algo = "sha256", file = TRUE)
  )
)
write_json(manifest,
           file.path(audit_dir, "GSE183904_ITGB2_scTenifoldKnk_matched_controls.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nITGB2 versus network-matched controls:\n")
print(comparison[order(target_empirical_upper_tail_p)])
