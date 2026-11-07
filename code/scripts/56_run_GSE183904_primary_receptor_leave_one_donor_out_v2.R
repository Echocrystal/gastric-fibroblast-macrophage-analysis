#!/usr/bin/env Rscript

# Leave-one-donor-out sensitivity for the four communication-derived receptors
# at one prespecified algorithm seed. Optional START:END selects donor indices
# for independent fit-only workers; a no-argument run resumes all fits and
# writes the combined summaries.

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
lodo_dir <- file.path(result_dir, "primary_receptor_leave_one_donor_out_seed2026072739")
resource_dir <- file.path(project, "resources", "gene_sets")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(lodo_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(input_dir, "GSE183904_C1QC_APOE_sanity_counts.rds")
metadata_file <- file.path(input_dir, "GSE183904_C1QC_APOE_sanity_metadata.tsv")
program_file <- file.path(resource_dir,
  "GSE183904_C1QC_APOE_receptor_screen_v2_scTenifoldKnk_programs_axis_excluded.tsv")
x <- readRDS(matrix_file)
metadata <- fread(metadata_file)
program_manifest <- fread(program_file)
if (!identical(colnames(x), metadata$cell_id)) stop("Count/metadata order mismatch")

primary_targets <- c("CD74", "CCRL2", "CXCR4", "ITGB2")
all_screen_targets <- c(primary_targets, "LRP1", "ITGAM", "CMKLR1", "ANXA2", "TREM2")
axis_genes <- unique(c(all_screen_targets, "PLAT", "ILK"))
programs <- split(program_manifest$gene_symbol, program_manifest$program)
core_programs <- c(
  "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
  "GOBP_PHAGOCYTOSIS",
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX"
)
patients <- unique(metadata$Patient_ID)
if (length(patients) != 26L) stop("Patient count changed")
patient_indices <- seq_along(patients)
fit_only_worker <- FALSE
worker_label <- NULL
if (length(args) >= 1L && nzchar(args[[1L]])) {
  bounds <- as.integer(strsplit(args[[1L]], ":", fixed = TRUE)[[1L]])
  if (length(bounds) != 2L || anyNA(bounds) || bounds[[1L]] < 1L ||
      bounds[[2L]] > length(patients) || bounds[[1L]] > bounds[[2L]]) {
    stop("Worker range must be START:END within the frozen patient order")
  }
  patient_indices <- seq.int(bounds[[1L]], bounds[[2L]])
  fit_only_worker <- TRUE
  worker_label <- paste0(bounds[[1L]], "_", bounds[[2L]])
}

seed <- 2026072739L
parameters <- list(
  qc = FALSE, nc_lambda = 0, nc_nNet = 3L, nc_nCells = 350L,
  nc_nComp = 3L, nc_scaleScores = TRUE, nc_symmetric = FALSE,
  nc_q = 0.90, td_K = 3L, td_maxIter = 1000L, td_maxError = 1e-5,
  td_nDecimal = 3L, ma_nDim = 2L, nCores = 4L
)

enrich_one <- function(dr, target, held_out_patient) {
  excluded <- unique(c(axis_genes, target))
  dr <- dr[!gene %in% excluded]
  dr[, distance_percentile := frank(distance, ties.method = "average") / .N]
  rows <- lapply(names(programs), function(program_id) {
    genes <- setdiff(unique(programs[[program_id]]), excluded)
    in_set <- dr$gene %in% genes
    n_in <- sum(in_set)
    n_out <- sum(!in_set)
    if (n_in < 5L || n_out < 5L) stop("Degenerate LODO program")
    wt <- suppressWarnings(wilcox.test(dr$distance[in_set], dr$distance[!in_set],
                                       alternative = "greater", exact = FALSE))
    data.table(
      held_out_patient = held_out_patient, target = target,
      program = program_id, n_in_universe = n_in,
      enrichment_auc = as.numeric(wt$statistic) / (n_in * n_out),
      wilcoxon_p = wt$p.value
    )
  })
  out <- rbindlist(rows)
  out[, wilcoxon_FDR_within_target_LODO := p.adjust(wilcoxon_p, "BH")]
  out
}

program_rows <- list()
run_rows <- list()
for (patient_index in patient_indices) {
  held_out <- patients[[patient_index]]
  keep <- metadata$Patient_ID != held_out
  x_lodo <- x[, keep, drop = FALSE]
  safe_patient <- gsub("[^A-Za-z0-9_.-]", "_", held_out)
  wt_reference <- NULL
  for (target in primary_targets) {
    message("LODO patient ", held_out, " target ", target)
    fit_file <- file.path(lodo_dir,
      paste0("scTenifoldKnk_LODO__", safe_patient, "__", target, ".rds"))
    if (file.exists(fit_file)) {
      fit <- readRDS(fit_file)
      elapsed <- NA_real_
    } else {
      set.seed(seed)
      started <- Sys.time()
      fit <- scTenifoldKnk::scTenifoldKnk(
        countMatrix = x_lodo, qc = parameters$qc, gKO = target,
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
    wt <- as.matrix(fit$tensorNetworks$WT)
    if (is.null(wt_reference)) {
      wt_reference <- wt
      wt_difference <- 0
    } else {
      wt_difference <- max(abs(wt - wt_reference))
      if (!is.finite(wt_difference) || wt_difference > 1e-12) {
        stop("Same-LODO WT network changed across targets: ", held_out)
      }
    }
    key <- paste(held_out, target, sep = "__")
    program_rows[[key]] <- enrich_one(as.data.table(fit$diffRegulation), target, held_out)
    run_rows[[key]] <- data.table(
      held_out_patient = held_out, target = target,
      cells_after_holdout = ncol(x_lodo), held_out_cells = sum(!keep),
      elapsed_seconds = elapsed, wt_max_abs_difference_within_LODO = wt_difference
    )
    if (!file.exists(fit_file)) saveRDS(fit, fit_file, compress = "xz")
    rm(fit, wt)
    gc(verbose = FALSE)
  }
}

program_table <- rbindlist(program_rows)
run_table <- rbindlist(run_rows)
if (fit_only_worker) {
  fwrite(run_table, file.path(lodo_dir,
    paste0("worker_", worker_label, "_run_summary.tsv")), sep = "\t")
  cat("LODO worker", worker_label, "completed", nrow(run_table), "fits\n")
  quit(save = "no", status = 0L)
}

if (nrow(program_table) != 26L * length(primary_targets) * length(programs)) {
  stop("Combined LODO program table is incomplete")
}
program_summary <- program_table[, .(
  leave_one_donor_runs = .N,
  median_auc = median(enrichment_auc),
  minimum_auc = min(enrichment_auc),
  maximum_auc = max(enrichment_auc),
  fraction_auc_above_0_5 = mean(enrichment_auc > 0.5),
  runs_FDR_below_0_05 = sum(wilcoxon_FDR_within_target_LODO < 0.05),
  fraction_FDR_below_0_05 = mean(wilcoxon_FDR_within_target_LODO < 0.05)
), by = .(target, program)]
qualification <- program_summary[program %chin% core_programs,
  .(core_programs_above_0_5_in_at_least_80pct = sum(fraction_auc_above_0_5 >= 0.8),
    minimum_core_direction_fraction = min(fraction_auc_above_0_5),
    core_programs_FDR_below_0_05_in_at_least_80pct =
      sum(fraction_FDR_below_0_05 >= 0.8)), by = target]

fwrite(program_table, file.path(lodo_dir,
  "GSE183904_primary_receptor_LODO_program_enrichment.tsv"), sep = "\t")
fwrite(program_summary, file.path(lodo_dir,
  "GSE183904_primary_receptor_LODO_program_summary.tsv"), sep = "\t")
fwrite(qualification, file.path(lodo_dir,
  "GSE183904_primary_receptor_LODO_qualification.tsv"), sep = "\t")
fwrite(run_table, file.path(lodo_dir,
  "GSE183904_primary_receptor_LODO_run_summary.tsv"), sep = "\t")

manifest <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE183904", state = "C1QC_APOE_macrophage",
  targets = primary_targets, patients = patients, seed = seed,
  interpretation_unit = "leave-one-donor-out sensitivity; seed is algorithmic",
  parameters = parameters, core_programs = core_programs,
  direction_stability_gate = "AUC > 0.5 in at least 80% of 26 holdouts",
  same_holdout_common_WT = all(run_table$wt_max_abs_difference_within_LODO <= 1e-12),
  differential_regulation_semantics = "unsigned manifold distance",
  input_sha256 = list(
    matrix = digest(matrix_file, algo = "sha256", file = TRUE),
    metadata = digest(metadata_file, algo = "sha256", file = TRUE),
    programs = digest(program_file, algo = "sha256", file = TRUE)
  ),
  claim_boundary = paste(
    "LODO evaluates sensitivity to donor omission within one discovery cohort.",
    "It does not supply an independent cohort perturbation or causal evidence."
  )
)
write_json(manifest, file.path(audit_dir,
  "GSE183904_primary_receptors_v2_scTenifoldKnk_leave_one_donor_out.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nLODO core-program qualification:\n")
print(qualification[order(-core_programs_above_0_5_in_at_least_80pct,
                          -minimum_core_direction_fraction)])

