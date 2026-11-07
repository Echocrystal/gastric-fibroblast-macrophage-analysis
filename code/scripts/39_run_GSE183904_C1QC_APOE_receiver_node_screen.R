#!/usr/bin/env Rscript

# Three-seed exploratory scTenifoldKnk screen of receiver-associated nodes in
# the donor-capped recurrent C1QC+APOE+ macrophage state. Differential-
# regulation distances are unsigned throughout.

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
screen_version <- if (length(args)) args[[1L]] else "v1"
if (!screen_version %in% c("v1", "receptor_screen_v2")) {
  stop("Unknown screen version: ", screen_version)
}

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

analysis_label <- if (screen_version == "v1") {
  "C1QC_APOE_macrophage_v1"
} else {
  "C1QC_APOE_macrophage_receptor_screen_v2"
}
resource_prefix <- if (screen_version == "v1") {
  "GSE183904_C1QC_APOE_scTenifoldKnk"
} else {
  "GSE183904_C1QC_APOE_receptor_screen_v2_scTenifoldKnk"
}
input_dir <- file.path(project, "analysis", "processed", "virtual_knockout", "GSE183904",
                       analysis_label)
result_dir <- file.path(project, "analysis", "results", "virtual_knockout", "GSE183904",
                        analysis_label)
audit_dir <- file.path(project, "analysis", "audit")
program_file <- file.path(project, "resources", "gene_sets",
                          paste0(resource_prefix, "_programs_axis_excluded.tsv"))
matrix_file <- file.path(input_dir, "GSE183904_C1QC_APOE_sanity_counts.rds")
metadata_file <- file.path(input_dir, "GSE183904_C1QC_APOE_sanity_metadata.tsv")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(matrix_file)
metadata <- fread(metadata_file)
program_manifest <- fread(program_file)
expected_genes <- if (screen_version == "v1") 774L else 775L
if (nrow(x) != expected_genes || ncol(x) != 706L || uniqueN(metadata$Patient_ID) != 26L) {
  stop("Frozen receiver virtual-KO input changed")
}

frozen_targets <- if (screen_version == "v1") {
  c("ANXA2", "ITGAM", "ITGB2", "LRP1")
} else {
  c("CD74", "CCRL2", "CXCR4", "ITGB2", "LRP1", "ITGAM", "CMKLR1", "ANXA2", "TREM2")
}
targets <- frozen_targets
if (length(args) >= 2L && nzchar(args[[2L]])) {
  requested_targets <- unique(strsplit(args[[2L]], ",", fixed = TRUE)[[1L]])
  if (!all(requested_targets %in% frozen_targets)) {
    stop("Requested target is outside the frozen set: ",
         paste(setdiff(requested_targets, frozen_targets), collapse = ", "))
  }
  targets <- requested_targets
}
axis_genes <- c(frozen_targets, "PLAT", "ILK")
if (!all(frozen_targets %in% rownames(x))) stop("A frozen receiver screen target is missing")
if (any(program_manifest$gene_symbol %in% axis_genes)) stop("Axis-gene leakage into program scores")
programs <- split(program_manifest$gene_symbol, program_manifest$program)

seeds <- c(2026072739L, 2026072740L, 2026072741L)
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

enrich_one <- function(dr, target, seed) {
  excluded <- unique(c(axis_genes, target))
  dr <- dr[!gene %in% excluded]
  dr[, distance_percentile := frank(distance, ties.method = "average") / .N]
  rows <- lapply(names(programs), function(program_id) {
    genes <- setdiff(unique(programs[[program_id]]), excluded)
    in_set <- dr$gene %in% genes
    n_in <- sum(in_set)
    n_out <- sum(!in_set)
    if (n_in < 5L || n_out < 5L) {
      return(data.table(program = program_id, n_in_universe = n_in,
                        mean_distance_percentile = NA_real_, enrichment_auc = NA_real_,
                        wilcoxon_p = NA_real_))
    }
    wt <- suppressWarnings(wilcox.test(dr$distance[in_set], dr$distance[!in_set],
                                       alternative = "greater", exact = FALSE))
    data.table(
      program = program_id,
      n_in_universe = n_in,
      mean_distance_percentile = mean(dr$distance_percentile[in_set]),
      enrichment_auc = as.numeric(wt$statistic) / (n_in * n_out),
      wilcoxon_p = wt$p.value
    )
  })
  out <- rbindlist(rows)
  out[, `:=`(target = target, seed = seed)]
  out[, wilcoxon_FDR_within_target_seed := p.adjust(wilcoxon_p, "BH")]
  setcolorder(out, c("target", "seed", "program", "n_in_universe",
                    "mean_distance_percentile", "enrichment_auc", "wilcoxon_p",
                    "wilcoxon_FDR_within_target_seed"))
  out
}

dr_rows <- list()
program_rows <- list()
summary_rows <- list()
wt_checks <- list()

for (seed in seeds) {
  wt_reference <- NULL
  for (target in targets) {
    message("Running seed ", seed, " target ", target)
    fit_file <- file.path(result_dir,
                          paste0("scTenifoldKnk_seed", seed, "__", target, ".rds"))
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
        gKO = target,
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
    if (is.null(wt_reference)) {
      wt_reference <- wt
      wt_difference <- 0
    } else {
      wt_difference <- max(abs(wt - wt_reference))
      if (!is.finite(wt_difference) || wt_difference > 1e-12) {
        stop("Same-seed WT network changed across targets")
      }
    }

    dr <- as.data.table(fit$diffRegulation)
    dr[, `:=`(
      target = target,
      seed = seed,
      rank_by_distance = frank(-distance, ties.method = "min"),
      out_strength_WT = rowSums(abs(wt))[gene],
      in_strength_WT = colSums(abs(wt))[gene],
      mean_raw_count = Matrix::rowMeans(x)[gene],
      detection_fraction = Matrix::rowMeans(x > 0)[gene]
    )]
    target_row <- dr[gene == target]
    if (nrow(target_row) != 1L) stop("Target row missing from differential regulation")

    key <- paste(seed, target, sep = "__")
    dr_rows[[key]] <- dr
    program_rows[[key]] <- enrich_one(copy(dr), target, seed)
    summary_rows[[key]] <- data.table(
      target = target,
      seed = seed,
      elapsed_seconds = elapsed,
      wt_max_abs_difference_within_seed = wt_difference,
      target_rank = target_row$rank_by_distance,
      target_distance = target_row$distance,
      target_p_adj = target_row$p.adj,
      target_out_strength_WT = target_row$out_strength_WT,
      target_in_strength_WT = target_row$in_strength_WT,
      target_detection_fraction = target_row$detection_fraction,
      target_mean_raw_count = target_row$mean_raw_count,
      n_diffreg_FDR_0_05 = sum(dr$p.adj < 0.05, na.rm = TRUE)
    )
    wt_checks[[key]] <- wt_difference
    if (!file.exists(fit_file)) saveRDS(fit, fit_file, compress = "xz")
    rm(fit, wt)
    gc(verbose = FALSE)
  }
}

dr_table <- rbindlist(dr_rows)
program_table <- rbindlist(program_rows)
run_summary <- rbindlist(summary_rows)
program_summary <- program_table[, .(
  seeds = .N,
  median_enrichment_auc = median(enrichment_auc),
  minimum_enrichment_auc = min(enrichment_auc),
  maximum_enrichment_auc = max(enrichment_auc),
  median_mean_distance_percentile = median(mean_distance_percentile),
  maximum_wilcoxon_p = max(wilcoxon_p),
  maximum_wilcoxon_FDR = max(wilcoxon_FDR_within_target_seed),
  seeds_FDR_below_0_05 = sum(wilcoxon_FDR_within_target_seed < 0.05)
), by = .(target, program)]
program_summary[, target_rank_by_median_auc :=
                  frank(-median_enrichment_auc, ties.method = "min"), by = program]

fwrite(dr_table, file.path(result_dir, "GSE183904_C1QC_APOE_scTenifoldKnk_all_diffreg.tsv.gz"),
       sep = "\t")
fwrite(program_table,
       file.path(result_dir, "GSE183904_C1QC_APOE_scTenifoldKnk_program_enrichment_by_seed.tsv"),
       sep = "\t")
fwrite(program_summary,
       file.path(result_dir, "GSE183904_C1QC_APOE_scTenifoldKnk_program_enrichment_summary.tsv"),
       sep = "\t")
fwrite(run_summary,
       file.path(result_dir, "GSE183904_C1QC_APOE_scTenifoldKnk_run_summary.tsv"), sep = "\t")

manifest <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE183904",
  state = "C1QC_APOE_macrophage",
  screen_version = screen_version,
  analysis_label = analysis_label,
  analysis = "three-seed exploratory receiver-node scTenifoldKnk screen",
  frozen_targets = frozen_targets,
  executed_targets = targets,
  seeds = seeds,
  parameters = parameters,
  same_seed_common_WT = all(unlist(wt_checks) <= 1e-12),
  differential_regulation_semantics = "unsigned manifold distance; expression direction is not inferred",
  program_axis_exclusion = axis_genes,
  multiplicity = "BH across six programs within each target and seed",
  interpretation = paste(
    "Target ranking is an exploratory screen across four receiver-associated nodes.",
    "Matched expression/degree/centrality controls are required before retaining a model-predicted dependency.",
    "The calculation does not establish extracellular binding, signaling direction or biological causality."
  ),
  input_sha256 = list(
    matrix = digest(matrix_file, algo = "sha256", file = TRUE),
    metadata = digest(metadata_file, algo = "sha256", file = TRUE),
    programs = digest(program_file, algo = "sha256", file = TRUE)
  ),
  software = list(
    R = R.version.string,
    scTenifoldKnk = as.character(packageVersion("scTenifoldKnk")),
    scTenifoldNet = as.character(packageVersion("scTenifoldNet"))
  )
)
audit_file <- if (screen_version == "v1") {
  "GSE183904_C1QC_APOE_scTenifoldKnk_screen.json"
} else {
  "GSE183904_C1QC_APOE_receptor_screen_v2_scTenifoldKnk_screen.json"
}
write_json(manifest,
           file.path(audit_dir, audit_file),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nProgram enrichment summary:\n")
print(program_summary[order(program, target_rank_by_median_auc)])
