#!/usr/bin/env Rscript

# Run the frozen ITGAM matched-control set at the representative seed and
# compute directional WT-versus-KO programme connectivity changes.

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
start_end <- if (length(args)) as.integer(strsplit(args[[1L]], ":", fixed = TRUE)[[1L]]) else c(1L, 20L)
if (length(start_end) != 2L || anyNA(start_end)) stop("Use START:END")

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(.libPaths(), file.path(project, "runtime/R_virtual_knockout"),
            "I:/ECM-new/scRNA/20.rebuild_GSE163558/R_libs_4.4",
            file.path(source_project, "runtime/R_library")))
suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(scTenifoldKnk); library(scTenifoldNet)
})

analysis_dir <- file.path(project, "analysis/results/virtual_knockout/GSE183904",
                          "C1QC_APOE_macrophage_receptor_screen_v2")
input_dir <- file.path(project, "analysis/processed/virtual_knockout/GSE183904",
                       "C1QC_APOE_macrophage_receptor_screen_v2")
control_dir <- file.path(analysis_dir, "ITGAM_directional_matched_controls_seed2026072739")
resource_dir <- file.path(project, "resources/gene_sets")
out_dir <- file.path(analysis_dir, "ITGAM_vs_ITGB2_audit")
dir.create(control_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- readRDS(file.path(input_dir, "GSE183904_C1QC_APOE_sanity_counts.rds"))
program_file <- file.path(resource_dir,
  "GSE183904_C1QC_APOE_receptor_screen_v2_scTenifoldKnk_programs_axis_excluded.tsv")
program_manifest <- fread(program_file)
programs <- split(program_manifest$gene_symbol, program_manifest$program)
control_manifest <- fread(file.path(resource_dir,
  "GSE183904_ITGAM_directional_network_matched_controls.tsv"))
controls <- control_manifest$control_gene
controls_to_run <- controls[seq.int(start_end[[1L]], start_end[[2L]])]
seed <- 2026072739L
axis_genes <- c("CD74", "CCRL2", "CXCR4", "ITGB2", "LRP1", "ITGAM", "CMKLR1", "ANXA2", "TREM2", "PLAT", "ILK")
parameters <- list(qc=FALSE, nc_lambda=0, nc_nNet=3L, nc_nCells=350L,
                   nc_nComp=3L, nc_scaleScores=TRUE, nc_symmetric=FALSE,
                   nc_q=0.90, td_K=3L, td_maxIter=1000L, td_maxError=1e-5,
                   td_nDecimal=3L, ma_nDim=2L, nCores=4L)

strength <- function(M) rowSums(abs(M)) + colSums(abs(M))
enrich_directional <- function(fit, ko_gene) {
  WT <- as.matrix(fit$tensorNetworks$WT); KO <- as.matrix(fit$tensorNetworks$KO)
  genes <- rownames(WT); sw <- strength(WT); sk <- strength(KO)
  rbindlist(lapply(names(programs), function(program_id) {
    pg <- intersect(setdiff(programs[[program_id]], axis_genes), genes)
    ix <- match(pg, genes); wt_mean <- mean(sw[ix]); ko_mean <- mean(sk[ix])
    data.table(ko_gene=ko_gene, role="matched_control", program=program_id,
               n_program_genes=length(pg), WT_mean_total_abs_strength=wt_mean,
               KO_mean_total_abs_strength=ko_mean,
               delta_total_abs_strength=ko_mean-wt_mean,
               percent_change_total_abs_strength=100*(ko_mean-wt_mean)/wt_mean)
  }))
}

rows <- list()
for (control in controls_to_run) {
  message("Running ITGAM matched control ", control)
  fit_file <- file.path(control_dir, paste0("scTenifoldKnk_control__", control, ".rds"))
  if (file.exists(fit_file)) {
    fit <- readRDS(fit_file)
  } else {
    set.seed(seed)
    fit <- scTenifoldKnk::scTenifoldKnk(
      countMatrix=x, qc=parameters$qc, gKO=control,
      nc_lambda=parameters$nc_lambda, nc_nNet=parameters$nc_nNet,
      nc_nCells=parameters$nc_nCells, nc_nComp=parameters$nc_nComp,
      nc_scaleScores=parameters$nc_scaleScores, nc_symmetric=parameters$nc_symmetric,
      nc_q=parameters$nc_q, td_K=parameters$td_K, td_maxIter=parameters$td_maxIter,
      td_maxError=parameters$td_maxError, td_nDecimal=parameters$td_nDecimal,
      ma_nDim=parameters$ma_nDim, nCores=parameters$nCores)
    saveRDS(fit, fit_file, compress="xz")
  }
  rows[[control]] <- enrich_directional(fit, control)
  rm(fit); gc(verbose=FALSE)
}
values <- rbindlist(rows)
if (start_end[[1L]] == 1L && start_end[[2L]] == 20L) {
  target_fit <- readRDS(file.path(analysis_dir, "scTenifoldKnk_seed2026072739__ITGAM.rds"))
  target_rows <- enrich_directional(target_fit, "ITGAM")[, role := "target"]
  all_values <- rbind(target_rows, values)
  summary <- all_values[, {
    tval <- percent_change_total_abs_strength[role == "target"]
    null <- percent_change_total_abs_strength[role == "matched_control"]
    data.table(target_percent_change=tval, control_median_percent_change=median(null),
               controls_at_least_as_negative=sum(null <= tval),
               empirical_lower_tail_p=(1+sum(null <= tval))/(1+length(null)))
  }, by=program]
  summary[, empirical_FDR_BH := p.adjust(empirical_lower_tail_p, "BH")]
  fwrite(all_values, file.path(out_dir, "ITGAM_directional_connectivity_all_nodes.tsv"), sep="\t")
  fwrite(summary, file.path(out_dir, "ITGAM_directional_connectivity_summary.tsv"), sep="\t")
  rm(target_fit)
} else {
  fwrite(values, file.path(out_dir, paste0("ITGAM_directional_controls_", start_end[[1L]], "_", start_end[[2L]], ".tsv")), sep="\t")
}
cat("Completed controls", start_end[[1L]], "to", start_end[[2L]], "\n")
