#!/usr/bin/env Rscript

# Representative-seed matched-control KOs for the sender-side SPON2 screen.
# Controls were frozen without reading SPON2 perturbation outcomes.

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
start_end <- if (length(args)) as.integer(strsplit(args[[1L]], ":", fixed = TRUE)[[1L]]) else c(1L, 20L)
if (length(start_end) != 2L || anyNA(start_end)) stop("Use START:END")
project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(project, "runtime", "R_virtual_knockout"),
            "I:/ECM-new/scRNA/20.rebuild_GSE163558/R_libs_4.4",
            file.path(source_project, "runtime", "R_library"), .libPaths()))
suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(scTenifoldKnk); library(scTenifoldNet)
  library(jsonlite)
})

input_dir <- file.path(project, "analysis", "processed", "virtual_knockout", "GSE206785")
result_dir <- file.path(project, "analysis", "results", "virtual_knockout", "GSE206785", "SPON2_sender_screen_v1")
control_dir <- file.path(result_dir, "matched_controls_seed2026081011")
resource_dir <- file.path(project, "resources", "gene_sets")
dir.create(control_dir, recursive = TRUE, showWarnings = FALSE)

genes <- fread(file.path(result_dir, "SPON2_sender_all_diffreg.tsv"))[seed == 2026081011L, gene]
counts_file <- file.path(input_dir, "GSE206785_sender_CTHRC1_CAF_full_raw_counts.mtx.gz")
all_genes <- fread(file.path(input_dir, "GSE206785_sender_CTHRC1_CAF_full_raw_counts_genes.tsv"))$gene_symbol
cells <- fread(file.path(input_dir, "GSE206785_sender_CTHRC1_CAF_full_raw_counts_metadata.tsv"))
x_full <- readMM(gzfile(counts_file)); rownames(x_full) <- all_genes
x <- x_full[genes, , drop = FALSE]
program_table <- rbind(
  fread(file.path(resource_dir, "virtual_knockout_sanity_programs_axis_excluded.tsv"))[
    program %chin% c("GOBP_CELL_MATRIX_ADHESION", "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY",
                     "REACTOME_INTEGRIN_SIGNALING"), .(program, gene_symbol)],
  fread(file.path(resource_dir, "GSE183904_fine_state_modules_frozen.tsv"))[
    lineage == "Fibroblast" & module %chin% c("Activated_matrix", "Inflammatory_CAF",
                                               "Adventitial_quiescent", "Contractile"),
    .(program = module, gene_symbol)]
)
programs <- split(program_table$gene_symbol[program_table$gene_symbol %chin% genes],
                  program_table$program[program_table$gene_symbol %chin% genes])
programs <- lapply(programs, unique)
controls <- fread(file.path(resource_dir, "GSE206785_SPON2_sender_matched_controls.tsv"))$control_gene
controls_to_run <- controls[seq.int(start_end[[1L]], start_end[[2L]])]
seed <- 2026081011L
parameters <- list(qc=FALSE, nc_lambda=0, nc_nNet=3L, nc_nCells=200L,
                   nc_nComp=3L, nc_scaleScores=TRUE, nc_symmetric=FALSE,
                   nc_q=0.90, td_K=3L, td_maxIter=1000L, td_maxError=1e-5,
                   td_nDecimal=3L, ma_nDim=2L, nCores=4L)

enrich_one <- function(fit, ko_gene) {
  wt <- as.matrix(fit$tensorNetworks$WT); ko <- as.matrix(fit$tensorNetworks$KO)
  sw <- rowSums(abs(wt)) + colSums(abs(wt)); sk <- rowSums(abs(ko)) + colSums(abs(ko))
  dr <- as.data.table(fit$diffRegulation)[gene != ko_gene]
  dr[, distance_percentile := frank(distance, ties.method = "average") / .N]
  rows <- lapply(names(programs), function(pid) {
    g <- setdiff(unique(programs[[pid]]), ko_gene); in_set <- dr$gene %in% g
    n_in <- sum(in_set); n_out <- sum(!in_set)
    wt_test <- suppressWarnings(wilcox.test(dr$distance[in_set], dr$distance[!in_set],
                                            alternative = "greater", exact = FALSE))
    ix <- match(g, rownames(wt)); wt_mean <- mean(sw[ix], na.rm = TRUE); ko_mean <- mean(sk[ix], na.rm = TRUE)
    data.table(ko_gene=ko_gene, program=pid, n_program_genes=n_in,
               enrichment_auc=as.numeric(wt_test$statistic)/(n_in*n_out),
               wilcoxon_p=wt_test$p.value,
               WT_mean_total_abs_strength=wt_mean, KO_mean_total_abs_strength=ko_mean,
               percent_change_total_abs_strength=100*(ko_mean-wt_mean)/wt_mean)
  })
  rbindlist(rows)
}

rows <- list()
for (control in controls_to_run) {
  message("Running sender matched control ", control)
  f <- file.path(control_dir, paste0("scTenifoldKnk_control__", control, ".rds"))
  if (file.exists(f)) fit <- readRDS(f) else {
    set.seed(seed)
    fit <- scTenifoldKnk::scTenifoldKnk(
      countMatrix=x, qc=parameters$qc, gKO=control,
      nc_lambda=parameters$nc_lambda, nc_nNet=parameters$nc_nNet,
      nc_nCells=parameters$nc_nCells, nc_nComp=parameters$nc_nComp,
      nc_scaleScores=parameters$nc_scaleScores, nc_symmetric=parameters$nc_symmetric,
      nc_q=parameters$nc_q, td_K=parameters$td_K, td_maxIter=parameters$td_maxIter,
      td_maxError=parameters$td_maxError, td_nDecimal=parameters$td_nDecimal,
      ma_nDim=parameters$ma_nDim, nCores=parameters$nCores)
    saveRDS(fit, f, compress="xz")
  }
  rows[[control]] <- enrich_one(fit, control)
  rm(fit); gc(verbose=FALSE)
}
controls_long <- rbindlist(rows)
fwrite(controls_long, file.path(result_dir, paste0("SPON2_sender_controls_", start_end[[1L]], "_", start_end[[2L]], ".tsv")), sep="\t")
if (start_end[[1L]] == 1L && start_end[[2L]] == 20L) {
  target_fit <- readRDS(file.path(result_dir, "scTenifoldKnk_seed2026081011__SPON2.rds"))
  target_long <- enrich_one(target_fit, "SPON2")
  target_long[, role := "target"]; controls_long[, role := "matched_control"]
  all_long <- rbind(target_long, controls_long)
  summary <- all_long[, {
    tv <- enrichment_auc[role == "target"]; cv <- enrichment_auc[role == "matched_control"]
    td <- percent_change_total_abs_strength[role == "target"]; cd <- percent_change_total_abs_strength[role == "matched_control"]
    data.table(target_auc=tv, control_median_auc=median(cv), controls_at_least_as_high_auc=sum(cv >= tv),
               auc_lower_tail_p=(1+sum(cv >= tv))/(1+length(cv)),
               target_percent_change=td, control_median_percent_change=median(cd),
               controls_at_least_as_negative=sum(cd <= td),
               directional_lower_tail_p=(1+sum(cd <= td))/(1+length(cd)))
  }, by=program]
  summary[, auc_FDR_BH := p.adjust(auc_lower_tail_p, "BH")]
  summary[, directional_FDR_BH := p.adjust(directional_lower_tail_p, "BH")]
  fwrite(all_long, file.path(result_dir, "SPON2_sender_matched_controls_all_nodes.tsv"), sep="\t")
  fwrite(summary, file.path(result_dir, "SPON2_sender_matched_controls_summary.tsv"), sep="\t")
  write_json(list(target="SPON2", seed=seed, controls=controls, controls_frozen_without_KO_outcomes=TRUE,
                  parameters=parameters, interpretation="sender-side specificity calibration"),
             file.path(result_dir, "SPON2_sender_matched_controls_manifest.json"), pretty=TRUE, auto_unbox=TRUE)
  print(summary)
} else {
  cat("Completed controls", start_end[[1L]], "to", start_end[[2L]], "\n")
}
