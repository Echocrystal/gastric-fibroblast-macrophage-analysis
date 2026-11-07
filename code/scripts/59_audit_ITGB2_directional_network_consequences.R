#!/usr/bin/env Rscript

# Audit whether the frozen scTenifoldKnk objects support a directional
# interpretation through WT-versus-KO network connectivity.  This is a
# descriptive network-connectivity analysis; it is not an expression change.

suppressPackageStartupMessages(library(data.table))

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
base <- file.path(project, "analysis/results/virtual_knockout/GSE183904",
                  "C1QC_APOE_macrophage_receptor_screen_v2")
program_file <- file.path(project, "resources/gene_sets",
  "GSE183904_C1QC_APOE_receptor_screen_v2_scTenifoldKnk_programs_axis_excluded.tsv")
out_dir <- file.path(base, "ITGB2_directional_network_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

programs <- fread(program_file)
seeds <- c(2026072739L, 2026072740L, 2026072741L)

node_strength <- function(M) rowSums(abs(M)) + colSums(abs(M))
within_strength <- function(M, genes, universe) {
  idx <- which(universe %in% genes)
  if (length(idx) < 2L) return(NA_real_)
  sum(abs(M[idx, idx, drop = FALSE])) / (length(idx) * (length(idx) - 1L))
}

rows <- list()
for (seed in seeds) {
  fit <- readRDS(file.path(base, paste0("scTenifoldKnk_seed", seed, "__ITGB2.rds")))
  WT <- as.matrix(fit$tensorNetworks$WT)
  KO <- as.matrix(fit$tensorNetworks$KO)
  genes <- rownames(WT)
  delta <- node_strength(KO) - node_strength(WT)
  for (program_id in unique(programs$program)) {
    program_genes <- intersect(programs[program == program_id]$gene_symbol, genes)
    wt_strength <- node_strength(WT)[match(program_genes, genes)]
    ko_strength <- node_strength(KO)[match(program_genes, genes)]
    rows[[paste(seed, program_id)]] <- data.table(
      seed = seed,
      program = program_id,
      n_genes = length(program_genes),
      mean_WT_total_abs_strength = mean(wt_strength),
      mean_KO_total_abs_strength = mean(ko_strength),
      mean_delta_total_abs_strength = mean(delta[match(program_genes, genes)]),
      median_delta_total_abs_strength = median(delta[match(program_genes, genes)]),
      mean_percent_change_total_abs_strength =
        100 * (mean(ko_strength) - mean(wt_strength)) / mean(wt_strength),
      delta_mean_within_program_abs_weight =
        within_strength(KO, program_genes, genes) - within_strength(WT, program_genes, genes)
    )
  }
}

by_seed <- rbindlist(rows)
summary <- by_seed[, .(
  seeds = .N,
  median_of_seed_medians = median(median_delta_total_abs_strength),
  minimum_seed_median = min(median_delta_total_abs_strength),
  maximum_seed_median = max(median_delta_total_abs_strength),
  seeds_with_strength_loss = sum(median_delta_total_abs_strength < 0),
  median_percent_change_total_abs_strength = median(mean_percent_change_total_abs_strength),
  minimum_percent_change_total_abs_strength = min(mean_percent_change_total_abs_strength),
  maximum_percent_change_total_abs_strength = max(mean_percent_change_total_abs_strength),
  median_delta_within_program_abs_weight = median(delta_mean_within_program_abs_weight),
  seeds_with_within_program_loss = sum(delta_mean_within_program_abs_weight < 0)
), by = program]

fwrite(by_seed, file.path(out_dir, "ITGB2_WT_KO_program_connectivity_by_seed.tsv"), sep = "\t")
fwrite(summary, file.path(out_dir, "ITGB2_WT_KO_program_connectivity_summary.tsv"), sep = "\t")
print(by_seed)
print(summary)
