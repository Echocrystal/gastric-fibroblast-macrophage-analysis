#!/usr/bin/env Rscript

# Parallel audit of ITGAM and ITGB2 using the already frozen three-seed
# scTenifoldKnk fits.  This stage deliberately does not use target outcomes to
# choose controls; it reports programme-level perturbation and WT-KO network
# connectivity changes before any new ITGAM control run is frozen.

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
base <- file.path(project, "analysis/results/virtual_knockout/GSE183904",
                  "C1QC_APOE_macrophage_receptor_screen_v2")
program_file <- file.path(project, "resources/gene_sets",
  "GSE183904_CFD_APOE_receptor_screen_v2_scTenifoldKnk_programs_axis_excluded.tsv")
if (!file.exists(program_file)) {
  program_file <- file.path(project, "resources/gene_sets",
    "GSE183904_C1QC_APOE_receptor_screen_v2_scTenifoldKnk_programs_axis_excluded.tsv")
}
out_dir <- file.path(base, "ITGAM_vs_ITGB2_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

programs <- fread(program_file)
program_ids <- unique(programs$program)
seeds <- c(2026072739L, 2026072740L, 2026072741L)
targets <- c("ITGB2", "ITGAM")
axis_genes <- c("CD74", "CCRL2", "CXCR4", "ITGB2", "LRP1", "ITGAM",
                "CMKLR1", "ANXA2", "TREM2", "PLAT", "ILK")

strength <- function(M) rowSums(abs(M)) + colSums(abs(M))

rows <- list()
for (target in targets) {
  for (seed in seeds) {
    f <- file.path(base, paste0("scTenifoldKnk_seed", seed, "__", target, ".rds"))
    fit <- readRDS(f)
    WT <- as.matrix(fit$tensorNetworks$WT)
    KO <- as.matrix(fit$tensorNetworks$KO)
    genes <- rownames(WT)
    wt_strength <- strength(WT)
    ko_strength <- strength(KO)
    dr <- as.data.table(fit$diffRegulation)
    target_row <- dr[gene == target]
    for (program_id in program_ids) {
      pg <- intersect(setdiff(programs[program == program_id]$gene_symbol, axis_genes), genes)
      ix <- match(pg, genes)
      wt_mean <- mean(wt_strength[ix]); ko_mean <- mean(ko_strength[ix])
      in_set <- dr$gene %in% pg
      out_set <- !in_set & !dr$gene %in% c(axis_genes, target)
      p <- suppressWarnings(wilcox.test(dr$distance[in_set], dr$distance[out_set], alternative = "greater", exact = FALSE)$p.value)
      rows[[paste(target, seed, program_id)]] <- data.table(
        target = target, seed = seed, program = program_id,
        n_program_genes = length(pg),
        enrichment_auc = as.numeric(wilcox.test(dr$distance[in_set], dr$distance[out_set], alternative = "greater", exact = FALSE)$statistic) / (sum(in_set) * sum(out_set)),
        program_wilcoxon_p = p,
        program_wilcoxon_q_within_seed = NA_real_,
        target_diffreg_Z = target_row$Z,
        target_diffreg_q = target_row$p.adj,
        WT_mean_total_abs_strength = wt_mean,
        KO_mean_total_abs_strength = ko_mean,
        percent_change_total_abs_strength = 100 * (ko_mean - wt_mean) / wt_mean
      )
    }
    rm(fit, WT, KO, dr)
  }
}
by_seed <- rbindlist(rows)
by_seed[, program_wilcoxon_q_within_seed := p.adjust(program_wilcoxon_p, "BH"), by = .(target, seed)]
summary <- by_seed[, .(
  seeds = .N,
  median_AUC = median(enrichment_auc),
  seeds_q05 = sum(program_wilcoxon_q_within_seed < 0.05),
  median_percent_change = median(percent_change_total_abs_strength),
  min_percent_change = min(percent_change_total_abs_strength),
  max_percent_change = max(percent_change_total_abs_strength)
), by = .(target, program)]

fwrite(by_seed, file.path(out_dir, "ITGAM_ITGB2_directional_by_seed.tsv"), sep = "\t")
fwrite(summary, file.path(out_dir, "ITGAM_ITGB2_directional_summary.tsv"), sep = "\t")
print(summary[order(program, target)])
