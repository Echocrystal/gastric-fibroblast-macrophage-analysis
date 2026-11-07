#!/usr/bin/env Rscript

# Directional WT-versus-KO network consequence analysis for ITGB2 and the
# prespecified matched-node control set.  The outcome is programme-gene total
# absolute network-strength change, not RNA expression change.

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
analysis_dir <- file.path(project, "analysis/results/virtual_knockout/GSE183904",
                           "C1QC_APOE_macrophage_receptor_screen_v2")
control_dir <- file.path(analysis_dir, "primary_receptor_network_matched_controls_seed2026072739")
program_file <- file.path(project, "resources/gene_sets",
  "GSE183904_C1QC_APOE_receptor_screen_v2_scTenifoldKnk_programs_axis_excluded.tsv")
manifest_file <- file.path(project, "resources/gene_sets",
  "GSE183904_primary_receptors_v2_scTenifoldKnk_network_matched_controls.tsv")
out_dir <- file.path(analysis_dir, "ITGB2_directional_network_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

programs <- fread(program_file)
controls <- fread(manifest_file)[target == "ITGB2", control_gene]
if (length(controls) != 20L) stop("Expected 20 frozen ITGB2 matched controls")
program_ids <- unique(programs$program)
axis_genes <- c("CD74", "CCRL2", "CXCR4", "ITGB2", "LRP1", "ITGAM",
                "CMKLR1", "ANXA2", "TREM2", "PLAT", "ILK")

strength <- function(M) rowSums(abs(M)) + colSums(abs(M))

read_fit <- function(gene, is_target = FALSE) {
  if (is_target) {
    f <- file.path(analysis_dir, "scTenifoldKnk_seed2026072739__ITGB2.rds")
  } else {
    f <- file.path(control_dir, paste0("scTenifoldKnk_control__", gene, ".rds"))
  }
  if (!file.exists(f)) stop("Missing fit: ", f)
  readRDS(f)
}

one_gene <- function(gene, role, fit) {
  WT <- as.matrix(fit$tensorNetworks$WT)
  KO <- as.matrix(fit$tensorNetworks$KO)
  genes <- rownames(WT)
  wt_strength <- strength(WT)
  ko_strength <- strength(KO)
  rows <- lapply(program_ids, function(program_id) {
    pg <- intersect(setdiff(programs[program == program_id]$gene_symbol, axis_genes), genes)
    if (length(pg) < 5L) stop("Programme too small: ", program_id)
    ix <- match(pg, genes)
    wt_mean <- mean(wt_strength[ix])
    ko_mean <- mean(ko_strength[ix])
    # Direction convention: negative = loss of total absolute network
    # connectivity after removing the node.
    data.table(
      ko_gene = gene,
      role = role,
      program = program_id,
      n_program_genes = length(pg),
      WT_mean_total_abs_strength = wt_mean,
      KO_mean_total_abs_strength = ko_mean,
      delta_total_abs_strength = ko_mean - wt_mean,
      percent_change_total_abs_strength = 100 * (ko_mean - wt_mean) / wt_mean
    )
  })
  rbindlist(rows)
}

target_fit <- read_fit("ITGB2", is_target = TRUE)
rows <- list(one_gene("ITGB2", "target", target_fit))
rm(target_fit)
for (gene in controls) {
  message("Reading matched control: ", gene)
  rows[[length(rows) + 1L]] <- one_gene(gene, "matched_control", read_fit(gene))
}
all_values <- rbindlist(rows)

summary <- all_values[, {
  target_value <- percent_change_total_abs_strength[role == "target"]
  null <- percent_change_total_abs_strength[role == "matched_control"]
  data.table(
    target_percent_change = target_value,
    control_median_percent_change = median(null),
    control_min_percent_change = min(null),
    control_max_percent_change = max(null),
    controls_at_least_as_negative = sum(null <= target_value),
    empirical_lower_tail_p = (1 + sum(null <= target_value)) / (1 + length(null)),
    target_percentile_for_loss = rank(c(null, target_value), ties.method = "average")[length(null) + 1L] /
      (length(null) + 1L)
  )
}, by = program]
summary[, empirical_FDR_BH := p.adjust(empirical_lower_tail_p, method = "BH")]

fwrite(all_values, file.path(out_dir, "Figure7G_ITGB2_directional_connectivity_all_nodes.tsv"), sep = "\t")
fwrite(summary, file.path(out_dir, "Figure7G_ITGB2_directional_connectivity_summary.tsv"), sep = "\t")

manifest <- list(
  analysis = "WT-versus-KO programme-gene total absolute network-strength change",
  direction = "KO minus WT; negative values indicate connectivity loss",
  target = "ITGB2",
  matched_controls = controls,
  controls_frozen_before_directional_outcome = TRUE,
  n_controls = length(controls),
  programmes = program_ids,
  empirical_test = "lower-tail matched-node rank test",
  multiplicity = "BH across six programmes",
  interpretation_boundary = "network-connectivity change; not expression change, phenotypic suppression or therapeutic efficacy"
)
jsonlite::write_json(manifest, file.path(out_dir, "Figure7G_directional_connectivity_manifest.json"),
                     pretty = TRUE, auto_unbox = TRUE)
print(summary)
