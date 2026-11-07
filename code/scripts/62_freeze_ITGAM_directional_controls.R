#!/usr/bin/env Rscript

# Freeze 20 ITGAM controls using the same outcome-blind matching logic used for
# ITGB2. No ITGAM directional KO results are read by this script.

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(data.table))

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
feature_file <- file.path(project, "analysis/audit",
  "GSE183904_primary_receptors_v2_scTenifoldKnk_control_matching_universe.tsv")
out_dir <- file.path(project, "resources/gene_sets")
audit_dir <- file.path(project, "analysis/audit")
features <- fread(feature_file)
axis_genes <- c("CD74", "CCRL2", "CXCR4", "ITGB2", "LRP1", "ITGAM",
                "CMKLR1", "ANXA2", "TREM2", "PLAT", "ILK")
target <- "ITGAM"
target_row <- features[gene_symbol == target]
if (nrow(target_row) != 1L) stop("Missing ITGAM feature row")
match_variables <- c("mean_raw_count", "detection_fraction",
                     "mean_total_abs_strength", "mean_nonzero_degree")
pool <- features[
  !gene_symbol %chin% axis_genes & technical == FALSE & patient_count >= 3L &
    detection_fraction >= 0.05 &
    program_membership_count == target_row$program_membership_count
]
if (nrow(pool) < 20L) stop("Insufficient exact-membership ITGAM pool")

transform_feature <- function(name, value) {
  if (name == "detection_fraction") {
    value <- pmin(1 - 1e-4, pmax(1e-4, value))
    return(qlogis(value))
  }
  log1p(value)
}
for (name in match_variables) {
  combined <- c(target_row[[name]], pool[[name]])
  transformed <- transform_feature(name, combined)
  center <- median(transformed)
  scale <- mad(transformed, center = center, constant = 1)
  if (!is.finite(scale) || scale <= 0) scale <- sd(transformed)
  target_z <- (transformed[1] - center) / scale
  pool[, (paste0("z_", name)) := (transformed[-1] - center) / scale]
  pool[, (paste0("target_z_", name)) := target_z]
}
z_cols <- paste0("z_", match_variables)
target_z_cols <- paste0("target_z_", match_variables)
pool[, feature_distance := sqrt(rowSums(
  (as.matrix(.SD[, ..z_cols]) - as.matrix(.SD[, ..target_z_cols]))^2
))]
pool[, membership_distance := abs(program_membership_count - target_row$program_membership_count)]
pool[, total_match_distance := feature_distance + 2 * membership_distance]
setorder(pool, total_match_distance, gene_symbol)
controls <- pool[1:20]
controls[, `:=`(target = target, control_rank = .I,
                exact_program_membership_pool = TRUE)]
out_file <- file.path(out_dir,
  "GSE183904_ITGAM_directional_network_matched_controls.tsv")
fwrite(controls[, .(
  target, control_rank, control_gene = gene_symbol,
  mean_raw_count, detection_fraction, patient_count,
  mean_total_abs_strength, mean_nonzero_degree, strength_seed_CV,
  program_membership_count, feature_distance, membership_distance,
  total_match_distance, exact_program_membership_pool
)], out_file, sep = "\t")
audit <- list(
  target = target,
  eligible_exact_membership_pool = nrow(pool),
  controls_frozen = nrow(controls),
  selection_reads_directional_outcomes = FALSE,
  match_variables = match_variables,
  axis_genes_excluded = axis_genes,
  control_file = out_file
)
jsonlite::write_json(
  audit,
  file.path(audit_dir, "GSE183904_ITGAM_directional_controls.json"),
  pretty = TRUE, auto_unbox = TRUE
)
cat("Eligible pool:", nrow(pool), " Frozen controls:", nrow(controls), "\n")
print(controls[, .(control_rank, control_gene = gene_symbol, total_match_distance)])
