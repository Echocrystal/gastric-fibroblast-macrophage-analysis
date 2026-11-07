#!/usr/bin/env Rscript

# Freeze network/expression-matched controls for ITGB2 without reading any
# virtual-knockout differential-regulation or program-enrichment outcome.

options(stringsAsFactors = FALSE, warn = 1)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(digest)
})

input_dir <- file.path(project, "analysis", "processed", "virtual_knockout", "GSE183904",
                       "C1QC_APOE_macrophage_v1")
result_dir <- file.path(project, "analysis", "results", "virtual_knockout", "GSE183904",
                        "C1QC_APOE_macrophage_v1")
resource_dir <- file.path(project, "resources", "gene_sets")
audit_dir <- file.path(project, "analysis", "audit")
matrix_file <- file.path(input_dir, "GSE183904_C1QC_APOE_sanity_counts.rds")
universe_file <- file.path(resource_dir, "GSE183904_C1QC_APOE_scTenifoldKnk_universe.tsv")
program_file <- file.path(resource_dir, "myeloid_receiver_programs_frozen.tsv")

x <- readRDS(matrix_file)
universe <- fread(universe_file)
focus_programs <- c(
  "GOBP_PHAGOCYTOSIS",
  "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
  "ID_INFLAMMATORY_MONO",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "GOBP_CELL_MATRIX_ADHESION",
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX"
)
program_manifest <- fread(program_file)[program %chin% focus_programs]
target <- "ITGB2"
axis_genes <- c("ANXA2", "ITGAM", "ITGB2", "LRP1", "PLAT", "ILK")
seeds <- c(2026072739L, 2026072740L, 2026072741L)

wt_rows <- list()
for (seed in seeds) {
  fit_file <- file.path(result_dir,
                        paste0("scTenifoldKnk_seed", seed, "__ITGB2.rds"))
  if (!file.exists(fit_file)) stop("Missing completed ITGB2 fit for WT matching")
  fit <- readRDS(fit_file)
  wt <- as.matrix(fit$tensorNetworks$WT)
  if (!identical(rownames(wt), rownames(x))) stop("WT/count gene-order mismatch")
  wt_rows[[as.character(seed)]] <- data.table(
    gene_symbol = rownames(wt),
    seed = seed,
    total_abs_strength = rowSums(abs(wt)) + colSums(abs(wt)),
    nonzero_degree = rowSums(abs(wt) > 0) + colSums(abs(wt) > 0)
  )
  rm(fit, wt)
}
wt_features <- rbindlist(wt_rows)[, .(
  mean_total_abs_strength = mean(total_abs_strength),
  mean_nonzero_degree = mean(nonzero_degree),
  strength_seed_CV = sd(total_abs_strength) / mean(total_abs_strength)
), by = gene_symbol]

features <- data.table(
  gene_symbol = rownames(x),
  mean_raw_count = as.numeric(Matrix::rowMeans(x)),
  detection_fraction = as.numeric(Matrix::rowMeans(x > 0))
)
features <- merge(features, wt_features, by = "gene_symbol")
features <- merge(features, universe[, .(gene_symbol, patient_count, technical)],
                  by = "gene_symbol", all.x = TRUE)
membership <- program_manifest[, .(program_membership_count = uniqueN(program)),
                               by = gene_symbol]
features <- merge(features, membership, by = "gene_symbol", all.x = TRUE)
features[is.na(program_membership_count), program_membership_count := 0L]

target_row <- features[gene_symbol == target]
if (nrow(target_row) != 1L) stop("ITGB2 feature row missing")
target_membership <- target_row$program_membership_count

pool <- features[
  !gene_symbol %chin% axis_genes &
    technical == FALSE &
    patient_count >= 3L &
    detection_fraction >= 0.05 &
    program_membership_count == target_membership
]
if (nrow(pool) < 20L) {
  pool <- features[
    !gene_symbol %chin% axis_genes &
      technical == FALSE &
      patient_count >= 3L &
      detection_fraction >= 0.05 &
      abs(program_membership_count - target_membership) <= 1L
  ]
}
if (nrow(pool) < 20L) stop("Insufficient matched-control pool")

match_variables <- c("mean_raw_count", "detection_fraction",
                     "mean_total_abs_strength", "mean_nonzero_degree")
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
  if (!is.finite(scale) || scale <= 0) stop("Degenerate match feature: ", name)
  target_z <- (transformed[1] - center) / scale
  pool[, (paste0("z_", name)) := (transformed[-1] - center) / scale]
  pool[, (paste0("target_z_", name)) := target_z]
}
z_cols <- paste0("z_", match_variables)
target_z_cols <- paste0("target_z_", match_variables)
pool[, feature_distance := sqrt(rowSums(
  (as.matrix(.SD[, ..z_cols]) - as.matrix(.SD[, ..target_z_cols]))^2
))]
pool[, membership_distance := abs(program_membership_count - target_membership)]
pool[, total_match_distance := feature_distance + 2 * membership_distance]
setorder(pool, total_match_distance, gene_symbol)
controls <- pool[1:20]
controls[, `:=`(target = target, control_rank = .I)]

control_file <- file.path(resource_dir,
                          "GSE183904_ITGB2_scTenifoldKnk_network_matched_controls.tsv")
fwrite(controls[, .(
  target, control_rank, control_gene = gene_symbol,
  mean_raw_count, detection_fraction, patient_count,
  mean_total_abs_strength, mean_nonzero_degree, strength_seed_CV,
  program_membership_count, feature_distance, membership_distance,
  total_match_distance
)], control_file, sep = "\t")
fwrite(features,
       file.path(audit_dir, "GSE183904_ITGB2_scTenifoldKnk_control_matching_universe.tsv"),
       sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE183904",
  state = "C1QC_APOE_macrophage",
  target = target,
  selection_reads_KO_outcomes = FALSE,
  WT_sources = paste0("same-seed WT network from ITGB2 fit: ", seeds),
  match_variables = match_variables,
  transformations = c(
    mean_raw_count = "log1p",
    detection_fraction = "logit after [1e-4,1-1e-4] clipping",
    mean_total_abs_strength = "log1p",
    mean_nonzero_degree = "log1p"
  ),
  scaling = "median and unscaled MAD across target plus eligible pool",
  program_membership_matching = paste(
    "exact number of memberships among the six axis-excluded frozen programs when at least 20 genes are available;",
    "otherwise within one membership with distance penalty 2"
  ),
  target_features = as.list(target_row),
  eligible_pool = nrow(pool),
  controls = controls$gene_symbol,
  control_file = control_file,
  control_file_sha256 = digest(control_file, algo = "sha256", file = TRUE),
  claim_boundary = paste(
    "Controls match observed expression and inferred WT topology within the compact universe.",
    "They do not turn virtual knockout into biological causality."
  )
)
write_json(audit,
           file.path(audit_dir, "GSE183904_ITGB2_scTenifoldKnk_network_matched_controls.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("ITGB2 target features:\n")
print(target_row)
cat("\nFrozen controls:\n")
print(controls[, .(control_rank, gene_symbol, mean_raw_count, detection_fraction,
                   mean_total_abs_strength, mean_nonzero_degree,
                   program_membership_count, total_match_distance)])
