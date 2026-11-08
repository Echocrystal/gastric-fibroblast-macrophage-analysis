#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else
  "I:/ECM-new/ECM-public"
legacy_lib <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026/runtime/R_library"
project_lib <- file.path(project_dir, "environment", "R_library")
.libPaths(c(project_lib, legacy_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

frozen_file <- file.path(project_dir, "models", "ecotype", "tcga_ecm_centroid_gene_contrasts.tsv")
marker_file <- file.path(project_dir, "gene_sets", "GSE183904_coarse_annotation_markers_locked.tsv")
program_file <- file.path(project_dir, "gene_sets", "scrna_spatial_programs_locked.tsv")
manifest_file <- file.path(project_dir, "manifests", "v4", "GSE270678_spatial_sample_manifest.tsv")
required <- c(frozen_file, marker_file, program_file, manifest_file)
if (!all(file.exists(required))) {
  stop("Missing required locked inputs: ", paste(required[!file.exists(required)], collapse = ", "))
}

model_dir <- file.path(project_dir, "models", "spatial", "v4")
audit_dir <- file.path(project_dir, "audit", "v4")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

frozen <- fread(frozen_file)
if (nrow(frozen) != 420L || uniqueN(frozen$gene_symbol) != 420L) {
  stop("Frozen ECM model is not the expected 420-gene model")
}
frozen <- frozen[order(-centroid_contrast_C2_minus_C1)]
ecm_high50 <- head(frozen$gene_symbol, 50L)
ecm_low50 <- tail(frozen$gene_symbol, 50L)

markers <- fread(marker_file)
programs <- fread(program_file)
required_markers <- c("Epithelial", "Fibroblast", "Endothelial", "T_cell")
required_programs <- c("T_CELL_CYTOTOXIC", "T_CELL_EXHAUSTION")
if (!all(required_markers %chin% markers$coarse_label)) stop("Missing a required locked coarse marker set")
if (!all(required_programs %chin% programs$program)) stop("Missing a required locked T-cell program")

gene_sets <- rbindlist(list(
  data.table(score = "FROZEN_ECM_C2_HIGH50", gene_symbol = ecm_high50,
             role = "ECM_C2_side", source = "frozen_TCGA_420_centroid_contrast"),
  data.table(score = "FROZEN_ECM_C1_HIGH50", gene_symbol = ecm_low50,
             role = "ECM_C1_side", source = "frozen_TCGA_420_centroid_contrast"),
  markers[coarse_label %chin% required_markers,
          .(score = paste0("MARKER_", coarse_label), gene_symbol,
            role = "cell_compartment_marker", source = "GSE183904_coarse_annotation_markers_locked")],
  programs[program %chin% required_programs,
           .(score = paste0("PROGRAM_", program), gene_symbol,
             role = "T_cell_state_program", source = "scrna_spatial_programs_locked")]
), fill = TRUE)
if (anyDuplicated(gene_sets[, .(score, gene_symbol)])) stop("Duplicate gene within a locked score")

endpoints <- data.table(
  endpoint_id = c(
    "P1_epi_high_vs_low_ECM_same_T",
    "P2_epi_high_vs_low_ECM_neighbor_T",
    "P3_ECM_specific_T_confinement_DiD",
    "P4_ECM_specific_cytotoxic_confinement_DiD",
    "S1_epi_ECM_same_T_rho",
    "S2_epi_ECM_neighbor_T_rho",
    "S3_allspot_ECM_T_partial_rho",
    "S4_allspot_ECM_neighbor_T_partial_rho",
    "S5_fibroblast_T_partial_rho",
    "S6_fibroblast_neighbor_T_partial_rho",
    "S7_endothelial_T_partial_rho",
    "S8_endothelial_neighbor_T_partial_rho"
  ),
  tier = c(rep("primary", 4), rep("secondary", 8)),
  outcome = c(
    "MARKER_T_cell", "lag_MARKER_T_cell", "MARKER_T_cell",
    "PROGRAM_T_CELL_CYTOTOXIC", "MARKER_T_cell", "lag_MARKER_T_cell",
    "MARKER_T_cell", "lag_MARKER_T_cell", "MARKER_T_cell",
    "lag_MARKER_T_cell", "MARKER_T_cell", "lag_MARKER_T_cell"
  ),
  estimand = c(
    "Within epithelial-rich spots: median T-cell score in top ECM quartile minus bottom ECM quartile",
    "Within epithelial-rich spots: median first-ring-neighbour T-cell score in top ECM quartile minus bottom ECM quartile",
    "[stromal-rich minus epithelial-rich T-cell gap] in high-ECM context minus the same gap in low-ECM context",
    "[stromal-rich minus epithelial-rich cytotoxic gap] in high-ECM context minus the same gap in low-ECM context",
    "Within epithelial-rich spots: Spearman ECM-T-cell correlation",
    "Within epithelial-rich spots: Spearman ECM-neighbour-T-cell correlation",
    "All-spot ECM-T-cell partial rank correlation adjusted for epithelial, fibroblast, and endothelial scores",
    "All-spot ECM-neighbour-T-cell partial rank correlation adjusted for epithelial, fibroblast, and endothelial scores",
    "Fibroblast-T-cell partial rank correlation adjusted for ECM and epithelial scores",
    "Fibroblast-neighbour-T-cell partial rank correlation adjusted for ECM and epithelial scores",
    "Endothelial-T-cell partial rank correlation adjusted for ECM, epithelial, and fibroblast scores",
    "Endothelial-neighbour-T-cell partial rank correlation adjusted for ECM, epithelial, and fibroblast scores"
  ),
  supportive_direction = c("negative", "negative", "positive", "positive",
                           "negative", "negative", NA, NA, "positive", "positive",
                           NA, NA),
  multiplicity_family = c(rep("four_primary_endpoints_BH", 4), rep("eight_secondary_axes_BH", 8))
)

definitions <- list(
  public_region_status = "source-study barcode-level LAR/SMSR/ICR/NCR not evaluable",
  analysis_label = "project-reconstructed expression compartments",
  primary_cohort = "33 main ST slides from 19 donors",
  primary_unit = "donor",
  replicate_cohort = "39 ST_R slides; sensitivity only; D30 extension-only excluded from primary cohort",
  matrix = "GEO filtered_feature_bc_matrix.h5 intersected with in_tissue=1 positions",
  qc = list(
    counts_lower = "10^(median(log10(counts+1)) - 4*MAD) - 1 within slide",
    features_lower = "10^(median(log10(features+1)) - 4*MAD) - 1 within slide",
    mitochondrial_upper = "median(percent_mito) + 5*MAD within slide",
    fallback = "MAD zero/nonfinite replaced by IQR/1.349; no fixed 200-feature threshold"
  ),
  normalisation = "log1p(CP10K), followed by within-slide gene-wise z standardisation",
  zero_variance_gene = "z score set to zero; coverage and variance reported",
  frozen_ECM_score = "mean z of fixed C2 high50 minus mean z of fixed C1 high50",
  marker_score = "mean within-slide gene z score for the locked marker/program",
  epithelial_rich = "top within-slide quartile of MARKER_Epithelial minus MARKER_Fibroblast",
  stromal_rich = "top within-slide quartile of MARKER_Fibroblast minus MARKER_Epithelial",
  high_ECM_context = "at or above within-slide median frozen ECM score",
  low_ECM_context = "below within-slide median frozen ECM score",
  ECM_extremes = "top and bottom within-slide quartiles",
  spatial_neighbor = "exact Visium first ring: (-1,-1),(-1,+1),(0,-2),(0,+2),(+1,-1),(+1,+1)",
  minimum_contrast_group_spots = 20L,
  slide_to_patient = "median slide effect per donor; Fisher-z mean for correlations before returning to rho",
  inference = "two-sided Wilcoxon signed-rank against zero across 19 donor-level estimates",
  interval = "10,000 bootstrap resamples of donors",
  primary_multiplicity = "BH across P1-P4",
  stage_lauren = "within-slide contrasts are primary; stage and Lauren direction-stratified sensitivity only when each stratum has at least five donors",
  author_region_DEG_use = "prohibited for classifying spots in the independent compartment analysis"
)

manifest <- fread(manifest_file)
if (manifest[platform_role == "ST", .N] != 33L ||
    uniqueN(manifest[platform_role == "ST", donor_id]) != 19L) {
  stop("Unexpected GSE270678 main-cohort manifest")
}

gene_set_file <- file.path(model_dir, "GSE270678_project_compartment_gene_sets_v4.tsv")
endpoint_file <- file.path(model_dir, "GSE270678_project_compartment_endpoints_v4.tsv")
model_file <- file.path(model_dir, "GSE270678_project_compartment_analysis_lock_v4.rds")
fwrite(gene_sets, gene_set_file, sep = "\t", quote = FALSE)
fwrite(endpoints, endpoint_file, sep = "\t", quote = FALSE)
saveRDS(list(
  version = "V4.1_pre_expression_lock",
  locked_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  definitions = definitions,
  endpoints = endpoints,
  gene_sets = gene_sets,
  outcome_or_expression_loaded = FALSE
), model_file, compress = "xz")

artifacts <- c(gene_set_file, endpoint_file, model_file, frozen_file, marker_file, program_file, manifest_file)
hash_table <- data.table(
  artifact = gsub("\\\\", "/", artifacts),
  sha256 = vapply(artifacts, digest, character(1), algo = "sha256", file = TRUE)
)
hash_file <- file.path(model_dir, "GSE270678_project_compartment_lock_hashes_v4.tsv")
fwrite(hash_table, hash_file, sep = "\t", quote = FALSE)

audit_lines <- c(
  "# GSE270678 project-compartment pre-expression lock",
  "",
  paste0("Locked: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Status",
  "",
  "The source-study barcode-level LAR/SMSR/ICR/NCR endpoint is not evaluable from the retrievable public artefacts. This lock defines a distinct project-reconstructed expression-compartment analysis before the GSE270678 expression TAR is opened or scored.",
  "",
  "## Non-circular inputs",
  "",
  "- Frozen TCGA ECM top-50 positive and bottom-50 negative directions.",
  "- Previously locked epithelial, fibroblast, endothelial, and T-cell coarse markers.",
  "- Previously locked cytotoxic and exhaustion programs.",
  "- Source-study region DE genes are prohibited from spot classification.",
  "",
  "## Primary endpoints",
  "",
  "1. Epithelial-rich top-versus-bottom ECM-quartile same-spot T-cell difference (supportive direction: negative).",
  "2. Epithelial-rich top-versus-bottom ECM-quartile first-ring-neighbour T-cell difference (negative).",
  "3. ECM-specific stromal-minus-epithelial T-cell confinement difference-in-differences (positive).",
  "4. ECM-specific stromal-minus-epithelial cytotoxic confinement difference-in-differences (positive).",
  "",
  "BH correction is applied across these four primary endpoints. All other correlations are secondary and corrected as a separate eight-axis family.",
  "",
  "## Unit and replication",
  "",
  "The donor is the inferential unit. Slide-level effects are first aggregated within each of the 19 main-cohort donors. ST_R slides are sensitivity data; they are not additional patients, and extension-only donor D30 is not added to the primary cohort.",
  "",
  "## Adjustment logic",
  "",
  "Primary contrasts are defined within slide using within-slide ranks/quantiles, limiting confounding by stage, Lauren type, sampling location, library depth, and global tissue-composition differences between slides. Donors then receive equal inferential weight. Stage and Lauren sensitivity summaries are descriptive and require at least five donors per stratum; no 19-patient high-dimensional adjustment model will be fit.",
  "",
  "## Guardrails",
  "",
  "Thresholds, marker genes, endpoint directions, group definitions, and multiplicity families will not be altered after inspecting results. Failure of a primary endpoint will be reported rather than rescued by author-region DEGs, alternate quartiles, spot-level pseudo-replication, or relabelling the reconstructed compartments as LAR/SMSR/ICR/NCR.",
  "",
  "Artifact hashes are stored in `models/spatial/v4/GSE270678_project_compartment_lock_hashes_v4.tsv`."
)
writeLines(
  audit_lines,
  file.path(audit_dir, "GSE270678_project_compartment_pre_expression_lock_v4.md"),
  useBytes = TRUE
)

message("Locked GSE270678 project-compartment analysis before expression data were loaded.")
message("Primary endpoints: 4; main donors: 19; main slides: 33.")
