#!/usr/bin/env Rscript

# Exploratory sender-side scTenifoldKnk screen for SPON2 in the
# CTHRC1-activated fibroblast state of GSE206785.  This analysis tests the
# intracellular network sensitivity of the sender state; it does not simulate
# propagation of a ligand knockout into macrophage cells.

options(stringsAsFactors = FALSE, warn = 1)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(project, "runtime", "R_virtual_knockout"),
            "I:/ECM-new/scRNA/20.rebuild_GSE163558/R_libs_4.4",
            file.path(source_project, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(scTenifoldKnk)
  library(scTenifoldNet)
  library(jsonlite)
})

input_dir <- file.path(project, "analysis", "processed", "virtual_knockout", "GSE206785")
result_dir <- file.path(project, "analysis", "results", "virtual_knockout", "GSE206785", "SPON2_sender_screen_v1")
audit_dir <- file.path(project, "analysis", "audit")
resource_dir <- file.path(project, "resources", "gene_sets")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- file.path(input_dir, "GSE206785_sender_CTHRC1_CAF_full_raw_counts.mtx.gz")
genes_file <- file.path(input_dir, "GSE206785_sender_CTHRC1_CAF_full_raw_counts_genes.tsv")
cells_file <- file.path(input_dir, "GSE206785_sender_CTHRC1_CAF_full_raw_counts_cells.tsv")
metadata_file <- file.path(input_dir, "GSE206785_sender_CTHRC1_CAF_full_raw_counts_metadata.tsv")

genes <- fread(genes_file)$gene_symbol
cells <- fread(cells_file)$cell_id
metadata <- fread(metadata_file)
x_full <- readMM(gzfile(counts_file))
rownames(x_full) <- genes
colnames(x_full) <- cells
stopifnot(nrow(x_full) == length(genes), ncol(x_full) == nrow(metadata),
          "SPON2" %in% rownames(x_full), uniqueN(metadata$Patient) == 17L)

# Freeze the input gene universe before reading any SPON2 perturbation output.
base_universe <- fread(file.path(resource_dir, "virtual_knockout_sanity_gene_universe.tsv"))$gene_symbol
fibro_module_table <- fread(file.path(resource_dir, "GSE183904_fine_state_modules_frozen.tsv"))[
  lineage == "Fibroblast" & module %chin% c("Activated_matrix", "Inflammatory_CAF",
                                             "Adventitial_quiescent", "Contractile")]
fibro_module_genes <- fibro_module_table$gene_symbol
universe <- sort(unique(c(base_universe, fibro_module_genes, "SPON2")))
universe <- intersect(universe, rownames(x_full))
stopifnot("SPON2" %in% universe, length(universe) >= 700L)
x <- x_full[universe, , drop = FALSE]

programs <- list(
  Activated_matrix = intersect(fibro_module_table[module == "Activated_matrix", gene_symbol], universe),
  Inflammatory_CAF = intersect(fibro_module_table[module == "Inflammatory_CAF", gene_symbol], universe),
  Adventitial_quiescent = intersect(fibro_module_table[module == "Adventitial_quiescent", gene_symbol], universe),
  Contractile = intersect(fibro_module_table[module == "Contractile", gene_symbol], universe),
  GOBP_CELL_MATRIX_ADHESION = intersect(
    fread(file.path(resource_dir, "virtual_knockout_sanity_programs_axis_excluded.tsv"))[
      program == "GOBP_CELL_MATRIX_ADHESION", gene_symbol], universe),
  GOBP_INTEGRIN_MEDIATED_SIGNALING = intersect(
    fread(file.path(resource_dir, "virtual_knockout_sanity_programs_axis_excluded.tsv"))[
      program == "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY", gene_symbol], universe),
  REACTOME_INTEGRIN_SIGNALING = intersect(
    fread(file.path(resource_dir, "virtual_knockout_sanity_programs_axis_excluded.tsv"))[
      program == "REACTOME_INTEGRIN_SIGNALING", gene_symbol], universe)
)
programs <- lapply(programs, setdiff, y = "SPON2")
if (any(vapply(programs, length, integer(1)) < 5L)) stop("A sender programme has fewer than five genes")

seeds <- c(2026081011L, 2026081012L, 2026081013L)
parameters <- list(
  qc = FALSE, nc_lambda = 0, nc_nNet = 3L, nc_nCells = 200L,
  nc_nComp = 3L, nc_scaleScores = TRUE, nc_symmetric = FALSE,
  nc_q = 0.90, td_K = 3L, td_maxIter = 1000L, td_maxError = 1e-5,
  td_nDecimal = 3L, ma_nDim = 2L, nCores = 4L
)

enrich_one <- function(dr, seed) {
  dr <- as.data.table(dr)[gene != "SPON2"]
  dr[, distance_percentile := frank(distance, ties.method = "average") / .N]
  rows <- lapply(names(programs), function(program_id) {
    g <- programs[[program_id]]
    in_set <- dr$gene %in% g
    n_in <- sum(in_set); n_out <- sum(!in_set)
    if (n_in < 5L || n_out < 5L) {
      return(data.table(seed = seed, program = program_id, n_in_universe = n_in,
                        enrichment_auc = NA_real_, wilcoxon_p = NA_real_))
    }
    wt <- suppressWarnings(wilcox.test(dr$distance[in_set], dr$distance[!in_set],
                                       alternative = "greater", exact = FALSE))
    data.table(seed = seed, program = program_id, n_in_universe = n_in,
               enrichment_auc = as.numeric(wt$statistic) / (n_in * n_out),
               wilcoxon_p = wt$p.value)
  })
  out <- rbindlist(rows)
  out[, wilcoxon_fdr_within_seed := p.adjust(wilcoxon_p, "BH")]
  out[]
}

dr_all <- list(); program_all <- list(); summary_all <- list()
for (seed in seeds) {
  message("Running SPON2 sender KO, seed ", seed)
  set.seed(seed)
  started <- Sys.time()
  fit <- scTenifoldKnk::scTenifoldKnk(
    countMatrix = x, qc = parameters$qc, gKO = "SPON2",
    nc_lambda = parameters$nc_lambda, nc_nNet = parameters$nc_nNet,
    nc_nCells = parameters$nc_nCells, nc_nComp = parameters$nc_nComp,
    nc_scaleScores = parameters$nc_scaleScores, nc_symmetric = parameters$nc_symmetric,
    nc_q = parameters$nc_q, td_K = parameters$td_K,
    td_maxIter = parameters$td_maxIter, td_maxError = parameters$td_maxError,
    td_nDecimal = parameters$td_nDecimal, ma_nDim = parameters$ma_nDim,
    nCores = parameters$nCores)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  wt <- as.matrix(fit$tensorNetworks$WT)
  dr <- as.data.table(fit$diffRegulation)
  dr[, `:=`(target = "SPON2", seed = seed,
            rank_by_distance = frank(-distance, ties.method = "min"),
            out_strength_WT = rowSums(abs(wt))[gene],
            in_strength_WT = colSums(abs(wt))[gene],
            mean_raw_count = Matrix::rowMeans(x)[gene],
            detection_fraction = Matrix::rowMeans(x > 0)[gene])]
  dr_all[[as.character(seed)]] <- dr
  program_all[[as.character(seed)]] <- enrich_one(dr, seed)
  tr <- dr[gene == "SPON2"]
  summary_all[[as.character(seed)]] <- data.table(
    target = "SPON2", seed = seed, cells = ncol(x), genes = nrow(x),
    elapsed_seconds = elapsed, target_distance = tr$distance,
    target_rank = tr$rank_by_distance,
    target_detection_fraction = tr$detection_fraction,
    target_mean_raw_count = tr$mean_raw_count,
    n_diffreg_fdr_0_05 = sum(dr$p.adj < 0.05, na.rm = TRUE),
    n_diffreg_fdr_0_10 = sum(dr$p.adj < 0.10, na.rm = TRUE))
  saveRDS(fit, file.path(result_dir, paste0("scTenifoldKnk_seed", seed, "__SPON2.rds")), compress = "xz")
  rm(fit, wt); invisible(gc())
}

fwrite(rbindlist(dr_all), file.path(result_dir, "SPON2_sender_all_diffreg.tsv"), sep = "\t")
fwrite(rbindlist(program_all), file.path(result_dir, "SPON2_sender_program_enrichment_by_seed.tsv"), sep = "\t")
fwrite(rbindlist(summary_all), file.path(result_dir, "SPON2_sender_run_summary.tsv"), sep = "\t")
manifest <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE206785", state = "F13_Activated_CTHRC1_fibroblast",
  target = "SPON2", patients = uniqueN(metadata$Patient), cells = ncol(x),
  universe_genes = nrow(x), seeds = seeds, parameters = parameters,
  interpretation = "sender-side intracellular network sensitivity; no receiver propagation",
  controls_not_run = TRUE)
write_json(manifest, file.path(audit_dir, "GSE206785_SPON2_sender_screen_manifest.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)
print(rbindlist(summary_all))
print(rbindlist(program_all)[order(wilcoxon_p)])
