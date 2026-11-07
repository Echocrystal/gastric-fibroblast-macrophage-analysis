#!/usr/bin/env Rscript

# Freeze outcome-blind sender controls for the SPON2 screen.  No SPON2 KO
# output is read; matching uses expression, detection, donor coverage and the
# WT network only.

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({ library(Matrix); library(data.table); library(jsonlite) })

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
input_dir <- file.path(project, "analysis", "processed", "virtual_knockout", "GSE206785")
result_dir <- file.path(project, "analysis", "results", "virtual_knockout", "GSE206785", "SPON2_sender_screen_v1")
resource_dir <- file.path(project, "resources", "gene_sets")
audit_dir <- file.path(project, "analysis", "audit")

matrix_file <- file.path(result_dir, "scTenifoldKnk_seed2026081011__SPON2.rds")
counts_file <- file.path(input_dir, "GSE206785_sender_CTHRC1_CAF_full_raw_counts.mtx.gz")
genes_file <- file.path(input_dir, "GSE206785_sender_CTHRC1_CAF_full_raw_counts_genes.tsv")
metadata_file <- file.path(input_dir, "GSE206785_sender_CTHRC1_CAF_full_raw_counts_metadata.tsv")
fit <- readRDS(matrix_file)
wt <- as.matrix(fit$tensorNetworks$WT)
genes <- rownames(wt)
counts <- readMM(gzfile(counts_file))
all_genes <- fread(genes_file)$gene_symbol
meta <- fread(metadata_file)
rownames(counts) <- all_genes
stopifnot(identical(genes, intersect(genes, all_genes)), ncol(counts) == nrow(meta))

program_table <- rbind(
  fread(file.path(resource_dir, "virtual_knockout_sanity_programs_axis_excluded.tsv"))[
    program %chin% c("GOBP_CELL_MATRIX_ADHESION", "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY",
                     "REACTOME_INTEGRIN_SIGNALING"), .(program, gene_symbol)],
  fread(file.path(resource_dir, "GSE183904_fine_state_modules_frozen.tsv"))[
    lineage == "Fibroblast" & module %chin% c("Activated_matrix", "Inflammatory_CAF",
                                               "Adventitial_quiescent", "Contractile"),
    .(program = module, gene_symbol)]
)
program_table <- unique(program_table[gene_symbol %chin% genes & gene_symbol != "SPON2"])
program_count <- program_table[, .(program_membership_count = uniqueN(program)), by = gene_symbol]

detect <- counts[genes, , drop = FALSE] > 0
metrics <- data.table(
  gene_symbol = genes,
  mean_raw_count = as.numeric(Matrix::rowMeans(counts[genes, , drop = FALSE])),
  detection_fraction = as.numeric(Matrix::rowMeans(detect)),
  mean_total_abs_strength = rowMeans(abs(wt)) + colMeans(abs(wt)),
  mean_nonzero_degree = rowMeans((abs(wt) > 0)) + colMeans((abs(wt) > 0))
)
patient_count <- vapply(seq_len(nrow(metrics)), function(i) {
  sum(vapply(unique(meta$Patient), function(p) {
    idx <- which(meta$Patient == p)
    any(detect[i, idx])
  }, logical(1)))
}, integer(1))
metrics[, patient_count := patient_count]
metrics <- merge(metrics, program_count, by = "gene_symbol", all.x = TRUE)
metrics[is.na(program_membership_count), program_membership_count := 0L]
metrics[, technical := grepl("^(MT-|RPL[0-9]|RPS[0-9]|HB[ABDEGQZ])", gene_symbol, ignore.case = TRUE)]
target <- metrics[gene_symbol == "SPON2"]
pool <- metrics[gene_symbol != "SPON2" & !technical & patient_count >= 3L &
                detection_fraction >= 0.05 &
                program_membership_count == target$program_membership_count]
stopifnot(nrow(pool) >= 20L)

transform_feature <- function(name, value) {
  if (name == "detection_fraction") {
    value <- pmin(1 - 1e-4, pmax(1e-4, value)); return(qlogis(value))
  }
  log1p(value)
}
match_vars <- c("mean_raw_count", "detection_fraction", "mean_total_abs_strength", "mean_nonzero_degree")
for (nm in match_vars) {
  combined <- c(target[[nm]], pool[[nm]])
  z <- transform_feature(nm, combined)
  center <- median(z); scale <- mad(z, center = center, constant = 1)
  if (!is.finite(scale) || scale <= 0) scale <- sd(z)
  target_z <- (z[1] - center) / scale
  pool[, (paste0("z_", nm)) := (z[-1] - center) / scale]
  pool[, (paste0("target_z_", nm)) := target_z]
}
zcols <- paste0("z_", match_vars); tzcols <- paste0("target_z_", match_vars)
pool[, feature_distance := sqrt(rowSums((as.matrix(.SD[, ..zcols]) - as.matrix(.SD[, ..tzcols]))^2))]
setorder(pool, feature_distance, gene_symbol)
controls <- pool[1:20]
controls[, control_rank := .I]
out <- file.path(resource_dir, "GSE206785_SPON2_sender_matched_controls.tsv")
fwrite(controls[, .(target = "SPON2", control_rank = .I, control_gene = gene_symbol,
                    mean_raw_count, detection_fraction, patient_count,
                    mean_total_abs_strength, mean_nonzero_degree,
                    program_membership_count, feature_distance)], out, sep = "\t")
write_json(list(target = "SPON2", eligible_pool = nrow(pool), controls_frozen = nrow(controls),
                selection_reads_KO_outcomes = FALSE, match_variables = match_vars,
                control_file = out),
           file.path(audit_dir, "GSE206785_SPON2_sender_controls.json"),
           pretty = TRUE, auto_unbox = TRUE)
cat("eligible pool:", nrow(pool), " controls:", nrow(controls), "\n")
print(controls[, .(control_rank, control_gene = gene_symbol, feature_distance)])
