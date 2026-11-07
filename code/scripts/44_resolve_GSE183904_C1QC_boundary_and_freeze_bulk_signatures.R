#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
options(error = function() {
  traceback(20L)
  quit(save = "no", status = 1L, runLast = FALSE)
})
set.seed(2026072944L)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(edgeR)
  library(jsonlite)
  library(digest)
})

manifest_file <- file.path(source_project, "data", "processed", "scrna", "GSE183904",
                           "GSE183904_sample_patient_clinical_manifest.tsv")
sce_dir <- file.path(source_project, "data", "intermediate", "scrna", "GSE183904", "sample_sce")
fib_file <- file.path(project, "analysis", "processed", "GSE183904_fine_states",
                      "GSE183904_Fibroblast_fine_state_cell_metadata.tsv.gz")
my_file <- file.path(project, "analysis", "processed", "GSE183904_fine_states",
                     "GSE183904_Myeloid_fine_state_cell_metadata.tsv.gz")
feature_file <- file.path(project, "analysis", "results", "PLAT_axis", "GSE183904",
                          "fine_state_axis", "GSE183904_fine_state_axis_patient_features.tsv")
count_file <- file.path(project, "analysis", "results", "PLAT_axis", "GSE183904",
                        "fine_state_axis", "GSE183904_fine_state_counts_by_patient.tsv")
program_state_file <- file.path(project, "analysis", "results", "PLAT_axis", "GSE183904",
                                "receiver_state_characterization",
                                "GSE183904_myeloid_program_scores_by_patient_state.tsv")
tcga_file <- "I:/ECM-new/ECM-public/data/processed/bulk/tcga_stad_full_symbol_log2cpm_412.rds"
prereg_file <- file.path(project, "analysis", "audit", "TCGA_PLAT_ecosystem_analysis_preregistration.json")
out_dir <- file.path(project, "analysis", "results", "PLAT_axis", "GSE183904",
                     "state_boundary_and_signatures")
audit_dir <- file.path(project, "analysis", "audit")
signature_file <- file.path(project, "resources", "gene_sets",
                            "GSE183904_PLAT_ecosystem_bulk_signatures.tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(manifest_file, fib_file, my_file, feature_file, count_file,
            program_state_file, tcga_file, prereg_file)
if (!all(file.exists(inputs))) stop("Missing boundary/signature input")

manifest <- fread(manifest_file)[include_primary_replication == TRUE &
                                   specimen_class == "Primary_Tumor"]
fib <- fread(fib_file)
my <- fread(my_file)
features <- fread(feature_file)
state_counts <- fread(count_file)
program_state <- fread(program_state_file)
if (nrow(manifest) != 26L || uniqueN(manifest$Patient_ID) != 26L) {
  stop("Expected 26 independent primary tumors")
}

tcga_obj <- readRDS(tcga_file)
tcga_genes <- rownames(tcga_obj$expression)
if (length(tcga_genes) < 10000L || ncol(tcga_obj$expression) != 412L) {
  stop("Unexpected TCGA expression object")
}

group_definitions <- list(
  fib_c1 = list(lineage = "fib", clusters = 1L),
  fib_c4 = list(lineage = "fib", clusters = 4L),
  my_c3 = list(lineage = "my", clusters = 3L),
  my_c5 = list(lineage = "my", clusters = 5L),
  my_c35 = list(lineage = "my", clusters = c(3L, 5L)),
  my_c145 = list(lineage = "my", clusters = c(1L, 4L, 5L)),
  my_c14 = list(lineage = "my", clusters = c(1L, 4L))
)

group_counts <- setNames(vector("list", length(group_definitions)), names(group_definitions))
group_meta <- list()
reference_genes <- NULL

for (i in seq_len(nrow(manifest))) {
  info <- manifest[i]
  message(sprintf("[%d/26] aggregate state pseudobulks: %s", i, info$sample_id))
  sce <- readRDS(file.path(sce_dir, paste0(info$sample_id, "_secondary_qc_sce.rds")))
  counts_mat <- counts(sce)
  if (is.null(reference_genes)) reference_genes <- rownames(counts_mat)
  if (!identical(reference_genes, rownames(counts_mat))) stop("Gene-order mismatch")
  fib_i <- fib[sample_id == info$sample_id]
  my_i <- my[sample_id == info$sample_id]
  for (group_id in names(group_definitions)) {
    def <- group_definitions[[group_id]]
    cm <- if (def$lineage == "fib") fib_i else my_i
    cells <- cm[cluster %in% def$clusters, cell_id]
    idx <- match(cells, colnames(counts_mat))
    if (anyNA(idx)) stop("Cell/count alignment failure: ", info$sample_id, " ", group_id)
    sums <- if (length(idx)) Matrix::rowSums(counts_mat[, idx, drop = FALSE]) else
      rep(0, nrow(counts_mat))
    group_counts[[group_id]][[info$sample_id]] <- as.numeric(sums)
    group_meta[[length(group_meta) + 1L]] <- data.table(
      sample_id = info$sample_id,
      Patient_ID = info$Patient_ID,
      group_id = group_id,
      n_cells = length(idx)
    )
  }
  rm(sce, counts_mat)
  gc(verbose = FALSE)
}

group_meta <- rbindlist(group_meta)
for (group_id in names(group_counts)) {
  mat <- do.call(cbind, group_counts[[group_id]])
  rownames(mat) <- reference_genes
  colnames(mat) <- names(group_counts[[group_id]])
  group_counts[[group_id]] <- mat
}
fwrite(group_meta, file.path(out_dir, "GSE183904_boundary_group_cell_counts.tsv"), sep = "\t")

technical_exclude <- function(g) {
  grepl("^(MT-|RPL|RPS|HB[ABDEGQZ]|MKI67$|TOP2A$|TYMS$|UBE2C$|CENPF$)", g)
}

paired_de <- function(target_group, comparator_group, contrast_id,
                      min_target = 5L, min_comparator = 5L) {
  mt <- group_meta[group_id == target_group]
  mc <- group_meta[group_id == comparator_group]
  paired <- merge(mt[, .(sample_id, Patient_ID, n_target = n_cells)],
                  mc[, .(sample_id, n_comparator = n_cells)], by = "sample_id")
  paired <- paired[n_target >= min_target & n_comparator >= min_comparator]
  if (nrow(paired) < 10L) stop("Insufficient paired patients for ", contrast_id)
  samples <- paired$sample_id
  target <- group_counts[[target_group]][, samples, drop = FALSE]
  comparator <- group_counts[[comparator_group]][, samples, drop = FALSE]
  colnames(comparator) <- paste0("comparator_", samples)
  colnames(target) <- paste0("target_", samples)
  counts_all <- cbind(comparator, target)
  coldata <- data.table(
    Patient_ID = factor(rep(paired$Patient_ID, 2L)),
    state = factor(rep(c("comparator", "target"), each = nrow(paired)),
                   levels = c("comparator", "target"))
  )
  design <- model.matrix(~ Patient_ID + state, data = coldata)
  dge <- DGEList(counts = counts_all)
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design, robust = TRUE)
  fit <- glmQLFit(dge, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "statetarget")
  tab <- as.data.table(topTags(qlf, n = Inf, sort.by = "none")$table, keep.rownames = "gene_symbol")
  setnames(tab, c("logFC", "logCPM", "F", "PValue", "FDR"),
           c("logFC_target_minus_comparator", "mean_logCPM", "QLF", "p_value", "FDR"))

  lcpm <- cpm(dge, log = TRUE, prior.count = 0.5)
  n <- nrow(paired)
  comp_expr <- lcpm[, seq_len(n), drop = FALSE]
  target_expr <- lcpm[, n + seq_len(n), drop = FALSE]
  extra <- data.table(
    gene_symbol = rownames(lcpm),
    paired_direction_consistency = rowMeans(target_expr > comp_expr),
    target_patient_detection_fraction = rowMeans(target[rownames(lcpm), , drop = FALSE] > 0),
    comparator_patient_detection_fraction = rowMeans(comparator[rownames(lcpm), , drop = FALSE] > 0),
    target_median_logCPM = apply(target_expr, 1L, median),
    comparator_median_logCPM = apply(comp_expr, 1L, median)
  )
  tab <- merge(tab, extra, by = "gene_symbol", all.x = TRUE)
  tab[, `:=`(
    contrast_id = contrast_id,
    target_group = target_group,
    comparator_group = comparator_group,
    paired_patients = nrow(paired),
    available_in_TCGA = gene_symbol %in% tcga_genes,
    technical_excluded = technical_exclude(gene_symbol)
  )]
  setorder(tab, FDR, -logFC_target_minus_comparator)
  fwrite(tab, file.path(out_dir, paste0("GSE183904_", contrast_id, "_paired_pseudobulk_DE.tsv")), sep = "\t")
  list(table = tab, paired = paired, lcpm = lcpm)
}

boundary <- paired_de("my_c3", "my_c5", "C1QC_APOE_vs_HLAII_C1QC", 5L, 5L)
sender <- paired_de("fib_c1", "fib_c4", "CXCL14_POSTN_PLAT_vs_activated_matrix", 5L, 5L)
receiver_specific <- paired_de("my_c3", "my_c145", "C1QC_APOE_vs_other_recurrent", 10L, 10L)
receiver_merged <- paired_de("my_c35", "my_c14", "merged_C1QC_vs_inflammatory_recurrent", 10L, 10L)

axis_genes <- c("PLAT", "ANXA2", "ITGAM", "ITGB2", "LRP1", "ILK")
freeze_signature <- function(de, signature_id, allow_axis = FALSE, n_max = 30L) {
  eligible <- de[
    FDR <= 0.01 & logFC_target_minus_comparator >= 0.5 &
      target_patient_detection_fraction >= 0.5 &
      paired_direction_consistency >= 0.65 &
      available_in_TCGA & !technical_excluded
  ]
  if (!allow_axis) eligible <- eligible[!gene_symbol %in% axis_genes]
  setorder(eligible, FDR, -paired_direction_consistency, -logFC_target_minus_comparator)
  selected <- head(eligible, n_max)
  if (nrow(selected) < 15L) stop("Fewer than 15 eligible genes for ", signature_id)
  selected[, `:=`(
    signature_id = signature_id,
    signature_rank = seq_len(.N),
    axis_genes_allowed = allow_axis,
    derivation_status = "frozen_before_TCGA_outcome_analysis"
  )]
  selected
}

sig_sender_primary <- freeze_signature(sender$table, "sender_CXCL14_POSTN_state_axis_excluded", FALSE)
sig_sender_inclusive <- freeze_signature(sender$table, "sender_CXCL14_POSTN_PLAT_state_inclusive", TRUE)
if ("PLAT" %in% sender$table$gene_symbol &&
    !"PLAT" %in% sig_sender_inclusive$gene_symbol) {
  plat_row <- sender$table[gene_symbol == "PLAT"]
  if (nrow(plat_row) == 1L && plat_row$FDR <= 0.01 &&
      plat_row$logFC_target_minus_comparator >= 0.5 &&
      plat_row$target_patient_detection_fraction >= 0.5 &&
      plat_row$paired_direction_consistency >= 0.65 &&
      plat_row$available_in_TCGA && !plat_row$technical_excluded) {
    retained_n <- min(29L, nrow(sig_sender_inclusive))
    sig_sender_inclusive <- rbind(head(sig_sender_inclusive, retained_n),
                                  plat_row, fill = TRUE)
    setorder(sig_sender_inclusive, FDR, -paired_direction_consistency,
             -logFC_target_minus_comparator)
    sig_sender_inclusive[, `:=`(
      signature_id = "sender_CXCL14_POSTN_PLAT_state_inclusive",
      signature_rank = seq_len(.N), axis_genes_allowed = TRUE,
      derivation_status = "frozen_before_TCGA_outcome_analysis"
    )]
  }
}
sig_receiver_specific <- freeze_signature(receiver_specific$table,
                                          "receiver_C1QC_APOE_specific_axis_excluded", FALSE)
sig_receiver_merged <- freeze_signature(receiver_merged$table,
                                        "receiver_merged_C1QC_axis_excluded", FALSE)
signatures <- rbindlist(list(sig_sender_primary, sig_sender_inclusive,
                            sig_receiver_specific, sig_receiver_merged), fill = TRUE)
fwrite(signatures, signature_file, sep = "\t")
fwrite(signatures, file.path(out_dir, "GSE183904_PLAT_ecosystem_bulk_signatures_frozen.tsv"), sep = "\t")

# Quantify transcriptomic similarity of clusters 3 and 5 using paired TMM logCPM.
n_boundary <- nrow(boundary$paired)
boundary_cor <- vapply(seq_len(n_boundary), function(i) {
  suppressWarnings(cor(boundary$lcpm[, i], boundary$lcpm[, n_boundary + i],
                       method = "spearman"))
}, numeric(1))
boundary_similarity <- data.table(
  sample_id = boundary$paired$sample_id,
  Patient_ID = boundary$paired$Patient_ID,
  spearman_all_expressed_genes = boundary_cor
)
fwrite(boundary_similarity,
       file.path(out_dir, "GSE183904_C1QC_cluster3_cluster5_patient_expression_similarity.tsv"), sep = "\t")

# Program comparison for cluster 3, cluster 5 and their merged sensitivity state.
c3_label <- "C1QC_APOE_macrophage"
c5_label <- "HLAII_C1QC_macrophage"
other_labels <- c("FCN1_AREG_inflammatory_monocyte", "IL1B_CXCL8_inflammatory_myeloid")
c3 <- program_state[refined_label == c3_label,
                    .(sample_id, Patient_ID, program, c3_cells = n_cells, c3_mean = mean_program_score)]
c5 <- program_state[refined_label == c5_label,
                    .(sample_id, Patient_ID, program, c5_cells = n_cells, c5_mean = mean_program_score)]
direct_program <- merge(c3, c5, by = c("sample_id", "Patient_ID", "program"))
direct_program[, delta_c3_minus_c5 := c3_mean - c5_mean]
direct_summary <- direct_program[, {
  d <- delta_c3_minus_c5
  wt <- suppressWarnings(wilcox.test(d, mu = 0, exact = FALSE))
  .(n_patients = .N, median_delta = median(d), positive_patients = sum(d > 0),
    negative_patients = sum(d < 0), paired_wilcoxon_p = wt$p.value)
}, by = program]
direct_summary[, FDR_BH := p.adjust(paired_wilcoxon_p, "BH")]
fwrite(direct_program, file.path(out_dir, "GSE183904_cluster3_vs_cluster5_program_paired.tsv"), sep = "\t")
fwrite(direct_summary, file.path(out_dir, "GSE183904_cluster3_vs_cluster5_program_summary.tsv"), sep = "\t")

merged <- merge(c3, c5, by = c("sample_id", "Patient_ID", "program"))
merged[, `:=`(
  target_cells = c3_cells + c5_cells,
  target_mean = (c3_cells * c3_mean + c5_cells * c5_mean) / (c3_cells + c5_cells)
)]
others <- program_state[refined_label %in% other_labels, .(
  comparator_cells = sum(n_cells),
  comparator_mean = weighted.mean(mean_program_score, n_cells)
), by = .(sample_id, Patient_ID, program)]
merged_program <- merge(merged[, .(sample_id, Patient_ID, program, target_cells, target_mean)],
                        others, by = c("sample_id", "Patient_ID", "program"))
merged_program[, delta_merged_minus_inflammatory := target_mean - comparator_mean]
merged_program_summary <- merged_program[, {
  d <- delta_merged_minus_inflammatory
  wt <- suppressWarnings(wilcox.test(d, mu = 0, exact = FALSE))
  nz <- sum(d != 0)
  sp <- if (nz) binom.test(sum(d > 0), nz, p = 0.5)$p.value else NA_real_
  .(n_patients = .N, median_delta = median(d), positive_patients = sum(d > 0),
    negative_patients = sum(d < 0), paired_wilcoxon_p = wt$p.value, exact_sign_p = sp)
}, by = program]
merged_program_summary[, `:=`(
  paired_wilcoxon_FDR_BH = p.adjust(paired_wilcoxon_p, "BH"),
  exact_sign_FDR_BH = p.adjust(exact_sign_p, "BH")
)]
fwrite(merged_program, file.path(out_dir, "GSE183904_merged_C1QC_program_paired.tsv"), sep = "\t")
fwrite(merged_program_summary, file.path(out_dir, "GSE183904_merged_C1QC_program_summary.tsv"), sep = "\t")

# Corrected patient-level abundance association using merged clusters 3+5.
sender_counts <- state_counts[lineage == "Fibroblast" &
                                refined_label == "CXCL14_POSTN_PLAT_rich_fibroblast",
                              .(sample_id, sender_cells = state_cells,
                                fibroblast_cells_check = lineage_cells)]
receiver_counts <- state_counts[lineage == "Myeloid" &
                                  refined_label %in% c(c3_label, c5_label),
                                .(receiver_cells = sum(state_cells),
                                  myeloid_cells_check = unique(lineage_cells)[1L]),
                                by = sample_id]
dt <- merge(features, sender_counts, by = "sample_id", all.x = TRUE)
dt <- merge(dt, receiver_counts, by = "sample_id", all.x = TRUE)
dt[is.na(sender_cells), sender_cells := 0]
dt[is.na(receiver_cells), receiver_cells := 0]
pc <- 0.5
dt[, sender_logratio := log((sender_cells + pc) /
                              (fibroblast_cells - sender_cells + pc))]
dt[, receiver_logratio_merged := log((receiver_cells + pc) /
                                       (myeloid_cells - receiver_cells + pc))]
dt[, sender_asin := asin(sqrt(sender_cells / fibroblast_cells))]
dt[, receiver_asin_merged := asin(sqrt(receiver_cells / myeloid_cells))]
covariates <- cbind(ECM_direction = dt$ECM_direction,
                    log_fibroblast_cells = log1p(dt$fibroblast_cells),
                    log_myeloid_cells = log1p(dt$myeloid_cells))

hc3_fit <- function(y, X) {
  fit <- lm.fit(X, y)
  xtx_inv <- chol2inv(chol(crossprod(X)))
  h <- rowSums((X %*% xtx_inv) * X)
  meat <- crossprod(X, X * (fit$residuals / pmax(1e-8, 1 - h))^2)
  vc <- xtx_inv %*% meat %*% xtx_inv
  se <- sqrt(diag(vc))[2]
  tv <- fit$coefficients[2] / se
  df <- nrow(X) - ncol(X)
  list(beta = fit$coefficients[2], se = se, t = tv, df = df,
       p = 2 * pt(-abs(tv), df = df), leverage = h)
}

analyze_assoc <- function(x, y, representation, B = 9999L, seed = 1L) {
  ok <- is.finite(x) & is.finite(y) & complete.cases(covariates)
  xr <- rank(x[ok], ties.method = "average")
  yr <- rank(y[ok], ties.method = "average")
  C <- covariates[ok, , drop = FALSE]
  X0 <- cbind(intercept = 1, C)
  X1 <- cbind(intercept = 1, exposure_rank = xr, C)
  xres <- lm.fit(X0, xr)$residuals
  yres <- lm.fit(X0, yr)$residuals
  rho <- unname(cor(xres, yres))
  q <- ncol(C)
  pdf <- length(xr) - q - 2L
  tt <- rho * sqrt(pdf / (1 - rho^2))
  pp <- 2 * pt(-abs(tt), df = pdf)
  robust <- hc3_fit(yr, X1)
  reduced <- lm.fit(X0, yr)
  h0 <- rowSums((X0 %*% chol2inv(chol(crossprod(X0)))) * X0)
  null_resid <- reduced$residuals / pmax(1e-8, 1 - h0)
  set.seed(seed)
  null_t <- vapply(seq_len(B), function(i) {
    ystar <- reduced$fitted.values + null_resid * sample(c(-1, 1), length(yr), TRUE)
    hc3_fit(ystar, X1)$t
  }, numeric(1))
  wild_p <- (1 + sum(abs(null_t) >= abs(robust$t))) / (B + 1)
  loo <- vapply(seq_along(xr), function(drop) {
    keep <- seq_along(xr) != drop
    cor(lm.fit(X0[keep, , drop = FALSE], xr[keep])$residuals,
        lm.fit(X0[keep, , drop = FALSE], yr[keep])$residuals)
  }, numeric(1))
  data.table(representation = representation, n = length(xr),
             adjusted_rank_residual_rho = rho, correct_partial_df = pdf,
             correct_partial_p = pp, HC3_p = robust$p, wild_bootstrap_t_p = wild_p,
             loo_rho_min = min(loo), loo_rho_max = max(loo),
             loo_positive_fraction = mean(loo > 0))
}

merged_assoc <- rbindlist(list(
  analyze_assoc(dt$sender_logratio, dt$receiver_logratio_merged,
                "primary_sender_vs_merged_receiver_logratio", 9999L, 2026072944L),
  analyze_assoc(dt$sender_asin, dt$receiver_asin_merged,
                "sensitivity_sender_vs_merged_receiver_asin", 9999L, 2026072945L)
))
fwrite(dt[, .(sample_id, Patient_ID, fibroblast_cells, sender_cells, sender_logratio,
              myeloid_cells, receiver_cells, receiver_logratio_merged,
              sender_asin, receiver_asin_merged, ECM_direction)],
       file.path(out_dir, "GSE183904_merged_receiver_abundance_patient_data.tsv"), sep = "\t")
fwrite(merged_assoc, file.path(out_dir, "GSE183904_merged_receiver_abundance_inference.tsv"), sep = "\t")

boundary_summary <- data.table(
  metric = c("paired_patients", "cluster3_cells", "cluster5_cells",
             "median_patient_expression_spearman", "minimum_patient_expression_spearman",
             "DE_FDR_0.05_abs_logFC_0.5", "DE_cluster3_up", "DE_cluster5_up"),
  value = c(nrow(boundary$paired), sum(group_meta[group_id == "my_c3", n_cells]),
            sum(group_meta[group_id == "my_c5", n_cells]), median(boundary_cor), min(boundary_cor),
            boundary$table[FDR < 0.05 & abs(logFC_target_minus_comparator) >= 0.5, .N],
            boundary$table[FDR < 0.05 & logFC_target_minus_comparator >= 0.5, .N],
            boundary$table[FDR < 0.05 & logFC_target_minus_comparator <= -0.5, .N])
)
fwrite(boundary_summary, file.path(out_dir, "GSE183904_C1QC_boundary_summary.tsv"), sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE183904",
  inference_unit = "patient pseudobulk or patient-level paired state score",
  boundary_comparison = "cluster 3 C1QC+APOE+ versus cluster 5 HLA-II+C1QC+",
  signature_method = "patient-blocked edgeR quasi-likelihood; outcome-blind TCGA-mappable top genes",
  signatures = split(signatures$gene_symbol, signatures$signature_id),
  merged_receiver_association = as.list(merged_assoc),
  interpretation = paste(
    "Cluster-3/cluster-5 differential expression assesses distinguishability.",
    "Merged-state abundance and program analyses test whether their boundary changes the ecosystem-level result.",
    "TCGA outcomes were not used to select signature genes."
  ),
  input_sha256 = as.list(setNames(vapply(inputs, digest, character(1),
                                         algo = "sha256", file = TRUE), basename(inputs)))
)
write_json(audit, file.path(audit_dir, "GSE183904_C1QC_boundary_and_bulk_signatures.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nBoundary summary:\n")
print(boundary_summary)
cat("\nFrozen signatures:\n")
print(signatures[, .(genes = .N, first_genes = paste(head(gene_symbol, 8), collapse = ";")),
                 by = signature_id])
cat("\nMerged receiver abundance association:\n")
print(merged_assoc)
cat("\nMerged receiver program summary:\n")
print(merged_program_summary[order(paired_wilcoxon_p)])
