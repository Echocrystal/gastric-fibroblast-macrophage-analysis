#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else
  "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(project_dir, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(svglite)
  library(ragg)
})

set.seed(183904L)

processed_dir <- file.path(project_dir, "data", "processed", "scrna", "GSE183904")
result_dir <- file.path(project_dir, "results", "scrna", "GSE183904")
figure_dir <- file.path(project_dir, "figures", "main")
source_dir <- file.path(project_dir, "source_data")
model_dir <- file.path(project_dir, "models", "ecology")
audit_dir <- file.path(project_dir, "audit")
for (d in c(result_dir, figure_dir, source_dir, model_dir, audit_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

pb <- fread(file.path(processed_dir, "GSE183904_target_gene_pseudobulk_counts.tsv.gz"))
composition <- fread(file.path(result_dir, "GSE183904_coarse_cell_composition.tsv"))
manifest <- fread(file.path(processed_dir, "GSE183904_sample_patient_clinical_manifest.tsv"))
frozen <- fread(file.path(project_dir, "models", "ecotype", "tcga_ecm_centroid_gene_contrasts.tsv"))
naba <- fread(file.path(project_dir, "data", "raw", "gene_sets", "naba_ecm_analysis_universe.tsv"))$gene_symbol
program_manifest <- fread(file.path(project_dir, "data", "raw", "gene_sets", "scrna_spatial_programs_locked.tsv"))
programs <- split(program_manifest$gene_symbol, program_manifest$program)

primary_samples <- manifest[include_primary_replication == TRUE]$sample_id
if (length(primary_samples) != 26L || uniqueN(manifest[sample_id %chin% primary_samples]$Patient_ID) != 26L)
  stop("Expected 26 independent primary gastric tumors")

z_columns <- function(mat) {
  z <- matrix(NA_real_, nrow(mat), ncol(mat), dimnames = dimnames(mat))
  for (j in seq_len(ncol(mat))) {
    v <- mat[, j]
    mu <- mean(v, na.rm = TRUE); s <- sd(v, na.rm = TRUE)
    if (is.finite(s) && s > 0) z[, j] <- (v - mu) / s else z[is.finite(v), j] <- 0
  }
  z
}

target_tmm <- function(long_counts, eligible_samples) {
  x <- long_counts[sample_id %chin% eligible_samples, .(count = sum(count)),
                   by = .(sample_id, gene_symbol)]
  wide <- dcast(x, gene_symbol ~ sample_id, value.var = "count", fill = 0)
  genes <- wide$gene_symbol
  mat <- as.matrix(wide[, -"gene_symbol"])
  rownames(mat) <- genes
  storage.mode(mat) <- "integer"
  y <- calcNormFactors(DGEList(mat), method = "TMM")
  lcpm <- cpm(y, log = TRUE, prior.count = 1)
  long <- melt(as.data.table(t(lcpm), keep.rownames = "sample_id"), id.vars = "sample_id",
               variable.name = "gene_symbol", value.name = "logCPM")
  list(
    expression = long,
    normalization = data.table(
      sample_id = colnames(y$counts), target_panel_library_size = y$samples$lib.size,
      TMM_normalization_factor = y$samples$norm.factors,
      effective_library_size = y$samples$lib.size * y$samples$norm.factors
    )
  )
}

score_sets <- function(long_expr, eligible_samples, gene_sets) {
  wide <- dcast(long_expr[sample_id %chin% eligible_samples], sample_id ~ gene_symbol,
                value.var = "logCPM")
  ids <- wide$sample_id
  mat <- as.matrix(wide[, -"sample_id"])
  rownames(mat) <- ids
  z <- z_columns(mat)
  rbindlist(lapply(names(gene_sets), function(nm) {
    present <- intersect(unique(gene_sets[[nm]]), colnames(z))
    if (length(present) < 2L) return(NULL)
    data.table(sample_id = rownames(z), program = nm,
               score = rowMeans(z[, present, drop = FALSE]),
               available_genes = length(present), requested_genes = length(unique(gene_sets[[nm]])))
  }), fill = TRUE)
}

contrast <- frozen[order(-centroid_contrast_C2_minus_C1)]
high50 <- head(contrast$gene_symbol, 50L)
low50 <- tail(contrast$gene_symbol, 50L)
ecm_sets <- list(ECM_C2_HIGH50 = high50, ECM_C1_HIGH50 = low50, PAN_ECM = naba)

score_ecm <- function(analysis_set, coarse_label, min_cells = 0L) {
  aset <- analysis_set
  clabel <- coarse_label
  selected_pb <- pb[analysis_set == aset & coarse_label == clabel]
  cell_counts <- unique(selected_pb[, .(sample_id, n_cells)])
  eligible <- intersect(primary_samples, cell_counts[n_cells >= min_cells]$sample_id)
  norm <- target_tmm(selected_pb, eligible)
  scores_long <- score_sets(norm$expression, eligible, ecm_sets)
  scores <- dcast(scores_long, sample_id ~ program, value.var = "score")
  scores[, ECM_direction := ECM_C2_HIGH50 - ECM_C1_HIGH50]
  list(scores = scores, normalization = norm$normalization, eligibility = cell_counts)
}

all_primary <- score_ecm("secondary_qc_singlets", "All_cells", 0L)
all_submitted <- score_ecm("submitted_cells", "All_cells", 0L)
fib20 <- score_ecm("secondary_qc_singlets", "Fibroblast", 20L)
fib10 <- score_ecm("secondary_qc_singlets", "Fibroblast", 10L)
fwrite(all_primary$scores, file.path(result_dir, "GSE183904_primary_tumor_frozen_ECM_scores.tsv"), sep = "\t")
fwrite(all_submitted$scores, file.path(result_dir, "GSE183904_primary_tumor_frozen_ECM_scores_submitted_cell_sensitivity.tsv"), sep = "\t")
fwrite(fib20$scores, file.path(result_dir, "GSE183904_primary_tumor_fibroblast_ECM_scores_min20.tsv"), sep = "\t")
fwrite(fib10$scores, file.path(result_dir, "GSE183904_primary_tumor_fibroblast_ECM_scores_min10.tsv"), sep = "\t")
fwrite(all_primary$normalization, file.path(result_dir, "GSE183904_primary_tumor_target_panel_TMM_factors.tsv"), sep = "\t")

# Coarse composition features.
comp <- composition[analysis_set == "secondary_qc_singlets" & sample_id %chin% primary_samples]
features <- comp[coarse_label %chin% c("Endothelial", "Fibroblast", "Mural", "Myeloid", "Mast",
                                      "T_cell", "NK", "B_cell", "Plasma"),
                 .(sample_id, feature = paste0("Abundance: ", coarse_label), value = fraction,
                   feature_class = "coarse_composition")]
bplasma <- comp[coarse_label %chin% c("B_cell", "Plasma"),
                .(value = sum(n_cells) / unique(total_cells)), by = sample_id]
bplasma[, `:=`(feature = "Abundance: B plus plasma", feature_class = "coarse_composition")]
features <- rbind(features, bplasma, fill = TRUE)

lineage_programs <- list(
  Fibroblast = c("CAF_ACTIVATED_MATRIX", "CAF_CONTRACTILE", "CAF_INFLAMMATORY_CHEMOKINE"),
  Myeloid = c("MYELOID_SPP1_TREM2", "MYELOID_INFLAMMATORY"),
  Endothelial = "ENDOTHELIAL_ANGIOGENIC",
  T_cell = c("T_CELL_EXHAUSTION", "T_CELL_CYTOTOXIC", "TREG", "TLS_B_SUPPORT"),
  B_cell = "TLS_B_SUPPORT",
  Plasma = "TLS_B_SUPPORT"
)
lineage_norm <- list()
program_coverage <- list()
for (lin in names(lineage_programs)) {
  dt <- pb[analysis_set == "secondary_qc_singlets" & coarse_label == lin]
  cells <- unique(dt[, .(sample_id, n_cells)])
  eligible <- intersect(primary_samples, cells[n_cells >= 20L]$sample_id)
  if (length(eligible) < 6L) next
  norm <- target_tmm(dt, eligible)
  scored <- score_sets(norm$expression, eligible, programs[lineage_programs[[lin]]])
  scored[, `:=`(feature = paste0("Expression in ", lin, ": ", program),
                feature_class = "within_lineage_expression")]
  features <- rbind(features, scored[, .(sample_id, feature, value = score, feature_class)], fill = TRUE)
  lineage_norm[[lin]] <- copy(norm$normalization)[, lineage := lin]
  program_coverage[[lin]] <- unique(scored[, .(lineage = lin, program, available_genes, requested_genes)])
}
fwrite(rbindlist(lineage_norm, fill = TRUE), file.path(result_dir, "GSE183904_lineage_target_panel_TMM_factors.tsv"), sep = "\t")
fwrite(rbindlist(program_coverage, fill = TRUE), file.path(result_dir, "GSE183904_lineage_program_gene_coverage.tsv"), sep = "\t")

clinical_cols <- c("sample_id", "Patient_ID", "Age", "Gender", "Stage*", "Laurens", "MMR", "TCGA")
clinical <- manifest[sample_id %chin% primary_samples, ..clinical_cols]
clinical[, stage_ordinal := fcase(
  grepl("^I$", `Stage*`), 1, grepl("^II$", `Stage*`), 2,
  grepl("^III$", `Stage*`), 3, grepl("^IV$", `Stage*`), 4, default = NA_real_
)]
clinical[, diffuse_component := as.integer(grepl("diffuse", Laurens, ignore.case = TRUE))]
features <- merge(features, clinical, by = "sample_id", all.x = TRUE)
fwrite(features, file.path(result_dir, "GSE183904_primary_tumor_ecology_features.tsv"), sep = "\t", na = "")

bootstrap_spearman <- function(x, y, B = 2000L) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 6L) return(c(NA_real_, NA_real_))
  vals <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    if (sd(x[idx]) == 0 || sd(y[idx]) == 0) NA_real_ else cor(x[idx], y[idx], method = "spearman")
  })
  quantile(vals, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

partial_clinical <- function(x, y, dt) {
  ok <- is.finite(x) & is.finite(y) & is.finite(dt$stage_ordinal) & !is.na(dt$TCGA) &
    !is.na(dt$diffuse_component)
  if (sum(ok) < 12L) return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
  design <- model.matrix(~ stage_ordinal + diffuse_component + factor(TCGA), data = dt[ok])
  xr <- lm.fit(design, rank(x[ok]))$residuals
  yr <- lm.fit(design, rank(y[ok]))$residuals
  test <- suppressWarnings(cor.test(xr, yr, method = "pearson"))
  c(rho = unname(test$estimate), p = test$p.value, n = sum(ok))
}

associate <- function(score_dt, score_col, name) {
  m <- merge(features, score_dt[, .(sample_id, ecm_score = get(score_col))], by = "sample_id")
  out <- m[, {
    ok <- is.finite(value) & is.finite(ecm_score)
    n <- sum(ok)
    if (n < 10L || uniqueN(value[ok]) < 3L) {
      list(n = n, rho = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
           adjusted_n = NA_integer_, adjusted_rho = NA_real_, adjusted_p = NA_real_)
    } else {
      test <- suppressWarnings(cor.test(ecm_score[ok], value[ok], method = "spearman", exact = FALSE))
      ci <- bootstrap_spearman(ecm_score[ok], value[ok])
      adj <- partial_clinical(ecm_score, value, .SD)
      list(n = n, rho = unname(test$estimate), ci_low = ci[1], ci_high = ci[2], p_value = test$p.value,
           adjusted_n = as.integer(adj[["n"]]), adjusted_rho = adj[["rho"]], adjusted_p = adj[["p"]])
    }
  }, by = .(feature, feature_class)]
  out[, `:=`(analysis = name, fdr = p.adjust(p_value, method = "BH"),
             adjusted_fdr = p.adjust(adjusted_p, method = "BH"))]
  out[]
}

assoc_primary <- associate(all_primary$scores, "ECM_direction", "secondary_qc_primary")
assoc_submitted <- associate(all_submitted$scores, "ECM_direction", "submitted_cell_sensitivity")
assoc_fib20 <- associate(fib20$scores, "ECM_direction", "fibroblast_intrinsic_min20")
assoc_fib10 <- associate(fib10$scores, "ECM_direction", "fibroblast_intrinsic_min10_sensitivity")
associations <- rbindlist(list(assoc_primary, assoc_submitted, assoc_fib20, assoc_fib10), fill = TRUE)
setcolorder(associations, c("analysis", "feature", "feature_class", "n", "rho", "ci_low", "ci_high",
                           "p_value", "fdr", "adjusted_n", "adjusted_rho", "adjusted_p", "adjusted_fdr"))
fwrite(associations, file.path(result_dir, "GSE183904_ECM_ecology_donor_associations.tsv"), sep = "\t", na = "")

# Clinical context of the frozen score, without selecting or redefining it.
score_clin <- merge(all_primary$scores, clinical, by = "sample_id")
stage_test <- suppressWarnings(cor.test(score_clin$ECM_direction, score_clin$stage_ordinal,
                                       method = "spearman", exact = FALSE))
lauren_test <- suppressWarnings(wilcox.test(ECM_direction ~ diffuse_component, score_clin, exact = FALSE))
tcga_test <- kruskal.test(ECM_direction ~ factor(TCGA), score_clin)
clinical_tests <- data.table(
  variable = c("Stage ordinal", "Lauren diffuse component", "TCGA subtype"),
  test = c("Spearman", "Wilcoxon rank-sum", "Kruskal-Wallis"),
  effect = c(unname(stage_test$estimate),
             median(score_clin[diffuse_component == 1]$ECM_direction) -
               median(score_clin[diffuse_component == 0]$ECM_direction), NA_real_),
  p_value = c(stage_test$p.value, lauren_test$p.value, tcga_test$p.value)
)
fwrite(clinical_tests, file.path(result_dir, "GSE183904_frozen_ECM_clinical_context.tsv"), sep = "\t", na = "")

# Harmonized two-study direction audit. This is a replication display rather
# than a new pooled training step.
map206 <- data.table(
  harmonized_axis = c("Endothelial abundance", "Endothelial angiogenic program",
                      "T-cell exhaustion program", "T-cell Treg program",
                      "T-cell cytotoxic program", "Myeloid SPP1/TREM2 program",
                      "Myeloid inflammatory program", "Plasma TLS-support program",
                      "B-cell TLS-support program", "B plus plasma abundance"),
  feature206 = c("Coarse abundance: endothelial cells",
                 "Expression in Endothelial: ENDOTHELIAL_ANGIOGENIC",
                 "Expression in CD8+ T: T_CELL_EXHAUSTION",
                 "Expression in CD4+ T: TREG",
                 "Expression in CD8+ T: T_CELL_CYTOTOXIC",
                 "Expression in Myeloid: MYELOID_SPP1_TREM2",
                 "Expression in Myeloid: MYELOID_INFLAMMATORY",
                 "Expression in Plasma: TLS_B_SUPPORT",
                 "Expression in B: TLS_B_SUPPORT",
                 "TLS/B state: B plus plasma abundance"),
  feature183 = c("Abundance: Endothelial",
                 "Expression in Endothelial: ENDOTHELIAL_ANGIOGENIC",
                 "Expression in T_cell: T_CELL_EXHAUSTION",
                 "Expression in T_cell: TREG",
                 "Expression in T_cell: T_CELL_CYTOTOXIC",
                 "Expression in Myeloid: MYELOID_SPP1_TREM2",
                 "Expression in Myeloid: MYELOID_INFLAMMATORY",
                 "Expression in Plasma: TLS_B_SUPPORT",
                 "Expression in B_cell: TLS_B_SUPPORT",
                 "Abundance: B plus plasma")
)
assoc206 <- fread(file.path(project_dir, "results", "scrna", "GSE206785",
                            "GSE206785_ECM_ecology_donor_associations.tsv"))
assoc206 <- assoc206[analysis == "bulk_like_all_cell_ECM"]
replication <- merge(map206, assoc206[, .(feature206 = feature, n206 = n, rho206 = rho,
                                          p206 = p_value, fdr206 = fdr)], by = "feature206", all.x = TRUE)
replication <- merge(replication, assoc_primary[, .(feature183 = feature, n183 = n, rho183 = rho,
                                                    p183 = p_value, fdr183 = fdr,
                                                    adjusted_rho183 = adjusted_rho,
                                                    adjusted_p183 = adjusted_p)],
                     by = "feature183", all.x = TRUE)
replication[, sign_concordant := sign(rho206) == sign(rho183)]
replication[, independent_nominal_replication := sign_concordant & p183 < 0.05]
fwrite(replication, file.path(result_dir, "GSE183904_GSE206785_harmonized_replication.tsv"), sep = "\t", na = "")

save_all <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".svg")), plot, width = width, height = height, device = svglite)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         dpi = 300, device = agg_png, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         dpi = 600, device = agg_tiff, compression = "lzw", bg = "white")
}
theme_manuscript <- theme_classic(base_size = 9) +
  theme(axis.title = element_text(face = "bold"), plot.title = element_text(face = "bold"))

forest <- assoc_primary[is.finite(rho)][order(rho)]
forest[, feature := factor(feature, levels = feature)]
forest[, fdr_status := factor(ifelse(fdr < 0.05, "FDR < 0.05", "FDR ≥ 0.05"),
                              levels = c("FDR < 0.05", "FDR ≥ 0.05"))]
p3e <- ggplot(forest, aes(rho, feature, color = feature_class)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55", linewidth = 0.35) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0, linewidth = 0.45) +
  geom_point(aes(shape = fdr_status), size = 2) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR ≥ 0.05" = 1)) +
  scale_color_manual(
    values = c(coarse_composition = "#577590", within_lineage_expression = "#D1495B"),
    labels = c(coarse_composition = "Cell abundance",
               within_lineage_expression = "Within-lineage expression")
  ) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = "Spearman rho with frozen ECM direction", y = NULL, color = NULL, shape = NULL,
       title = "Independent donor-level replication in 26 primary gastric tumors") +
  theme_manuscript + theme(legend.position = "top", axis.text.y = element_text(size = 7))
save_all(p3e, "Fig3E_GSE183904_ECM_ecology_replication_forest", 9.2,
         max(5.2, 0.27 * nrow(forest) + 1.5))
fwrite(forest, file.path(source_dir, "Fig3E_GSE183904_ECM_ecology_replication_forest.tsv"), sep = "\t")

rep_long <- melt(replication,
                 id.vars = c("harmonized_axis", "sign_concordant", "independent_nominal_replication"),
                 measure.vars = list(rho = c("rho206", "rho183"), n = c("n206", "n183"),
                                     p = c("p206", "p183")),
                 variable.name = "study_index")
rep_long[, study := factor(study_index, levels = 1:2, labels = c("GSE206785", "GSE183904"))]
rep_long[, harmonized_axis := factor(harmonized_axis, levels = rev(map206$harmonized_axis))]
rep_long[, nominal_status := factor(ifelse(p < 0.05, "P < 0.05", "P ≥ 0.05"),
                                    levels = c("P < 0.05", "P ≥ 0.05"))]
p3f <- ggplot(rep_long, aes(rho, harmonized_axis, color = study, group = harmonized_axis)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.35) +
  geom_line(color = "grey75", linewidth = 0.45) +
  geom_point(aes(shape = nominal_status), size = 2.3) +
  scale_color_manual(values = c(GSE206785 = "#4C78A8", GSE183904 = "#E45756")) +
  scale_shape_manual(values = c("P < 0.05" = 16, "P ≥ 0.05" = 1)) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = "Donor-level Spearman rho", y = NULL, color = NULL, shape = NULL,
       title = "Cross-study direction audit of prespecified ECM–immune axes") +
  theme_manuscript + theme(legend.position = "top", axis.text.y = element_text(size = 7.5))
save_all(p3f, "Fig3F_cross_scRNA_replication", 8.6, 5.8)
fwrite(rep_long, file.path(source_dir, "Fig3F_cross_scRNA_replication.tsv"), sep = "\t")

saveRDS(list(
  primary_samples = primary_samples,
  score_definition = "secondary-QC target-panel pseudobulk; TMM logCPM; frozen high50 minus low50 gene-wise z score",
  high50 = high50, low50 = low50, programs = programs,
  all_cell_scores = all_primary$scores,
  fibroblast_scores_min20 = fib20$scores,
  associations = associations,
  harmonized_replication = replication
), file.path(model_dir, "GSE183904_frozen_ECM_donor_replication_object.rds"), compress = "xz")

top <- assoc_primary[is.finite(rho)][order(fdr, -abs(rho))]
audit <- c(
  "# GSE183904 donor-level ECM ecology replication audit",
  "",
  "- Primary cohort: 26 independent primary gastric tumors; peritoneal tumors are excluded from the main replication.",
  "- Primary cells: secondary-QC scDblFinder singlets with conservative coarse annotation; submitted-cell results are a sensitivity analysis.",
  "- Expression: reconstructed counts summed by sample x coarse lineage, TMM/logCPM within the pre-locked target panel.",
  "- Frozen ECM direction: identical TCGA C2-high50 minus C1-high50 definition; no scRNA outcome or state was used to redefine genes.",
  "- Within-lineage programs require at least 20 cells per tumor; all tests use tumors/patients as units.",
  "- Confidence intervals use 2,000 patient-level bootstrap resamples. BH FDR is applied within each analysis.",
  "- Clinical-adjusted sensitivity residualizes ranks for stage, a diffuse-component indicator, and TCGA subtype.",
  "",
  "## Strongest associations",
  "",
  "| Feature | n | rho | 95% bootstrap CI | FDR | clinically adjusted rho |",
  "|---|---:|---:|---:|---:|---:|"
)
for (i in seq_len(min(12L, nrow(top)))) {
  x <- top[i]
  audit <- c(audit, sprintf("| %s | %d | %.3f | %.3f to %.3f | %.3g | %.3f |",
                            x$feature, x$n, x$rho, x$ci_low, x$ci_high, x$fdr, x$adjusted_rho))
}
audit <- c(audit, "", "Cross-study direction concordance is descriptive replication evidence and is not a new trained model.")
writeLines(audit, file.path(audit_dir, "GSE183904_donor_replication.md"), useBytes = TRUE)

cat(paste(audit, collapse = "\n"), "\n")
cat("\nHarmonized replication:\n")
print(replication[, .(harmonized_axis, n206, rho206, p206, n183, rho183, p183,
                      adjusted_rho183, adjusted_p183, sign_concordant,
                      independent_nominal_replication)])
