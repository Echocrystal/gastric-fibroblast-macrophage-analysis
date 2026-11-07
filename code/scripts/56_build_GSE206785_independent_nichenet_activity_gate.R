#!/usr/bin/env Rscript

# Build an outcome-blind, receiver-nonprespecified NicheNet activity gate for
# Figure 1.  Receiver responses are donor-blocked tumour-versus-normal
# pseudobulk contrasts within broad immune lineages.  No fibroblast subtype,
# myeloid subtype, PLAT status, survival outcome, treatment outcome or spatial
# result is used to define the response genes or the activity threshold.

options(stringsAsFactors = FALSE, warn = 1)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
public <- "I:/ECM-new/ECM-public"
figure_lib <- file.path(project, "runtime", "R_figure")
vendor_lib <- "I:/ECM-new/scRNA/26.matrix_feedback_multicohort/R_libs"
seurat_lib <- "I:/ECM-new/scRNA/20.rebuild_GSE163558/R_libs_4.4"
.libPaths(c(figure_lib, seurat_lib, vendor_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(nichenetr)
})

processed <- file.path(
  project, "analysis", "processed",
  "GSE206785_tumor_normal_fine_state_lr_expression.tsv.gz"
)
metadata_path <- file.path(
  public, "data", "raw", "scrna", "GSE206785",
  "GSE206785_metadata.txt.gz"
)
old_membership_path <- file.path(
  project, "figures", "data", "Figure1_GSE206785_unbiased",
  "panel_b_four_set_membership.tsv.gz"
)
ligand_target_path <- file.path(
  public, "data", "public_gene_sets", "nichenet_v2",
  "ligand_target_matrix_nsga2r_final.rds"
)
lr_network_path <- file.path(
  project, "resources", "ligand_receptor",
  "lr_network_human_21122021.rds"
)
out_dir <- file.path(project, "analysis", "results", "figure1_gse206785_unbiased")
fig_dir <- file.path(project, "figures", "data", "Figure1_GSE206785_unbiased")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  file.exists(processed), file.exists(metadata_path),
  file.exists(old_membership_path), file.exists(ligand_target_path),
  file.exists(lr_network_path)
)

# Frozen, PLAT-independent eligibility and response rules.
min_cells_each_tissue <- 20L
min_paired_donors <- 5L
min_cpm <- 1
min_expr_samples <- 3L
response_fdr <- 0.05
response_logfc <- 0.25
response_fraction <- 0.10
min_response_genes <- 50L
max_response_genes <- 250L
min_positive_lineages <- 2L
activity_percentile_threshold <- 0.75

lineage_map <- c(
  "B" = "B cells",
  "CD4 T" = "CD4+ T",
  "CD8 T" = "CD8+ T",
  "Innate lymphoid" = "NK/ILC",
  "Mast" = "Mast",
  "Myeloid" = "Myeloid",
  "Plasma" = "Plasma"
)

message("Reading metadata and auditing tumour/normal coverage")
meta <- fread(metadata_path)
meta <- meta[
  Type %chin% names(lineage_map) &
    !grepl("lowQual|UnIdent", Annotation, ignore.case = TRUE)
]
meta[, lineage := unname(lineage_map[Type])]
cell_counts <- meta[, .N, by = .(lineage, Patient, Tissue)]
coverage_wide <- dcast(
  cell_counts, lineage + Patient ~ Tissue, value.var = "N", fill = 0
)
if (!"Normal" %in% names(coverage_wide)) coverage_wide[, Normal := 0L]
if (!"Tumor" %in% names(coverage_wide)) coverage_wide[, Tumor := 0L]
coverage <- coverage_wide[, .(
  tumour_cells = sum(Tumor),
  normal_cells = sum(Normal),
  tumour_patients = sum(Tumor > 0),
  normal_patients = sum(Normal > 0),
  paired_patients = sum(Tumor > 0 & Normal > 0),
  paired_ge10 = sum(Tumor >= 10 & Normal >= 10),
  paired_ge20 = sum(Tumor >= 20 & Normal >= 20),
  minimum_cells_in_eligible_pairs = if (any(Tumor >= 20 & Normal >= 20)) {
    min(c(Tumor[Tumor >= 20 & Normal >= 20], Normal[Tumor >= 20 & Normal >= 20]))
  } else 0L
), by = lineage]
setorder(coverage, lineage)
fwrite(coverage, file.path(out_dir, "GSE206785_immune_tumour_normal_coverage.tsv"), sep = "\t")
fwrite(coverage, file.path(fig_dir, "panel_a_receiver_coverage.tsv"), sep = "\t")

eligible_pairs <- coverage_wide[Tumor >= min_cells_each_tissue & Normal >= min_cells_each_tissue]
eligible_lineages <- coverage[paired_ge20 >= min_paired_donors, lineage]
if (!length(eligible_lineages)) stop("No immune lineage meets the paired-donor coverage gate")

message("Reading patient-by-state expression summaries")
x <- fread(processed)
x <- x[
  Type %chin% names(lineage_map) &
    !grepl("lowQual|UnIdent", Annotation, ignore.case = TRUE)
]
x[, lineage := unname(lineage_map[Type])]
x <- merge(
  x,
  eligible_pairs[, .(lineage, Patient)],
  by = c("lineage", "Patient"),
  all = FALSE
)

# Sum reconstructed counts and cells across fine annotations to broad receiver
# lineages. n_cells is constant across genes within an annotation, therefore
# max(n_cells) is used before summing annotation cell counts.
annotation_cells <- unique(x[, .(lineage, Patient, Tissue, Annotation, n_cells)])
sample_cells <- annotation_cells[, .(n_cells = sum(n_cells)), by = .(lineage, Patient, Tissue)]
counts <- x[, .(count = sum(reconstructed_count_sum)),
            by = .(lineage, Patient, Tissue, gene_symbol)]
counts <- merge(counts, sample_cells, by = c("lineage", "Patient", "Tissue"), all.x = TRUE)

old_membership <- fread(old_membership_path)
matrisome_set <- old_membership[Matrisome == TRUE, unique(gene_symbol)]
lr_ligand_set <- old_membership[`LR prior` == TRUE, unique(gene_symbol)]
fibroblast_expression_set <- old_membership[
  `Tumour fibroblast expression` == TRUE, unique(gene_symbol)
]

ligand_target_matrix <- readRDS(ligand_target_path)
lr_network <- readRDS(lr_network_path)
lr_network <- as.data.table(lr_network)
if (!all(c("from", "to") %in% names(lr_network))) {
  stop("NicheNet LR network lacks from/to columns")
}

sender_ligands <- intersect(lr_ligand_set, fibroblast_expression_set)
sender_ligands <- intersect(sender_ligands, colnames(ligand_target_matrix))
if (!length(sender_ligands)) stop("No sender ligands overlap the NicheNet matrix")

activity_list <- list()
de_list <- list()
response_summary <- list()

for (lin in eligible_lineages) {
  message("Receiver lineage: ", lin)
  pair_ids <- eligible_pairs[lineage == lin, Patient]
  z <- counts[lineage == lin & Patient %chin% pair_ids]
  z[, sample_id := paste(Patient, Tissue, sep = "__")]
  sample_meta <- unique(z[, .(sample_id, Patient, Tissue, n_cells)])
  sample_meta[, Tissue := factor(Tissue, levels = c("Normal", "Tumor"))]
  sample_meta[, Patient := factor(Patient)]
  setorder(sample_meta, Patient, Tissue)

  count_wide <- dcast(z, gene_symbol ~ sample_id, value.var = "count", fill = 0)
  genes <- count_wide$gene_symbol
  ymat <- as.matrix(count_wide[, -"gene_symbol"])
  rownames(ymat) <- genes
  ymat <- ymat[, sample_meta$sample_id, drop = FALSE]

  y <- DGEList(ymat)
  y <- calcNormFactors(y)
  keep <- rowSums(cpm(y) >= min_cpm) >= min_expr_samples
  y <- y[keep, , keep.lib.sizes = FALSE]
  design <- model.matrix(~ Patient + Tissue, data = sample_meta)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  tissue_coef <- grep("^TissueTumor$", colnames(design))
  if (length(tissue_coef) != 1L) stop("Tumour coefficient not found for ", lin)
  tab <- as.data.table(topTags(glmQLFTest(fit, coef = tissue_coef), n = Inf, sort.by = "none")$table,
                       keep.rownames = "gene_symbol")
  tab[, lineage := lin]
  tab[, strict_fdr_response := FDR <= response_fdr & logFC >= response_logfc]
  tab[, signed_ql_strength := sign(logFC) * sqrt(pmax(F, 0))]
  positive_ranked <- tab[logFC > 0][order(-signed_ql_strength, gene_symbol), gene_symbol]
  n_ranked_response <- min(
    max_response_genes,
    max(min_response_genes, as.integer(ceiling(response_fraction * nrow(tab))))
  )
  n_ranked_response <- min(n_ranked_response, length(positive_ranked))
  geneset <- head(positive_ranked, n_ranked_response)
  tab[, activity_response_gene := gene_symbol %chin% geneset]
  de_list[[lin]] <- tab

  background <- tab$gene_symbol
  expressed_receptors <- intersect(unique(lr_network$to), background)
  potential_ligands <- lr_network[
    from %chin% sender_ligands & to %chin% expressed_receptors,
    unique(from)
  ]
  potential_ligands <- intersect(potential_ligands, colnames(ligand_target_matrix))

  response_summary[[lin]] <- data.table(
    lineage = lin,
    paired_donors = uniqueN(sample_meta$Patient),
    strict_fdr_response_genes = sum(tab$strict_fdr_response),
    activity_response_genes = length(geneset),
    background_genes = length(background),
    expressed_receptors = length(expressed_receptors),
    potential_ligands = length(potential_ligands),
    eligible_for_activity = length(geneset) >= min_response_genes
  )

  if (length(geneset) < min_response_genes) next
  act <- as.data.table(predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  ))
  setnames(act, "test_ligand", "ligand")
  act[, lineage := lin]
  activity_list[[lin]] <- act
}

de_all <- rbindlist(de_list, use.names = TRUE, fill = TRUE)
response_all <- rbindlist(response_summary, use.names = TRUE, fill = TRUE)
activity_all <- rbindlist(activity_list, use.names = TRUE, fill = TRUE)
fwrite(de_all, file.path(out_dir, "GSE206785_immune_tumour_vs_normal_pseudobulk_DE.tsv.gz"), sep = "\t")
fwrite(response_all, file.path(out_dir, "GSE206785_immune_receiver_response_summary.tsv"), sep = "\t")
fwrite(activity_all, file.path(out_dir, "GSE206785_immune_NicheNet_ligand_activity.tsv.gz"), sep = "\t")
fwrite(response_all, file.path(fig_dir, "panel_a_receiver_response_summary.tsv"), sep = "\t")
fwrite(activity_all, file.path(fig_dir, "panel_a_receiver_nichenet_activity.tsv.gz"), sep = "\t")

eligible_activity_lineages <- response_all[eligible_for_activity == TRUE, lineage]
if (length(eligible_activity_lineages) < min_positive_lineages) {
  stop("Fewer than two immune lineages have an eligible receiver response")
}

activity_all[, activity_percentile := frank(aupr_corrected, ties.method = "average") / .N,
             by = lineage]
activity_all[, positive_activity := activity_percentile >= activity_percentile_threshold]
activity_summary <- activity_all[, .(
  eligible_lineages_tested = uniqueN(lineage),
  positive_lineages = sum(positive_activity),
  median_aupr_corrected = median(aupr_corrected, na.rm = TRUE),
  median_activity_percentile = median(activity_percentile, na.rm = TRUE),
  maximum_aupr_corrected = max(aupr_corrected, na.rm = TRUE),
  mean_aupr_corrected = mean(aupr_corrected, na.rm = TRUE)
), by = ligand]
activity_summary[, activity_gate := positive_lineages >= min_positive_lineages]
setorder(activity_summary, -positive_lineages, -median_aupr_corrected, -maximum_aupr_corrected, ligand)
activity_summary[, activity_rank := .I]
fwrite(activity_summary, file.path(fig_dir, "panel_a_independent_nichenet_activity_gate.tsv"), sep = "\t")

activity_set <- activity_summary[activity_gate == TRUE, ligand]
all_genes <- sort(unique(c(matrisome_set, lr_ligand_set, fibroblast_expression_set, activity_set)))
membership <- data.table(gene_symbol = all_genes)
membership[, Matrisome := gene_symbol %chin% matrisome_set]
membership[, `LR prior` := gene_symbol %chin% lr_ligand_set]
membership[, `Tumour fibroblast expression` := gene_symbol %chin% fibroblast_expression_set]
membership[, `NicheNet activity` := gene_symbol %chin% activity_set]
set_cols <- c("Matrisome", "LR prior", "Tumour fibroblast expression", "NicheNet activity")
membership[, membership_code := do.call(paste0, lapply(.SD, as.integer)), .SDcols = set_cols]
membership[, exclusive_region := apply(.SD, 1, function(v) paste(set_cols[as.logical(v)], collapse = " & ")),
           .SDcols = set_cols]
fwrite(membership, file.path(fig_dir, "panel_a_four_set_activity_membership.tsv.gz"), sep = "\t")

# intToBits is little-endian; use an explicit binary formatter for readability.
all_codes <- vapply(
  1:15,
  function(i) paste0(as.integer(intToBits(i))[4:1], collapse = ""),
  character(1)
)
region_counts <- membership[, .N, by = membership_code]
regions <- data.table(membership_code = all_codes)
regions <- merge(regions, region_counts, by = "membership_code", all.x = TRUE)
regions[is.na(N), N := 0L]
setnames(regions, "N", "exclusive_gene_count")
regions[, exclusive_region := vapply(membership_code, function(code) {
  paste(set_cols[strsplit(code, "", fixed = TRUE)[[1]] == "1"], collapse = " & ")
}, character(1))]
setorder(regions, membership_code)
fwrite(regions, file.path(fig_dir, "panel_a_four_set_activity_exclusive_regions.tsv"), sep = "\t")

final_candidates <- membership[
  Matrisome == TRUE & `LR prior` == TRUE &
    `Tumour fibroblast expression` == TRUE & `NicheNet activity` == TRUE,
  gene_symbol
]

# Update candidate ranking using independent activity evidence plus the frozen
# sender-side expression/specificity metrics.  All metrics are percentile
# transformed within the final candidate universe before equal-weight averaging.
old_ranking <- fread(file.path(fig_dir, "panel_c_receiver_unbiased_ligand_ranking.tsv"))
ranking <- merge(
  old_ranking[ligand %chin% final_candidates],
  activity_summary,
  by = "ligand",
  all.x = TRUE
)
rank_metrics <- c(
  "patient_equal_mean_log_expression",
  "patient_equal_detection",
  "patients_detected",
  "fibroblast_lineage_specificity_log2_ratio",
  "positive_lineages",
  "median_aupr_corrected"
)
for (metric in rank_metrics) {
  ranking[, paste0("percentile_", metric) := frank(get(metric), ties.method = "average") / .N]
}
ranking[, evidence_score_activity := rowMeans(.SD),
        .SDcols = paste0("percentile_", rank_metrics)]
setorder(ranking, -evidence_score_activity, -positive_lineages, -median_aupr_corrected, ligand)
ranking[, scan_rank_activity := .I]
fwrite(ranking, file.path(fig_dir, "panel_b_activity_gated_ligand_ranking.tsv"), sep = "\t")

manifest <- data.table(
  field = c(
    "dataset", "receiver_contrast", "eligible_lineages", "minimum_cells_each_tissue",
    "minimum_paired_donors", "response_rule", "minimum_response_genes",
    "activity_metric", "activity_gate", "forbidden_selection_inputs",
    "final_candidate_count", "PLAT_retained", "PLAT_rank"
  ),
  value = c(
    "GSE206785", "paired tumour versus normal within immune lineage",
    paste(eligible_activity_lineages, collapse = "; "), min_cells_each_tissue,
    min_paired_donors,
    sprintf("top positive donor-blocked QL-ranked %.0f%%, bounded at %d-%d genes; strict FDR counts reported separately",
            100 * response_fraction, min_response_genes, max_response_genes),
    min_response_genes, "NicheNet aupr_corrected",
    sprintf("within-lineage aupr_corrected percentile >= %.2f in >= %d eligible immune lineages",
            activity_percentile_threshold, min_positive_lineages),
    "PLAT; fibroblast subtype; myeloid subtype; survival; treatment response; spatial result",
    length(final_candidates), "PLAT" %chin% final_candidates,
    if ("PLAT" %chin% ranking$ligand) ranking[ligand == "PLAT", scan_rank_activity] else NA_integer_
  )
)
fwrite(manifest, file.path(audit_dir, "Figure1_GSE206785_independent_NicheNet_gate_manifest.tsv"), sep = "\t")

print(coverage)
print(response_all)
print(activity_summary[1:min(20L, .N)])
cat("Final ECM candidates:", length(final_candidates), "\n")
cat("PLAT retained:", "PLAT" %chin% final_candidates, "\n")
if ("PLAT" %chin% ranking$ligand) {
  cat("PLAT activity-gated rank:", ranking[ligand == "PLAT", scan_rank_activity],
      "/", nrow(ranking), "\n")
}
