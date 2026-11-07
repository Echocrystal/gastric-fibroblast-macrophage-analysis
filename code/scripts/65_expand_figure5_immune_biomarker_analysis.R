options(stringsAsFactors = FALSE, warn = 1)
options(error = function() {
  traceback(20L)
  quit(save = "no", status = 1L, runLast = FALSE)
})
set.seed(2026080601L)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
public_project <- "I:/ECM-new/ECM-public"
.libPaths(c(
  file.path(source_project, "runtime", "R_library"),
  file.path(public_project, "environment", "R_library"),
  "I:/ECM-new/scRNA/20.rebuild_GSE163558/R_libs_4.4",
  .libPaths()
))

suppressPackageStartupMessages({
  library(data.table)
  library(MCPcounter)
  library(sandwich)
  library(lmtest)
  library(jsonlite)
  library(digest)
})

score_dir <- file.path(project, "analysis", "results", "PLAT_axis")
bulk_dir <- file.path(source_project, "data", "processed", "bulk")
out_dir <- file.path(score_dir, "external_validation", "immune_biomarker_extension_v1")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expression_files <- c(
  `TCGA-STAD` = file.path(public_project, "data", "processed", "bulk",
                          "tcga_stad_full_symbol_log2cpm_412.rds"),
  GSE62254 = file.path(bulk_dir, "GSE62254_projectable_expression.rds"),
  GSE84437 = file.path(bulk_dir, "GSE84437_projectable_expression.rds"),
  GSE15459 = file.path(bulk_dir, "GSE15459_projectable_expression.rds")
)
score_files <- c(
  `TCGA-STAD` = file.path(score_dir, "TCGA_STAD", "frozen_state_projection",
                          "TCGA_STAD_PLAT_ecosystem_scores_clinical.tsv.gz"),
  GSE62254 = file.path(score_dir, "external_validation", "frozen_state_projection",
                       "GSE62254_frozen_PLAT_ecosystem_scores.tsv.gz"),
  GSE84437 = file.path(score_dir, "external_validation", "frozen_state_projection",
                       "GSE84437_frozen_PLAT_ecosystem_scores.tsv.gz"),
  GSE15459 = file.path(score_dir, "external_validation", "frozen_state_projection",
                       "GSE15459_frozen_PLAT_ecosystem_scores.tsv.gz")
)

signature_file <- file.path(project, "resources", "gene_sets",
                            "GSE183904_PLAT_ecosystem_bulk_signatures.tsv")
checkpoint_file <- file.path(project, "resources", "gene_sets",
                             "immunotherapy_biomarkers_frozen_v1.tsv")
mcp_gene_file <- file.path(source_project, "runtime", "sources", "MCPcounter-master",
                           "Signatures", "genes.txt")
mcp_probe_file <- file.path(source_project, "runtime", "sources", "MCPcounter-master",
                            "Signatures", "probesets.txt")
lineage_file <- file.path(score_dir, "TCGA_STAD", "leave_signature_out_sensitivity",
                          "TCGA_STAD_leave_signature_out_endpoint_scores.tsv.gz")
tide_source_file <- "I:/ECM-new/scRNA/26.matrix_feedback_multicohort/results/source_data/CONSENSUS_ECM_subtype_TIDE_sample_level.tsv"
tide_model_file <- file.path(score_dir, "TCGA_STAD", "continuous_ecosystem_function_immune_v1",
                             "TCGA_ecosystem_TIDE_continuous_HC3_models.tsv")
camera_file <- file.path(score_dir, "TCGA_STAD", "continuous_ecosystem_function_immune_v1",
                         "TCGA_continuous_ecosystem_CAMERA_pathways.tsv.gz")

inputs <- c(expression_files, score_files, signature_file, checkpoint_file,
            mcp_gene_file, mcp_probe_file, lineage_file, tide_source_file,
            tide_model_file, camera_file)
if (!all(file.exists(inputs))) {
  stop("Missing immune-biomarker input(s): ", paste(inputs[!file.exists(inputs)], collapse = "; "))
}

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

read_expression <- function(cohort, path) {
  obj <- readRDS(path)
  expr <- if (cohort == "TCGA-STAD") obj$expression else obj$symbol_expression
  if (is.null(expr) || is.null(rownames(expr)) || is.null(colnames(expr))) {
    stop("Expression matrix or dimnames absent for ", cohort)
  }
  if (anyDuplicated(rownames(expr))) stop("Duplicated gene symbols in ", cohort)
  if (anyDuplicated(colnames(expr))) stop("Duplicated sample identifiers in ", cohort)
  storage.mode(expr) <- "double"
  expr
}

read_scores <- function(cohort, path) {
  d <- fread(path)
  id_col <- if (cohort == "TCGA-STAD") "sample_submitter_id" else "sample_id"
  if (!all(c(id_col, "ecosystem_z", "ecotype4") %in% names(d))) {
    stop("Required score fields absent in ", cohort)
  }
  setnames(d, id_col, "sample_id")
  d
}

hc3_linear <- function(d, formula, term = "ecosystem_z") {
  fit <- lm(formula, data = d)
  vc <- vcovHC(fit, type = "HC3")
  if (!term %in% names(coef(fit))) stop("Term absent from model: ", term)
  beta <- unname(coef(fit)[term])
  se <- unname(sqrt(diag(vc))[term])
  data.table(
    n = nobs(fit), beta = beta, se_HC3 = se,
    ci_low = beta - qnorm(0.975) * se,
    ci_high = beta + qnorm(0.975) * se,
    p_value = 2 * pnorm(abs(beta / se), lower.tail = FALSE)
  )
}

dl_meta <- function(d) {
  d <- d[is.finite(beta) & is.finite(se_HC3) & se_HC3 > 0]
  k <- nrow(d)
  if (k < 2L) return(NULL)
  yi <- d$beta
  vi <- d$se_HC3^2
  wi <- 1 / vi
  fixed <- sum(wi * yi) / sum(wi)
  Q <- sum(wi * (yi - fixed)^2)
  cval <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- max(0, (Q - (k - 1)) / cval)
  wr <- 1 / (vi + tau2)
  beta <- sum(wr * yi) / sum(wr)
  se <- sqrt(1 / sum(wr))
  data.table(
    k = k, cohorts = paste(d$cohort, collapse = ";"), total_n = sum(d$n),
    beta_random_DL = beta, se_random_DL = se,
    ci_low = beta - qnorm(0.975) * se,
    ci_high = beta + qnorm(0.975) * se,
    p_value = 2 * pnorm(abs(beta / se), lower.tail = FALSE),
    tau2_DL = tau2, Q = Q, Q_df = k - 1,
    Q_p = pchisq(Q, df = k - 1, lower.tail = FALSE),
    I2_percent = if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0,
    positive_cohorts = sum(yi > 0), negative_cohorts = sum(yi < 0)
  )
}

sig <- fread(signature_file)
primary_ids <- c("sender_CXCL14_POSTN_state_axis_excluded",
                 "receiver_C1QC_APOE_specific_axis_excluded")
primary_genes <- unique(sig[signature_id %in% primary_ids, gene_symbol])
if (length(primary_genes) != 60L) stop("Expected 60 frozen primary signature genes")

checkpoints <- fread(checkpoint_file)
if (anyDuplicated(checkpoints$gene_symbol)) stop("Checkpoint manifest has duplicate symbols")

mcp_genes <- fread(mcp_gene_file)
mcp_probes <- read.table(mcp_probe_file, sep = "\t", header = FALSE,
                         stringsAsFactors = FALSE, colClasses = "character")
mcp_population_col <- "Cell population"
mcp_symbol_col <- "HUGO symbols"
mcp_overlap <- mcp_genes[get(mcp_symbol_col) %in% primary_genes,
                         .(removed_genes = paste(sort(unique(get(mcp_symbol_col))), collapse = ";"),
                           removed_n = uniqueN(get(mcp_symbol_col))),
                         by = mcp_population_col]
mcp_genes_leaveout <- mcp_genes[!get(mcp_symbol_col) %in% primary_genes]

cohort_objects <- list()
mcp_patient <- list()
mcp_models <- list()
mcp_groups <- list()
checkpoint_patient <- list()
checkpoint_models <- list()
checkpoint_groups <- list()
coverage_rows <- list()

for (cohort in names(expression_files)) {
  message("Process ", cohort)
  expr <- read_expression(cohort, expression_files[[cohort]])
  scores <- read_scores(cohort, score_files[[cohort]])
  idx <- match(colnames(expr), scores$sample_id)
  if (anyNA(idx)) stop("Expression/score sample mismatch in ", cohort)
  scores <- scores[idx]
  if (!identical(scores$sample_id, colnames(expr))) stop("Failed ordered alignment in ", cohort)
  scores[, ecosystem_group := ifelse(ecosystem_z > median(ecosystem_z, na.rm = TRUE),
                                     "High ecosystem", "Low ecosystem")]
  scores[, ecosystem_group := factor(ecosystem_group,
                                     levels = c("Low ecosystem", "High ecosystem"))]

  # MCP-counter is recalculated after removing every defining sender/receiver gene
  # from its marker table. Its scores are abundance indices, not proportions.
  mcp <- MCPcounter.estimate(
    expr, featuresType = "HUGO_symbols",
    probesets = mcp_probes, genes = as.data.frame(mcp_genes_leaveout)
  )
  if (!identical(colnames(mcp), colnames(expr))) stop("MCP-counter output mismatch")
  mcp_dt <- as.data.table(t(mcp), keep.rownames = "sample_id")
  mcp_long <- melt(mcp_dt, id.vars = "sample_id",
                   variable.name = "population", value.name = "score_raw")
  mcp_long <- merge(mcp_long,
                    scores[, .(sample_id, ecosystem_z, ecotype4, ecosystem_group)],
                    by = "sample_id", all.x = TRUE, sort = FALSE)
  mcp_long[, score_z := zscore(score_raw), by = population]
  mcp_long[, cohort := cohort]
  mcp_patient[[cohort]] <- mcp_long

  mcp_models[[cohort]] <- mcp_long[, {
    d <- .SD[complete.cases(.SD[, .(score_z, ecosystem_z)])]
    cbind(hc3_linear(d, score_z ~ ecosystem_z), model_tier = "unadjusted")
  }, by = .(cohort, population)]
  mcp_groups[[cohort]] <- mcp_long[, {
    d <- .SD[complete.cases(.SD[, .(score_z, ecosystem_group)])]
    w <- wilcox.test(score_z ~ ecosystem_group, data = d, exact = FALSE)
    .(n = nrow(d), low_n = sum(d$ecosystem_group == "Low ecosystem"),
      high_n = sum(d$ecosystem_group == "High ecosystem"),
      low_median_z = median(d$score_z[d$ecosystem_group == "Low ecosystem"]),
      high_median_z = median(d$score_z[d$ecosystem_group == "High ecosystem"]),
      median_difference_high_minus_low =
        median(d$score_z[d$ecosystem_group == "High ecosystem"]) -
        median(d$score_z[d$ecosystem_group == "Low ecosystem"]),
      wilcoxon_p = w$p.value)
  }, by = .(cohort, population)]

  checkpoint_present <- checkpoints$gene_symbol %in% rownames(expr)
  coverage_rows[[cohort]] <- data.table(
    cohort = cohort,
    expression_genes = nrow(expr), samples = ncol(expr),
    checkpoint_genes_requested = nrow(checkpoints),
    checkpoint_genes_present = sum(checkpoint_present),
    checkpoint_genes_missing = paste(checkpoints$gene_symbol[!checkpoint_present], collapse = ";")
  )
  cp_expr <- expr[intersect(checkpoints$gene_symbol, rownames(expr)), , drop = FALSE]
  cp_dt <- as.data.table(t(cp_expr), keep.rownames = "sample_id")
  cp_long <- melt(cp_dt, id.vars = "sample_id",
                  variable.name = "gene_symbol", value.name = "expression_raw")
  cp_long <- merge(cp_long, checkpoints, by = "gene_symbol", all.x = TRUE)
  cp_long <- merge(cp_long,
                   scores[, .(sample_id, ecosystem_z, ecotype4, ecosystem_group)],
                   by = "sample_id", all.x = TRUE, sort = FALSE)
  cp_long[, expression_z := zscore(expression_raw), by = gene_symbol]
  cp_long[, cohort := cohort]
  checkpoint_patient[[cohort]] <- cp_long
  checkpoint_models[[cohort]] <- cp_long[, {
    d <- .SD[complete.cases(.SD[, .(expression_z, ecosystem_z)])]
    cbind(hc3_linear(d, expression_z ~ ecosystem_z), model_tier = "unadjusted")
  }, by = .(cohort, gene_symbol, display_label, marker_class)]
  checkpoint_groups[[cohort]] <- cp_long[, {
    d <- .SD[complete.cases(.SD[, .(expression_z, ecosystem_group)])]
    w <- wilcox.test(expression_z ~ ecosystem_group, data = d, exact = FALSE)
    .(n = nrow(d), low_n = sum(d$ecosystem_group == "Low ecosystem"),
      high_n = sum(d$ecosystem_group == "High ecosystem"),
      low_median_z = median(d$expression_z[d$ecosystem_group == "Low ecosystem"]),
      high_median_z = median(d$expression_z[d$ecosystem_group == "High ecosystem"]),
      median_difference_high_minus_low =
        median(d$expression_z[d$ecosystem_group == "High ecosystem"]) -
        median(d$expression_z[d$ecosystem_group == "Low ecosystem"]),
      wilcoxon_p = w$p.value)
  }, by = .(cohort, gene_symbol, display_label, marker_class)]

  cohort_objects[[cohort]] <- list(expression = expr, scores = scores)
}

mcp_patient_dt <- rbindlist(mcp_patient, use.names = TRUE, fill = TRUE)
mcp_model_dt <- rbindlist(mcp_models, use.names = TRUE, fill = TRUE)
mcp_group_dt <- rbindlist(mcp_groups, use.names = TRUE, fill = TRUE)
mcp_model_dt[, FDR_BH_within_cohort := p.adjust(p_value, "BH"), by = cohort]
mcp_group_dt[, FDR_BH_within_cohort := p.adjust(wilcoxon_p, "BH"), by = cohort]

mcp_meta_dt <- mcp_model_dt[, dl_meta(.SD), by = population]
mcp_meta_dt[, FDR_BH := p.adjust(p_value, "BH")]

checkpoint_patient_dt <- rbindlist(checkpoint_patient, use.names = TRUE, fill = TRUE)
checkpoint_model_dt <- rbindlist(checkpoint_models, use.names = TRUE, fill = TRUE)
checkpoint_group_dt <- rbindlist(checkpoint_groups, use.names = TRUE, fill = TRUE)
checkpoint_model_dt[, FDR_BH_within_cohort := p.adjust(p_value, "BH"), by = cohort]
checkpoint_group_dt[, FDR_BH_within_cohort := p.adjust(wilcoxon_p, "BH"), by = cohort]
checkpoint_meta_dt <- checkpoint_model_dt[, dl_meta(.SD),
                                          by = .(gene_symbol, display_label, marker_class)]
checkpoint_meta_dt[, FDR_BH := p.adjust(p_value, "BH")]

# TCGA-specific calibration of checkpoint associations. The leave-signature-out
# broad lineage scores are used so the adjustment does not reuse any of the 60
# defining genes.
tcga_cp <- copy(checkpoint_patient_dt[cohort == "TCGA-STAD"])
tcga_scores <- copy(cohort_objects[["TCGA-STAD"]]$scores)
lineage <- fread(lineage_file)
setnames(lineage, "sample_submitter_id", "sample_id")
needed_lineage <- c("sample_id", "Fibroblasts", "Monocytic lineage")
if (!all(needed_lineage %in% names(lineage))) stop("Leave-out lineage fields absent")
tcga_cov <- merge(
  tcga_scores[, .(sample_id, absolute_purity, tcga2014_subtype)],
  lineage[, ..needed_lineage], by = "sample_id", all.x = TRUE, sort = FALSE
)
setnames(tcga_cov, c("Fibroblasts", "Monocytic lineage"),
         c("fibroblasts_LO", "monocytic_LO"))
tcga_cov[, `:=`(
  purity_z = zscore(absolute_purity),
  fibroblasts_z = zscore(fibroblasts_LO),
  monocytic_z = zscore(monocytic_LO),
  tcga2014_subtype = factor(tcga2014_subtype, levels = c("CIN", "GS", "MSI", "EBV"))
)]
tcga_cp <- merge(tcga_cp, tcga_cov, by = "sample_id", all.x = TRUE, sort = FALSE)

tcga_checkpoint_adjusted <- tcga_cp[, {
  d1 <- .SD[complete.cases(.SD[, .(expression_z, ecosystem_z, purity_z)])]
  d2 <- .SD[complete.cases(.SD[, .(expression_z, ecosystem_z, purity_z,
                                   fibroblasts_z, monocytic_z, tcga2014_subtype)])]
  rbind(
    cbind(hc3_linear(d1, expression_z ~ ecosystem_z + purity_z),
          model_tier = "purity_adjusted"),
    cbind(hc3_linear(d2, expression_z ~ ecosystem_z + purity_z +
                       fibroblasts_z + monocytic_z + tcga2014_subtype),
          model_tier = "lineage_subtype_adjusted")
  )
}, by = .(gene_symbol, display_label, marker_class)]
tcga_checkpoint_adjusted[, FDR_BH := p.adjust(p_value, "BH"), by = model_tier]

# TIDE-predicted responder is retained as a model-derived phenotype and is kept
# separate from the observed Kim2018 anti-PD-1 response analysis.
tide <- fread(tide_source_file)
setnames(tide, "sample", "case_submitter_id")
tide_map <- merge(
  tcga_scores[, .(sample_id, case_submitter_id, ecosystem_z, ecotype4)],
  tide[, .(case_submitter_id, tide_predicted_responder = as.logical(Responder),
           TIDE, Dysfunction, Exclusion, MDSC, CAF, `TAM M2`, IFNG, CD274, CD8, CTL)],
  by = "case_submitter_id", all = FALSE, sort = FALSE
)
tide_map <- merge(tide_map, tcga_cov, by = "sample_id", all.x = TRUE, sort = FALSE)
if (anyDuplicated(tide_map$sample_id)) stop("Duplicated TIDE mapping")

logistic_hc3 <- function(d, formula, term = "ecosystem_z") {
  fit <- glm(formula, data = d, family = binomial())
  vc <- vcovHC(fit, type = "HC3")
  beta <- unname(coef(fit)[term])
  se <- unname(sqrt(diag(vc))[term])
  data.table(
    n = nobs(fit), events = sum(model.response(model.frame(fit)) == 1),
    OR_per_ecosystem_SD = exp(beta),
    ci_low = exp(beta - qnorm(0.975) * se),
    ci_high = exp(beta + qnorm(0.975) * se),
    p_value = 2 * pnorm(abs(beta / se), lower.tail = FALSE),
    convergence = fit$converged
  )
}

tide_d0 <- tide_map[complete.cases(tide_map[, .(tide_predicted_responder, ecosystem_z)])]
tide_d1 <- tide_map[complete.cases(tide_map[, .(tide_predicted_responder, ecosystem_z, purity_z)])]
tide_d2 <- tide_map[complete.cases(tide_map[, .(tide_predicted_responder, ecosystem_z,
                                               purity_z, fibroblasts_z, monocytic_z,
                                               tcga2014_subtype)])]
tide_response_models <- rbindlist(list(
  cbind(logistic_hc3(tide_d0, tide_predicted_responder ~ ecosystem_z),
        model_tier = "unadjusted"),
  cbind(logistic_hc3(tide_d1, tide_predicted_responder ~ ecosystem_z + purity_z),
        model_tier = "purity_adjusted"),
  cbind(logistic_hc3(tide_d2, tide_predicted_responder ~ ecosystem_z + purity_z +
                       fibroblasts_z + monocytic_z + tcga2014_subtype),
        model_tier = "lineage_subtype_adjusted")
), use.names = TRUE)

tide_config <- tide_map[, .(
  n = .N,
  predicted_responders = sum(tide_predicted_responder, na.rm = TRUE),
  predicted_nonresponders = sum(!tide_predicted_responder, na.rm = TRUE),
  predicted_response_fraction = mean(tide_predicted_responder, na.rm = TRUE)
), by = ecotype4]
tide_tab <- table(tide_map$ecotype4, tide_map$tide_predicted_responder)
tide_fisher <- fisher.test(tide_tab, simulate.p.value = TRUE, B = 99999)
tide_config[, `:=`(
  omnibus_test = "Fisher exact test with 99,999 Monte Carlo replicates",
  omnibus_p = tide_fisher$p.value
)]
tide_numeric_models <- fread(tide_model_file)

# Existing competitive pathway analysis is curated into an immune/ICI-facing
# table. No new high-versus-low pathway search is introduced.
camera <- fread(camera_file)
immune_regex <- paste(c(
  "IMMUNE", "INTERFERON", "T_CELL", "B_CELL", "LEUKOCYTE", "LYMPHOCYTE",
  "MYELOID", "MACROPHAGE", "PHAGOCYT", "CHEMOKINE", "CYTOKINE",
  "ANTIGEN", "PD_1", "PD1", "CHECKPOINT", "COMPLEMENT", "COAGULATION",
  "EXTRACELLULAR_MATRIX", "ECM", "CELL_MATRIX", "INTEGRIN", "FOCAL_ADHESION"
), collapse = "|")
camera_immune <- camera[
  leave_primary_signature_out == TRUE & grepl(immune_regex, pathway, ignore.case = TRUE)
]
setorder(camera_immune, model_tier, FDR_BH, p_value)

fwrite(rbindlist(coverage_rows), file.path(out_dir, "immune_biomarker_gene_coverage.tsv"), sep = "\t")
fwrite(mcp_overlap, file.path(out_dir, "MCPcounter_primary_signature_overlap_removed.tsv"), sep = "\t")
fwrite(mcp_patient_dt, file.path(out_dir, "four_cohort_MCPcounter_leaveout_sample_scores.tsv.gz"), sep = "\t")
fwrite(mcp_model_dt, file.path(out_dir, "four_cohort_MCPcounter_continuous_HC3_models.tsv"), sep = "\t")
fwrite(mcp_meta_dt, file.path(out_dir, "four_cohort_MCPcounter_random_effects_meta.tsv"), sep = "\t")
fwrite(mcp_group_dt, file.path(out_dir, "four_cohort_MCPcounter_high_low_descriptive.tsv"), sep = "\t")
fwrite(checkpoint_patient_dt, file.path(out_dir, "four_cohort_checkpoint_sample_expression.tsv.gz"), sep = "\t")
fwrite(checkpoint_model_dt, file.path(out_dir, "four_cohort_checkpoint_continuous_HC3_models.tsv"), sep = "\t")
fwrite(checkpoint_meta_dt, file.path(out_dir, "four_cohort_checkpoint_random_effects_meta.tsv"), sep = "\t")
fwrite(checkpoint_group_dt, file.path(out_dir, "four_cohort_checkpoint_high_low_descriptive.tsv"), sep = "\t")
fwrite(tcga_checkpoint_adjusted, file.path(out_dir, "TCGA_checkpoint_adjusted_HC3_models.tsv"), sep = "\t")
fwrite(tide_map, file.path(out_dir, "TCGA_TIDE_predicted_response_sample_data.tsv.gz"), sep = "\t")
fwrite(tide_response_models, file.path(out_dir, "TCGA_TIDE_predicted_response_logistic_models.tsv"), sep = "\t")
fwrite(tide_config, file.path(out_dir, "TCGA_TIDE_predicted_response_by_configuration.tsv"), sep = "\t")
fwrite(tide_numeric_models, file.path(out_dir, "TCGA_TIDE_numeric_endpoint_models_existing.tsv"), sep = "\t")
fwrite(camera_immune, file.path(out_dir, "TCGA_CAMERA_immune_ICI_relevant_pathways.tsv.gz"), sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  objective = "Figure 5 extension: immune infiltration, ICI biomarkers and TIDE-predicted response",
  primary_exposure = "frozen continuous ecosystem_z per within-cohort standard deviation",
  groups = "cohort-median high/low ecosystem groups are descriptive only",
  MCPcounter = list(
    interpretation = "abundance indices, not cell fractions",
    signature_decontamination = "all 60 defining sender/receiver genes removed from MCP-counter gene markers before scoring",
    primary_analysis = "cohort-specific HC3 linear model followed by DerSimonian-Laird random-effects meta-analysis"
  ),
  checkpoint_panel = as.data.frame(checkpoints),
  TIDE = list(
    mapped_n = nrow(tide_map),
    predicted_responders = sum(tide_map$tide_predicted_responder, na.rm = TRUE),
    boundary = "TIDE responder is model-derived and is not an observed treatment outcome; Kim2018 observed response remains separate"
  ),
  pathway = list(
    method = "existing leave-60-genes-out CAMERA competitive tests",
    boundary = "pathway association does not establish mechanism"
  ),
  input_sha256 = as.list(setNames(vapply(inputs, digest, character(1),
                                         algo = "sha256", file = TRUE), basename(inputs)))
)
write_json(audit, file.path(audit_dir, "Figure5_immune_biomarker_extension_v1.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nTop cross-cohort MCP-counter associations:\n")
print(mcp_meta_dt[order(FDR_BH)])
cat("\nTop cross-cohort checkpoint associations:\n")
print(checkpoint_meta_dt[order(FDR_BH)])
cat("\nTCGA fully adjusted checkpoint associations:\n")
print(tcga_checkpoint_adjusted[model_tier == "lineage_subtype_adjusted"][order(FDR_BH)])
cat("\nTIDE predicted response models:\n")
print(tide_response_models)
cat("\nTIDE predicted response by configuration:\n")
print(tide_config)
