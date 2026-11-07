options(stringsAsFactors = FALSE, warn = 1)
options(error = function() {
  traceback(20L)
  quit(save = "no", status = 1L, runLast = FALSE)
})
set.seed(2026072946L)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
public_project <- "I:/ECM-new/ECM-public"
.libPaths(c(file.path(source_project, "runtime", "R_library"),
            file.path(public_project, "environment", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(brglm2)
  library(pROC)
  library(jsonlite)
  library(digest)
  library(ggplot2)
  library(patchwork)
})

signature_file <- file.path(project, "resources", "gene_sets",
                            "GSE183904_PLAT_ecosystem_bulk_signatures.tsv")
tcga_cox_file <- file.path(project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
                           "frozen_state_projection", "TCGA_STAD_continuous_Cox_models.tsv")
bulk_dir <- file.path(source_project, "data", "processed", "bulk")
cohort_files <- c(
  GSE62254 = file.path(bulk_dir, "GSE62254_projectable_expression.rds"),
  GSE84437 = file.path(bulk_dir, "GSE84437_projectable_expression.rds"),
  GSE15459 = file.path(bulk_dir, "GSE15459_projectable_expression.rds")
)
gse62254_clinical_file <- file.path(bulk_dir, "GSE62254_ACRG_clinical_standardized.tsv")
kim_tpm_file <- file.path(public_project, "data", "processed", "ici", "Kim2018",
                          "Kim2018_gene_symbol_TPM.tsv.gz")
kim_clinical_file <- file.path(public_project, "manifests",
                               "Kim2018_Nature_Table2_clinical_45.tsv")
out_dir <- file.path(project, "analysis", "results", "PLAT_axis",
                     "external_validation", "frozen_state_projection")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(signature_file, tcga_cox_file, cohort_files,
            gse62254_clinical_file, kim_tpm_file, kim_clinical_file)
if (!all(file.exists(inputs))) stop("Missing external-validation input")

signatures_dt <- fread(signature_file)
sender_id <- "sender_CXCL14_POSTN_state_axis_excluded"
sender_inclusive_id <- "sender_CXCL14_POSTN_PLAT_state_inclusive"
receiver_specific_id <- "receiver_C1QC_APOE_specific_axis_excluded"
receiver_merged_id <- "receiver_merged_C1QC_axis_excluded"
required_ids <- c(sender_id, sender_inclusive_id, receiver_specific_id, receiver_merged_id)
sig_list <- split(signatures_dt$gene_symbol, signatures_dt$signature_id)[required_ids]

zscore <- function(x) as.numeric((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))

score_expression <- function(expression, cohort, sample_ids = colnames(expression)) {
  if (is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("Expression dimnames absent: ", cohort)
  }
  coverage <- rbindlist(lapply(names(sig_list), function(id) {
    genes <- unique(sig_list[[id]])
    data.table(cohort = cohort, signature_id = id, frozen_genes = length(genes),
               mapped_genes = sum(genes %in% rownames(expression)),
               coverage_fraction = mean(genes %in% rownames(expression)),
               missing_genes = paste(setdiff(genes, rownames(expression)), collapse = ";"))
  }))
  if (coverage[, any(mapped_genes < 15L)]) stop("Insufficient signature coverage: ", cohort)
  indices <- lapply(sig_list, function(g) match(intersect(g, rownames(expression)),
                                                    rownames(expression)))
  raw <- vapply(seq_len(ncol(expression)), function(j) {
    nr <- rank(expression[, j], ties.method = "average") / nrow(expression)
    vapply(indices, function(idx) mean(nr[idx]), numeric(1))
  }, numeric(length(indices)))
  rownames(raw) <- names(indices)
  colnames(raw) <- colnames(expression)
  z <- t(apply(raw, 1L, zscore))
  scores <- data.table(sample_id = sample_ids)
  scores[, `:=`(
    sender_primary_z = z[sender_id, ],
    sender_PLAT_inclusive_z = z[sender_inclusive_id, ],
    receiver_specific_z = z[receiver_specific_id, ],
    receiver_merged_z = z[receiver_merged_id, ]
  )]
  scores[, ecosystem_z := zscore((sender_primary_z + receiver_specific_z) / 2)]
  scores[, ecosystem_merged_sensitivity_z :=
           zscore((sender_primary_z + receiver_merged_z) / 2)]
  sender_cut <- median(scores$sender_primary_z)
  receiver_cut <- median(scores$receiver_specific_z)
  scores[, ecotype4 := fifelse(
    sender_primary_z >= sender_cut & receiver_specific_z >= receiver_cut, "Dual-high",
    fifelse(sender_primary_z >= sender_cut & receiver_specific_z < receiver_cut,
            "Fibroblast-dominant",
            fifelse(sender_primary_z < sender_cut & receiver_specific_z >= receiver_cut,
                    "Macrophage-dominant", "Dual-low")))
  ]
  scores[, ecotype4 := factor(ecotype4,
                              levels = c("Dual-low", "Fibroblast-dominant",
                                         "Macrophage-dominant", "Dual-high"))]
  list(scores = scores, coverage = coverage)
}

cox_one <- function(data, cohort, exposure, covariates) {
  vars <- c("os_time_months", "os_event", exposure, covariates)
  d <- as.data.frame(data[, ..vars])
  d <- d[complete.cases(d) & is.finite(d$os_time_months) & d$os_time_months > 0, , drop = FALSE]
  names(d)[match(exposure, names(d))] <- "exposure"
  rhs <- paste(c("exposure", covariates), collapse = " + ")
  fit <- coxph(as.formula(paste("Surv(os_time_months, os_event) ~", rhs)),
               data = d, ties = "efron", x = TRUE)
  sm <- summary(fit)$coefficients
  ci <- summary(fit)$conf.int
  zph <- cox.zph(fit, transform = "km")
  data.table(
    cohort = cohort, exposure = exposure,
    model = if (length(covariates)) "cohort_adjusted" else "univariable",
    covariates = paste(covariates, collapse = "+"), n = fit$n, events = fit$nevent,
    logHR = unname(sm["exposure", "coef"]),
    se_logHR = unname(sm["exposure", "se(coef)"]),
    HR_per_SD = unname(ci["exposure", "exp(coef)"]),
    ci_low = unname(ci["exposure", "lower .95"]),
    ci_high = unname(ci["exposure", "upper .95"]),
    p_value = unname(sm["exposure", "Pr(>|z|)"]),
    PH_exposure_p = unname(zph$table["exposure", "p"]),
    PH_global_p = unname(zph$table["GLOBAL", "p"])
  )
}

cohort_scores <- list()
coverage_all <- list()
ecology_all <- list()
survival_all <- list()

for (cohort in names(cohort_files)) {
  obj <- readRDS(cohort_files[[cohort]])
  expression <- obj$symbol_expression
  scored <- score_expression(expression, cohort, colnames(expression))
  s <- scored$scores
  coverage_all[[cohort]] <- scored$coverage
  meta <- copy(obj$metadata)
  setnames(meta, "geo_sample_id", "sample_id")
  d <- merge(meta, s, by = "sample_id", all.x = TRUE, sort = FALSE)
  if (cohort == "GSE62254") {
    clin <- fread(gse62254_clinical_file)
    d[, patient_id := as.character(patient_id)]
    clin[, patient_id := as.character(patient_id)]
    d <- merge(d, clin, by = "patient_id", all.x = TRUE, sort = FALSE)
    d <- d[inclusion_expression == TRUE]
    adjusted_covariates <- c("age", "sex", "stage_numeric")
  } else if (cohort == "GSE84437") {
    d <- d[inclusion_expression == TRUE]
    adjusted_covariates <- c("age", "sex", "pT", "pN")
  } else {
    d <- d[inclusion_expression == TRUE & explicit_qc_exclusion == FALSE]
    adjusted_covariates <- c("age", "sex", "stage")
  }
  if (nrow(d) != ncol(expression)) stop("External cohort merge failure: ", cohort)
  cohort_scores[[cohort]] <- copy(d)
  ecology_all[[cohort]] <- data.table(
    cohort = cohort, n = nrow(d),
    sender_receiver_rho = cor(d$sender_primary_z, d$receiver_specific_z,
                              method = "spearman", use = "complete.obs"),
    sender_receiver_p = cor.test(d$sender_primary_z, d$receiver_specific_z,
                                 method = "spearman", exact = FALSE)$p.value,
    sender_merged_receiver_rho = cor(d$sender_primary_z, d$receiver_merged_z,
                                     method = "spearman", use = "complete.obs"),
    sender_merged_receiver_p = cor.test(d$sender_primary_z, d$receiver_merged_z,
                                        method = "spearman", exact = FALSE)$p.value,
    PLAT_inclusive_receiver_rho = cor(d$sender_PLAT_inclusive_z,
                                      d$receiver_specific_z,
                                      method = "spearman", use = "complete.obs")
  )
  for (exposure in c("sender_primary_z", "receiver_specific_z", "ecosystem_z",
                     "ecosystem_merged_sensitivity_z")) {
    survival_all[[paste(cohort, exposure, "uni")]] <-
      cox_one(d, cohort, exposure, character())
    survival_all[[paste(cohort, exposure, "adj")]] <-
      cox_one(d, cohort, exposure, adjusted_covariates)
  }
  fwrite(d, file.path(out_dir, paste0(cohort, "_frozen_PLAT_ecosystem_scores.tsv.gz")),
         sep = "\t")
}

coverage_table <- rbindlist(coverage_all)
ecology_table <- rbindlist(ecology_all)
survival_table <- rbindlist(survival_all)
survival_table[, FDR_BH_within_model := p.adjust(p_value, "BH"), by = model]
fwrite(coverage_table, file.path(out_dir, "external_GC_signature_coverage.tsv"), sep = "\t")
fwrite(ecology_table, file.path(out_dir, "external_GC_ecology_correlations.tsv"), sep = "\t")
fwrite(survival_table, file.path(out_dir, "external_GC_continuous_Cox_models.tsv"), sep = "\t")

# Fixed-effect and DerSimonian-Laird random-effects summaries for the ecosystem score.
tcga_cox <- fread(tcga_cox_file)
tcga_row <- tcga_cox[exposure == "ecosystem_z" & tier == "clinical_purity"]
tcga_meta <- data.table(
  cohort = "TCGA-STAD", exposure = "ecosystem_z", model = "cohort_adjusted",
  n = tcga_row$n, events = tcga_row$events,
  logHR = log(tcga_row$HR_per_SD),
  se_logHR = (log(tcga_row$ci_high) - log(tcga_row$ci_low)) / (2 * qnorm(0.975)),
  HR_per_SD = tcga_row$HR_per_SD, ci_low = tcga_row$ci_low,
  ci_high = tcga_row$ci_high, p_value = tcga_row$p_value
)
meta_input <- rbindlist(list(
  tcga_meta,
  survival_table[exposure == "ecosystem_z" & model == "cohort_adjusted",
                 .(cohort, exposure, model, n, events, logHR, se_logHR,
                   HR_per_SD, ci_low, ci_high, p_value)]
), fill = TRUE)

meta_analyze <- function(d) {
  wi <- 1 / d$se_logHR^2
  fixed <- sum(wi * d$logHR) / sum(wi)
  se_fixed <- sqrt(1 / sum(wi))
  Q <- sum(wi * (d$logHR - fixed)^2)
  df <- nrow(d) - 1L
  C <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- max(0, (Q - df) / C)
  wr <- 1 / (d$se_logHR^2 + tau2)
  random <- sum(wr * d$logHR) / sum(wr)
  se_random <- sqrt(1 / sum(wr))
  rbindlist(list(
    data.table(model = "fixed_effect", k = nrow(d), logHR = fixed,
               se_logHR = se_fixed, HR = exp(fixed),
               ci_low = exp(fixed - qnorm(0.975) * se_fixed),
               ci_high = exp(fixed + qnorm(0.975) * se_fixed),
               p_value = 2 * pnorm(-abs(fixed / se_fixed)),
               Q = Q, Q_df = df, Q_p = pchisq(Q, df, lower.tail = FALSE),
               I2_percent = max(0, (Q - df) / Q) * 100, tau2 = tau2),
    data.table(model = "random_effect_DL", k = nrow(d), logHR = random,
               se_logHR = se_random, HR = exp(random),
               ci_low = exp(random - qnorm(0.975) * se_random),
               ci_high = exp(random + qnorm(0.975) * se_random),
               p_value = 2 * pnorm(-abs(random / se_random)),
               Q = Q, Q_df = df, Q_p = pchisq(Q, df, lower.tail = FALSE),
               I2_percent = max(0, (Q - df) / Q) * 100, tau2 = tau2)
  ))
}
meta_summary <- meta_analyze(meta_input)
fwrite(meta_input, file.path(out_dir, "GC_cohort_adjusted_ecosystem_OS_meta_input.tsv"), sep = "\t")
fwrite(meta_summary, file.path(out_dir, "GC_cohort_adjusted_ecosystem_OS_meta_summary.tsv"), sep = "\t")

# Kim2018 anti-PD-1 cohort: frozen continuous scoring and bias-reduced logistic regression.
kim_tpm <- fread(kim_tpm_file)
setnames(kim_tpm, 1L, "gene_symbol")
kim_tpm[, gene_symbol := toupper(trimws(gene_symbol))]
kim_tpm <- kim_tpm[!is.na(gene_symbol) & nzchar(gene_symbol)]
kim_mat <- rowsum(as.matrix(kim_tpm[, -1]), group = kim_tpm$gene_symbol,
                  reorder = FALSE, na.rm = TRUE)
kim_log <- log2(kim_mat + 1)
kim_scored <- score_expression(kim_log, "Kim2018", colnames(kim_log))
kim <- copy(kim_scored$scores)
setnames(kim, "sample_id", "patient_id")
kim_clinical <- fread(kim_clinical_file)
kim <- merge(kim_clinical, kim, by = "patient_id", all.x = TRUE, sort = FALSE)
if (nrow(kim) != 45L || kim[, sum(response_binary == 1L)] != 12L) {
  stop("Kim2018 response alignment failure")
}
kim[, MSI_binary := as.integer(msi_type_pentaplex == "MSI")]

kim_logistic <- function(exposure, adjusted = FALSE) {
  vars <- c("response_binary", exposure, if (adjusted) "MSI_binary")
  d <- as.data.frame(kim[, ..vars])
  d <- d[complete.cases(d), , drop = FALSE]
  names(d)[match(exposure, names(d))] <- "exposure"
  if (adjusted) {
    fit <- glm(response_binary ~ exposure + MSI_binary,
               data = d, family = binomial("logit"), method = brglm2::brglmFit,
               type = "AS_mean")
  } else {
    fit <- glm(response_binary ~ exposure,
               data = d, family = binomial("logit"), method = brglm2::brglmFit,
               type = "AS_mean")
  }
  sm <- summary(fit)$coefficients
  beta <- unname(sm["exposure", "Estimate"])
  se <- unname(sm["exposure", "Std. Error"])
  data.table(
    exposure = exposure, model = if (adjusted) "MSI_adjusted" else "univariable",
    n = nrow(d), responders = sum(d$response_binary == 1L),
    OR_per_SD_for_response = exp(beta),
    ci_low = exp(beta - qnorm(0.975) * se),
    ci_high = exp(beta + qnorm(0.975) * se),
    p_value = unname(sm["exposure", "Pr(>|z|)"])
  )
}

kim_exposures <- c("sender_primary_z", "receiver_specific_z", "ecosystem_z",
                   "ecosystem_merged_sensitivity_z")
kim_models <- rbindlist(lapply(kim_exposures, function(exposure) {
  rbindlist(list(kim_logistic(exposure, FALSE), kim_logistic(exposure, TRUE)))
}))
kim_models[, FDR_BH_within_model := p.adjust(p_value, "BH"), by = model]

extract_joint_model <- function(fit, model_id) {
  sm <- summary(fit)$coefficients
  terms <- setdiff(rownames(sm), "(Intercept)")
  rbindlist(lapply(terms, function(term) {
    beta <- unname(sm[term, "Estimate"])
    se <- unname(sm[term, "Std. Error"])
    data.table(model_id = model_id, term = term, n = nobs(fit),
               OR = exp(beta), ci_low = exp(beta - qnorm(0.975) * se),
               ci_high = exp(beta + qnorm(0.975) * se),
               p_value = unname(sm[term, "Pr(>|z|)"]))
  }))
}
kim_joint_fits <- list(
  sender_receiver = glm(
    response_binary ~ sender_primary_z + receiver_specific_z,
    data = kim, family = binomial("logit"), method = brglm2::brglmFit,
    type = "AS_mean"),
  sender_receiver_MSI = glm(
    response_binary ~ sender_primary_z + receiver_specific_z + MSI_binary,
    data = kim, family = binomial("logit"), method = brglm2::brglmFit,
    type = "AS_mean"),
  sender_receiver_interaction = glm(
    response_binary ~ sender_primary_z * receiver_specific_z,
    data = kim, family = binomial("logit"), method = brglm2::brglmFit,
    type = "AS_mean")
)
kim_joint_models <- rbindlist(lapply(names(kim_joint_fits), function(model_id) {
  extract_joint_model(kim_joint_fits[[model_id]], model_id)
}))
kim_joint_models[, FDR_BH_exploratory_family := p.adjust(p_value, "BH")]

kim_auc <- rbindlist(lapply(kim_exposures, function(exposure) {
  roc_obj <- pROC::roc(kim$response_binary, kim[[exposure]], levels = c(0, 1),
                       direction = "<", quiet = TRUE)
  ci <- pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = 2000,
                     stratified = TRUE, progress = "none")
  data.table(exposure = exposure, n = nrow(kim), responders = sum(kim$response_binary),
             signed_AUC_higher_predicts_response = as.numeric(pROC::auc(roc_obj)),
             ci_low = as.numeric(ci[1]), ci_high = as.numeric(ci[3]))
}))
kim_group <- kim[, .(n = .N, responders = sum(response_binary),
                     response_rate = mean(response_binary)), by = ecotype4]
kim_fisher <- fisher.test(table(kim$ecotype4, kim$response_binary))

fwrite(kim_scored$coverage, file.path(out_dir, "Kim2018_signature_coverage.tsv"), sep = "\t")
fwrite(kim, file.path(out_dir, "Kim2018_frozen_PLAT_ecosystem_scores.tsv"), sep = "\t")
fwrite(kim_models, file.path(out_dir, "Kim2018_bias_reduced_response_models.tsv"), sep = "\t")
fwrite(kim_joint_models,
       file.path(out_dir, "Kim2018_continuous_joint_state_response_models.tsv"), sep = "\t")
fwrite(kim_auc, file.path(out_dir, "Kim2018_signed_response_AUC.tsv"), sep = "\t")
fwrite(kim_group, file.path(out_dir, "Kim2018_ecotype4_response_summary.tsv"), sep = "\t")

palette4 <- c("Dual-low" = "#264653", "Fibroblast-dominant" = "#E9C46A",
              "Macrophage-dominant" = "#2A9D8F", "Dual-high" = "#E76F51")
plot_ecology <- rbindlist(lapply(names(cohort_scores), function(cohort) {
  cohort_scores[[cohort]][, .(cohort, sender_primary_z, receiver_specific_z)]
}))
p1 <- ggplot(plot_ecology, aes(sender_primary_z, receiver_specific_z)) +
  geom_point(size = 0.8, alpha = 0.35, color = "#2A9D8F") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.55, color = "#264653") +
  facet_wrap(~cohort, nrow = 1) +
  labs(x = "Fibroblast-state score (z)", y = "C1QC/APOE macrophage score (z)",
       title = "A  External gastric-cancer cohorts") +
  theme_classic(base_size = 9)

forest <- rbindlist(list(meta_input[, .(cohort, HR = HR_per_SD, ci_low, ci_high)],
                         meta_summary[model == "random_effect_DL",
                                      .(cohort = "Random-effects summary", HR, ci_low, ci_high)]))
forest[, cohort := factor(cohort, levels = rev(cohort))]
p2 <- ggplot(forest, aes(HR, cohort)) +
  geom_vline(xintercept = 1, color = "grey70", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.1,
                 color = "#264653") +
  geom_point(size = 2.2, color = "#E76F51") +
  scale_x_log10() +
  labs(x = "Adjusted HR per ecosystem SD", y = NULL,
       title = "B  Overall survival") +
  theme_classic(base_size = 9)

p3 <- ggplot(kim, aes(ecotype4, ecosystem_z, fill = factor(response_binary))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.72) +
  geom_jitter(aes(shape = factor(response_binary)), width = 0.13,
              size = 1.6, alpha = 0.75) +
  scale_fill_manual(values = c(`0` = "#B8B8B8", `1` = "#E76F51")) +
  labs(x = NULL, y = "Ecosystem score (z)", fill = "Response", shape = "Response",
       title = "C  Kim2018 anti-PD-1 cohort") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")

diagnostic <- p1 / (p2 | p3) + plot_layout(heights = c(1, 1.05))
ggsave(file.path(out_dir, "external_GC_PLAT_ecosystem_diagnostic_panels.png"),
       diagnostic, width = 11, height = 7.5, dpi = 320, bg = "white")
ggsave(file.path(out_dir, "external_GC_PLAT_ecosystem_diagnostic_panels.pdf"),
       diagnostic, width = 11, height = 7.5, device = cairo_pdf, bg = "white")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  signature_status = "unchanged frozen GSE183904 signatures",
  scoring = "within-sample normalized mean gene rank; within-cohort z standardization",
  cohorts = lapply(cohort_scores, nrow),
  meta_analysis = meta_summary,
  Kim2018 = list(n = nrow(kim), responders = sum(kim$response_binary),
                 ecotype4_fisher_p = kim_fisher$p.value),
  claim_boundary = paste(
    "External survival and response associations are observational.",
    "Kim2018 is a single-arm treatment cohort and does not estimate treatment interaction."
  ),
  input_sha256 = as.list(setNames(vapply(inputs, digest, character(1),
                                         algo = "sha256", file = TRUE), basename(inputs)))
)
write_json(audit, file.path(audit_dir, "external_GC_PLAT_ecosystem_validation.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nExternal signature coverage:\n")
print(coverage_table)
cat("\nExternal ecology correlations:\n")
print(ecology_table)
cat("\nExternal continuous Cox models:\n")
print(survival_table)
cat("\nOS meta-analysis:\n")
print(meta_input)
print(meta_summary)
cat("\nKim2018 response models:\n")
print(kim_models)
cat("\nKim2018 exploratory continuous joint-state models:\n")
print(kim_joint_models)
cat("\nKim2018 signed AUC:\n")
print(kim_auc)
cat("\nKim2018 four-group response:\n")
print(kim_group)
cat("Fisher exact P = ", signif(kim_fisher$p.value, 5), "\n", sep = "")
