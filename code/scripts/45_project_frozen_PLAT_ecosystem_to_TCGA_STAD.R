options(stringsAsFactors = FALSE, warn = 1)
options(error = function() {
  traceback(20L)
  quit(save = "no", status = 1L, runLast = FALSE)
})
set.seed(2026072945L)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(sandwich)
  library(lmtest)
  library(jsonlite)
  library(digest)
  library(ggplot2)
  library(patchwork)
})

signature_file <- file.path(project, "resources", "gene_sets",
                            "GSE183904_PLAT_ecosystem_bulk_signatures.tsv")
expression_file <- "I:/ECM-new/ECM-public/data/processed/bulk/tcga_stad_full_symbol_log2cpm_412.rds"
clinical_file <- "I:/ECM-new/ECM-public/results/known_subtypes_v2/TCGA_STAD_current412_with_official_TCGA2014_subtype.tsv.gz"
prereg_file <- file.path(project, "analysis", "audit",
                         "TCGA_PLAT_ecosystem_analysis_preregistration.json")
boundary_file <- file.path(project, "analysis", "results", "PLAT_axis", "GSE183904",
                           "state_boundary_and_signatures",
                           "GSE183904_merged_receiver_abundance_inference.tsv")
out_dir <- file.path(project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
                     "frozen_state_projection")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(signature_file, expression_file, clinical_file, prereg_file, boundary_file)
if (!all(file.exists(inputs))) stop("Missing TCGA projection input")

signatures_dt <- fread(signature_file)
expr_obj <- readRDS(expression_file)
expr <- expr_obj$expression
clinical <- fread(clinical_file)
if (!identical(colnames(expr), clinical$sample_submitter_id)) {
  idx <- match(colnames(expr), clinical$sample_submitter_id)
  if (anyNA(idx)) stop("TCGA expression/clinical alignment failure")
  clinical <- clinical[idx]
}
if (ncol(expr) != 412L || nrow(clinical) != 412L) stop("Unexpected TCGA cohort size")

sender_id <- "sender_CXCL14_POSTN_state_axis_excluded"
sender_inclusive_id <- "sender_CXCL14_POSTN_PLAT_state_inclusive"
receiver_specific_id <- "receiver_C1QC_APOE_specific_axis_excluded"
receiver_merged_id <- "receiver_merged_C1QC_axis_excluded"
required_ids <- c(sender_id, sender_inclusive_id, receiver_specific_id, receiver_merged_id)
if (!all(required_ids %in% signatures_dt$signature_id)) stop("Frozen signature missing")

sig_list <- split(signatures_dt$gene_symbol, signatures_dt$signature_id)
sig_list <- sig_list[required_ids]
coverage <- rbindlist(lapply(names(sig_list), function(id) {
  genes <- unique(sig_list[[id]])
  data.table(signature_id = id, frozen_genes = length(genes),
             mapped_genes = sum(genes %in% rownames(expr)),
             coverage_fraction = mean(genes %in% rownames(expr)),
             missing_genes = paste(setdiff(genes, rownames(expr)), collapse = ";"))
}))
if (coverage[, any(mapped_genes < 15L)]) stop("Insufficient TCGA signature coverage")
fwrite(coverage, file.path(out_dir, "TCGA_STAD_frozen_signature_coverage.tsv"), sep = "\t")

# One within-sample ranking pass supplies all frozen signature scores.
sig_indices <- lapply(sig_list, function(g) match(intersect(g, rownames(expr)), rownames(expr)))
score_raw <- vapply(seq_len(ncol(expr)), function(j) {
  normalized_rank <- rank(expr[, j], ties.method = "average") / nrow(expr)
  vapply(sig_indices, function(idx) mean(normalized_rank[idx]), numeric(1))
}, numeric(length(sig_indices)))
rownames(score_raw) <- names(sig_indices)
colnames(score_raw) <- colnames(expr)
zscore <- function(x) as.numeric((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
score_z <- t(apply(score_raw, 1L, zscore))

score_names <- c(
  sender_primary = sender_id,
  sender_PLAT_inclusive = sender_inclusive_id,
  receiver_specific = receiver_specific_id,
  receiver_merged = receiver_merged_id
)
scores <- data.table(sample_submitter_id = colnames(expr))
for (nm in names(score_names)) {
  id <- unname(score_names[nm])
  scores[[paste0(nm, "_rankmean")]] <- score_raw[id, ]
  scores[[paste0(nm, "_z")]] <- score_z[id, ]
}
if (!"PLAT" %in% rownames(expr)) stop("PLAT absent from TCGA matrix")
scores[, PLAT_log2CPM := as.numeric(expr["PLAT", ])]
scores[, PLAT_z := zscore(PLAT_log2CPM)]
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

dt <- merge(clinical, scores, by = "sample_submitter_id", all.x = TRUE, sort = FALSE)
if (nrow(dt) != 412L || dt[, anyNA(sender_primary_z)]) stop("TCGA score merge failure")
dt[tcga_subtype == "", tcga_subtype := NA_character_]
dt[, stage_group_frozen := fifelse(grepl("STAGE IV", stage), "IV",
                             fifelse(grepl("STAGE III", stage), "III",
                               fifelse(grepl("STAGE II", stage), "II",
                                 fifelse(grepl("STAGE I", stage), "I", NA_character_))))]
dt[, stage_ordinal := match(stage_group_frozen, c("I", "II", "III", "IV"))]
dt[, os_event_frozen := fifelse(!is.na(os_event), as.integer(os_event),
                                fifelse(grepl("^1", os_status), 1L,
                                        fifelse(grepl("^0", os_status), 0L, NA_integer_)))]

fwrite(dt, file.path(out_dir, "TCGA_STAD_PLAT_ecosystem_scores_clinical.tsv.gz"),
       sep = "\t")
saveRDS(list(scores = scores, linked_data = dt, signature_coverage = coverage,
             frozen_signatures = sig_list,
             scoring = "within-sample normalized mean rank; cohort z standardization"),
        file.path(out_dir, "TCGA_STAD_PLAT_ecosystem_projection.rds"))

rank_residual_association <- function(data, x, y, covariates, analysis_id) {
  vars <- c(x, y, covariates)
  d <- data[, ..vars]
  numeric_matrix <- as.matrix(d)
  storage.mode(numeric_matrix) <- "double"
  ok <- complete.cases(d) & apply(numeric_matrix, 1L, function(v) all(is.finite(v)))
  d <- d[ok]
  xr <- rank(d[[x]], ties.method = "average")
  yr <- rank(d[[y]], ties.method = "average")
  if (length(covariates)) {
    C <- as.data.frame(d[, ..covariates])
    xr_res <- residuals(lm(xr ~ ., data = C))
    yr_res <- residuals(lm(yr ~ ., data = C))
  } else {
    C <- data.frame(row_id = seq_along(xr))[, FALSE, drop = FALSE]
    xr_res <- xr
    yr_res <- yr
  }
  rho <- unname(cor(xr_res, yr_res))
  partial_df <- length(xr) - length(covariates) - 2L
  partial_t <- rho * sqrt(partial_df / pmax(1e-12, 1 - rho^2))
  partial_p <- 2 * pt(-abs(partial_t), df = partial_df)
  model_df <- data.frame(y_rank = yr, x_rank = xr, C, check.names = TRUE)
  fit <- lm(y_rank ~ ., data = model_df)
  vc <- sandwich::vcovHC(fit, type = "HC3")
  ct <- lmtest::coeftest(fit, vcov. = vc)
  data.table(
    analysis_id = analysis_id, x = x, y = y,
    adjustment = if (length(covariates)) paste(covariates, collapse = "+") else "none",
    n = nrow(model_df), rank_residual_rho = rho,
    partial_df = partial_df, partial_p = partial_p,
    HC3_beta_x_rank = unname(ct["x_rank", 1]),
    HC3_se = unname(ct["x_rank", 2]),
    HC3_p = unname(ct["x_rank", 4])
  )
}

pair_specs <- list(
  list("primary_raw", "sender_primary_z", "receiver_specific_z", character()),
  list("primary_purity_adjusted", "sender_primary_z", "receiver_specific_z", "absolute_purity"),
  list("primary_lineage_adjusted", "sender_primary_z", "receiver_specific_z",
       c("Fibroblasts", "Monocytic lineage")),
  list("primary_full_adjusted", "sender_primary_z", "receiver_specific_z",
       c("absolute_purity", "Fibroblasts", "Monocytic lineage")),
  list("PLAT_inclusive_raw", "sender_PLAT_inclusive_z", "receiver_specific_z", character()),
  list("PLAT_inclusive_full_adjusted", "sender_PLAT_inclusive_z", "receiver_specific_z",
       c("absolute_purity", "Fibroblasts", "Monocytic lineage")),
  list("merged_receiver_raw", "sender_primary_z", "receiver_merged_z", character()),
  list("merged_receiver_full_adjusted", "sender_primary_z", "receiver_merged_z",
       c("absolute_purity", "Fibroblasts", "Monocytic lineage")),
  list("PLAT_sender_raw", "PLAT_log2CPM", "sender_primary_z", character()),
  list("PLAT_sender_full_adjusted", "PLAT_log2CPM", "sender_primary_z",
       c("absolute_purity", "Fibroblasts", "Monocytic lineage")),
  list("receiver_boundary_raw", "receiver_specific_z", "receiver_merged_z", character()),
  list("receiver_boundary_full_adjusted", "receiver_specific_z", "receiver_merged_z",
       c("absolute_purity", "Fibroblasts", "Monocytic lineage"))
)
ecology_associations <- rbindlist(lapply(pair_specs, function(s) {
  rank_residual_association(dt, s[[2]], s[[3]], s[[4]], s[[1]])
}))
ecology_associations[, `:=`(
  partial_FDR_BH = p.adjust(partial_p, "BH"),
  HC3_FDR_BH = p.adjust(HC3_p, "BH")
)]
fwrite(ecology_associations,
       file.path(out_dir, "TCGA_STAD_continuous_ecology_associations.tsv"), sep = "\t")

hc3_endpoint_model <- function(data, endpoint, exposure, covariates, family, transform) {
  y <- transform(data[[endpoint]])
  model_df <- data.frame(y = zscore(y), exposure = data[[exposure]])
  for (v in covariates) model_df[[make.names(v)]] <- data[[v]]
  ok <- complete.cases(model_df)
  model_df <- model_df[ok, , drop = FALSE]
  fit <- lm(y ~ ., data = model_df)
  vc <- sandwich::vcovHC(fit, type = "HC3")
  ct <- lmtest::coeftest(fit, vcov. = vc)
  data.table(
    family = family, endpoint = endpoint, exposure = exposure,
    endpoint_transform = deparse(substitute(transform)),
    covariates = paste(covariates, collapse = "+"), n = nrow(model_df),
    beta_per_exposure_SD = unname(ct["exposure", 1]),
    HC3_se = unname(ct["exposure", 2]),
    HC3_p = unname(ct["exposure", 4]),
    ci_low = unname(ct["exposure", 1] - qnorm(0.975) * ct["exposure", 2]),
    ci_high = unname(ct["exposure", 1] + qnorm(0.975) * ct["exposure", 2])
  )
}

identity_transform <- function(x) x
log1p_transform <- function(x) log1p(pmax(0, x))
exposures <- c("sender_primary_z", "receiver_specific_z", "ecosystem_z")
immune_endpoints <- c(
  "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "B lineage", "NK cells",
  "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils", "Endothelial cells",
  "Fibroblasts", "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_ANGIOGENESIS",
  "HALLMARK_COAGULATION", "HALLMARK_COMPLEMENT",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA",
  "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "T_CELL_INFLAMED_18", "CYTOLYTIC_ACTIVITY", "TLS_B_CELL_CHEMOKINE"
)
immune_covariates <- c("absolute_purity", "tcga_subtype",
                       "TMB_NONSYNONYMOUS.x", "FRACTION_GENOME_ALTERED.x")
immune_models <- rbindlist(lapply(exposures, function(exposure) {
  rbindlist(lapply(immune_endpoints, function(endpoint) {
    hc3_endpoint_model(dt, endpoint, exposure, immune_covariates,
                       "immune_microenvironment", identity_transform)
  }))
}))

molecular_specs <- list(
  ANEUPLOIDY_SCORE = identity_transform,
  FRACTION_GENOME_ALTERED.x = identity_transform,
  MSI_SCORE_MANTIS = log1p_transform,
  MSI_SENSOR_SCORE = log1p_transform,
  MUTATION_COUNT = log1p_transform,
  TMB_NONSYNONYMOUS.x = log1p_transform
)
molecular_models <- rbindlist(lapply(exposures, function(exposure) {
  rbindlist(lapply(names(molecular_specs), function(endpoint) {
    hc3_endpoint_model(dt, endpoint, exposure, "absolute_purity",
                       "molecular", molecular_specs[[endpoint]])
  }))
}))

clinical_models <- rbindlist(lapply(exposures, function(exposure) {
  hc3_endpoint_model(dt, "stage_ordinal", exposure,
                     c("age", "sex", "absolute_purity"),
                     "clinical_stage", identity_transform)
}))
endpoint_models <- rbindlist(list(immune_models, molecular_models, clinical_models), fill = TRUE)
endpoint_models[, FDR_BH_within_family := p.adjust(HC3_p, "BH"), by = family]
fwrite(endpoint_models, file.path(out_dir, "TCGA_STAD_endpoint_associations_HC3.tsv"),
       sep = "\t")

subtype_tests <- rbindlist(lapply(exposures, function(exposure) {
  d <- dt[!is.na(tcga_subtype) & is.finite(get(exposure))]
  kt <- kruskal.test(d[[exposure]] ~ factor(d$tcga_subtype))
  data.table(exposure = exposure, n = nrow(d), subtypes = uniqueN(d$tcga_subtype),
             kruskal_chisq = unname(kt$statistic), df = unname(kt$parameter),
             p_value = kt$p.value)
}))
subtype_tests[, FDR_BH := p.adjust(p_value, "BH")]
subtype_summary <- rbindlist(lapply(exposures, function(exposure) {
  dt[!is.na(tcga_subtype), .(
    n = sum(is.finite(get(exposure))),
    median = median(get(exposure), na.rm = TRUE),
    q1 = quantile(get(exposure), 0.25, na.rm = TRUE),
    q3 = quantile(get(exposure), 0.75, na.rm = TRUE)
  ), by = tcga_subtype][, exposure := exposure]
}))
fwrite(subtype_tests, file.path(out_dir, "TCGA_STAD_subtype_omnibus_tests.tsv"), sep = "\t")
fwrite(subtype_summary, file.path(out_dir, "TCGA_STAD_subtype_score_summary.tsv"), sep = "\t")

ecotype_counts <- dt[, .N, by = ecotype4][order(ecotype4)]
ecotype_subtype <- dcast(dt[!is.na(tcga_subtype)], ecotype4 ~ tcga_subtype,
                         fun.aggregate = length, value.var = "sample_submitter_id")
ecotype_stage <- dcast(dt[!is.na(stage_group_frozen)], ecotype4 ~ stage_group_frozen,
                       fun.aggregate = length, value.var = "sample_submitter_id")
ecotype_subtype_test <- chisq.test(table(dt$ecotype4, dt$tcga_subtype))
ecotype_stage_test <- chisq.test(table(dt$ecotype4, dt$stage_group_frozen))
ecotype_continuous_summary <- dt[, lapply(.SD, function(x) {
  paste(c(n = sum(is.finite(x)), median = median(x, na.rm = TRUE),
          q1 = quantile(x, 0.25, na.rm = TRUE), q3 = quantile(x, 0.75, na.rm = TRUE)),
        collapse = ";")
}), by = ecotype4, .SDcols = c("ecosystem_z", "TMB_NONSYNONYMOUS.x",
                              "MSI_SCORE_MANTIS", "FRACTION_GENOME_ALTERED.x",
                              immune_endpoints)]
fwrite(ecotype_counts, file.path(out_dir, "TCGA_STAD_ecotype4_counts.tsv"), sep = "\t")
fwrite(ecotype_subtype, file.path(out_dir, "TCGA_STAD_ecotype4_by_subtype.tsv"), sep = "\t")
fwrite(ecotype_stage, file.path(out_dir, "TCGA_STAD_ecotype4_by_stage.tsv"), sep = "\t")
fwrite(ecotype_continuous_summary,
       file.path(out_dir, "TCGA_STAD_ecotype4_continuous_summary.tsv"), sep = "\t")

cox_one <- function(data, exposure, tier) {
  covariates <- switch(tier,
                       univariable = character(),
                       clinical = c("age", "sex", "stage_group_frozen"),
                       clinical_purity = c("age", "sex", "stage_group_frozen", "absolute_purity"))
  vars <- c("os_months", "os_event_frozen", exposure, covariates)
  d <- as.data.frame(data[, ..vars])
  d <- d[complete.cases(d) & d$os_months > 0, , drop = FALSE]
  names(d)[match(exposure, names(d))] <- "exposure"
  rhs <- paste(c("exposure", covariates), collapse = " + ")
  fit <- coxph(as.formula(paste("Surv(os_months, os_event_frozen) ~", rhs)),
               data = d, ties = "efron", x = TRUE)
  sm <- summary(fit)$coefficients
  ci <- summary(fit)$conf.int
  zph <- cox.zph(fit, transform = "km")
  data.table(
    exposure = exposure, tier = tier, covariates = paste(covariates, collapse = "+"),
    n = fit$n, events = fit$nevent,
    HR_per_SD = unname(ci["exposure", "exp(coef)"]),
    ci_low = unname(ci["exposure", "lower .95"]),
    ci_high = unname(ci["exposure", "upper .95"]),
    p_value = unname(sm["exposure", "Pr(>|z|)"]),
    PH_exposure_p = unname(zph$table["exposure", "p"]),
    PH_global_p = unname(zph$table["GLOBAL", "p"]), AIC = AIC(fit)
  )
}

survival_exposures <- c("sender_primary_z", "receiver_specific_z", "ecosystem_z",
                        "ecosystem_merged_sensitivity_z")
cox_models <- rbindlist(lapply(survival_exposures, function(exposure) {
  rbindlist(lapply(c("univariable", "clinical", "clinical_purity"), function(tier) {
    cox_one(dt, exposure, tier)
  }))
}))
cox_models[, FDR_BH_within_tier := p.adjust(p_value, "BH"), by = tier]
fwrite(cox_models, file.path(out_dir, "TCGA_STAD_continuous_Cox_models.tsv"), sep = "\t")

km_data <- dt[is.finite(os_months) & os_months > 0 & !is.na(os_event_frozen) & !is.na(ecotype4)]
km_fit <- survfit(Surv(os_months, os_event_frozen) ~ ecotype4, data = km_data)
km_diff <- survdiff(Surv(os_months, os_event_frozen) ~ ecotype4, data = km_data)
km_logrank_p <- pchisq(km_diff$chisq, df = length(km_diff$n) - 1L, lower.tail = FALSE)
km_table <- summary(km_fit)$table
km_summary <- data.table(
  group = sub("^ecotype4=", "", rownames(km_table)),
  n = as.integer(km_table[, "records"]),
  median_os_months = unname(km_table[, "median"]),
  lower_95 = unname(km_table[, "0.95LCL"]),
  upper_95 = unname(km_table[, "0.95UCL"]),
  global_logrank_p = km_logrank_p
)
fwrite(km_summary, file.path(out_dir, "TCGA_STAD_ecotype4_KM_summary.tsv"), sep = "\t")

# Diagnostic panels are analysis outputs; manuscript panels will be redrawn after validation.
palette4 <- c("Dual-low" = "#264653", "Fibroblast-dominant" = "#E9C46A",
              "Macrophage-dominant" = "#2A9D8F", "Dual-high" = "#E76F51")
p1 <- ggplot(dt, aes(sender_primary_z, receiver_specific_z, color = ecotype4)) +
  geom_hline(yintercept = receiver_cut, linewidth = 0.35, color = "grey70") +
  geom_vline(xintercept = sender_cut, linewidth = 0.35, color = "grey70") +
  geom_point(size = 1.8, alpha = 0.78) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black",
              fill = "grey75", linewidth = 0.55) +
  scale_color_manual(values = palette4) +
  labs(x = "Fibroblast-state score (z)", y = "C1QC/APOE macrophage score (z)",
       color = NULL, title = "A  Frozen-state projection") +
  theme_classic(base_size = 10) + theme(legend.position = "bottom")

subtype_plot <- dt[!is.na(tcga_subtype)]
p2 <- ggplot(subtype_plot, aes(tcga_subtype, ecosystem_z, fill = tcga_subtype)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.78) +
  geom_jitter(width = 0.15, size = 0.65, alpha = 0.35) +
  labs(x = NULL, y = "Ecosystem score (z)", title = "B  TCGA molecular subtype") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))

immune_plot_dt <- endpoint_models[
  family == "immune_microenvironment" & exposure == "ecosystem_z"
][order(beta_per_exposure_SD)]
immune_plot_dt[, endpoint_label := factor(endpoint, levels = endpoint)]
p3 <- ggplot(immune_plot_dt, aes(beta_per_exposure_SD, endpoint_label,
                                 color = FDR_BH_within_family < 0.05)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.4) +
  geom_point(size = 1.8) +
  scale_color_manual(values = c(`TRUE` = "#E76F51", `FALSE` = "#7A7A7A")) +
  labs(x = "Adjusted beta per ecosystem SD", y = NULL,
       title = "C  Immune and ECM context", color = "BH FDR < 0.05") +
  theme_classic(base_size = 8) + theme(legend.position = "bottom")

diagnostic <- (p1 | p2) / p3 + plot_layout(heights = c(1, 1.55))
ggsave(file.path(out_dir, "TCGA_STAD_PLAT_ecosystem_diagnostic_panels.png"),
       diagnostic, width = 10.5, height = 10.5, dpi = 320, bg = "white")
ggsave(file.path(out_dir, "TCGA_STAD_PLAT_ecosystem_diagnostic_panels.pdf"),
       diagnostic, width = 10.5, height = 10.5, device = cairo_pdf, bg = "white")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  cohort = "TCGA-STAD primary tumors", n = nrow(dt),
  receiver_decision = paste(
    "Cluster-3-specific receiver is primary because the preregistered merged cluster 3+5",
    "single-cell abundance sensitivity did not preserve the sender-receiver ecology."
  ),
  score_method = "within-sample normalized mean gene rank; cohort z standardization",
  ecotype4_cutpoints = list(sender_median_z = sender_cut, receiver_median_z = receiver_cut),
  signature_coverage = coverage,
  ecotype4_counts = ecotype_counts,
  ecotype_subtype_chisq = list(statistic = unname(ecotype_subtype_test$statistic),
                               df = unname(ecotype_subtype_test$parameter),
                               p_value = ecotype_subtype_test$p.value),
  ecotype_stage_chisq = list(statistic = unname(ecotype_stage_test$statistic),
                             df = unname(ecotype_stage_test$parameter),
                             p_value = ecotype_stage_test$p.value),
  km_global_logrank_p = km_logrank_p,
  claim_boundary = paste(
    "Observational projection of transcript-defined cell states.",
    "No direct ligand-receptor mechanism or cell-cell causality is inferred."
  ),
  input_sha256 = as.list(setNames(vapply(inputs, digest, character(1),
                                         algo = "sha256", file = TRUE), basename(inputs)))
)
write_json(audit, file.path(audit_dir, "TCGA_STAD_PLAT_ecosystem_projection.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nSignature coverage:\n")
print(coverage)
cat("\nEcology associations:\n")
print(ecology_associations)
cat("\nFour-group counts:\n")
print(ecotype_counts)
cat("\nSubtype tests:\n")
print(subtype_tests)
cat("\nTop ecosystem immune associations:\n")
print(endpoint_models[family == "immune_microenvironment" & exposure == "ecosystem_z"]
      [order(HC3_p)][1:12])
cat("\nContinuous Cox models:\n")
print(cox_models)
cat("\nKM summary:\n")
print(km_summary)
