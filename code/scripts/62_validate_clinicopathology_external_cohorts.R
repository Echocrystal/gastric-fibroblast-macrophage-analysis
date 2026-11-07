suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
})

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
ext_dir <- file.path(
  project, "analysis", "results", "PLAT_axis", "external_validation",
  "frozen_state_projection"
)
tcga_file <- file.path(
  project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
  "ecosystem_clinicopathology_v2", "TCGA_clinicopathology_patient_level.tsv.gz"
)
out_dir <- file.path(
  project, "analysis", "results", "PLAT_axis", "external_validation",
  "clinicopathology_validation_v1"
)
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  TCGA_STAD = tcga_file,
  GSE62254 = file.path(ext_dir, "GSE62254_frozen_PLAT_ecosystem_scores.tsv.gz"),
  GSE84437 = file.path(ext_dir, "GSE84437_frozen_PLAT_ecosystem_scores.tsv.gz"),
  GSE15459 = file.path(ext_dir, "GSE15459_frozen_PLAT_ecosystem_scores.tsv.gz")
)
if (!all(file.exists(files))) stop("Missing frozen external clinicopathology input")

norm_sex <- function(x) {
  z <- tolower(trimws(as.character(x)))
  fifelse(z %chin% c("m", "male"), 1,
          fifelse(z %chin% c("f", "female"), 0, NA_real_))
}
norm_lauren <- function(x) {
  z <- tolower(trimws(as.character(x)))
  fifelse(z == "diffuse", 1, fifelse(z == "intestinal", 0, NA_real_))
}
norm_t <- function(x) {
  z <- toupper(trimws(as.character(x)))
  fifelse(z %chin% c("T1", "T1A", "T1B"), 1,
          fifelse(z == "T2", 2,
                  fifelse(z == "T3", 3,
                          fifelse(z %chin% c("T4", "T4A", "T4B"), 4, NA_real_))))
}
norm_n <- function(x) {
  z <- toupper(trimws(as.character(x)))
  fifelse(z == "N0", 0, fifelse(z == "N1", 1,
                                 fifelse(z == "N2", 2,
                                         fifelse(z == "N3", 3, NA_real_))))
}

tcga <- fread(files[["TCGA_STAD"]])
tcga_h <- tcga[, .(
  cohort = "TCGA-STAD", sample_id, ecosystem_z, ecological_configuration,
  age10 = as.numeric(age) / 10, male = norm_sex(sex),
  stage = match(as.character(pathological_stage), c("I", "II", "III", "IV")),
  t_stage = norm_t(pathological_T_stage), n_stage = norm_n(pathological_N_stage),
  lauren_diffuse = norm_lauren(Lauren_class)
)]

g62254 <- fread(files[["GSE62254"]])
g62254_h <- g62254[, .(
  cohort = "GSE62254", sample_id, ecosystem_z,
  ecological_configuration = ecotype4,
  age10 = as.numeric(age) / 10, male = norm_sex(sex),
  stage = as.numeric(stage_numeric), t_stage = NA_real_, n_stage = NA_real_,
  lauren_diffuse = norm_lauren(lauren)
)]

g84437 <- fread(files[["GSE84437"]])
g84437_h <- g84437[, .(
  cohort = "GSE84437", sample_id, ecosystem_z,
  ecological_configuration = ecotype4,
  age10 = as.numeric(age) / 10, male = norm_sex(sex),
  stage = NA_real_, t_stage = norm_t(pT), n_stage = norm_n(pN),
  lauren_diffuse = NA_real_
)]

g15459 <- fread(files[["GSE15459"]])
g15459_h <- g15459[, .(
  cohort = "GSE15459", sample_id, ecosystem_z,
  ecological_configuration = ecotype4,
  age10 = as.numeric(age) / 10, male = norm_sex(sex),
  stage = as.numeric(stage), t_stage = NA_real_, n_stage = NA_real_,
  lauren_diffuse = norm_lauren(lauren)
)]

d <- rbindlist(list(tcga_h, g62254_h, g84437_h, g15459_h), fill = TRUE)
d[, ecological_configuration := factor(
  ecological_configuration,
  levels = c("Dual-low", "Fibroblast-dominant", "Macrophage-dominant", "Dual-high")
)]
expected_n <- c(`TCGA-STAD` = 412L, GSE62254 = 300L, GSE84437 = 483L, GSE15459 = 192L)
observed_n <- d[, .N, by = cohort]
if (!identical(setNames(observed_n$N, observed_n$cohort)[names(expected_n)], expected_n)) {
  stop("Unexpected harmonized cohort size")
}

feature_labels <- c(
  age10 = "Age (per 10 years)", male = "Sex (male vs female)",
  stage = "Pathological stage (per level)",
  t_stage = "Pathological T stage (per level)",
  n_stage = "Pathological N stage (per level)",
  lauren_diffuse = "Lauren class (diffuse vs intestinal)"
)

coverage <- rbindlist(lapply(names(expected_n), function(cohort_name) {
  z <- d[cohort == cohort_name]
  rbindlist(lapply(names(feature_labels), function(feature) {
    x <- z[[feature]]
    vals <- x[is.finite(x)]
    data.table(
      cohort = cohort_name, feature, label = feature_labels[[feature]],
      cohort_n = nrow(z), nonmissing_n = length(vals),
      missing_n = nrow(z) - length(vals), missing_fraction = 1 - length(vals) / nrow(z),
      unique_values = uniqueN(vals), minimum = if (length(vals)) min(vals) else NA_real_,
      maximum = if (length(vals)) max(vals) else NA_real_
    )
  }))
}))

hc3_fit <- function(z, exposure, covariates, tier) {
  vars <- c("ecosystem_z", exposure, covariates)
  q <- z[complete.cases(z[, ..vars])]
  if (nrow(q) < 30L || uniqueN(q[[exposure]]) < 2L) return(NULL)
  X <- cbind(Intercept = 1, as.matrix(q[, c(exposure, covariates), with = FALSE]))
  storage.mode(X) <- "double"
  y <- q$ecosystem_z
  xtx <- crossprod(X)
  if (qr(xtx)$rank < ncol(X)) return(NULL)
  xtx_inv <- solve(xtx)
  beta <- drop(xtx_inv %*% crossprod(X, y))
  resid <- y - drop(X %*% beta)
  hat <- rowSums((X %*% xtx_inv) * X)
  u <- resid / pmax(1 - hat, 1e-8)
  meat <- crossprod(X * u)
  vc <- xtx_inv %*% meat %*% xtx_inv
  se <- sqrt(diag(vc))[match(exposure, colnames(X))]
  estimate <- beta[match(exposure, colnames(X))]
  stat <- estimate / se
  data.table(
    cohort = q$cohort[1], feature = exposure, label = feature_labels[[exposure]],
    model_tier = tier, covariates = if (length(covariates)) paste(covariates, collapse = "+") else "none",
    n = nrow(q), beta = estimate, se_HC3 = se,
    ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se,
    p_value = 2 * pnorm(abs(stat), lower.tail = FALSE)
  )
}

models <- list()
for (cohort_name in names(expected_n)) {
  z <- d[cohort == cohort_name]
  available <- names(feature_labels)[vapply(names(feature_labels), function(x) {
    sum(is.finite(z[[x]])) >= 30L && uniqueN(z[[x]][is.finite(z[[x]])]) >= 2L
  }, logical(1))]
  for (feature in available) {
    models[[length(models) + 1L]] <- hc3_fit(z, feature, character(), "univariable")
    basic_cov <- if (feature == "age10") "male" else if (feature == "male") "age10" else c("age10", "male")
    if (all(basic_cov %in% available)) {
      models[[length(models) + 1L]] <- hc3_fit(z, feature, basic_cov, "age_sex_adjusted")
    }
    full_cov <- NULL
    if (feature == "stage" && "lauren_diffuse" %in% available) {
      full_cov <- c("age10", "male", "lauren_diffuse")
    } else if (feature == "lauren_diffuse" && "stage" %in% available) {
      full_cov <- c("age10", "male", "stage")
    } else if (feature == "t_stage" && "n_stage" %in% available) {
      full_cov <- c("age10", "male", "n_stage")
    } else if (feature == "n_stage" && "t_stage" %in% available) {
      full_cov <- c("age10", "male", "t_stage")
    }
    if (!is.null(full_cov) && all(full_cov %in% available)) {
      models[[length(models) + 1L]] <- hc3_fit(z, feature, full_cov, "clinicopath_adjusted")
    }
  }
}
model_results <- rbindlist(models, fill = TRUE)
model_results[, p_BH_within_tier := p.adjust(p_value, method = "BH"), by = model_tier]

meta_dl <- function(q, analysis_set) {
  q <- q[is.finite(beta) & is.finite(se_HC3) & se_HC3 > 0]
  if (nrow(q) < 2L) return(NULL)
  w <- 1 / q$se_HC3^2
  fixed_beta <- sum(w * q$beta) / sum(w)
  Q <- sum(w * (q$beta - fixed_beta)^2)
  df_q <- nrow(q) - 1L
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - df_q) / C)
  wr <- 1 / (q$se_HC3^2 + tau2)
  pooled <- sum(wr * q$beta) / sum(wr)
  pooled_se <- sqrt(1 / sum(wr))
  z <- pooled / pooled_se
  data.table(
    feature = q$feature[1], label = q$label[1], model_tier = q$model_tier[1],
    analysis_set, k = nrow(q), cohorts = paste(q$cohort, collapse = ";"),
    total_n = sum(q$n), beta_random_DL = pooled, se_random_DL = pooled_se,
    ci_low = pooled - 1.96 * pooled_se, ci_high = pooled + 1.96 * pooled_se,
    p_value = 2 * pnorm(abs(z), lower.tail = FALSE), tau2_DL = tau2,
    Q = Q, Q_df = df_q, Q_p = pchisq(Q, df = df_q, lower.tail = FALSE),
    I2_percent = if (Q > 0) max(0, (Q - df_q) / Q) * 100 else 0,
    positive_cohorts = sum(q$beta > 0), negative_cohorts = sum(q$beta < 0)
  )
}

meta_rows <- list()
for (tier in unique(model_results$model_tier)) {
  for (feat in unique(model_results$feature)) {
    q <- model_results[model_tier == tier & feature == feat]
    ext <- q[cohort != "TCGA-STAD"]
    if (nrow(ext) >= 2L) {
      meta_rows[[length(meta_rows) + 1L]] <- meta_dl(ext, "external_only")
    }
    if (nrow(q) >= 2L) {
      meta_rows[[length(meta_rows) + 1L]] <- meta_dl(q, "all_available_including_TCGA")
    }
  }
}
meta_results <- rbindlist(meta_rows, fill = TRUE)
meta_results[, p_BH_within_analysis_tier := p.adjust(p_value, method = "BH"),
             by = .(analysis_set, model_tier)]

# Configuration-level tests remain descriptive and secondary.
kw_config <- function(z, feature) {
  q <- z[is.finite(get(feature)) & !is.na(ecological_configuration)]
  if (nrow(q) < 30L) return(NULL)
  test <- kruskal.test(q[[feature]], q$ecological_configuration)
  k <- nlevels(droplevels(q$ecological_configuration))
  eps <- max(0, (unname(test$statistic) - k + 1) / (nrow(q) - k))
  data.table(n = nrow(q), test = "Kruskal-Wallis", statistic = unname(test$statistic),
             df = unname(test$parameter), p_value = test$p.value,
             effect_size_name = "epsilon-squared", effect_size = eps,
             sparse_table = FALSE, min_expected = NA_real_)
}
cat_config <- function(z, feature, seed) {
  q <- z[is.finite(get(feature)) & !is.na(ecological_configuration)]
  if (nrow(q) < 30L) return(NULL)
  tab <- table(q$ecological_configuration, q[[feature]])
  if (ncol(tab) < 2L) return(NULL)
  asym <- suppressWarnings(chisq.test(tab, correct = FALSE))
  sparse <- any(asym$expected < 1) || mean(asym$expected < 5) > .20
  if (sparse) {
    set.seed(seed)
    test <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 99999))
    method <- "Monte Carlo Pearson chi-squared (99,999 replicates)"
  } else {
    test <- asym
    method <- "Pearson chi-squared"
  }
  v <- sqrt(unname(asym$statistic) /
              (sum(tab) * min(nrow(tab) - 1, ncol(tab) - 1)))
  data.table(n = sum(tab), test = method, statistic = unname(asym$statistic),
             df = unname(asym$parameter), p_value = test$p.value,
             effect_size_name = "Cramer's V", effect_size = v,
             sparse_table = sparse, min_expected = min(asym$expected))
}
config_rows <- list()
seed <- 20260805L
for (cohort_name in names(expected_n)) {
  z <- d[cohort == cohort_name]
  available <- names(feature_labels)[vapply(names(feature_labels), function(x) {
    sum(is.finite(z[[x]])) >= 30L && uniqueN(z[[x]][is.finite(z[[x]])]) >= 2L
  }, logical(1))]
  for (feature in available) {
    seed <- seed + 1L
    test <- if (feature == "age10") kw_config(z, feature) else cat_config(z, feature, seed)
    if (!is.null(test)) {
      config_rows[[length(config_rows) + 1L]] <- cbind(
        data.table(cohort = cohort_name, feature, label = feature_labels[[feature]]), test
      )
    }
  }
}
configuration_tests <- rbindlist(config_rows, fill = TRUE)
configuration_tests[, p_BH_all_tests := p.adjust(p_value, method = "BH")]

configuration_counts <- d[, .N, by = .(cohort, ecological_configuration)]
setorder(configuration_counts, cohort, ecological_configuration)

fwrite(d, file.path(out_dir, "harmonized_four_cohort_clinicopathology.tsv.gz"), sep = "\t")
fwrite(coverage, file.path(out_dir, "four_cohort_clinicopathology_coverage.tsv"), sep = "\t")
fwrite(model_results, file.path(out_dir, "four_cohort_clinicopathology_HC3_models.tsv"), sep = "\t")
fwrite(meta_results, file.path(out_dir, "four_cohort_clinicopathology_random_effects_meta.tsv"), sep = "\t")
fwrite(configuration_tests, file.path(out_dir, "four_cohort_configuration_clinicopathology_tests.tsv"), sep = "\t")
fwrite(configuration_counts, file.path(out_dir, "four_cohort_configuration_counts.tsv"), sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  objective = "External validation of clinicopathological associations for the frozen continuous PLAT ecosystem score",
  cohorts = as.list(expected_n),
  primary_exposure = "within-cohort standardized continuous ecosystem_z",
  model_tiers = list(
    univariable = "simple clinicopathological correspondence",
    age_sex_adjusted = "mutual age/sex or age+sex adjustment",
    clinicopath_adjusted = "stage and Lauren mutually adjusted; T and N mutually adjusted"
  ),
  meta_analysis = "DerSimonian-Laird random effects on HC3 coefficients; external-only and all-available summaries are separated",
  configuration_boundary = "median-defined descriptive configurations; secondary only",
  multiplicity = "BH within each model tier for cohort-specific models and within analysis-set/model-tier for meta summaries",
  source_files_sha256 = as.list(setNames(vapply(files, digest, character(1),
                                                 algo = "sha256", file = TRUE), names(files)))
)
write_json(audit, file.path(audit_dir, "external_clinicopathology_validation_v1.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("External clinicopathology validation complete\n")
print(coverage)
print(model_results[model_tier == "univariable"])
print(meta_results[model_tier == "univariable"])
