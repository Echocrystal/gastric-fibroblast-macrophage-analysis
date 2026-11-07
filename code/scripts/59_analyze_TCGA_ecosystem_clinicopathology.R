suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
})

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
input_file <- file.path(
  project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
  "frozen_state_projection", "TCGA_STAD_PLAT_ecosystem_scores_clinical.tsv.gz"
)
out_dir <- file.path(
  project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
  "ecosystem_clinicopathology_v1"
)
source_dir <- file.path(project, "figures", "source_data", "Figure5_revised")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Missing frozen TCGA ecosystem table")
d <- fread(input_file)
required <- c(
  "sample_submitter_id", "ecosystem_z", "ecotype4", "age", "sex",
  "stage_group_frozen", "stage_ordinal"
)
if (!all(required %in% names(d))) stop("Missing required clinicopathological fields")
if (nrow(d) != 412L) stop("Unexpected TCGA cohort size")

d[, sex := fifelse(sex %chin% c("MALE", "FEMALE"), sex, NA_character_)]
d[, stage_group_frozen := fifelse(
  stage_group_frozen %chin% c("I", "II", "III", "IV"),
  stage_group_frozen, NA_character_
)]
d[, stage_ordinal := match(stage_group_frozen, c("I", "II", "III", "IV"))]
d[, ecotype4 := factor(
  ecotype4,
  levels = c("Dual-low", "Fibroblast-dominant", "Macrophage-dominant", "Dual-high")
)]

# Coverage is prespecified and is not filtered using association P values.
coverage <- rbindlist(list(
  data.table(
    variable = "Age", nonmissing_n = sum(!is.na(d$age)),
    missing_n = sum(is.na(d$age)), missing_fraction = mean(is.na(d$age)),
    n_levels = uniqueN(d$age, na.rm = TRUE),
    smallest_display_group_n = NA_integer_, display_role = "continuous annotation"
  ),
  data.table(
    variable = "Sex", nonmissing_n = sum(!is.na(d$sex)),
    missing_n = sum(is.na(d$sex)), missing_fraction = mean(is.na(d$sex)),
    n_levels = uniqueN(d$sex, na.rm = TRUE),
    smallest_display_group_n = min(table(d$sex)), display_role = "categorical annotation"
  ),
  data.table(
    variable = "Pathological stage", nonmissing_n = sum(!is.na(d$stage_group_frozen)),
    missing_n = sum(is.na(d$stage_group_frozen)),
    missing_fraction = mean(is.na(d$stage_group_frozen)),
    n_levels = uniqueN(d$stage_group_frozen, na.rm = TRUE),
    smallest_display_group_n = min(table(d$stage_group_frozen)),
    display_role = "ordinal categorical annotation"
  )
), fill = TRUE)
coverage[, passes_main_figure_gate :=
           nonmissing_n >= 250L & missing_fraction <= 0.35 &
           (is.na(smallest_display_group_n) | smallest_display_group_n >= 20L)]

# Primary continuous-score associations.
age_complete <- d[!is.na(age) & !is.na(ecosystem_z)]
age_test <- cor.test(age_complete$age, age_complete$ecosystem_z,
                     method = "spearman", exact = FALSE)

sex_complete <- d[!is.na(sex) & !is.na(ecosystem_z)]
sex_test <- wilcox.test(ecosystem_z ~ sex, data = sex_complete,
                        exact = FALSE, conf.int = FALSE)
n_female <- sum(sex_complete$sex == "FEMALE")
n_male <- sum(sex_complete$sex == "MALE")
# For the two-sample interface, R reports the Mann-Whitney U statistic for the
# first lexicographic group (FEMALE). Convert it to rank-biserial correlation,
# positive when FEMALE values exceed MALE values.
u_female <- unname(sex_test$statistic)
rank_biserial <- 2 * u_female / (n_female * n_male) - 1

stage_complete <- d[!is.na(stage_group_frozen) & !is.na(ecosystem_z)]
stage_kw <- kruskal.test(ecosystem_z ~ factor(stage_group_frozen,
                                              levels = c("I", "II", "III", "IV")),
                         data = stage_complete)
k_stage <- uniqueN(stage_complete$stage_group_frozen)
epsilon_sq <- max(0, (unname(stage_kw$statistic) - k_stage + 1) /
                        (nrow(stage_complete) - k_stage))
stage_trend <- cor.test(stage_complete$stage_ordinal, stage_complete$ecosystem_z,
                        method = "spearman", exact = FALSE)

continuous_results <- rbindlist(list(
  data.table(
    endpoint = "Age", test = "Spearman rank correlation",
    estimand = "rho per increasing year", estimate = unname(age_test$estimate),
    n = nrow(age_complete), statistic = unname(age_test$statistic),
    df = NA_real_, p_value = age_test$p.value,
    group_summary = sprintf("median age %.1f y (IQR %.1f-%.1f)",
                            median(age_complete$age), quantile(age_complete$age, .25),
                            quantile(age_complete$age, .75))
  ),
  data.table(
    endpoint = "Sex", test = "Wilcoxon rank-sum",
    estimand = "rank-biserial r (female vs male)", estimate = rank_biserial,
    n = nrow(sex_complete), statistic = unname(sex_test$statistic),
    df = NA_real_, p_value = sex_test$p.value,
    group_summary = sprintf("female n=%d, median %.3f; male n=%d, median %.3f",
                            n_female, median(sex_complete[sex == "FEMALE", ecosystem_z]),
                            n_male, median(sex_complete[sex == "MALE", ecosystem_z]))
  ),
  data.table(
    endpoint = "Pathological stage (global)", test = "Kruskal-Wallis",
    estimand = "epsilon-squared", estimate = epsilon_sq,
    n = nrow(stage_complete), statistic = unname(stage_kw$statistic),
    df = unname(stage_kw$parameter), p_value = stage_kw$p.value,
    group_summary = paste(
      stage_complete[, .(n = .N, median = median(ecosystem_z)),
                     by = stage_group_frozen][
                       match(c("I", "II", "III", "IV"), stage_group_frozen),
                       sprintf("%s n=%d, median %.3f", stage_group_frozen, n, median)
                     ],
      collapse = "; "
    )
  ),
  data.table(
    endpoint = "Pathological stage (trend)", test = "Spearman ordinal trend",
    estimand = "rho per increasing stage", estimate = unname(stage_trend$estimate),
    n = nrow(stage_complete), statistic = unname(stage_trend$statistic),
    df = NA_real_, p_value = stage_trend$p.value,
    group_summary = "stage I-IV coded 1-4"
  )
), fill = TRUE)
continuous_results[, p_BH := p.adjust(p_value, method = "BH")]

# Secondary descriptive-configuration calibration. Sparse RxC tables use a
# fixed-seed Monte Carlo chi-squared P value; effect size is Cramer's V.
categorical_test <- function(variable, label) {
  z <- d[!is.na(ecotype4) & !is.na(get(variable))]
  tab <- table(z$ecotype4, z[[variable]])
  asym <- suppressWarnings(chisq.test(tab, correct = FALSE))
  sparse <- any(asym$expected < 1) || mean(asym$expected < 5) > 0.20
  if (sparse) {
    set.seed(20260804)
    test <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 99999))
    method <- "Monte Carlo Pearson chi-squared (99,999 replicates)"
  } else {
    test <- asym
    method <- "Pearson chi-squared"
  }
  n <- sum(tab)
  v <- sqrt(unname(asym$statistic) / (n * min(nrow(tab) - 1, ncol(tab) - 1)))
  list(
    result = data.table(
      endpoint = label, n = n, rows = nrow(tab), columns = ncol(tab),
      min_expected = min(asym$expected), fraction_expected_below_5 = mean(asym$expected < 5),
      sparse_table = sparse, test = method,
      statistic = unname(asym$statistic), df = unname(asym$parameter),
      p_value = test$p.value, cramers_v = v
    ),
    counts = as.data.table(tab)[, endpoint := label][]
  )
}
sex_cat <- categorical_test("sex", "Sex")
stage_cat <- categorical_test("stage_group_frozen", "Pathological stage")
categorical_results <- rbindlist(list(sex_cat$result, stage_cat$result), fill = TRUE)
categorical_results[, p_BH := p.adjust(p_value, method = "BH")]
contingency <- rbindlist(list(sex_cat$counts, stage_cat$counts), fill = TRUE,
                         use.names = TRUE)

patient_data <- d[, .(
  sample_id = sample_submitter_id, ecosystem_z,
  ecological_configuration = as.character(ecotype4), age, sex,
  pathological_stage = stage_group_frozen, stage_ordinal
)]

fwrite(coverage, file.path(out_dir, "TCGA_ecosystem_clinicopathology_missingness.tsv"), sep = "\t")
fwrite(continuous_results,
       file.path(out_dir, "TCGA_ecosystem_clinicopathology_continuous_tests.tsv"), sep = "\t")
fwrite(categorical_results,
       file.path(out_dir, "TCGA_ecotype4_clinicopathology_categorical_tests.tsv"), sep = "\t")
fwrite(contingency,
       file.path(out_dir, "TCGA_ecotype4_clinicopathology_contingency.tsv"), sep = "\t")
fwrite(patient_data,
       file.path(out_dir, "TCGA_ecosystem_clinicopathology_patient_data.tsv.gz"), sep = "\t")

# Copy the inferential tables into the Figure 5 source-data package.
fwrite(coverage, file.path(source_dir, "Figure5b_clinicopathology_coverage.tsv"), sep = "\t")
fwrite(continuous_results,
       file.path(source_dir, "Figure5b_continuous_clinicopathology_tests.tsv"), sep = "\t")
fwrite(categorical_results,
       file.path(source_dir, "Figure5b_configuration_clinicopathology_tests.tsv"), sep = "\t")
fwrite(contingency,
       file.path(source_dir, "Figure5b_configuration_clinicopathology_contingency.tsv"), sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  cohort = "TCGA-STAD frozen 412-tumour projection",
  main_inference = "continuous ecosystem score; clinicopathological endpoints were not used to define signatures or configurations",
  multiplicity = "Benjamini-Hochberg within the four continuous-score tests and within the two configuration tests",
  display_gate = "nonmissing n >= 250, missing fraction <= 0.35, and each categorical level n >= 20",
  displayed_variables = coverage[passes_main_figure_gate == TRUE, variable],
  input_file = input_file,
  input_sha256 = digest(input_file, algo = "sha256", file = TRUE),
  note = paste(
    "The frozen table does not contain reliable Lauren subtype, histological grade or T/N/M fields.",
    "These variables were not imputed or borrowed from external cohorts."
  )
)
write_json(audit, file.path(audit_dir, "Figure5_TCGA_clinicopathology_v1.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("Clinicopathology analysis complete\n")
print(coverage)
print(continuous_results)
print(categorical_results)
