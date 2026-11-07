suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
})

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
score_file <- file.path(
  project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
  "frozen_state_projection", "TCGA_STAD_PLAT_ecosystem_scores_clinical.tsv.gz"
)
cbio_file <- "I:/ECM-new/ECM-public/references/cbioportal_stad_tcga_pub_patient_clinical_2026-07-19.json"
out_dir <- file.path(
  project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
  "ecosystem_clinicopathology_v2"
)
source_dir <- file.path(project, "figures", "source_data", "Figure5_revised")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

if (!all(file.exists(c(score_file, cbio_file)))) stop("Missing clinicopathology input")
scores <- fread(score_file)
if (nrow(scores) != 412L || uniqueN(scores$case_submitter_id) != 412L) {
  stop("Unexpected frozen TCGA patient universe")
}

raw <- as.data.table(fromJSON(cbio_file))
raw <- raw[patientId %chin% scores$case_submitter_id,
           .(case_submitter_id = patientId, clinicalAttributeId, value)]
if (raw[, anyDuplicated(paste(case_submitter_id, clinicalAttributeId))] != 0L) {
  stop("Duplicated cBioPortal patient attribute")
}
wide <- dcast(raw, case_submitter_id ~ clinicalAttributeId, value.var = "value")
d <- merge(scores, wide, by = "case_submitter_id", all.x = TRUE, sort = FALSE)
d <- d[match(scores$case_submitter_id, case_submitter_id)]

blank_to_na <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "NA", "N/A", "[Not Available]", "[Not Applicable]", "Unknown")] <- NA_character_
  x
}
as_num <- function(x) suppressWarnings(as.numeric(blank_to_na(x)))

# The frozen table is authoritative for age, sex and overall pathological stage
# because it covers the full 412-tumour expression universe. cBioPortal adds
# detailed pathology for its curated TCGA publication subset.
d[, clinic_age := as.numeric(age)]
d[, clinic_sex := fifelse(sex == "MALE", "Male",
                          fifelse(sex == "FEMALE", "Female", NA_character_))]
d[, clinic_stage := factor(stage_group_frozen, levels = c("I", "II", "III", "IV"), ordered = TRUE)]

d[, path_t_raw := blank_to_na(PATH_T_STAGE)]
d[path_t_raw == "TX", path_t_raw := NA_character_]
d[, path_t_ordinal := fifelse(path_t_raw %chin% c("T1", "T1a", "T1b"), 1,
                              fifelse(path_t_raw == "T2", 2,
                                      fifelse(path_t_raw == "T3", 3,
                                              fifelse(path_t_raw %chin% c("T4", "T4a", "T4b"), 4, NA_real_))))]
d[, path_t_group := factor(path_t_ordinal, levels = 1:4, labels = paste0("T", 1:4), ordered = TRUE)]

d[, path_n_raw := blank_to_na(PATH_N_STAGE)]
d[path_n_raw == "NX", path_n_raw := NA_character_]
d[, path_n_ordinal := fifelse(path_n_raw == "N0", 0,
                              fifelse(path_n_raw == "N1", 1,
                                      fifelse(path_n_raw == "N2", 2,
                                              fifelse(path_n_raw == "N3", 3, NA_real_))))]
d[, path_n_group := factor(path_n_ordinal, levels = 0:3, labels = paste0("N", 0:3), ordered = TRUE)]

d[, path_m_raw := blank_to_na(PATH_M_STAGE)]
d[path_m_raw == "MX", path_m_raw := NA_character_]
d[, path_m_binary := fifelse(path_m_raw == "M1", 1,
                             fifelse(path_m_raw == "M0", 0, NA_real_))]
d[, path_m_group := factor(path_m_binary, levels = 0:1, labels = c("M0", "M1"))]

d[, lauren_class := factor(blank_to_na(LAUREN_CLASS),
                           levels = c("Intestinal", "Diffuse", "Mixed"))]
d[, anatomic_region := factor(blank_to_na(ANATOMIC_REGION),
                              levels = c("ANTRUM", "FUNDUS_BODY", "GEJ_CARDIA"),
                              labels = c("Antrum", "Fundus/body", "GEJ/cardia"))]
d[, who_class := factor(blank_to_na(WHO_CLASS),
                        levels = c("Tubular", "Poorly_Cohesive", "Papillary", "Mucinous", "Mixed"),
                        labels = c("Tubular", "Poorly cohesive", "Papillary", "Mucinous", "Mixed"))]
d[, intestinal_subclass := factor(blank_to_na(INTESTINAL_TYPE_SUBCLASS),
                                  levels = c("Tubular", "Papillary", "Mucinous"))]
d[, residual_tumour := factor(blank_to_na(RESIDUAL_TUMOR),
                              levels = c("R0", "R1", "R2", "RX"))]
d[, signet_ring := factor(blank_to_na(SIGNET_RING), levels = c("0", "1"),
                          labels = c("Absent", "Present"))]
d[, race_group := factor(blank_to_na(RACE),
                         levels = c("WHITE", "ASIAN", "BLACK_OR_AFRICAN_AMERICAN"),
                         labels = c("White", "Asian", "Black/African American"))]
d[, lymphocyte_infiltration_pct := as_num(PERCENT_LYMPHOCYTE_INFILTRATION)]
d[, tumour_cells_pct := as_num(PERCENT_TUMOR_CELLS)]
d[, tumour_nuclei_pct := as_num(PERCENT_TUMOR_NUCLEI)]
d[, ecological_configuration := factor(
  ecotype4,
  levels = c("Dual-low", "Fibroblast-dominant", "Macrophage-dominant", "Dual-high")
)]

# Verify the high-coverage fields agree where both sources are present.
sex_overlap <- d[!is.na(SEX) & !is.na(clinic_sex)]
sex_concordance <- mean(tolower(sex_overlap$SEX) == tolower(sex_overlap$clinic_sex))
age_overlap <- d[!is.na(AGE) & !is.na(clinic_age)]
age_abs_difference <- abs(as_num(age_overlap$AGE) - age_overlap$clinic_age)

variables <- data.table(
  variable_id = c(
    "age", "sex", "overall_stage", "path_t", "path_n", "path_m", "lauren",
    "anatomic_region", "who_class", "intestinal_subclass", "residual_tumour",
    "signet_ring", "race", "lymphocyte_infiltration", "tumour_cells", "tumour_nuclei"
  ),
  display_name = c(
    "Age at diagnosis, years", "Sex", "Pathological stage", "Pathological T stage",
    "Pathological N stage", "Pathological M stage", "Lauren class", "Anatomic region",
    "WHO histological class", "Intestinal-type subclass", "Residual tumour",
    "Signet-ring component", "Race", "Lymphocyte infiltration, %",
    "Tumour cells, %", "Tumour nuclei, %"
  ),
  column = c(
    "clinic_age", "clinic_sex", "clinic_stage", "path_t_group", "path_n_group",
    "path_m_group", "lauren_class", "anatomic_region", "who_class",
    "intestinal_subclass", "residual_tumour", "signet_ring", "race_group",
    "lymphocyte_infiltration_pct", "tumour_cells_pct", "tumour_nuclei_pct"
  ),
  type = c(
    "continuous", "categorical", rep("categorical", 11),
    "continuous", "continuous", "continuous"
  ),
  main_forest = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, rep(FALSE, 9))
)

fmt_cont <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  sprintf("%.1f [%.1f, %.1f]", median(x), quantile(x, .25), quantile(x, .75))
}
fmt_cat <- function(x, level) {
  den <- sum(!is.na(x))
  if (!den) return("NA")
  sprintf("%d (%.1f%%)", sum(x == level, na.rm = TRUE), 100 * mean(x == level, na.rm = TRUE))
}
kw_effect <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(g[ok])
  if (length(x) < 10L || nlevels(g) < 2L) return(NULL)
  z <- kruskal.test(x ~ g)
  k <- nlevels(g)
  eps <- max(0, (unname(z$statistic) - k + 1) / (length(x) - k))
  list(n = length(x), test = "Kruskal-Wallis", statistic = unname(z$statistic),
       df = unname(z$parameter), p = z$p.value, effect_name = "epsilon-squared",
       effect = eps, sparse = FALSE, min_expected = NA_real_)
}
cat_effect <- function(x, g, seed_offset = 0L) {
  ok <- !is.na(x) & !is.na(g)
  x <- droplevels(as.factor(x[ok])); g <- droplevels(g[ok])
  if (length(x) < 10L || nlevels(x) < 2L || nlevels(g) < 2L) return(NULL)
  tab <- table(g, x)
  asym <- suppressWarnings(chisq.test(tab, correct = FALSE))
  sparse <- any(asym$expected < 1) || mean(asym$expected < 5) > .20
  if (sparse) {
    set.seed(20260804 + seed_offset)
    z <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 99999))
    method <- "Monte Carlo Pearson chi-squared (99,999 replicates)"
  } else {
    z <- asym
    method <- "Pearson chi-squared"
  }
  v <- sqrt(unname(asym$statistic) /
              (sum(tab) * min(nrow(tab) - 1, ncol(tab) - 1)))
  list(n = sum(tab), test = method, statistic = unname(asym$statistic),
       df = unname(asym$parameter), p = z$p.value, effect_name = "Cramer's V",
       effect = v, sparse = sparse, min_expected = min(asym$expected))
}

group_levels <- levels(d$ecological_configuration)
table_rows <- list()
global_rows <- list()
for (i in seq_len(nrow(variables))) {
  spec <- variables[i]
  x <- d[[spec$column]]
  nonmissing_n <- sum(!is.na(x))
  missing_n <- nrow(d) - nonmissing_n
  missing_fraction <- missing_n / nrow(d)
  if (spec$type == "continuous") {
    x <- as.numeric(x)
    test <- kw_effect(x, d$ecological_configuration)
    row <- data.table(
      variable_id = spec$variable_id, characteristic = spec$display_name,
      level_or_summary = "Median [IQR]", overall = fmt_cont(x),
      `Dual-low` = fmt_cont(x[d$ecological_configuration == "Dual-low"]),
      `Fibroblast-dominant` = fmt_cont(x[d$ecological_configuration == "Fibroblast-dominant"]),
      `Macrophage-dominant` = fmt_cont(x[d$ecological_configuration == "Macrophage-dominant"]),
      `Dual-high` = fmt_cont(x[d$ecological_configuration == "Dual-high"]),
      nonmissing_n, missing_n, missing_fraction
    )
    table_rows[[length(table_rows) + 1L]] <- row
  } else {
    x <- droplevels(as.factor(x))
    test <- cat_effect(x, d$ecological_configuration, i)
    levs <- levels(x)
    for (lev in levs) {
      row <- data.table(
        variable_id = spec$variable_id, characteristic = spec$display_name,
        level_or_summary = lev, overall = fmt_cat(x, lev),
        `Dual-low` = fmt_cat(x[d$ecological_configuration == "Dual-low"], lev),
        `Fibroblast-dominant` = fmt_cat(x[d$ecological_configuration == "Fibroblast-dominant"], lev),
        `Macrophage-dominant` = fmt_cat(x[d$ecological_configuration == "Macrophage-dominant"], lev),
        `Dual-high` = fmt_cat(x[d$ecological_configuration == "Dual-high"], lev),
        nonmissing_n, missing_n, missing_fraction
      )
      table_rows[[length(table_rows) + 1L]] <- row
    }
  }
  if (!is.null(test)) {
    global_rows[[length(global_rows) + 1L]] <- data.table(
      variable_id = spec$variable_id, characteristic = spec$display_name,
      n = test$n, test = test$test, statistic = test$statistic, df = test$df,
      p_value = test$p, effect_size_name = test$effect_name,
      effect_size = test$effect, sparse_table = test$sparse,
      min_expected_count = test$min_expected, main_forest = spec$main_forest
    )
  }
}
table1 <- rbindlist(table_rows, fill = TRUE)
global_tests <- rbindlist(global_rows, fill = TRUE)
global_tests[, p_BH := p.adjust(p_value, method = "BH")]
table1 <- merge(table1, global_tests[, .(
  variable_id, association_test = test, global_p_value = p_value,
  global_p_BH = p_BH, effect_size_name, effect_size
)], by = "variable_id", all.x = TRUE, sort = FALSE)
table1 <- table1[order(match(variable_id, variables$variable_id))]
table1[, show_test_on_this_row := sequence(.N), by = variable_id]
table1[show_test_on_this_row > 1L,
       c("association_test", "global_p_value", "global_p_BH",
         "effect_size_name", "effect_size") := list(NA_character_, NA_real_, NA_real_,
                                                     NA_character_, NA_real_)]
table1[, show_test_on_this_row := NULL]

# HC3 coefficient and confidence interval for prespecified main-text features.
hc3_single <- function(y, x, label, estimand, reference, n_min = 30L) {
  ok <- is.finite(y) & is.finite(x)
  y <- y[ok]; x <- x[ok]
  if (length(y) < n_min || uniqueN(x) < 2L) return(NULL)
  X <- cbind(Intercept = 1, exposure = x)
  xtx_inv <- solve(crossprod(X))
  beta <- drop(xtx_inv %*% crossprod(X, y))
  residual <- y - drop(X %*% beta)
  hat <- rowSums((X %*% xtx_inv) * X)
  adjusted_residual <- residual / pmax(1 - hat, 1e-8)
  meat <- crossprod(X * adjusted_residual)
  vcov_hc3 <- xtx_inv %*% meat %*% xtx_inv
  se <- sqrt(diag(vcov_hc3))[2]
  stat <- beta[2] / se
  p <- 2 * pnorm(abs(stat), lower.tail = FALSE)
  data.table(
    feature = label, estimand, reference, n = length(y), beta = beta[2],
    se_HC3 = se, ci_low = beta[2] - 1.96 * se, ci_high = beta[2] + 1.96 * se,
    p_value = p
  )
}

forest <- rbindlist(list(
  hc3_single(d$ecosystem_z, d$clinic_age / 10,
             "Age", "per 10 years", "continuous"),
  hc3_single(d$ecosystem_z, as.numeric(d$clinic_sex == "Male"),
             "Sex", "male vs female", "female"),
  hc3_single(d$ecosystem_z, as.numeric(d$clinic_stage),
             "Pathological stage", "per I-to-IV level", "lower stage"),
  hc3_single(d$ecosystem_z, d$path_t_ordinal,
             "Pathological T stage", "per T1-to-T4 level", "lower T stage"),
  hc3_single(d$ecosystem_z, d$path_n_ordinal,
             "Pathological N stage", "per N0-to-N3 level", "lower N stage"),
  hc3_single(d$ecosystem_z, d$path_m_binary,
             "Pathological M stage", "M1 vs M0", "M0"),
  hc3_single(d$ecosystem_z,
             fifelse(d$lauren_class == "Diffuse", 1,
                     fifelse(d$lauren_class == "Intestinal", 0, NA_real_)),
             "Lauren class", "diffuse vs intestinal", "intestinal")
), fill = TRUE)
forest[, p_BH := p.adjust(p_value, method = "BH")]
forest[, display_label := paste0(feature, " (", estimand, ")")]
forest[, selection_rule := paste(
  "Prespecified for main text by clinical importance and interpretable coding;",
  "not selected using association P values"
)]

coverage <- variables[, .(
  variable_id, characteristic = display_name, source_column = column,
  variable_type = type, prespecified_main_forest = main_forest
)]
coverage[, nonmissing_n := vapply(source_column, function(nm) sum(!is.na(d[[nm]])), integer(1))]
coverage[, missing_n := nrow(d) - nonmissing_n]
coverage[, missing_fraction := missing_n / nrow(d)]
coverage[, main_figure_annotation_eligible :=
           characteristic %chin% c("Age at diagnosis, years", "Sex", "Pathological stage")]

patient_export <- d[, .(
  sample_id = sample_submitter_id, case_submitter_id, ecosystem_z,
  ecological_configuration = as.character(ecological_configuration),
  age = clinic_age, sex = clinic_sex, pathological_stage = clinic_stage,
  pathological_T_stage = path_t_group, pathological_N_stage = path_n_group,
  pathological_M_stage = path_m_group, Lauren_class = lauren_class,
  anatomic_region, WHO_histological_class = who_class, intestinal_subclass,
  residual_tumour, signet_ring_component = signet_ring, race = race_group,
  lymphocyte_infiltration_pct, tumour_cells_pct, tumour_nuclei_pct
)]

fwrite(table1, file.path(out_dir, "TableS_TCGA_clinicopathology_by_configuration.tsv"), sep = "\t")
fwrite(global_tests, file.path(out_dir, "TCGA_clinicopathology_global_tests.tsv"), sep = "\t")
fwrite(forest, file.path(out_dir, "TCGA_clinicopathology_main_forest.tsv"), sep = "\t")
fwrite(coverage, file.path(out_dir, "TCGA_clinicopathology_coverage_audit.tsv"), sep = "\t")
fwrite(patient_export, file.path(out_dir, "TCGA_clinicopathology_patient_level.tsv.gz"), sep = "\t")

fwrite(table1, file.path(source_dir, "TableS_TCGA_clinicopathology_by_configuration.tsv"), sep = "\t")
fwrite(global_tests, file.path(source_dir, "TableS_TCGA_clinicopathology_global_tests.tsv"), sep = "\t")
fwrite(forest, file.path(source_dir, "Figure5f_TCGA_clinicopathology_forest.tsv"), sep = "\t")
fwrite(coverage, file.path(source_dir, "TableS_TCGA_clinicopathology_coverage.tsv"), sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  tcga_expression_universe_n = nrow(d),
  cbioportal_matched_cases_n = uniqueN(raw$case_submitter_id),
  source_hierarchy = list(
    age_sex_overall_stage = "frozen 412-tumour TCGA table",
    detailed_pathology = "cBioPortal stad_tcga_pub patient-level clinical export dated 2026-07-19"
  ),
  source_concordance = list(
    sex_overlap_n = nrow(sex_overlap), sex_exact_concordance = sex_concordance,
    age_overlap_n = nrow(age_overlap),
    age_median_absolute_difference_years = median(age_abs_difference, na.rm = TRUE)
  ),
  supplementary_table = paste(
    "All prespecified clinicopathological variables are retained regardless of P value.",
    "Table cells are median [IQR] or n (% among nonmissing within column)."
  ),
  main_forest = paste(
    "Seven clinically central, interpretable variables were prespecified without outcome-based selection.",
    "Coefficients are univariable HC3 estimates in ecosystem-score SD units."
  ),
  multiplicity = list(
    table_global_tests = "BH across all clinicopathological characteristics",
    forest = "BH across seven prespecified forest estimates"
  ),
  limitations = c(
    "Detailed pathology is available only for the curated cBioPortal TCGA publication subset.",
    "No histological grade field is present in the audited source.",
    "Sparse categories remain visible in the supplementary table but are not promoted to the main forest."
  ),
  input_sha256 = list(
    score_file = digest(score_file, algo = "sha256", file = TRUE),
    cbioportal_file = digest(cbio_file, algo = "sha256", file = TRUE)
  )
)
write_json(audit, file.path(audit_dir, "Figure5_TCGA_clinicopathology_v2.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("Built clinicopathology table and main-text forest sources\n")
print(coverage)
print(forest)
