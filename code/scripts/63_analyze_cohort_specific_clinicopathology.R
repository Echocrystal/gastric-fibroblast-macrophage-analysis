suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
})

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
tcga_patient_file <- file.path(
  project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
  "ecosystem_clinicopathology_v2", "TCGA_clinicopathology_patient_level.tsv.gz"
)
tcga_score_file <- file.path(
  project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
  "frozen_state_projection", "TCGA_STAD_PLAT_ecosystem_scores_clinical.tsv.gz"
)
ext_dir <- file.path(
  project, "analysis", "results", "PLAT_axis", "external_validation",
  "frozen_state_projection"
)
g62254_file <- file.path(ext_dir, "GSE62254_frozen_PLAT_ecosystem_scores.tsv.gz")
g84437_file <- file.path(ext_dir, "GSE84437_frozen_PLAT_ecosystem_scores.tsv.gz")
g15459_file <- file.path(ext_dir, "GSE15459_frozen_PLAT_ecosystem_scores.tsv.gz")
files <- c(tcga_patient_file, tcga_score_file, g62254_file, g84437_file, g15459_file)
if (!all(file.exists(files))) stop("Missing cohort-specific clinical input")

out_dir <- file.path(
  project, "analysis", "results", "PLAT_axis", "external_validation",
  "cohort_specific_clinicopathology_v1"
)
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

clean_char <- function(x) {
  z <- trimws(as.character(x))
  z[z %chin% c("", "NA", "N/A", "Unknown", "Indeterminate")] <- NA_character_
  z
}

long_parts <- list()
add_categorical <- function(cohort, sample_id, ecosystem_z, config, variable_id,
                            label, values, source_field, family = "clinicopathology") {
  data.table(
    cohort, sample_id, ecosystem_z, ecological_configuration = config,
    variable_id, variable_label = label, variable_type = "categorical",
    family, source_field, category = clean_char(values), value = NA_real_
  )
}
add_continuous <- function(cohort, sample_id, ecosystem_z, config, variable_id,
                           label, values, source_field, family = "clinicopathology") {
  data.table(
    cohort, sample_id, ecosystem_z, ecological_configuration = config,
    variable_id, variable_label = label, variable_type = "continuous",
    family, source_field, category = NA_character_,
    value = suppressWarnings(as.numeric(values))
  )
}

# TCGA-specific clinical, histological and molecular context.
tcga_p <- fread(tcga_patient_file)
tcga_s <- fread(tcga_score_file)[, .(
  sample_id = sample_submitter_id, tcga2014_subtype, MSI_STATUS, EBV_PRESENT
)]
tcga <- merge(tcga_p, tcga_s, by = "sample_id", all.x = TRUE, sort = FALSE)
tcga_cfg <- tcga$ecological_configuration
long_parts <- c(long_parts, list(
  add_categorical("TCGA-STAD", tcga$sample_id, tcga$ecosystem_z, tcga_cfg,
                  "tcga2014_subtype", "TCGA 2014 molecular subtype",
                  tcga$tcga2014_subtype, "tcga2014_subtype", "source_molecular_subtype"),
  add_categorical("TCGA-STAD", tcga$sample_id, tcga$ecosystem_z, tcga_cfg,
                  "who_histology", "WHO histological class",
                  tcga$WHO_histological_class, "WHO_CLASS"),
  add_categorical("TCGA-STAD", tcga$sample_id, tcga$ecosystem_z, tcga_cfg,
                  "anatomic_region", "Anatomic region",
                  tcga$anatomic_region, "ANATOMIC_REGION"),
  add_categorical("TCGA-STAD", tcga$sample_id, tcga$ecosystem_z, tcga_cfg,
                  "intestinal_subclass", "Intestinal-type subclass",
                  tcga$intestinal_subclass, "INTESTINAL_TYPE_SUBCLASS"),
  add_categorical("TCGA-STAD", tcga$sample_id, tcga$ecosystem_z, tcga_cfg,
                  "residual_tumour", "Residual tumour",
                  tcga$residual_tumour, "RESIDUAL_TUMOR"),
  add_categorical("TCGA-STAD", tcga$sample_id, tcga$ecosystem_z, tcga_cfg,
                  "signet_ring", "Signet-ring component",
                  tcga$signet_ring_component, "SIGNET_RING"),
  add_categorical("TCGA-STAD", tcga$sample_id, tcga$ecosystem_z, tcga_cfg,
                  "race", "Race", tcga$race, "RACE", "demographic_context"),
  add_continuous("TCGA-STAD", tcga$sample_id, tcga$ecosystem_z, tcga_cfg,
                 "lymphocyte_infiltration", "Lymphocyte infiltration (%)",
                 tcga$lymphocyte_infiltration_pct, "PERCENT_LYMPHOCYTE_INFILTRATION",
                 "specimen_pathology"),
  add_continuous("TCGA-STAD", tcga$sample_id, tcga$ecosystem_z, tcga_cfg,
                 "tumour_cells", "Tumour cells (%)", tcga$tumour_cells_pct,
                 "PERCENT_TUMOR_CELLS", "specimen_pathology"),
  add_continuous("TCGA-STAD", tcga$sample_id, tcga$ecosystem_z, tcga_cfg,
                 "tumour_nuclei", "Tumour nuclei (%)", tcga$tumour_nuclei_pct,
                 "PERCENT_TUMOR_NUCLEI", "specimen_pathology")
))

# GSE62254/ACRG-specific clinical and molecular context.
g62254 <- fread(g62254_file)
g62254[, mlh1_status := fifelse(
  grepl("negative|partial loss", mlh1_ihc, ignore.case = TRUE), "Loss",
  fifelse(grepl("positive", mlh1_ihc, ignore.case = TRUE), "Retained", NA_character_)
)]
g62254[, ebv_status := fifelse(ebv_ish == 1, "Positive",
                               fifelse(ebv_ish == 0, "Negative", NA_character_))]
long_parts <- c(long_parts, list(
  add_categorical("GSE62254", g62254$sample_id, g62254$ecosystem_z, g62254$ecotype4,
                  "acrg_subtype", "ACRG molecular subtype", g62254$acrg_subtype,
                  "acrg_subtype", "source_molecular_subtype"),
  add_categorical("GSE62254", g62254$sample_id, g62254$ecosystem_z, g62254$ecotype4,
                  "ebv_ish", "EBV-ISH status", g62254$ebv_status,
                  "ebv_ish", "molecular_pathology"),
  add_categorical("GSE62254", g62254$sample_id, g62254$ecosystem_z, g62254$ecotype4,
                  "mlh1_ihc", "MLH1-IHC status", g62254$mlh1_status,
                  "mlh1_ihc", "molecular_pathology")
))

# GSE15459-specific source molecular subtype.
g15459 <- fread(g15459_file)
long_parts <- c(long_parts, list(
  add_categorical("GSE15459", g15459$sample_id, g15459$ecosystem_z, g15459$ecotype4,
                  "lei_molecular_subtype", "Lei molecular subtype",
                  g15459$molecular_subtype, "molecular_subtype",
                  "source_molecular_subtype")
))

# GSE84437 has no additional audited clinicopathological field beyond age, sex,
# pT and pN in the frozen projection table. This absence is retained in coverage.
g84437 <- fread(g84437_file)

long <- rbindlist(long_parts, fill = TRUE)
long[, ecological_configuration := factor(
  ecological_configuration,
  levels = c("Dual-low", "Fibroblast-dominant", "Macrophage-dominant", "Dual-high")
)]

cohort_sizes <- c(`TCGA-STAD` = 412L, GSE62254 = 300L, GSE84437 = 483L, GSE15459 = 192L)
coverage <- unique(long[, .(
  cohort, variable_id, variable_label, variable_type, family, source_field
)])
coverage[, cohort_n := cohort_sizes[cohort]]
coverage[, nonmissing_n := vapply(seq_len(.N), function(i) {
  q <- long[cohort == coverage$cohort[i] & variable_id == coverage$variable_id[i]]
  if (coverage$variable_type[i] == "categorical") sum(!is.na(q$category)) else sum(is.finite(q$value))
}, integer(1))]
coverage[, missing_n := cohort_n - nonmissing_n]
coverage[, missing_fraction := missing_n / cohort_n]
coverage[, usable_levels := vapply(seq_len(.N), function(i) {
  q <- long[cohort == coverage$cohort[i] & variable_id == coverage$variable_id[i]]
  if (coverage$variable_type[i] == "categorical") uniqueN(q$category, na.rm = TRUE) else uniqueN(q$value[is.finite(q$value)])
}, integer(1))]

global_rows <- list()
category_rows <- list()
pair_rows <- list()
configuration_rows <- list()
seed <- 20260805L

for (i in seq_len(nrow(coverage))) {
  spec <- coverage[i]
  q <- long[cohort == spec$cohort & variable_id == spec$variable_id]
  if (spec$variable_type == "continuous") {
    z <- q[is.finite(value) & is.finite(ecosystem_z)]
    if (nrow(z) >= 30L && uniqueN(z$value) >= 3L) {
      test <- cor.test(z$value, z$ecosystem_z, method = "spearman", exact = FALSE)
      global_rows[[length(global_rows) + 1L]] <- data.table(
        cohort = spec$cohort, variable_id = spec$variable_id,
        variable_label = spec$variable_label, variable_type = spec$variable_type,
        family = spec$family, n = nrow(z), test = "Spearman rank correlation",
        effect_size_name = "rho", effect_size = unname(test$estimate),
        statistic = unname(test$statistic), df = NA_real_, p_value = test$p.value
      )
      cfg_test <- kruskal.test(value ~ ecological_configuration,
                               data = z[!is.na(ecological_configuration)])
      k <- nlevels(droplevels(z[!is.na(ecological_configuration)]$ecological_configuration))
      eps <- max(0, (unname(cfg_test$statistic) - k + 1) /
                   (sum(!is.na(z$ecological_configuration)) - k))
      configuration_rows[[length(configuration_rows) + 1L]] <- data.table(
        cohort = spec$cohort, variable_id = spec$variable_id,
        variable_label = spec$variable_label, n = sum(!is.na(z$ecological_configuration)),
        test = "Kruskal-Wallis", effect_size_name = "epsilon-squared",
        effect_size = eps, p_value = cfg_test$p.value,
        sparse_table = FALSE, min_expected = NA_real_
      )
    }
  } else {
    z <- q[!is.na(category) & is.finite(ecosystem_z)]
    z[, category := factor(category)]
    if (nrow(z) >= 30L && nlevels(z$category) >= 2L) {
      test <- kruskal.test(ecosystem_z ~ category, data = z)
      k <- nlevels(z$category)
      eps <- max(0, (unname(test$statistic) - k + 1) / (nrow(z) - k))
      global_rows[[length(global_rows) + 1L]] <- data.table(
        cohort = spec$cohort, variable_id = spec$variable_id,
        variable_label = spec$variable_label, variable_type = spec$variable_type,
        family = spec$family, n = nrow(z), test = "Kruskal-Wallis",
        effect_size_name = "epsilon-squared", effect_size = eps,
        statistic = unname(test$statistic), df = unname(test$parameter), p_value = test$p.value
      )
      category_rows[[length(category_rows) + 1L]] <- z[, .(
        n = .N, median_ecosystem_z = median(ecosystem_z),
        q1 = quantile(ecosystem_z, .25), q3 = quantile(ecosystem_z, .75),
        mean_ecosystem_z = mean(ecosystem_z), sd_ecosystem_z = sd(ecosystem_z)
      ), by = .(cohort, variable_id, variable_label, category)]
      levs <- levels(z$category)
      if (length(levs) >= 2L) {
        pairs <- combn(levs, 2, simplify = FALSE)
        pair_results <- rbindlist(lapply(pairs, function(pair) {
          xa <- z[category == pair[1], ecosystem_z]
          xb <- z[category == pair[2], ecosystem_z]
          wt <- wilcox.test(xb, xa, exact = FALSE)
          rb <- 2 * unname(wt$statistic) / (length(xa) * length(xb)) - 1
          data.table(
            cohort = spec$cohort, variable_id = spec$variable_id,
            variable_label = spec$variable_label,
            group_a = pair[1], group_b = pair[2], n_a = length(xa), n_b = length(xb),
            median_a = median(xa), median_b = median(xb),
            median_difference_b_minus_a = median(xb) - median(xa),
            rank_biserial_b_vs_a = rb, statistic = unname(wt$statistic),
            p_value = wt$p.value
          )
        }))
        pair_results[, p_BH_within_variable := p.adjust(p_value, method = "BH")]
        pair_rows[[length(pair_rows) + 1L]] <- pair_results
      }

      tab <- table(z$ecological_configuration, z$category)
      asym <- suppressWarnings(chisq.test(tab, correct = FALSE))
      sparse <- any(asym$expected < 1) || mean(asym$expected < 5) > .20
      if (sparse) {
        seed <- seed + 1L; set.seed(seed)
        ct <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 99999))
        method <- "Monte Carlo Pearson chi-squared (99,999 replicates)"
      } else {
        ct <- asym; method <- "Pearson chi-squared"
      }
      v <- sqrt(unname(asym$statistic) /
                (sum(tab) * min(nrow(tab) - 1, ncol(tab) - 1)))
      configuration_rows[[length(configuration_rows) + 1L]] <- data.table(
        cohort = spec$cohort, variable_id = spec$variable_id,
        variable_label = spec$variable_label, n = sum(tab), test = method,
        effect_size_name = "Cramer's V", effect_size = v, p_value = ct$p.value,
        sparse_table = sparse, min_expected = min(asym$expected)
      )
    }
  }
}

global_tests <- rbindlist(global_rows, fill = TRUE)
global_tests[, p_BH_within_cohort := p.adjust(p_value, method = "BH"), by = cohort]
global_tests[, p_BH_all_specific_tests := p.adjust(p_value, method = "BH")]
category_summary <- rbindlist(category_rows, fill = TRUE)
pairwise_tests <- rbindlist(pair_rows, fill = TRUE)
configuration_tests <- rbindlist(configuration_rows, fill = TRUE)
configuration_tests[, p_BH_within_cohort := p.adjust(p_value, method = "BH"), by = cohort]
configuration_tests[, p_BH_all_specific_tests := p.adjust(p_value, method = "BH")]

# Explicitly record the absence of extra fields in GSE84437.
coverage <- rbind(
  coverage,
  data.table(
    cohort = "GSE84437", variable_id = "no_additional_fields",
    variable_label = "No additional audited fields beyond age, sex, pT and pN",
    variable_type = "audit_note", family = "audit",
    source_field = "frozen projection table", cohort_n = nrow(g84437),
    nonmissing_n = NA_integer_, missing_n = NA_integer_, missing_fraction = NA_real_,
    usable_levels = NA_integer_
  ), fill = TRUE
)

fwrite(long, file.path(out_dir, "cohort_specific_patient_variable_long.tsv.gz"), sep = "\t")
fwrite(coverage, file.path(out_dir, "cohort_specific_variable_coverage.tsv"), sep = "\t")
fwrite(global_tests, file.path(out_dir, "cohort_specific_ecosystem_global_tests.tsv"), sep = "\t")
fwrite(category_summary, file.path(out_dir, "cohort_specific_category_ecosystem_summary.tsv"), sep = "\t")
fwrite(pairwise_tests, file.path(out_dir, "cohort_specific_pairwise_tests.tsv"), sep = "\t")
fwrite(configuration_tests, file.path(out_dir, "cohort_specific_configuration_tests.tsv"), sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  objective = "Audit all clinically or pathologically relevant fields unique to each frozen cohort",
  common_variables = "Analyzed separately in script 62; this package retains cohort-specific information",
  GSE84437_note = "No additional audited clinicopathological field beyond age, sex, pT and pN was found in the frozen or existing processed clinical table",
  exclusions = list(
    GSE62254_fu_status = "follow-up/outcome status, not a baseline clinicopathological characteristic",
    survival_fields = "outcomes are not used for clinical-context screening"
  ),
  tests = list(
    categorical_vs_continuous_ecosystem = "Kruskal-Wallis with epsilon-squared; pairwise Wilcoxon with rank-biserial effect",
    continuous_vs_ecosystem = "Spearman rho",
    descriptive_configuration = "Pearson chi-square or 99,999-replicate Monte Carlo chi-square after expected-count audit"
  ),
  source_sha256 = as.list(setNames(vapply(files, digest, character(1),
                                           algo = "sha256", file = TRUE), basename(files)))
)
write_json(audit, file.path(audit_dir, "cohort_specific_clinicopathology_v1.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("Cohort-specific clinicopathology analysis complete\n")
print(coverage)
print(global_tests)
