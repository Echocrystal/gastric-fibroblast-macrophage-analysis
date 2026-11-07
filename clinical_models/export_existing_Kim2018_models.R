#!/usr/bin/env Rscript
# Recover complete coefficient tables for existing models only.
# The fixed input scores, response labels and model specifications are unchanged.
args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1L])
runtime_args <- commandArgs(trailingOnly = TRUE)
# Optional first argument is this clinical_models package root; default is script directory.
task_dir <- if (length(runtime_args)) {
  normalizePath(runtime_args[1L], winslash = "/", mustWork = TRUE)
} else dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
extra_libs <- strsplit(Sys.getenv("KIM2018_R_LIBRARIES", ""), ";", fixed = TRUE)[[1L]]
.libPaths(c(extra_libs[dir.exists(extra_libs)], .libPaths()))
suppressPackageStartupMessages({library(data.table); library(brglm2); library(jsonlite)})
out_dir <- file.path(task_dir, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dat <- fread(file.path(task_dir, "inputs", "Kim2018_model_inputs_deidentified.tsv"),
             na.strings = c("NA", "N/A", ""))
stopifnot(nrow(dat) == 45L, uniqueN(dat$patient_id) == 45L,
          sum(dat$response_binary == 1L) == 12L,
          identical(as.integer(dat$MSI_binary), as.integer(dat$msi_type_pentaplex == "MSI")))

exposures <- c("sender_primary_z", "receiver_specific_z", "ecosystem_z",
               "ecosystem_merged_sensitivity_z")
specs <- list()
for (exposure in exposures) for (adjusted in c(FALSE, TRUE)) {
  adjustment <- if (adjusted) "MSI_adjusted" else "univariable"
  model_id <- paste(exposure, adjustment, sep = "__")
  specs[[model_id]] <- list(
    kind = "single_score", exposure = exposure, original_model = adjustment,
    formula = as.formula(paste("response_binary ~", exposure,
                              if (adjusted) "+ MSI_binary" else "")),
    original_source_lines = "46_validate_frozen_PLAT_ecosystem_external_GC_cohorts.R:262-294"
  )
}
specs$sender_receiver <- list(kind = "joint", exposure = NA_character_,
  original_model = "sender_receiver",
  formula = response_binary ~ sender_primary_z + receiver_specific_z,
  original_source_lines = "46_validate_frozen_PLAT_ecosystem_external_GC_cohorts.R:308-325")
specs$sender_receiver_MSI <- list(kind = "joint", exposure = NA_character_,
  original_model = "sender_receiver_MSI",
  formula = response_binary ~ sender_primary_z + receiver_specific_z + MSI_binary,
  original_source_lines = "46_validate_frozen_PLAT_ecosystem_external_GC_cohorts.R:308-325")
specs$sender_receiver_interaction <- list(kind = "joint", exposure = NA_character_,
  original_model = "sender_receiver_interaction",
  formula = response_binary ~ sender_primary_z * receiver_specific_z,
  original_source_lines = "46_validate_frozen_PLAT_ecosystem_external_GC_cohorts.R:308-325")

fits <- list(); all_terms <- list(); summaries <- list(); inclusion <- list()
for (model_id in names(specs)) {
  s <- specs[[model_id]]
  vars <- all.vars(s$formula)
  cc <- complete.cases(dat[, ..vars])
  d <- as.data.frame(dat[cc, ..vars])
  fit_formula <- s$formula
  # Match the original helper, which named each single-score predictor 'exposure'.
  if (s$kind == "single_score") {
    names(d)[names(d) == s$exposure] <- "exposure"
    fit_formula <- if (s$original_model == "MSI_adjusted")
      response_binary ~ exposure + MSI_binary else response_binary ~ exposure
  }
  warnings <- character()
  fit <- withCallingHandlers(
    glm(fit_formula, data = d, family = binomial("logit"),
        method = brglm2::brglmFit, type = "AS_mean"),
    warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}
  )
  fits[[model_id]] <- fit
  sm <- as.data.table(summary(fit)$coefficients, keep.rownames = "fit_term")
  setnames(sm, c("Estimate", "Std. Error", "z value", "Pr(>|z|)"),
           c("coefficient_log_odds", "standard_error", "z_statistic", "p_value"))
  sm[, term := fit_term]
  if (s$kind == "single_score") sm[fit_term == "exposure", term := s$exposure]
  converged <- isTRUE(fit$converged)
  status <- if (!converged) "not_converged" else if (length(warnings)) "converged_with_warnings" else "converged"
  sm[, `:=`(
    model_id = model_id, formula = paste(deparse(s$formula), collapse = " "),
    n = nrow(d), events = sum(d$response_binary == 1L),
    non_events = sum(d$response_binary == 0L), input_n = nrow(dat),
    missing_or_excluded_n = sum(!cc), OR_or_baseline_odds = exp(coefficient_log_odds),
    ci_low = exp(coefficient_log_odds - qnorm(0.975) * standard_error),
    ci_high = exp(coefficient_log_odds + qnorm(0.975) * standard_error),
    fit_status = status, converged = converged,
    original_kind = s$kind, original_model = s$original_model,
    original_source_lines = s$original_source_lines
  )]
  sm[, effect_definition := fifelse(term == "(Intercept)", "Baseline odds (not an odds ratio)",
    fifelse(term == "MSI_binary", "Odds ratio for MSI versus MSS",
    fifelse(grepl(":", term, fixed = TRUE), "Interaction odds ratio per unit product of sender and receiver Z scores",
            "Odds ratio per one within-cohort standard deviation")))]
  if (model_id == "sender_receiver_interaction") {
    sm[term == "sender_primary_z", effect_definition := "Sender odds ratio per one SD at receiver Z=0"]
    sm[term == "receiver_specific_z", effect_definition := "Receiver odds ratio per one SD at sender Z=0"]
  }
  summaries[[model_id]] <- data.table(
    model_id = model_id, formula = paste(deparse(s$formula), collapse = " "),
    input_n = nrow(dat), n = nrow(d), events = sum(d$response_binary == 1L),
    non_events = sum(d$response_binary == 0L), missing_or_excluded_n = sum(!cc),
    fit_status = status, converged = converged, iterations = fit$iter,
    design_rank = fit$rank, coefficient_n = length(coef(fit)),
    residual_degrees_of_freedom = fit$df.residual,
    boundary = if (is.null(fit$boundary)) NA else fit$boundary,
    warning_text = paste(unique(warnings), collapse = " | "),
    estimator = "brglm2::brglmFit; type=AS_mean; binomial(logit)",
    original_source_lines = s$original_source_lines
  )
  all_terms[[model_id]] <- sm
  inclusion[[model_id]] <- data.table(
    model_id = model_id, patient_id = dat$patient_id, included = cc,
    response_binary = dat$response_binary,
    missing_model_variables = vapply(seq_len(nrow(dat)), function(i) {
      paste(vars[vapply(vars, function(v) is.na(dat[[v]][i]), logical(1))], collapse = ";")
    }, character(1))
  )
}
terms <- rbindlist(all_terms, fill = TRUE)
model_summary <- rbindlist(summaries, fill = TRUE)
inclusion <- rbindlist(inclusion)

# Preserve the source correction families; do not adjust a new all-coefficient family.
old_single <- fread(file.path(task_dir, "inputs", "Kim2018_bias_reduced_response_models_original.tsv"))
old_joint <- fread(file.path(task_dir, "inputs", "Kim2018_continuous_joint_state_response_models_original.tsv"))
terms[, `:=`(original_BH_q = NA_real_, original_BH_family = NA_character_,
             originally_reported = FALSE)]
comparison <- list()
for (i in seq_len(nrow(old_single))) {
  o <- old_single[i]
  id <- paste(o$exposure, o$model, sep = "__")
  idx <- which(terms$model_id == id & terms$term == o$exposure)
  stopifnot(length(idx) == 1L)
  set(terms, idx, "original_BH_q", o$FDR_BH_within_model)
  set(terms, idx, "original_BH_family", paste0("Original four single-score exposure tests: ", o$model))
  set(terms, idx, "originally_reported", TRUE)
  for (metric in c("OR_or_baseline_odds", "ci_low", "ci_high", "p_value")) {
    old_column <- if (metric == "OR_or_baseline_odds") "OR_per_SD_for_response" else metric
    comparison[[length(comparison)+1L]] <- data.table(model_id = id, term = o$exposure,
      metric = metric, original_value = o[[old_column]], reproduced_value = terms[[metric]][idx],
      absolute_difference = abs(o[[old_column]] - terms[[metric]][idx]))
  }
}
for (i in seq_len(nrow(old_joint))) {
  o <- old_joint[i]
  idx <- which(terms$model_id == o$model_id & terms$term == o$term)
  stopifnot(length(idx) == 1L)
  set(terms, idx, "original_BH_q", o$FDR_BH_exploratory_family)
  set(terms, idx, "original_BH_family", "Original eight non-intercept terms across three joint/interaction models")
  set(terms, idx, "originally_reported", TRUE)
  for (metric in c("OR_or_baseline_odds", "ci_low", "ci_high", "p_value")) {
    old_column <- if (metric == "OR_or_baseline_odds") "OR" else metric
    comparison[[length(comparison)+1L]] <- data.table(model_id = o$model_id, term = o$term,
      metric = metric, original_value = o[[old_column]], reproduced_value = terms[[metric]][idx],
      absolute_difference = abs(o[[old_column]] - terms[[metric]][idx]))
  }
}
comparison <- rbindlist(comparison)
comparison[, pass := absolute_difference <= 1e-9 * pmax(1, abs(original_value))]
stopifnot(all(comparison$pass), all(model_summary$converged))

plot_map <- data.table(
  figure_term = c("Average ecosystem", "Fibroblast sender", "Macrophage receiver", "Sender x receiver"),
  model_id = c("ecosystem_z__MSI_adjusted", "sender_receiver_MSI", "sender_receiver_MSI", "sender_receiver_interaction"),
  term = c("ecosystem_z", "sender_primary_z", "receiver_specific_z", "sender_primary_z:receiver_specific_z")
)
plotted <- merge(plot_map, terms, by = c("model_id", "term"), sort = FALSE)
old_plot <- fread(file.path(task_dir, "inputs", "Figure6h_Kim2018_continuous_models_original.tsv"))
plot_check <- rbindlist(lapply(seq_len(nrow(plot_map)), function(i) {
  m <- plot_map[i]; t <- terms[model_id == m$model_id & term == m$term]
  old_values <- unlist(old_plot[i, .(OR, low, high, p)], use.names = FALSE)
  new_values <- unlist(t[, .(OR_or_baseline_odds, ci_low, ci_high, p_value)], use.names = FALSE)
  data.table(figure_term = m$figure_term, model_id = m$model_id, term = m$term,
    metric = c("OR", "CI_low", "CI_high", "P"), original_value = old_values,
    reproduced_value = new_values, absolute_difference = abs(old_values - new_values))
}))
plot_check[, pass := absolute_difference <= 1e-9 * pmax(1, abs(original_value))]
stopifnot(all(plot_check$pass), nrow(terms) == 31L, nrow(model_summary) == 11L)

setcolorder(terms, c("model_id", "formula", "term", "effect_definition", "coefficient_log_odds",
  "standard_error", "z_statistic", "OR_or_baseline_odds", "ci_low", "ci_high", "p_value",
  "original_BH_q", "original_BH_family", "originally_reported", "n", "events", "non_events",
  "input_n", "missing_or_excluded_n", "fit_status", "converged", "fit_term", "original_kind",
  "original_model", "original_source_lines"))
fwrite(terms, file.path(out_dir, "Kim2018_existing_models_ALL_TERMS.tsv"), sep = "\t", na = "")
fwrite(model_summary, file.path(out_dir, "Kim2018_existing_models_FIT_SUMMARY.tsv"), sep = "\t", na = "")
fwrite(inclusion, file.path(out_dir, "Kim2018_model_patient_inclusion.tsv"), sep = "\t", na = "")
fwrite(comparison, file.path(out_dir, "Kim2018_original16_terms_REPRODUCTION_QA.tsv"), sep = "\t")
fwrite(plot_check, file.path(out_dir, "Figure6H_four_terms_REPRODUCTION_QA.tsv"), sep = "\t")
fwrite(plotted, file.path(out_dir, "Figure6H_four_terms_complete.tsv"), sep = "\t", na = "")
write_json(list(
  scope = "Complete coefficient recovery for the 11 existing models in original script 46; no new model specifications",
  input_n = nrow(dat), events = sum(dat$response_binary), term_n = nrow(terms), model_n = nrow(model_summary),
  all_converged = all(model_summary$converged), all_original16_terms_match = all(comparison$pass),
  all_Figure6H_four_terms_match = all(plot_check$pass),
  maximum_original_output_absolute_difference = max(comparison$absolute_difference),
  maximum_Figure6H_absolute_difference = max(plot_check$absolute_difference),
  correction_policy = "Original BH q values preserved only for the originally tested terms and their existing families; new intercept/covariate reporting did not create a new correction family.",
  model_summary = model_summary, coefficients = terms,
  software = list(R = R.version.string, brglm2 = as.character(packageVersion("brglm2")),
                  data_table = as.character(packageVersion("data.table")), jsonlite = as.character(packageVersion("jsonlite")))
), file.path(out_dir, "Kim2018_existing_models_complete.json"), pretty = TRUE, auto_unbox = TRUE, na = "null", digits = 16)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
cat("Recovered", nrow(terms), "terms in", nrow(model_summary), "existing models.\n")
cat("All original 16 terms and Figure6H 4 terms reproduce within tolerance.\n")
print(model_summary[, .(model_id, n, events, missing_or_excluded_n, fit_status, iterations)])
