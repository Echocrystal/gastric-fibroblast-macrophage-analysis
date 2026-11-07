#!/usr/bin/env Rscript

# Corrective patient-level inference for the fine-state abundance association.
# Uses correct partial-correlation degrees of freedom, HC3 robust regression,
# restricted-model Rademacher wild bootstrap, and an explicit target-versus-
# remainder log-ratio representation for both lineages.

options(stringsAsFactors = FALSE, warn = 1)
set.seed(2026072742L)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(digest)
})

feature_file <- file.path(project, "analysis", "results", "PLAT_axis", "GSE183904",
                          "fine_state_axis", "GSE183904_fine_state_axis_patient_features.tsv")
state_count_file <- file.path(project, "analysis", "results", "PLAT_axis", "GSE183904",
                              "fine_state_axis", "GSE183904_fine_state_counts_by_patient.tsv")
amendment_file <- file.path(project, "analysis", "audit",
                            "GSE183904_fine_state_axis_association_amendment_01.json")
out_dir <- file.path(project, "analysis", "results", "PLAT_axis", "GSE183904",
                     "corrected_state_abundance_inference")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

features <- fread(feature_file)
state_counts <- fread(state_count_file)
if (nrow(features) != 26L || uniqueN(features$Patient_ID) != 26L) {
  stop("Expected 26 patient-independent tumors")
}

sender_state <- "CXCL14_POSTN_PLAT_rich_fibroblast"
receiver_state <- "C1QC_APOE_macrophage"

sender <- state_counts[lineage == "Fibroblast" & refined_label == sender_state,
                       .(sample_id, sender_cells = state_cells,
                         fibroblast_cells_check = lineage_cells)]
receiver <- state_counts[lineage == "Myeloid" & refined_label == receiver_state,
                         .(sample_id, receiver_cells = state_cells,
                           myeloid_cells_check = lineage_cells)]
dt <- merge(features, sender, by = "sample_id", all.x = TRUE)
dt <- merge(dt, receiver, by = "sample_id", all.x = TRUE)
dt[is.na(sender_cells), sender_cells := 0L]
dt[is.na(receiver_cells), receiver_cells := 0L]
if (any(dt[!is.na(fibroblast_cells_check), fibroblast_cells != fibroblast_cells_check]) ||
    any(dt[!is.na(myeloid_cells_check), myeloid_cells != myeloid_cells_check])) {
  stop("Lineage denominator mismatch")
}

pseudocount <- 0.5
dt[, sender_logratio := log((sender_cells + pseudocount) /
                              (fibroblast_cells - sender_cells + pseudocount))]
dt[, receiver_logratio := log((receiver_cells + pseudocount) /
                                (myeloid_cells - receiver_cells + pseudocount))]
dt[, sender_asin_fraction := asin(sqrt(sender_cells / fibroblast_cells))]
dt[, receiver_asin_fraction := asin(sqrt(receiver_cells / myeloid_cells))]

covariates <- cbind(
  ECM_direction = dt$ECM_direction,
  log_fibroblast_cells = log1p(dt$fibroblast_cells),
  log_myeloid_cells = log1p(dt$myeloid_cells)
)

hc3_fit <- function(y, X) {
  fit <- lm.fit(X, y)
  if (fit$rank != ncol(X)) stop("Rank-deficient corrective model")
  xtx_inv <- chol2inv(chol(crossprod(X)))
  h <- rowSums((X %*% xtx_inv) * X)
  adjusted_sq <- (fit$residuals / pmax(1e-8, 1 - h))^2
  meat <- crossprod(X, X * adjusted_sq)
  vcov_hc3 <- xtx_inv %*% meat %*% xtx_inv
  se <- sqrt(diag(vcov_hc3))
  df <- nrow(X) - ncol(X)
  t_value <- fit$coefficients[2] / se[2]
  list(
    fit = fit,
    beta = fit$coefficients[2],
    se = se[2],
    t = t_value,
    df = df,
    p = 2 * pt(-abs(t_value), df = df),
    leverage = h
  )
}

analyze_one <- function(x, y, covariates, representation, B = 9999L, seed = 1L) {
  ok <- is.finite(x) & is.finite(y) & complete.cases(covariates)
  x <- x[ok]
  y <- y[ok]
  C <- as.matrix(covariates[ok, , drop = FALSE])
  xr <- rank(x, ties.method = "average")
  yr <- rank(y, ties.method = "average")
  X0 <- cbind(intercept = 1, C)
  X1 <- cbind(intercept = 1, exposure_rank = xr, C)

  x_resid <- lm.fit(X0, xr)$residuals
  y_resid <- lm.fit(X0, yr)$residuals
  adjusted_rho <- unname(cor(x_resid, y_resid))
  q <- ncol(C)
  partial_df <- length(xr) - q - 2L
  partial_t <- adjusted_rho * sqrt(partial_df / (1 - adjusted_rho^2))
  partial_p <- 2 * pt(-abs(partial_t), df = partial_df)

  robust <- hc3_fit(yr, X1)
  reduced <- lm.fit(X0, yr)
  xtx0_inv <- chol2inv(chol(crossprod(X0)))
  h0 <- rowSums((X0 %*% xtx0_inv) * X0)
  null_resid_hc3 <- reduced$residuals / pmax(1e-8, 1 - h0)

  set.seed(seed)
  null_t <- vapply(seq_len(B), function(i) {
    weights <- sample(c(-1, 1), length(yr), replace = TRUE)
    y_star <- reduced$fitted.values + null_resid_hc3 * weights
    hc3_fit(y_star, X1)$t
  }, numeric(1))
  wild_p <- (1 + sum(abs(null_t) >= abs(robust$t))) / (B + 1)

  loo <- vapply(seq_along(xr), function(drop) {
    keep <- seq_along(xr) != drop
    xr_i <- lm.fit(X0[keep, , drop = FALSE], xr[keep])$residuals
    yr_i <- lm.fit(X0[keep, , drop = FALSE], yr[keep])$residuals
    cor(xr_i, yr_i)
  }, numeric(1))

  data.table(
    representation = representation,
    n = length(xr),
    covariates = q,
    adjusted_rank_residual_rho = adjusted_rho,
    correct_partial_df = partial_df,
    correct_partial_t = partial_t,
    correct_partial_asymptotic_p = partial_p,
    exposure_rank_beta = robust$beta,
    HC3_SE = robust$se,
    HC3_t = robust$t,
    HC3_df = robust$df,
    HC3_p = robust$p,
    wild_bootstrap_draws = B,
    wild_bootstrap_t_p = wild_p,
    full_model_max_leverage = max(robust$leverage),
    reduced_model_max_leverage = max(h0),
    loo_rho_min = min(loo),
    loo_rho_max = max(loo),
    loo_positive_fraction = mean(loo > 0)
  )
}

results <- rbindlist(list(
  analyze_one(dt$sender_logratio, dt$receiver_logratio, covariates,
              "primary_target_vs_remainder_logratio", 9999L, 2026072742L),
  analyze_one(dt$sender_asin_fraction, dt$receiver_asin_fraction, covariates,
              "sensitivity_arcsin_sqrt_fraction", 9999L, 2026072743L)
))

fwrite(dt[, .(
  sample_id, Patient_ID, ECM_direction,
  fibroblast_cells, sender_cells, sender_logratio, sender_asin_fraction,
  myeloid_cells, receiver_cells, receiver_logratio, receiver_asin_fraction
)], file.path(out_dir, "GSE183904_corrected_state_abundance_patient_data.tsv"), sep = "\t")
fwrite(results,
       file.path(out_dir, "GSE183904_corrected_state_abundance_inference.tsv"), sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE183904",
  purpose = "post-audit corrective patient-level state-abundance inference",
  sender_state = sender_state,
  receiver_state = receiver_state,
  inference_unit = "patient/sample",
  primary_representation = paste(
    "log((target state cells + 0.5)/(remaining lineage cells + 0.5))",
    "for both fibroblast sender and myeloid receiver"
  ),
  sensitivity_representation = "arcsin-square-root lineage fraction",
  covariates = c("ECM direction", "log1p fibroblast cells", "log1p myeloid cells"),
  methods = c(
    "partial rank-residual correlation with n-q-2 degrees of freedom",
    "rank regression with HC3 standard error",
    "restricted-model Rademacher wild bootstrap-t with HC3 studentization",
    "leave-one-patient-out adjusted rank-residual correlation"
  ),
  multiplicity = paste(
    "The target-versus-remainder log-ratio is the primary corrective representation;",
    "the arcsin-square-root fraction repeats the original representation as sensitivity."
  ),
  interpretation = paste(
    "Post-audit corrective exploratory association in the same cohort.",
    "The result remains observational and may reflect unmeasured tumor-state composition or other confounding."
  ),
  input_sha256 = list(
    patient_features = digest(feature_file, algo = "sha256", file = TRUE),
    state_counts = digest(state_count_file, algo = "sha256", file = TRUE),
    amendment = digest(amendment_file, algo = "sha256", file = TRUE)
  )
)
write_json(audit,
           file.path(audit_dir, "GSE183904_corrected_state_abundance_inference.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nCorrected patient-level state-abundance inference:\n")
print(results)
