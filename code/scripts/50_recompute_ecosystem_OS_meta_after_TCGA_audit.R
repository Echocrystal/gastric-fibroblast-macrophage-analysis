options(stringsAsFactors = FALSE, warn = 1)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(digest)
})

external_file <- file.path(project, "analysis", "results", "PLAT_axis",
                           "external_validation", "frozen_state_projection",
                           "external_GC_continuous_Cox_models.tsv")
tcga_file <- file.path(project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
                       "leave_signature_out_sensitivity",
                       "TCGA_STAD_lineage_adjusted_Cox_sensitivity.tsv")
out_dir <- file.path(project, "analysis", "results", "PLAT_axis",
                     "external_validation", "frozen_state_projection")
audit_dir <- file.path(project, "analysis", "audit")
if (!all(file.exists(c(external_file, tcga_file)))) stop("Meta-analysis input absent")

external <- fread(external_file)[exposure == "ecosystem_z" & model == "cohort_adjusted",
  .(cohort, n, events, logHR, se_logHR, HR_per_SD, ci_low, ci_high, p_value)]
tcga_all <- fread(tcga_file)[exposure == "ecosystem_z"]
tcga_rows <- tcga_all[, .(
  cohort = paste0("TCGA-STAD_", tier), n, events,
  logHR = log(HR_per_SD),
  se_logHR = (log(ci_high) - log(ci_low)) / (2 * qnorm(0.975)),
  HR_per_SD, ci_low, ci_high, p_value, tier
)]

meta_one <- function(d, analysis_id) {
  wi <- 1 / d$se_logHR^2
  fixed <- sum(wi * d$logHR) / sum(wi)
  se_fixed <- sqrt(1 / sum(wi))
  Q <- sum(wi * (d$logHR - fixed)^2)
  Q_df <- nrow(d) - 1L
  C <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- max(0, (Q - Q_df) / C)
  wr <- 1 / (d$se_logHR^2 + tau2)
  random <- sum(wr * d$logHR) / sum(wr)
  se_random <- sqrt(1 / sum(wr))
  rbindlist(list(
    data.table(analysis_id, meta_model = "fixed_effect", k = nrow(d),
               HR = exp(fixed), ci_low = exp(fixed - qnorm(0.975) * se_fixed),
               ci_high = exp(fixed + qnorm(0.975) * se_fixed),
               p_value = 2 * pnorm(-abs(fixed / se_fixed))),
    data.table(analysis_id, meta_model = "random_effect_DL", k = nrow(d),
               HR = exp(random), ci_low = exp(random - qnorm(0.975) * se_random),
               ci_high = exp(random + qnorm(0.975) * se_random),
               p_value = 2 * pnorm(-abs(random / se_random)))
  ))[, `:=`(Q = Q, Q_df = Q_df, Q_p = pchisq(Q, Q_df, lower.tail = FALSE),
            I2_percent = if (Q > 0) max(0, (Q - Q_df) / Q) * 100 else 0,
            tau2 = tau2)]
}

inputs <- list(
  external_only = external,
  external_plus_TCGA_preregistered = rbindlist(list(
    external, tcga_rows[tier == "preregistered",
                        .(cohort, n, events, logHR, se_logHR,
                          HR_per_SD, ci_low, ci_high, p_value)])),
  external_plus_TCGA_broad_lineage = rbindlist(list(
    external, tcga_rows[tier == "broad_lineage",
                        .(cohort, n, events, logHR, se_logHR,
                          HR_per_SD, ci_low, ci_high, p_value)])),
  external_plus_TCGA_broad_lineage_subtype = rbindlist(list(
    external, tcga_rows[tier == "broad_lineage_subtype",
                        .(cohort, n, events, logHR, se_logHR,
                          HR_per_SD, ci_low, ci_high, p_value)]))
)
meta_input <- rbindlist(lapply(names(inputs), function(id) {
  copy(inputs[[id]])[, analysis_id := id]
}), fill = TRUE)
meta_summary <- rbindlist(lapply(names(inputs), function(id) meta_one(inputs[[id]], id)))

fwrite(meta_input,
       file.path(out_dir, "GC_ecosystem_OS_meta_input_audit_corrected.tsv"), sep = "\t")
fwrite(meta_summary,
       file.path(out_dir, "GC_ecosystem_OS_meta_summary_audit_corrected.tsv"), sep = "\t")
write_json(
  list(generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
       rationale = paste(
         "Recomputed after the TCGA audit using the corrected 372-patient preregistered model",
         "and explicit broad-lineage TCGA sensitivity models."
       ),
       summary = meta_summary,
       input_sha256 = list(
         external = digest(external_file, algo = "sha256", file = TRUE),
         TCGA = digest(tcga_file, algo = "sha256", file = TRUE))),
  file.path(audit_dir, "GC_ecosystem_OS_meta_audit_corrected.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA
)
print(meta_summary)

