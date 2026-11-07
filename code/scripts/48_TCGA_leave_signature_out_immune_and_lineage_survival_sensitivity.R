options(stringsAsFactors = FALSE, warn = 1)
options(error = function() {
  traceback(20L)
  quit(save = "no", status = 1L, runLast = FALSE)
})
set.seed(2026072948L)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(MCPcounter)
  library(survival)
  library(sandwich)
  library(lmtest)
  library(jsonlite)
  library(digest)
})

signature_file <- file.path(project, "resources", "gene_sets",
                            "GSE183904_PLAT_ecosystem_bulk_signatures.tsv")
expression_file <- "I:/ECM-new/ECM-public/data/processed/bulk/tcga_stad_full_symbol_log2cpm_412.rds"
score_file <- file.path(project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
                        "frozen_state_projection",
                        "TCGA_STAD_PLAT_ecosystem_scores_clinical.tsv.gz")
program_file <- file.path(source_project, "data", "raw", "gene_sets",
                          "bulk_immune_programs_locked.tsv")
mcp_gene_file <- file.path(source_project, "runtime", "sources", "MCPcounter-master",
                           "Signatures", "genes.txt")
mcp_probe_file <- file.path(source_project, "runtime", "sources", "MCPcounter-master",
                            "Signatures", "probesets.txt")
out_dir <- file.path(project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
                     "leave_signature_out_sensitivity")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(signature_file, expression_file, score_file,
            program_file, mcp_gene_file, mcp_probe_file)
if (!all(file.exists(inputs))) stop("Missing leave-signature-out input")

sig <- fread(signature_file)
primary_ids <- c("sender_CXCL14_POSTN_state_axis_excluded",
                 "receiver_C1QC_APOE_specific_axis_excluded")
signature_union <- sort(unique(sig[signature_id %in% primary_ids, gene_symbol]))
if (length(signature_union) != 60L) stop("Expected 60 non-overlapping primary signature genes")

expr <- readRDS(expression_file)$expression
dt <- fread(score_file)
if (!identical(colnames(expr), dt$sample_submitter_id)) {
  idx <- match(colnames(expr), dt$sample_submitter_id)
  if (anyNA(idx)) stop("TCGA expression/score alignment failure")
  dt <- dt[idx]
}

program_manifest <- fread(program_file)
programs <- split(program_manifest$gene_symbol, program_manifest$program)
mcp_genes <- read.table(mcp_gene_file, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, colClasses = "character",
                        check.names = FALSE)
mcp_probesets <- read.table(mcp_probe_file, sep = "\t", header = FALSE,
                            stringsAsFactors = FALSE, colClasses = "character")

overlap_program <- rbindlist(lapply(names(programs), function(program) {
  genes <- unique(programs[[program]])
  overlap <- intersect(genes, signature_union)
  data.table(endpoint_type = "program", endpoint = program,
             original_signature_genes = length(genes), overlap_n = length(overlap),
             overlap_fraction = length(overlap) / length(genes),
             overlap_genes = paste(overlap, collapse = ";"))
}))
mcp_split <- split(mcp_genes[["HUGO symbols"]], mcp_genes[["Cell population"]])
overlap_mcp <- rbindlist(lapply(names(mcp_split), function(endpoint) {
  genes <- unique(mcp_split[[endpoint]])
  overlap <- intersect(genes, signature_union)
  data.table(endpoint_type = "MCPcounter", endpoint = endpoint,
             original_signature_genes = length(genes), overlap_n = length(overlap),
             overlap_fraction = length(overlap) / length(genes),
             overlap_genes = paste(overlap, collapse = ";"))
}))
overlap_table <- rbindlist(list(overlap_program, overlap_mcp), fill = TRUE)
fwrite(overlap_table,
       file.path(out_dir, "TCGA_STAD_endpoint_primary_signature_overlap.tsv"), sep = "\t")

expr_leaveout <- expr[!rownames(expr) %in% signature_union, , drop = FALSE]
module_scores <- function(expression, gene_sets) {
  z <- t(scale(t(expression), center = TRUE, scale = TRUE))
  z[!is.finite(z)] <- 0
  out <- vapply(gene_sets, function(gs) {
    present <- intersect(setdiff(unique(gs), signature_union), rownames(z))
    if (length(present) < 2L) return(rep(NA_real_, ncol(z)))
    colMeans(z[present, , drop = FALSE])
  }, numeric(ncol(z)))
  out <- t(out)
  colnames(out) <- colnames(expression)
  out
}
program_leaveout <- module_scores(expr_leaveout, programs)
mcp_leaveout <- MCPcounter.estimate(
  expr_leaveout, featuresType = "HUGO_symbols",
  probesets = mcp_probesets, genes = mcp_genes
)
leaveout_scores <- cbind(t(mcp_leaveout), t(program_leaveout))
leaveout_dt <- as.data.table(leaveout_scores, keep.rownames = "sample_submitter_id")
dt <- merge(dt, leaveout_dt, by = "sample_submitter_id", all.x = TRUE,
            sort = FALSE, suffixes = c("_original", "_leaveout"))
if (nrow(dt) != 412L) stop("Leaveout score merge failure")
dt[!is.na(stage_group_frozen) & !nzchar(stage_group_frozen),
   stage_group_frozen := NA_character_]
dt[!is.na(sex) & !nzchar(sex), sex := NA_character_]
dt[!is.na(tcga_subtype) & !nzchar(tcga_subtype), tcga_subtype := NA_character_]

endpoint_names <- c(rownames(mcp_leaveout), rownames(program_leaveout))
# Colliding endpoint names were suffixed by merge; select the recomputed fields explicitly.
leaveout_name <- function(endpoint) {
  if (paste0(endpoint, "_leaveout") %in% names(dt)) paste0(endpoint, "_leaveout") else endpoint
}

zscore <- function(x) as.numeric((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
hc3_model <- function(endpoint, exposure) {
  endpoint_col <- leaveout_name(endpoint)
  d <- data.frame(
    y = zscore(dt[[endpoint_col]]), exposure = dt[[exposure]],
    purity = dt$absolute_purity,
    subtype = factor(dt$tcga_subtype),
    logTMB = log1p(pmax(0, dt$TMB_NONSYNONYMOUS.x)),
    FGA = dt$FRACTION_GENOME_ALTERED.x
  )
  d <- d[complete.cases(d), , drop = FALSE]
  fit <- lm(y ~ exposure + purity + subtype + logTMB + FGA, data = d)
  ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
  data.table(
    endpoint = endpoint, leaveout_column = endpoint_col, exposure = exposure,
    n = nrow(d), covariates = "ABSOLUTE_purity+TCGA_subtype+log1p_TMB+FGA",
    beta_per_exposure_SD = unname(ct["exposure", 1]),
    HC3_se = unname(ct["exposure", 2]),
    ci_low = unname(ct["exposure", 1] - qnorm(0.975) * ct["exposure", 2]),
    ci_high = unname(ct["exposure", 1] + qnorm(0.975) * ct["exposure", 2]),
    HC3_p = unname(ct["exposure", 4])
  )
}

exposures <- c("sender_primary_z", "receiver_specific_z", "ecosystem_z")
immune_models <- rbindlist(lapply(exposures, function(exposure) {
  rbindlist(lapply(endpoint_names, function(endpoint) hc3_model(endpoint, exposure)))
}))
immune_models[, FDR_BH_family := p.adjust(HC3_p, "BH")]
fwrite(immune_models,
       file.path(out_dir, "TCGA_STAD_leave_signature_out_immune_models_HC3.tsv"), sep = "\t")
fwrite(leaveout_dt,
       file.path(out_dir, "TCGA_STAD_leave_signature_out_endpoint_scores.tsv.gz"), sep = "\t")

cox_sensitivity <- function(exposure, tier) {
  covariates <- switch(
    tier,
    preregistered = c("age", "sex", "stage_group_frozen", "absolute_purity"),
    broad_lineage = c("age", "sex", "stage_group_frozen", "absolute_purity",
                      leaveout_name("Fibroblasts"), leaveout_name("Monocytic lineage")),
    broad_lineage_subtype = c("age", "sex", "stage_group_frozen", "absolute_purity",
                              leaveout_name("Fibroblasts"), leaveout_name("Monocytic lineage"),
                              "tcga_subtype")
  )
  vars <- c("os_months", "os_event_frozen", exposure, covariates)
  d <- as.data.frame(dt[, ..vars])
  d <- d[complete.cases(d) & is.finite(d$os_months) & d$os_months > 0, , drop = FALSE]
  names(d)[match(exposure, names(d))] <- "exposure"
  names(d) <- make.names(names(d), unique = TRUE)
  rhs <- paste(setdiff(names(d), c("os_months", "os_event_frozen", "exposure")),
               collapse = " + ")
  formula <- as.formula(paste("Surv(os_months, os_event_frozen) ~ exposure +", rhs))
  fit <- coxph(formula, data = d, ties = "efron", x = TRUE)
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
    PH_global_p = unname(zph$table["GLOBAL", "p"])
  )
}

cox_models <- rbindlist(lapply(exposures, function(exposure) {
  rbindlist(lapply(c("preregistered", "broad_lineage", "broad_lineage_subtype"),
                   function(tier) cox_sensitivity(exposure, tier)))
}))
cox_models[, FDR_BH_within_tier := p.adjust(p_value, "BH"), by = tier]
fwrite(cox_models,
       file.path(out_dir, "TCGA_STAD_lineage_adjusted_Cox_sensitivity.tsv"), sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  correction = c(
    "TMB covariate uses the existing TMB_NONSYNONYMOUS.x column.",
    "All 60 primary signature genes are removed before recomputing MCPcounter and immune programs.",
    "Broad-lineage survival sensitivity uses leave-signature-out Fibroblasts and Monocytic lineage scores."
  ),
  excluded_signature_genes = signature_union,
  overlapping_endpoints = overlap_table[overlap_n > 0],
  input_sha256 = as.list(setNames(vapply(inputs, digest, character(1),
                                         algo = "sha256", file = TRUE), basename(inputs)))
)
write_json(audit,
           file.path(audit_dir, "TCGA_STAD_leave_signature_out_sensitivity.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nEndpoint/signature overlap:\n")
print(overlap_table[overlap_n > 0][order(-overlap_n)])
cat("\nTop leave-signature-out ecosystem associations:\n")
print(immune_models[exposure == "ecosystem_z"][order(HC3_p)][1:20])
cat("\nLineage-adjusted Cox sensitivity:\n")
print(cox_models)
