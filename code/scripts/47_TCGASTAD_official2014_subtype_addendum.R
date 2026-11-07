options(stringsAsFactors = FALSE, warn = 1)
options(error = function() {
  traceback(20L)
  quit(save = "no", status = 1L, runLast = FALSE)
})

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(digest)
})

score_file <- file.path(project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
                        "frozen_state_projection",
                        "TCGA_STAD_PLAT_ecosystem_scores_clinical.tsv.gz")
out_dir <- file.path(project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
                     "frozen_state_projection", "official_TCGA2014_subtype_addendum")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(score_file)) stop("TCGA score table absent")

dt <- fread(score_file)
official_levels <- c("CIN", "GS", "MSI", "EBV")
official <- dt[official_2014_available == TRUE & tcga2014_subtype %in% official_levels]
official[, tcga2014_subtype := factor(tcga2014_subtype, levels = official_levels)]
if (nrow(official) != 274L) stop("Unexpected official TCGA2014 subtype subset")

exposures <- c("sender_primary_z", "receiver_specific_z", "ecosystem_z")
tests <- rbindlist(lapply(exposures, function(exposure) {
  kt <- kruskal.test(official[[exposure]] ~ official$tcga2014_subtype)
  data.table(classification = "official_TCGA2014_four_subtypes",
             exposure = exposure, n = nrow(official), groups = 4L,
             kruskal_chisq = unname(kt$statistic), df = unname(kt$parameter),
             p_value = kt$p.value)
}))
tests[, FDR_BH := p.adjust(p_value, "BH")]

summary_table <- rbindlist(lapply(exposures, function(exposure) {
  official[, .(
    n = .N, median = median(get(exposure)),
    q1 = quantile(get(exposure), 0.25), q3 = quantile(get(exposure), 0.75)
  ), by = tcga2014_subtype][, exposure := exposure]
}))
ecotype_table <- dcast(official, ecotype4 ~ tcga2014_subtype,
                       fun.aggregate = length, value.var = "sample_submitter_id")
ecotype_test <- chisq.test(table(official$ecotype4, official$tcga2014_subtype))
set.seed(2026072947L)
ecotype_fisher <- fisher.test(table(official$ecotype4, official$tcga2014_subtype),
                              simulate.p.value = TRUE, B = 99999L)

fwrite(tests, file.path(out_dir, "TCGA_STAD_official2014_subtype_omnibus_tests.tsv"), sep = "\t")
fwrite(summary_table,
       file.path(out_dir, "TCGA_STAD_official2014_subtype_score_summary.tsv"), sep = "\t")
fwrite(ecotype_table,
       file.path(out_dir, "TCGA_STAD_ecotype4_by_official2014_subtype.tsv"), sep = "\t")

palette_subtype <- c(CIN = "#E76F51", GS = "#2A9D8F", MSI = "#457B9D", EBV = "#E9C46A")
p <- ggplot(official, aes(tcga2014_subtype, ecosystem_z, fill = tcga2014_subtype)) +
  geom_boxplot(width = 0.64, outlier.shape = NA, alpha = 0.82) +
  geom_jitter(width = 0.14, size = 0.75, alpha = 0.38) +
  scale_fill_manual(values = palette_subtype) +
  labs(x = NULL, y = "Ecosystem score (z)",
       title = "Official TCGA 2014 gastric-cancer molecular subtypes",
       subtitle = "274 tumors with an official four-subtype annotation") +
  theme_classic(base_size = 10) + theme(legend.position = "none")
ggsave(file.path(out_dir, "TCGA_STAD_official2014_subtype_ecosystem.png"),
       p, width = 5.4, height = 4.5, dpi = 320, bg = "white")
ggsave(file.path(out_dir, "TCGA_STAD_official2014_subtype_ecosystem.pdf"),
       p, width = 5.4, height = 4.5, device = cairo_pdf, bg = "white")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  reason = paste(
    "The existing tcga_subtype field is a broader PanCan classifier containing STAD_POLE.",
    "The official TCGA2014 four-subtype field is used here for gastric-cancer subtype claims."
  ),
  n_official = nrow(official),
  counts = official[, .N, by = tcga2014_subtype],
  ecotype_subtype_chisq = list(statistic = unname(ecotype_test$statistic),
                               df = unname(ecotype_test$parameter),
                               p_value = ecotype_test$p.value),
  ecotype_subtype_fisher_monte_carlo = list(
    replicates = 99999L, p_value = ecotype_fisher$p.value),
  input_sha256 = digest(score_file, algo = "sha256", file = TRUE)
)
write_json(audit,
           file.path(audit_dir, "TCGA_STAD_official2014_subtype_addendum.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nOfficial TCGA2014 subtype counts:\n")
print(official[, .N, by = tcga2014_subtype])
cat("\nOfficial TCGA2014 subtype tests:\n")
print(tests)
cat("\nOfficial TCGA2014 subtype ecosystem summary:\n")
print(summary_table[exposure == "ecosystem_z"][order(-median)])
cat("\nEcotype by official subtype:\n")
print(ecotype_table)
cat("Ecotype/subtype chi-square P = ", signif(ecotype_test$p.value, 5), "\n", sep = "")
cat("Ecotype/subtype Fisher Monte Carlo P = ",
    signif(ecotype_fisher$p.value, 5), "\n", sep = "")
