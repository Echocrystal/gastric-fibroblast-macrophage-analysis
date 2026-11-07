suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(msigdbr)
  library(sandwich)
  library(lmtest)
  library(jsonlite)
  library(digest)
})

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
tcga_file <- "I:/ECM-new/ECM-public/data/processed/bulk/tcga_stad_full_symbol_log2cpm_412.rds"
score_file <- file.path(project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
                        "frozen_state_projection",
                        "TCGA_STAD_PLAT_ecosystem_scores_clinical.tsv.gz")
signature_file <- file.path(project, "resources", "gene_sets",
                            "GSE183904_PLAT_ecosystem_bulk_signatures.tsv")
lineage_file <- file.path(project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
                          "leave_signature_out_sensitivity",
                          "TCGA_STAD_leave_signature_out_endpoint_scores.tsv.gz")
tide_file <- "I:/ECM-new/scRNA/26.matrix_feedback_multicohort/results/source_data/CONSENSUS_ECM_subtype_TIDE_sample_level.tsv"
out_dir <- file.path(project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
                     "continuous_ecosystem_function_immune_v1")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(tcga_file, score_file, signature_file, lineage_file, tide_file)
if (!all(file.exists(inputs))) stop("Missing continuous-ecosystem input")

tcga <- readRDS(tcga_file)
expr <- tcga$expression
scores <- fread(score_file)
lineage <- fread(lineage_file)
sig <- fread(signature_file)
primary_ids <- c("sender_CXCL14_POSTN_state_axis_excluded",
                 "receiver_C1QC_APOE_specific_axis_excluded")
primary_genes <- unique(sig[signature_id %in% primary_ids, gene_symbol])
if (length(primary_genes) != 60L) stop("Expected 60 primary signature genes")

meta <- merge(scores, lineage[, .(sample_submitter_id,
                                   broad_fibroblasts_LO = Fibroblasts,
                                   broad_monocytic_LO = `Monocytic lineage`)],
              by = "sample_submitter_id", all.x = TRUE)
meta <- meta[match(colnames(expr), sample_submitter_id)]
if (!identical(meta$sample_submitter_id, colnames(expr))) stop("TCGA expression/metadata mismatch")
meta[, tcga2014_subtype := factor(tcga2014_subtype,
                                  levels = c("CIN", "GS", "MSI", "EBV"))]
meta[, purity_z_model := as.numeric(scale(absolute_purity))]
meta[, fibroblast_z_model := as.numeric(scale(broad_fibroblasts_LO))]
meta[, monocytic_z_model := as.numeric(scale(broad_monocytic_LO))]

run_limma_tier <- function(tier) {
  if (tier == "purity_adjusted") {
    vars <- c("ecosystem_z", "purity_z_model")
    design_formula <- ~ ecosystem_z + purity_z_model
  } else if (tier == "lineage_subtype_adjusted") {
    vars <- c("ecosystem_z", "purity_z_model", "fibroblast_z_model",
              "monocytic_z_model", "tcga2014_subtype")
    design_formula <- ~ ecosystem_z + purity_z_model + fibroblast_z_model +
      monocytic_z_model + tcga2014_subtype
  } else stop("Unknown tier")
  keep <- complete.cases(meta[, ..vars])
  design <- model.matrix(design_formula, data = meta[keep])
  y <- expr[, keep, drop = FALSE]
  fit <- eBayes(lmFit(y, design), trend = TRUE, robust = TRUE)
  coef_name <- "ecosystem_z"
  tab <- as.data.table(topTable(fit, coef = coef_name, number = Inf,
                               sort.by = "none"), keep.rownames = "gene_symbol")
  setnames(tab, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"),
           c("beta_per_ecosystem_SD", "mean_expression", "moderated_t",
             "p_value", "FDR_BH", "B_statistic"))
  tab[, `:=`(model_tier = tier, n_samples = sum(keep),
             primary_signature_gene = gene_symbol %in% primary_genes)]
  list(table = tab, expression = y, design = design, keep = keep,
       coefficient = which(colnames(design) == coef_name))
}

limma_primary <- run_limma_tier("purity_adjusted")
limma_full <- run_limma_tier("lineage_subtype_adjusted")
gene_stats <- rbindlist(list(limma_primary$table, limma_full$table), fill = TRUE)
fwrite(gene_stats, file.path(out_dir, "TCGA_continuous_ecosystem_limma_gene_statistics.tsv.gz"),
       sep = "\t")

collections <- list(
  Hallmark = list(collection = "H", subcollection = NULL),
  Reactome = list(collection = "C2", subcollection = "CP:REACTOME"),
  KEGG = list(collection = "C2", subcollection = "CP:KEGG_LEGACY"),
  GO_BP = list(collection = "C5", subcollection = "GO:BP")
)

load_sets <- function(spec) {
  if (is.null(spec$subcollection)) {
    m <- msigdbr(species = "Homo sapiens", collection = spec$collection)
  } else {
    m <- msigdbr(species = "Homo sapiens", collection = spec$collection,
                 subcollection = spec$subcollection)
  }
  m <- as.data.table(m)[, .(gene_symbol = unique(gene_symbol)), by = gs_name]
  split(m, by = "gs_name", keep.by = FALSE)
}

run_camera_collection <- function(limma_obj, collection_name, sets, leave_signature_out) {
  genes <- rownames(limma_obj$expression)
  index <- lapply(sets, function(d) {
    g <- unique(d$gene_symbol)
    if (leave_signature_out) g <- setdiff(g, primary_genes)
    which(genes %in% g)
  })
  index <- index[lengths(index) >= 10L & lengths(index) <= 500L]
  ans <- camera(limma_obj$expression, index = index,
                design = limma_obj$design, contrast = limma_obj$coefficient,
                inter.gene.cor = NA, trend.var = TRUE, use.ranks = FALSE,
                sort = FALSE)
  dt <- as.data.table(ans, keep.rownames = "pathway")
  setnames(dt, c("NGenes", "Direction", "PValue", "FDR"),
           c("n_genes", "direction", "p_value", "FDR_BH"))
  dt[, `:=`(collection = collection_name,
            model_tier = unique(limma_obj$table$model_tier),
            leave_primary_signature_out = leave_signature_out)]
  dt
}

pathways <- list()
for (collection_name in names(collections)) {
  message("Load and test ", collection_name)
  sets <- load_sets(collections[[collection_name]])
  pathways[[paste(collection_name, "primary_LO", sep = "_")]] <-
    run_camera_collection(limma_primary, collection_name, sets, TRUE)
  pathways[[paste(collection_name, "full_LO", sep = "_")]] <-
    run_camera_collection(limma_full, collection_name, sets, TRUE)
  pathways[[paste(collection_name, "primary_all", sep = "_")]] <-
    run_camera_collection(limma_primary, collection_name, sets, FALSE)
}
pathway_dt <- rbindlist(pathways, fill = TRUE)
fwrite(pathway_dt, file.path(out_dir, "TCGA_continuous_ecosystem_CAMERA_pathways.tsv.gz"),
       sep = "\t")

# TIDE is a pre-existing model-derived TCGA phenotype. It is explicitly not
# treated as observed immunotherapy response. Map by TCGA case identifier.
tide <- fread(tide_file)
setnames(tide, "sample", "case_submitter_id")
tide_cols <- c("TIDE", "Dysfunction", "Exclusion", "MDSC", "CAF", "TAM M2",
               "IFNG", "CD274", "CD8", "CTL", "MSI Score")
tide_map <- merge(meta, tide[, c("case_submitter_id", tide_cols), with = FALSE],
                  by = "case_submitter_id", all = FALSE)
if (uniqueN(tide_map$sample_submitter_id) != nrow(tide_map)) stop("TIDE duplicate mapping")
fwrite(tide_map[, c("sample_submitter_id", "case_submitter_id", "ecotype4",
                    "sender_primary_z", "receiver_specific_z", "ecosystem_z",
                    tide_cols), with = FALSE],
       file.path(out_dir, "TCGA_ecosystem_TIDE_sample_level.tsv.gz"), sep = "\t")

hc3_term <- function(data, endpoint, tier) {
  needed <- c(endpoint, "ecosystem_z", "purity_z_model")
  if (tier == "lineage_subtype_adjusted") {
    needed <- c(needed, "fibroblast_z_model", "monocytic_z_model", "tcga2014_subtype")
    formula <- as.formula(sprintf("`%s` ~ ecosystem_z + purity_z_model + fibroblast_z_model + monocytic_z_model + tcga2014_subtype", endpoint))
  } else {
    formula <- as.formula(sprintf("`%s` ~ ecosystem_z + purity_z_model", endpoint))
  }
  d <- data[complete.cases(data[, ..needed])]
  fit <- lm(formula, data = d)
  ct <- coeftest(fit, vcov. = vcovHC(fit, type = "HC3"))
  beta <- unname(coef(fit)["ecosystem_z"])
  se <- unname(sqrt(diag(vcovHC(fit, type = "HC3")))["ecosystem_z"])
  data.table(endpoint = endpoint, model_tier = tier, n = nrow(d), beta = beta,
             se_HC3 = se, ci_low = beta - qnorm(0.975) * se,
             ci_high = beta + qnorm(0.975) * se,
             p_value = unname(ct["ecosystem_z", "Pr(>|t|)"]))
}
tide_models <- rbindlist(lapply(tide_cols, function(endpoint) {
  rbindlist(list(hc3_term(tide_map, endpoint, "purity_adjusted"),
                 hc3_term(tide_map, endpoint, "lineage_subtype_adjusted")))
}))
tide_models[, FDR_BH := p.adjust(p_value, "BH"), by = model_tier]
fwrite(tide_models, file.path(out_dir, "TCGA_ecosystem_TIDE_continuous_HC3_models.tsv"),
       sep = "\t")

tide_group_tests <- rbindlist(lapply(tide_cols, function(endpoint) {
  d <- tide_map[!is.na(get(endpoint)) & !is.na(ecotype4)]
  kw <- kruskal.test(d[[endpoint]] ~ d$ecotype4)
  sums <- d[, .(n = .N, median = median(get(endpoint)),
                q1 = quantile(get(endpoint), 0.25), q3 = quantile(get(endpoint), 0.75)),
            by = ecotype4]
  sums[, `:=`(endpoint = endpoint, omnibus_test = "Kruskal-Wallis",
              omnibus_p = kw$p.value)]
  sums
}))
tide_group_tests[, omnibus_FDR_BH := p.adjust(omnibus_p, "BH"), by = ecotype4]
fwrite(tide_group_tests, file.path(out_dir, "TCGA_ecotype4_TIDE_descriptive_tests.tsv"),
       sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  pathway_analysis = list(
    exposure = "continuous frozen ecosystem score per within-cohort SD",
    primary_model = "limma robust trend: expression ~ ecosystem_z + ABSOLUTE purity",
    sensitivity_model = paste("plus leave-signature-out broad fibroblast and monocytic scores",
                              "and official TCGA 2014 subtype"),
    method = "CAMERA competitive gene-set test with estimated inter-gene correlation",
    collections = names(collections),
    primary_display = "all 60 primary sender/receiver signature genes removed from gene sets"
  ),
  TIDE_analysis = list(
    mapped_n = nrow(tide_map),
    provenance = tide_file,
    interpretation = paste(
      "TIDE, Dysfunction, Exclusion and related values are model-derived immune phenotypes.",
      "They are not observed treatment response and do not estimate ICI treatment benefit."
    )
  ),
  input_sha256 = as.list(setNames(vapply(inputs, digest, character(1),
                                         algo = "sha256", file = TRUE), basename(inputs)))
)
write_json(audit, file.path(audit_dir, "TCGA_continuous_ecosystem_pathway_TIDE_v1.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nLimma samples:\n")
print(gene_stats[, .(n_samples = unique(n_samples)), by = model_tier])
cat("\nTop leave-signature-out purity-adjusted CAMERA results:\n")
print(pathway_dt[model_tier == "purity_adjusted" & leave_primary_signature_out == TRUE]
      [order(FDR_BH)][1:30])
cat("\nTIDE mapped samples:", nrow(tide_map), "\n")
print(tide_models[order(model_tier, FDR_BH)][1:22])
