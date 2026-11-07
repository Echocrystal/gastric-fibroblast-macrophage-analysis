#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260726L)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(jsonlite)
  library(digest)
})

processed <- file.path(project, "analysis", "processed", "receiver_pathways")
out_dir <- file.path(project, "analysis", "results", "PLAT_axis", "GSE183904")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expression_file <- file.path(
  processed,
  "GSE183904_primary_fibroblast_myeloid_receiver_program_expression_plat_axis.tsv.gz"
)
codetection_file <- file.path(
  processed,
  "GSE183904_primary_state_receptor_codetection_plat_axis.tsv"
)
program_file <- file.path(project, "resources", "gene_sets", "myeloid_receiver_programs_frozen.tsv")
ecm_file <- file.path(
  source_project, "results", "scrna", "GSE183904",
  "GSE183904_primary_tumor_frozen_ECM_scores.tsv"
)
inputs <- c(expression_file, codetection_file, program_file, ecm_file)
if (!all(file.exists(inputs))) stop("Missing PLAT coarse-replication inputs")

expr <- fread(expression_file)
codetect <- fread(codetection_file)[coarse_label == "Myeloid"]
programs <- fread(program_file)
ecm <- fread(ecm_file)[, .(sample_id, ECM_direction)]

focus_programs <- c(
  "ID_INFLAMMATORY_MONO",
  "ID_MONO_FCN1_VCAN",
  "GOBP_CELL_MATRIX_ADHESION",
  "GOBP_FOCAL_ADHESION_ASSEMBLY",
  "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
  "GOBP_PHAGOCYTOSIS",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
  "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
  "REACTOME_INTEGRIN_SIGNALING"
)
axis_genes <- c("PLAT", "ANXA2", "ITGAM", "ITGB2", "LRP1", "ILK")
programs <- programs[program %chin% focus_programs & !gene_symbol %chin% axis_genes]

safe_z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

normalize_panel <- function(dt, samples) {
  x <- dt[sample_id %chin% samples, .(count_sum = sum(count_sum)), by = .(gene_symbol, sample_id)]
  wide <- dcast(x, gene_symbol ~ sample_id, value.var = "count_sum", fill = 0)
  genes <- wide$gene_symbol
  mat <- as.matrix(wide[, -1])
  rownames(mat) <- genes
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  dge <- calcNormFactors(DGEList(counts = mat), method = "TMM")
  lcpm <- cpm(dge, log = TRUE, prior.count = 0.5)
  melt(as.data.table(lcpm, keep.rownames = "gene_symbol"), id.vars = "gene_symbol",
       variable.name = "sample_id", value.name = "logCPM")
}

partial_test <- function(x, y, covariates, permutations = 9999L) {
  keep <- complete.cases(x, y, covariates)
  if (sum(keep) < 8L) {
    return(list(n = sum(keep), rho = NA_real_, p = NA_real_, partial_rho = NA_real_,
                partial_p = NA_real_, permutation_p = NA_real_))
  }
  xx <- x[keep]; yy <- y[keep]; cc <- as.data.frame(covariates[keep, , drop = FALSE])
  raw <- suppressWarnings(cor.test(xx, yy, method = "spearman", exact = FALSE))
  rx <- resid(lm(rank(xx, ties.method = "average") ~ ., data = cc))
  ry <- resid(lm(rank(yy, ties.method = "average") ~ ., data = cc))
  observed <- cor(rx, ry)
  asymptotic <- suppressWarnings(cor.test(rx, ry, method = "pearson"))$p.value
  null <- replicate(permutations, cor(rx, sample(ry, replace = FALSE)))
  perm_p <- (1 + sum(abs(null) >= abs(observed))) / (permutations + 1)
  list(
    n = sum(keep), rho = unname(raw$estimate), p = raw$p.value,
    partial_rho = observed, partial_p = asymptotic, permutation_p = perm_p
  )
}

loo_summary <- function(x, y) {
  keep <- complete.cases(x, y)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 9L) return(c(min = NA, max = NA, positive_fraction = NA))
  vals <- vapply(seq_along(x), function(i) {
    suppressWarnings(cor(x[-i], y[-i], method = "spearman"))
  }, numeric(1))
  c(min = min(vals), max = max(vals), positive_fraction = mean(vals > 0))
}

cell_counts <- unique(expr[, .(sample_id, Patient_ID, coarse_label, n_cells)])
cell_counts <- dcast(cell_counts, sample_id + Patient_ID ~ coarse_label,
                     value.var = "n_cells", fill = 0)
setnames(cell_counts, c("Fibroblast", "Myeloid"), c("fibroblast_cells", "myeloid_cells"))

variants <- list(
  permissive_5_5 = c(fibroblast = 5L, myeloid = 5L),
  primary_20_20 = c(fibroblast = 20L, myeloid = 20L),
  stringent_30_30 = c(fibroblast = 30L, myeloid = 30L)
)

score_rows <- list()
association_rows <- list()
component_rows <- list()

for (variant in names(variants)) {
  cut <- variants[[variant]]
  eligible <- cell_counts[
    fibroblast_cells >= cut[["fibroblast"]] & myeloid_cells >= cut[["myeloid"]], sample_id
  ]
  fibro <- expr[coarse_label == "Fibroblast"]
  myeloid <- expr[coarse_label == "Myeloid"]
  f_norm <- normalize_panel(fibro, eligible)
  m_norm <- normalize_panel(myeloid, eligible)
  f_wide <- dcast(f_norm, sample_id ~ gene_symbol, value.var = "logCPM")
  m_wide <- dcast(m_norm, sample_id ~ gene_symbol, value.var = "logCPM")
  common <- Reduce(intersect, list(eligible, f_wide$sample_id, m_wide$sample_id))
  f_wide <- f_wide[sample_id %chin% common]; setorder(f_wide, sample_id)
  m_wide <- m_wide[sample_id %chin% common]; setorder(m_wide, sample_id)
  if (!identical(f_wide$sample_id, m_wide$sample_id)) stop("Sample alignment failure")
  if (!"PLAT" %chin% names(f_wide)) stop("PLAT absent from fibroblast pseudobulk")
  if (!all(c("ANXA2", "ITGAM", "ITGB2", "LRP1") %chin% names(m_wide))) {
    stop("Receiver component absent from myeloid pseudobulk")
  }

  receiver_z <- copy(m_wide)
  for (gene in setdiff(names(receiver_z), "sample_id")) {
    set(receiver_z, j = gene, value = safe_z(receiver_z[[gene]]))
  }
  sender_plat_z <- safe_z(f_wide$PLAT)
  receptor_expr <- rowMeans(receiver_z[, .(ANXA2, ITGAM, ITGB2)], na.rm = TRUE)

  cd <- codetect[sample_id %chin% common]
  setorder(cd, sample_id)
  cd <- cd[match(common, sample_id)]
  triple <- asin(sqrt(pmin(pmax(cd$ANXA2_ITGAM_ITGB2_triple_positive_fraction, 0), 1)))
  lrp_context <- asin(sqrt(pmin(pmax(cd$ANXA2_LRP1_double_positive_fraction, 0), 1)))
  triple_z <- safe_z(triple)
  plat_triple_axis <- rowMeans(cbind(sender_plat_z, triple_z), na.rm = TRUE)

  cc <- cell_counts[sample_id %chin% common]
  setorder(cc, sample_id)
  cc <- cc[match(common, sample_id)]
  cov <- merge(data.table(sample_id = common), ecm, by = "sample_id", all.x = TRUE)
  setorder(cov, sample_id)
  cov[, `:=`(
    log_fibroblast_cells = log1p(cc$fibroblast_cells),
    log_myeloid_cells = log1p(cc$myeloid_cells)
  )]

  program_score_list <- list()
  for (program_id in unique(programs$program)) {
    genes <- intersect(programs[program == program_id, gene_symbol], names(receiver_z))
    informative <- genes[vapply(genes, function(g) any(is.finite(receiver_z[[g]])), logical(1))]
    if (length(informative) < 5L) next
    program_score_list[[program_id]] <- rowMeans(receiver_z[, ..informative], na.rm = TRUE)
    score_rows[[paste(variant, program_id, sep = "::")]] <- data.table(
      variant = variant,
      sample_id = common,
      Patient_ID = cc$Patient_ID,
      fibroblast_cells = cc$fibroblast_cells,
      myeloid_cells = cc$myeloid_cells,
      fibroblast_PLAT_logCPM = f_wide$PLAT,
      fibroblast_PLAT_z = sender_plat_z,
      myeloid_receptor_expression_score = receptor_expr,
      myeloid_ANXA2_ITGAM_ITGB2_fraction = cd$ANXA2_ITGAM_ITGB2_triple_positive_fraction,
      myeloid_ANXA2_LRP1_fraction = cd$ANXA2_LRP1_double_positive_fraction,
      PLAT_triple_axis_score = plat_triple_axis,
      program = program_id,
      n_program_genes = length(informative),
      program_score = program_score_list[[program_id]],
      ECM_direction = cov$ECM_direction
    )
  }

  exposures <- list(
    fibroblast_PLAT = sender_plat_z,
    PLAT_triple_axis = plat_triple_axis,
    myeloid_receptor_expression = receptor_expr,
    myeloid_triple_codetection = triple_z
  )
  cov_mat <- as.matrix(cov[, .(ECM_direction, log_fibroblast_cells, log_myeloid_cells)])
  for (program_id in names(program_score_list)) {
    y <- program_score_list[[program_id]]
    for (exposure_id in names(exposures)) {
      x <- exposures[[exposure_id]]
      tst <- partial_test(x, y, cov_mat)
      loo <- loo_summary(x, y)
      association_rows[[paste(variant, program_id, exposure_id, sep = "::")]] <- data.table(
        variant = variant, exposure = exposure_id, program = program_id,
        n = tst$n, rho = tst$rho, p_value = tst$p,
        adjusted_rho = tst$partial_rho, adjusted_p = tst$partial_p,
        adjusted_permutation_p = tst$permutation_p,
        loo_rho_min = loo[["min"]], loo_rho_max = loo[["max"]],
        loo_positive_fraction = loo[["positive_fraction"]]
      )
    }
  }
  for (component_id in c("triple", "ANXA2_LRP1")) {
    y <- if (component_id == "triple") triple else lrp_context
    tst <- partial_test(sender_plat_z, y, cov_mat)
    component_rows[[paste(variant, component_id, sep = "::")]] <- data.table(
      variant = variant, component = component_id, n = tst$n,
      rho = tst$rho, p_value = tst$p, adjusted_rho = tst$partial_rho,
      adjusted_p = tst$partial_p, adjusted_permutation_p = tst$permutation_p
    )
  }
}

scores <- rbindlist(score_rows)
associations <- rbindlist(association_rows)
components <- rbindlist(component_rows)
associations[, raw_FDR_BH := p.adjust(p_value, "BH"), by = .(variant, exposure)]
associations[, adjusted_permutation_FDR_BH := p.adjust(adjusted_permutation_p, "BH"),
             by = .(variant, exposure)]

fwrite(scores, file.path(out_dir, "GSE183904_PLAT_coarse_replication_scores.tsv.gz"), sep = "\t")
fwrite(associations, file.path(out_dir, "GSE183904_PLAT_coarse_replication_associations.tsv"), sep = "\t")
fwrite(components, file.path(out_dir, "GSE183904_PLAT_coarse_component_associations.tsv"), sep = "\t")
fwrite(codetect, file.path(out_dir, "GSE183904_PLAT_receiver_codetection_by_sample.tsv"), sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE183904",
  inference_unit = "one primary-tumor sample per patient",
  variants = variants,
  adjustment = c("ECM_direction", "log1p fibroblast cells", "log1p myeloid cells"),
  partial_test = "Pearson correlation of rank residuals; 9,999 response-residual permutations",
  interpretation = paste(
    "Coarse-lineage replication only. Association does not establish CTHRC1/VSTM2A sender state,",
    "surface receptor localization, PLAT secretion, protease activity, direction, or causality."
  ),
  input_sha256 = as.list(setNames(vapply(inputs, digest, character(1), algo = "sha256", file = TRUE), basename(inputs)))
)
write_json(audit, file.path(audit_dir, "GSE183904_PLAT_coarse_replication.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("Eligible samples by variant:\n")
print(unique(scores[, .(variant, sample_id)])[, .N, by = variant])
cat("\nPrimary 20/20 fibroblast-PLAT associations:\n")
print(associations[variant == "primary_20_20" & exposure == "fibroblast_PLAT"][
  order(adjusted_permutation_p)
])
cat("\nReceiver co-detection summary:\n")
print(codetect[, .(
  samples = .N,
  median_triple_fraction = median(ANXA2_ITGAM_ITGB2_triple_positive_fraction),
  samples_with_triple_positive_cells = sum(ANXA2_ITGAM_ITGB2_triple_positive_cells > 0),
  median_ANXA2_LRP1_fraction = median(ANXA2_LRP1_double_positive_fraction)
)])
