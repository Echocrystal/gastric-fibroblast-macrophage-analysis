#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
set.seed(2026072635L)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(scran)
  library(igraph)
  library(irlba)
  library(uwot)
  library(ggplot2)
  library(jsonlite)
  library(digest)
})

manifest_file <- file.path(
  source_project, "data", "processed", "scrna", "GSE183904",
  "GSE183904_sample_patient_clinical_manifest.tsv"
)
cell_meta_file <- file.path(
  source_project, "data", "processed", "scrna", "GSE183904",
  "GSE183904_secondary_qc_coarse_annotation_metadata.tsv.gz"
)
sce_dir <- file.path(
  source_project, "data", "intermediate", "scrna", "GSE183904", "sample_sce"
)
module_file <- file.path(
  project, "resources", "gene_sets", "GSE183904_fine_state_modules_frozen.tsv"
)
prereg_file <- file.path(
  project, "analysis", "audit", "GSE183904_fine_state_analysis_preregistration.json"
)
out_dir <- file.path(project, "analysis", "results", "PLAT_axis", "GSE183904", "fine_states")
processed_dir <- file.path(project, "analysis", "processed", "GSE183904_fine_states")
figure_dir <- file.path(project, "figures", "PLAT_axis", "GSE183904_fine_states")
audit_dir <- file.path(project, "analysis", "audit")
for (d in c(out_dir, processed_dir, figure_dir, audit_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

manifest <- fread(manifest_file)[
  include_primary_replication == TRUE & specimen_class == "Primary_Tumor"
]
cell_meta <- fread(cell_meta_file)
modules <- fread(module_file)
if (nrow(manifest) != 26L || uniqueN(manifest$Patient_ID) != 26L) {
  stop("Expected 26 patient-independent primary gastric tumors")
}
if (!file.exists(prereg_file)) stop("Missing pre-run fine-state analysis manifest")

axis_genes <- c("PLAT", "ANXA2", "ITGAM", "ITGB2", "LRP1", "ILK")
cell_cycle_genes <- c("MKI67", "TOP2A", "TYMS", "UBE2C", "CENPF")
target_genes <- c(
  axis_genes, "CTHRC1", "VSTM2A", "FAP", "POSTN", "FN1", "COL1A1", "COL1A2",
  "COL3A1", "CXCL12", "CCL2", "CFD", "SPP1", "TREM2", "FCN1", "VCAN",
  "C1QA", "C1QB", "C1QC", "IL1B", "S100A8", "S100A9", "RELA", "NFKB1"
)

safe_scale <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
  z <- (x - mean(x, na.rm = TRUE)) / s
  z[!is.finite(z)] <- 0
  z
}

load_lineage <- function(lineage) {
  mats <- list()
  metas <- list()
  reference_genes <- NULL
  total_alignment_mismatches <- 0L
  for (i in seq_len(nrow(manifest))) {
    info <- manifest[i]
    message(sprintf("[%s %d/26] load %s", lineage, i, info$sample_id))
    sce <- readRDS(file.path(sce_dir, paste0(info$sample_id, "_secondary_qc_sce.rds")))
    cm <- cell_meta[sample_id == info$sample_id]
    idx <- match(colnames(sce), cm$cell_id)
    total_alignment_mismatches <- total_alignment_mismatches + sum(is.na(idx))
    if (anyNA(idx)) stop("Cell alignment failure: ", info$sample_id)
    cm <- cm[idx]
    keep <- which(
      cm$passes_secondary_qc_singlet &
        cm$coarse_label == lineage
    )
    if (!length(keep)) next
    mat <- counts(sce)[, keep, drop = FALSE]
    if (is.null(reference_genes)) {
      reference_genes <- rownames(mat)
    } else if (!identical(reference_genes, rownames(mat))) {
      stop("Gene-order mismatch in ", info$sample_id)
    }
    mats[[info$sample_id]] <- mat
    metas[[info$sample_id]] <- data.table(
      cell_id = colnames(mat),
      sample_id = info$sample_id,
      Patient_ID = info$Patient_ID,
      lineage = lineage,
      sum_counts = cm$sum_counts[keep],
      detected_genes = cm$detected_genes[keep],
      mito_percent = cm$mito_percent[keep],
      secondary_scDblFinder_score = cm$secondary_scDblFinder_score[keep],
      annotation_margin = cm$annotation_margin[keep]
    )
    rm(sce, mat, cm)
    gc(verbose = FALSE)
  }
  if (total_alignment_mismatches != 0L) stop("Non-zero cell alignment mismatches")
  list(counts = do.call(cbind, mats), metadata = rbindlist(metas))
}

within_patient_variance <- function(lognorm, patient) {
  patients <- unique(patient)
  variance_sum <- numeric(nrow(lognorm))
  variance_n <- integer(nrow(lognorm))
  patient_detected <- integer(nrow(lognorm))
  for (id in patients) {
    idx <- which(patient == id)
    if (length(idx) < 2L) next
    x <- lognorm[, idx, drop = FALSE]
    mu <- Matrix::rowMeans(x)
    x2 <- x
    x2@x <- x2@x^2
    v <- (Matrix::rowSums(x2) - length(idx) * mu^2) / (length(idx) - 1L)
    v[v < 0 & v > -1e-10] <- 0
    ok <- is.finite(v)
    variance_sum[ok] <- variance_sum[ok] + v[ok]
    variance_n[ok] <- variance_n[ok] + 1L
    patient_detected <- patient_detected + as.integer(Matrix::rowSums(x > 0) > 0)
  }
  data.table(
    gene_symbol = rownames(lognorm),
    mean_within_patient_variance = variance_sum / pmax(variance_n, 1L),
    contributing_patients = variance_n,
    detected_patients = patient_detected,
    detected_fraction = Matrix::rowMeans(lognorm > 0)
  )
}

module_scores <- function(lognorm, meta, lineage_id) {
  defs <- modules[modules[["lineage"]] == lineage_id]
  split_defs <- split(defs$gene_symbol, defs$module)
  score <- matrix(
    NA_real_, nrow = ncol(lognorm), ncol = length(split_defs),
    dimnames = list(colnames(lognorm), names(split_defs))
  )
  coverage <- list()
  for (module_id in names(split_defs)) {
    present <- intersect(unique(split_defs[[module_id]]), rownames(lognorm))
    coverage[[module_id]] <- data.table(
      lineage = lineage_id, module = module_id,
      requested_genes = length(unique(split_defs[[module_id]])),
      available_genes = length(present)
    )
    if (!length(present)) next
    score[, module_id] <- Matrix::colMeans(lognorm[present, , drop = FALSE])
  }
  centered <- score
  for (id in unique(meta$Patient_ID)) {
    idx <- which(meta$Patient_ID == id)
    centered[idx, ] <- sweep(score[idx, , drop = FALSE], 2L,
                             colMeans(score[idx, , drop = FALSE], na.rm = TRUE), "-")
  }
  centered_z <- apply(centered, 2L, safe_scale)
  if (is.null(dim(centered_z))) {
    centered_z <- matrix(centered_z, ncol = 1L, dimnames = dimnames(centered))
  }
  list(raw = score, centered_z = centered_z, coverage = rbindlist(coverage))
}

cluster_markers <- function(lognorm, cluster, lineage, top_n = 30L) {
  expressed <- Matrix::rowMeans(lognorm > 0) >= 0.01
  x <- lognorm[expressed, , drop = FALSE]
  genes <- rownames(x)
  out <- list()
  for (cl in sort(unique(cluster))) {
    idx <- which(cluster == cl)
    other <- which(cluster != cl)
    mean_in <- Matrix::rowMeans(x[, idx, drop = FALSE])
    mean_out <- Matrix::rowMeans(x[, other, drop = FALSE])
    det_in <- Matrix::rowMeans(x[, idx, drop = FALSE] > 0)
    det_out <- Matrix::rowMeans(x[, other, drop = FALSE] > 0)
    dt <- data.table(
      lineage = lineage, cluster = cl, gene_symbol = genes,
      mean_log_in = mean_in, mean_log_out = mean_out,
      mean_log_difference = mean_in - mean_out,
      detected_fraction_in = det_in, detected_fraction_out = det_out,
      detected_fraction_difference = det_in - det_out
    )
    setorder(dt, -mean_log_difference, -detected_fraction_difference)
    out[[as.character(cl)]] <- head(dt, top_n)
  }
  rbindlist(out)
}

run_lineage <- function(lineage) {
  lineage_id <- lineage
  obj <- load_lineage(lineage)
  counts_mat <- obj$counts
  meta <- obj$metadata
  if (!identical(colnames(counts_mat), meta$cell_id)) stop("Metadata order mismatch")
  lib <- Matrix::colSums(counts_mat)
  lib[!is.finite(lib) | lib <= 0] <- 1
  lognorm <- counts_mat %*% Diagonal(x = 1e4 / lib)
  lognorm@x <- log1p(lognorm@x)

  var_table <- within_patient_variance(lognorm, meta$Patient_ID)
  excluded <- grepl("^(MT-|RPL|RPS|HB)", var_table$gene_symbol) |
    var_table$gene_symbol %chin% c(axis_genes, cell_cycle_genes)
  eligible <- var_table[
    !excluded & detected_fraction >= 0.01 & detected_patients >= 5L &
      is.finite(mean_within_patient_variance)
  ]
  setorder(eligible, -mean_within_patient_variance, gene_symbol)
  hvg <- head(eligible$gene_symbol, 2000L)
  if (length(hvg) < 1000L) stop("Insufficient patient-aware HVGs for ", lineage)

  centered <- as.matrix(lognorm[hvg, , drop = FALSE])
  for (id in unique(meta$Patient_ID)) {
    idx <- which(meta$Patient_ID == id)
    centered[, idx] <- centered[, idx, drop = FALSE] -
      rowMeans(centered[, idx, drop = FALSE])
  }
  gene_sd <- sqrt(rowSums(centered^2) / pmax(ncol(centered) - 1L, 1L))
  gene_sd[!is.finite(gene_sd) | gene_sd <= 0] <- 1
  centered <- centered / gene_sd
  centered[centered > 10] <- 10
  centered[centered < -10] <- -10

  set.seed(2026072635L + match(lineage, c("Fibroblast", "Myeloid")))
  pca <- irlba::prcomp_irlba(t(centered), n = 30L, center = FALSE, scale. = FALSE)$x
  rownames(pca) <- meta$cell_id
  snn <- scran::buildSNNGraph(t(pca), k = 30L, type = "rank")
  set.seed(2026072735L + match(lineage, c("Fibroblast", "Myeloid")))
  community <- igraph::cluster_leiden(
    snn, objective_function = "modularity", weights = igraph::E(snn)$weight,
    resolution = 0.6, n_iterations = 20L
  )
  cluster <- as.integer(community$membership)
  names(cluster) <- igraph::V(snn)$name
  if (!all(meta$cell_id %chin% names(cluster))) {
    # buildSNNGraph may use numeric vertex names; vertex order still follows PCA rows.
    cluster <- as.integer(community$membership)
  } else {
    cluster <- cluster[meta$cell_id]
  }

  set.seed(2026072835L + match(lineage, c("Fibroblast", "Myeloid")))
  umap <- uwot::umap(
    pca, n_neighbors = 30L, min_dist = 0.3, metric = "euclidean",
    n_components = 2L, n_threads = 4L, verbose = TRUE, ret_model = FALSE
  )

  ms <- module_scores(lognorm, meta, lineage)
  module_wide <- data.table(cell_id = meta$cell_id, as.data.table(ms$centered_z))
  module_long <- melt(
    module_wide,
    id.vars = "cell_id", variable.name = "module", value.name = "module_centered_z"
  )
  module_long[, cluster := cluster[match(cell_id, meta$cell_id)]]
  module_cluster <- module_long[, .(
    mean_module_centered_z = mean(module_centered_z, na.rm = TRUE),
    median_module_centered_z = median(module_centered_z, na.rm = TRUE)
  ), by = .(cluster, module)]
  module_cluster <- merge(
    module_cluster,
    unique(modules[modules[["lineage"]] == lineage_id, .(module, role)]),
    by = "module", all.x = TRUE
  )

  cluster_patient <- data.table(
    cluster = cluster, Patient_ID = meta$Patient_ID
  )[, .N, by = .(cluster, Patient_ID)]
  cluster_recurrence <- cluster_patient[, .(
    n_cells = sum(N),
    n_patients = uniqueN(Patient_ID),
    patients_with_at_least_5_cells = sum(N >= 5L),
    max_patient_fraction = max(N) / sum(N)
  ), by = cluster]
  cluster_recurrence[, recurrent_cluster := (
    n_cells >= 50L & patients_with_at_least_5_cells >= 5L &
      max_patient_fraction <= 0.5
  )]

  annotation_scores <- module_cluster[role == "annotation"]
  state_call <- annotation_scores[
    order(cluster, -mean_module_centered_z),
    .SD[1L], by = cluster
  ][, .(cluster, leading_annotation_module = module,
        leading_annotation_score = mean_module_centered_z)]
  second_score <- annotation_scores[
    order(cluster, -mean_module_centered_z),
    .SD[min(2L, .N)], by = cluster
  ][, .(cluster, second_annotation_score = mean_module_centered_z)]
  state_call <- merge(state_call, second_score, by = "cluster", all.x = TRUE)
  state_call[, annotation_margin := leading_annotation_score - second_annotation_score]

  quality_scores <- module_cluster[role == "quality_review"][
    order(cluster, -mean_module_centered_z),
    .SD[1L], by = cluster
  ][, .(cluster, leading_quality_module = module,
        leading_quality_score = mean_module_centered_z)]
  state_call <- merge(state_call, quality_scores, by = "cluster", all.x = TRUE)
  state_call <- merge(state_call, cluster_recurrence, by = "cluster", all.x = TRUE)
  state_call[, state_label := leading_annotation_module]
  state_call[
    is.finite(leading_quality_score) & leading_quality_score > 0.5 &
      leading_quality_score > leading_annotation_score,
    state_label := leading_quality_module
  ]
  state_call[leading_annotation_score < 0 | annotation_margin < 0.05,
             state_label := "Unresolved_mixed"]
  state_call[recurrent_cluster == FALSE,
             state_label := paste0("Donor_restricted__", state_label)]

  target_present <- intersect(target_genes, rownames(lognorm))
  target_norm <- t(as.matrix(lognorm[target_present, , drop = FALSE]))
  colnames(target_norm) <- paste0("expr_", target_present)
  target_detect <- t(as.matrix(counts_mat[target_present, , drop = FALSE] > 0))
  colnames(target_detect) <- paste0("detect_", target_present)

  cell_out <- cbind(
    meta,
    data.table(
      cluster = cluster,
      state_label = state_call$state_label[match(cluster, state_call$cluster)],
      recurrent_cluster = state_call$recurrent_cluster[match(cluster, state_call$cluster)],
      UMAP1 = umap[, 1L], UMAP2 = umap[, 2L]
    ),
    as.data.table(target_norm),
    as.data.table(target_detect)
  )
  for (module_id in colnames(ms$centered_z)) {
    cell_out[[paste0("module_", module_id, "_centered_z")]] <- ms$centered_z[, module_id]
  }
  if (lineage == "Myeloid") {
    cell_out[, ANXA2_ITGAM_ITGB2_triple_detected := (
      detect_ANXA2 & detect_ITGAM & detect_ITGB2
    )]
    cell_out[, ANXA2_LRP1_double_detected := detect_ANXA2 & detect_LRP1]
  }

  summary_rows <- list()
  for (cl in sort(unique(cluster))) {
    idx <- which(cluster == cl)
    base <- state_call[cluster == cl]
    extra <- data.table()
    if (lineage == "Fibroblast") {
      extra <- data.table(
        PLAT_mean_log1p_CP10K = mean(cell_out$expr_PLAT[idx]),
        PLAT_detected_fraction = mean(cell_out$detect_PLAT[idx]),
        CTHRC1_mean_log1p_CP10K = mean(cell_out$expr_CTHRC1[idx]),
        VSTM2A_mean_log1p_CP10K = mean(cell_out$expr_VSTM2A[idx])
      )
    } else {
      extra <- data.table(
        ANXA2_mean_log1p_CP10K = mean(cell_out$expr_ANXA2[idx]),
        ITGAM_detected_fraction = mean(cell_out$detect_ITGAM[idx]),
        ITGB2_detected_fraction = mean(cell_out$detect_ITGB2[idx]),
        ANXA2_ITGAM_ITGB2_triple_fraction = mean(
          cell_out$ANXA2_ITGAM_ITGB2_triple_detected[idx]
        ),
        ANXA2_LRP1_double_fraction = mean(cell_out$ANXA2_LRP1_double_detected[idx])
      )
    }
    summary_rows[[as.character(cl)]] <- cbind(
      data.table(lineage = lineage), base, extra
    )
  }
  cluster_summary <- rbindlist(summary_rows, fill = TRUE)
  markers <- cluster_markers(lognorm, cluster, lineage)

  hvg_out <- merge(var_table, data.table(gene_symbol = hvg, selected_HVG = TRUE),
                   by = "gene_symbol", all.x = TRUE)
  hvg_out[is.na(selected_HVG), selected_HVG := FALSE]
  pca_out <- data.table(cell_id = meta$cell_id, as.data.table(pca))

  fwrite(cell_out, file.path(processed_dir, paste0("GSE183904_", lineage,
                                                   "_fine_state_cell_metadata.tsv.gz")), sep = "\t")
  fwrite(cluster_summary, file.path(out_dir, paste0("GSE183904_", lineage,
                                                    "_fine_state_cluster_summary.tsv")), sep = "\t")
  fwrite(module_cluster, file.path(out_dir, paste0("GSE183904_", lineage,
                                                   "_fine_state_module_by_cluster.tsv")), sep = "\t")
  fwrite(markers, file.path(out_dir, paste0("GSE183904_", lineage,
                                            "_fine_state_top_markers.tsv")), sep = "\t")
  fwrite(hvg_out, file.path(out_dir, paste0("GSE183904_", lineage,
                                            "_patient_aware_HVGs.tsv")), sep = "\t")
  fwrite(pca_out, file.path(processed_dir, paste0("GSE183904_", lineage,
                                                  "_patient_centered_PCA.tsv.gz")), sep = "\t")
  fwrite(ms$coverage, file.path(out_dir, paste0("GSE183904_", lineage,
                                                "_fine_state_module_coverage.tsv")), sep = "\t")

  plot_dt <- cell_out[, .(UMAP1, UMAP2, cluster = factor(cluster), state_label,
                          Patient_ID, recurrent_cluster)]
  p_cluster <- ggplot(plot_dt, aes(UMAP1, UMAP2, color = cluster)) +
    geom_point(size = 0.22, alpha = 0.65) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 2)) +
    labs(title = paste("GSE183904", lineage, "patient-centered fine states"),
         color = "Cluster") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold"))
  p_state <- ggplot(plot_dt, aes(UMAP1, UMAP2, color = state_label)) +
    geom_point(size = 0.22, alpha = 0.65) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1)) +
    labs(title = paste("GSE183904", lineage, "conservative state labels"),
         color = "State") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(figure_dir, paste0("GSE183904_", lineage, "_fine_state_clusters.pdf")),
         p_cluster, width = 8.5, height = 6.5, device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0("GSE183904_", lineage, "_fine_state_clusters.png")),
         p_cluster, width = 8.5, height = 6.5, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, paste0("GSE183904_", lineage, "_fine_state_labels.pdf")),
         p_state, width = 9.5, height = 6.5, device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0("GSE183904_", lineage, "_fine_state_labels.png")),
         p_state, width = 9.5, height = 6.5, dpi = 300, bg = "white")

  list(
    cells = nrow(cell_out), patients = uniqueN(cell_out$Patient_ID),
    clusters = uniqueN(cell_out$cluster),
    recurrent_clusters = sum(cluster_summary$recurrent_cluster),
    cluster_summary = cluster_summary
  )
}

fib <- run_lineage("Fibroblast")
gc(verbose = FALSE)
myeloid <- run_lineage("Myeloid")
gc(verbose = FALSE)

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  dataset = "GSE183904",
  samples = 26L,
  patients = 26L,
  fibroblast_cells = fib$cells,
  myeloid_cells = myeloid$cells,
  fibroblast_clusters = fib$clusters,
  myeloid_clusters = myeloid$clusters,
  fibroblast_recurrent_clusters = fib$recurrent_clusters,
  myeloid_recurrent_clusters = myeloid$recurrent_clusters,
  clustering_features_exclude_axis_genes = axis_genes,
  inference_note = "UMAP and clusters are cell-level descriptive outputs; all downstream inference uses patients.",
  ambient_RNA_note = paste(
    "Empty droplets were not deposited. Ambient-RNA estimation cannot be reconstructed from",
    "the author-filtered matrices and remains a source-data limitation."
  ),
  input_sha256 = as.list(setNames(
    vapply(c(manifest_file, cell_meta_file, module_file, prereg_file),
           digest, character(1), algo = "sha256", file = TRUE),
    basename(c(manifest_file, cell_meta_file, module_file, prereg_file))
  ))
)
write_json(
  audit, file.path(audit_dir, "GSE183904_fine_state_subclustering.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA
)

cat("\nFibroblast cluster summary:\n")
print(fib$cluster_summary[order(-PLAT_mean_log1p_CP10K)])
cat("\nMyeloid cluster summary:\n")
print(myeloid$cluster_summary[order(-ANXA2_ITGAM_ITGB2_triple_fraction)])
