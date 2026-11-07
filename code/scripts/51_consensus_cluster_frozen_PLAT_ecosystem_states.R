suppressPackageStartupMessages({
  library(data.table)
  library(ConsensusClusterPlus)
  library(cluster)
  library(jsonlite)
  library(digest)
})

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
tcga_file <- "I:/ECM-new/ECM-public/data/processed/bulk/tcga_stad_full_symbol_log2cpm_412.rds"
signature_file <- file.path(project, "resources", "gene_sets",
                            "GSE183904_PLAT_ecosystem_bulk_signatures.tsv")
score_file <- file.path(project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
                        "frozen_state_projection",
                        "TCGA_STAD_PLAT_ecosystem_scores_clinical.tsv.gz")
out_dir <- file.path(project, "analysis", "results", "PLAT_axis", "TCGA_STAD",
                     "bulk_ecological_subtypes_v1")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(tcga_file, signature_file, score_file)
if (!all(file.exists(inputs))) stop("Missing consensus-clustering input")

sender_id <- "sender_CXCL14_POSTN_state_axis_excluded"
receiver_id <- "receiver_C1QC_APOE_specific_axis_excluded"
sig <- fread(signature_file)
primary <- sig[signature_id %in% c(sender_id, receiver_id)]
if (primary[, uniqueN(signature_id)] != 2L) stop("Primary signatures absent")
if (primary[, .N, by = signature_id][, any(N != 30L)]) stop("Expected 30 genes per signature")
feature_genes <- unique(primary$gene_symbol)

tcga <- readRDS(tcga_file)
expr <- tcga$expression
if (!is.matrix(expr) || ncol(expr) != 412L) stop("Unexpected TCGA expression object")
available <- intersect(feature_genes, rownames(expr))
if (length(available) < 55L) stop("Insufficient frozen feature coverage")
expr_frozen <- expr[available, , drop = FALSE]
keep_var <- apply(expr_frozen, 1L, sd, na.rm = TRUE) > 0
expr_frozen <- expr_frozen[keep_var, , drop = FALSE]
z <- t(scale(t(expr_frozen)))
z[!is.finite(z)] <- 0

# ConsensusClusterPlus expects features in rows and samples in columns.
# Primary: PAM with 1-Spearman distance. Sensitivity: HC with Euclidean/Ward.D2.
consensus_reps <- as.integer(Sys.getenv("CC_REPS", "1000"))
if (!is.finite(consensus_reps) || consensus_reps < 20L) stop("Invalid CC_REPS")
set.seed(2026080151)
primary_dir <- file.path(out_dir, "consensus_primary_pam_spearman")
sensitivity_dir <- file.path(out_dir, "consensus_sensitivity_hc_wardD2")
cc_primary <- ConsensusClusterPlus(
  z, maxK = 6, reps = consensus_reps, pItem = 0.8, pFeature = 1,
  clusterAlg = "pam", distance = "spearman", seed = 2026080151,
  title = primary_dir, plot = "png", writeTable = FALSE, verbose = FALSE,
  corUse = "pairwise.complete.obs"
)
cc_sensitivity <- ConsensusClusterPlus(
  z, maxK = 6, reps = consensus_reps, pItem = 0.8, pFeature = 1,
  clusterAlg = "hc", distance = "euclidean",
  innerLinkage = "ward.D2", finalLinkage = "ward.D2", seed = 2026080152,
  title = sensitivity_dir, plot = "png", writeTable = FALSE, verbose = FALSE
)

adjusted_rand <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(v) v * (v - 1) / 2
  n <- sum(tab)
  if (n < 2L) return(NA_real_)
  index <- sum(choose2(tab))
  a <- sum(choose2(rowSums(tab)))
  b <- sum(choose2(colSums(tab)))
  expected <- a * b / choose2(n)
  max_index <- (a + b) / 2
  if (max_index == expected) return(1)
  (index - expected) / (max_index - expected)
}

cdf_area <- function(cm, grid = seq(0, 1, length.out = 1001L)) {
  v <- cm[upper.tri(cm)]
  f <- vapply(grid, function(t) mean(v <= t), numeric(1))
  sum(diff(grid) * (head(f, -1L) + tail(f, -1L)) / 2)
}

sample_stability <- function(cm, cls, method, k) {
  ids <- colnames(cm)
  if (is.null(ids)) ids <- names(cls)
  if (is.null(ids) || length(ids) != length(cls)) stop("Consensus sample IDs unavailable")
  rbindlist(lapply(seq_along(ids), function(i) {
    same <- which(cls == cls[i] & seq_along(cls) != i)
    other_groups <- setdiff(sort(unique(cls)), cls[i])
    within <- if (length(same)) mean(cm[i, same]) else NA_real_
    other_means <- vapply(other_groups, function(g) mean(cm[i, cls == g]), numeric(1))
    best_other <- if (length(other_means)) max(other_means) else NA_real_
    data.table(sample_id = ids[i], k = k, method = method,
               cluster = as.integer(cls[i]), within_consensus = within,
               best_other_consensus = best_other,
               consensus_margin = within - best_other)
  }))
}

sample_distance <- as.dist(1 - cor(z, method = "spearman", use = "pairwise.complete.obs"))
metrics <- list()
stability_all <- list()
areas <- numeric(6L)
for (k in 2:6) {
  cp <- cc_primary[[k]]
  cs <- cc_sensitivity[[k]]
  cls_p <- cp$consensusClass
  cls_s <- cs$consensusClass
  names(cls_p) <- colnames(z)
  names(cls_s) <- colnames(z)
  if (is.null(colnames(cp$consensusMatrix))) {
    dimnames(cp$consensusMatrix) <- list(colnames(z), colnames(z))
  }
  if (is.null(colnames(cs$consensusMatrix))) {
    dimnames(cs$consensusMatrix) <- list(colnames(z), colnames(z))
  }
  cls_s <- cls_s[names(cls_p)]
  sil <- silhouette(cls_p, sample_distance)
  cm <- cp$consensusMatrix
  v <- cm[upper.tri(cm)]
  areas[k] <- cdf_area(cm)
  st <- sample_stability(cm, cls_p, "primary_pam_spearman", k)
  stability_all[[paste0("primary_", k)]] <- st
  stability_all[[paste0("sensitivity_", k)]] <-
    sample_stability(cs$consensusMatrix, cc_sensitivity[[k]]$consensusClass,
                     "sensitivity_hc_wardD2", k)
  cluster_sizes <- sort(table(cls_p))
  within_pair <- outer(cls_p, cls_p, FUN = "==") & upper.tri(cm)
  metrics[[as.character(k)]] <- data.table(
    k = k,
    PAC_0.1_0.9 = mean(v > 0.1 & v < 0.9),
    CDF_area = areas[k],
    mean_silhouette = mean(sil[, "sil_width"]),
    median_silhouette = median(sil[, "sil_width"]),
    negative_silhouette_fraction = mean(sil[, "sil_width"] < 0),
    mean_within_cluster_consensus = mean(cm[within_pair]),
    minimum_cluster_n = min(cluster_sizes),
    minimum_cluster_fraction = min(cluster_sizes) / ncol(z),
    maximum_cluster_n = max(cluster_sizes),
    primary_vs_sensitivity_ARI = adjusted_rand(cls_p, cls_s),
    cluster_sizes = paste(names(cluster_sizes), cluster_sizes, sep = ":", collapse = ";")
  )
}
metrics_dt <- rbindlist(metrics)
metrics_dt[, CDF_delta_area := CDF_area - shift(CDF_area)]

# Predefined, outcome-blind selection gate and composite ranking.
metrics_dt[, gate_eligible := minimum_cluster_fraction >= 0.10 &
             PAC_0.1_0.9 <= 0.25 & mean_silhouette >= 0.25 &
             primary_vs_sensitivity_ARI >= 0.50]
metrics_dt[, composite_rank :=
             frank(PAC_0.1_0.9, ties.method = "average") +
             frank(-mean_silhouette, ties.method = "average") +
             frank(-mean_within_cluster_consensus, ties.method = "average") +
             frank(-primary_vs_sensitivity_ARI, ties.method = "average")]
eligible <- metrics_dt[gate_eligible == TRUE][order(composite_rank, k)]
gate_passed <- nrow(eligible) > 0L
selected_k <- if (gate_passed) eligible$k[1L] else NA_integer_

fwrite(metrics_dt, file.path(out_dir, "TCGA_frozen_state_consensus_k_metrics.tsv"), sep = "\t")
fwrite(rbindlist(stability_all),
       file.path(out_dir, "TCGA_frozen_state_consensus_sample_stability_all_k.tsv"), sep = "\t")

if (gate_passed) {
  cls <- cc_primary[[selected_k]]$consensusClass
  names(cls) <- colnames(z)
  cls <- cls[colnames(z)]
  if (is.null(colnames(cc_primary[[selected_k]]$consensusMatrix))) {
    dimnames(cc_primary[[selected_k]]$consensusMatrix) <- list(colnames(z), colnames(z))
  }
  score_dt <- fread(score_file)
  assignments <- data.table(sample_submitter_id = names(cls),
                            ecological_cluster_numeric = as.integer(cls))
  assignments[, ecological_cluster := sprintf("EC%d", ecological_cluster_numeric)]
  st_sel <- sample_stability(cc_primary[[selected_k]]$consensusMatrix,
                             cc_primary[[selected_k]]$consensusClass,
                             "primary_pam_spearman", selected_k)
  setnames(st_sel, "sample_id", "sample_submitter_id")
  assignments <- merge(assignments, st_sel[, .(sample_submitter_id,
                                               within_consensus,
                                               best_other_consensus,
                                               consensus_margin)],
                       by = "sample_submitter_id", all.x = TRUE)
  assignments <- merge(assignments,
                       score_dt[, .(sample_submitter_id, sender_primary_z,
                                    receiver_specific_z, ecosystem_z,
                                    absolute_purity, tcga2014_subtype,
                                    stage_group_frozen, os_months, os_event_frozen)],
                       by = "sample_submitter_id", all.x = TRUE)
  cluster_summary <- assignments[, .(
    n = .N,
    sender_mean = mean(sender_primary_z), sender_median = median(sender_primary_z),
    receiver_mean = mean(receiver_specific_z), receiver_median = median(receiver_specific_z),
    ecosystem_mean = mean(ecosystem_z), ecosystem_median = median(ecosystem_z),
    median_consensus_margin = median(consensus_margin),
    minimum_consensus_margin = min(consensus_margin)
  ), by = ecological_cluster][order(ecological_cluster)]

  # Neutral EC labels are frozen here. Biological names are added only after
  # inspecting sender/receiver profiles, never using outcome or treatment data.
  centroids <- sapply(sort(unique(assignments$ecological_cluster)), function(cl) {
    rowMeans(z[, assignments[ecological_cluster == cl,
                             match(sample_submitter_id, colnames(z))], drop = FALSE])
  })
  colnames(centroids) <- sort(unique(assignments$ecological_cluster))
  centroid_dt <- data.table(gene_symbol = rownames(centroids), centroids)
  fwrite(assignments, file.path(out_dir, "TCGA_frozen_state_ecological_cluster_assignments.tsv"), sep = "\t")
  fwrite(cluster_summary, file.path(out_dir, "TCGA_frozen_state_ecological_cluster_summary.tsv"), sep = "\t")
  fwrite(centroid_dt, file.path(out_dir, "TCGA_frozen_state_ecological_cluster_centroids.tsv"), sep = "\t")
  saveRDS(list(selected_k = selected_k, features = rownames(z), z_expression = z,
               primary = cc_primary, sensitivity = cc_sensitivity,
               assignments = assignments, centroids = centroids,
               metrics = metrics_dt),
          file.path(out_dir, "TCGA_frozen_state_consensus_clustering.rds"))
}

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  analysis = "Outcome-blind consensus clustering of frozen single-cell-derived state genes",
  feature_definition = list(sender = sender_id, receiver = receiver_id,
                            frozen_union_n = length(feature_genes),
                            available_n = nrow(z), genes = rownames(z)),
  primary_method = list(algorithm = "PAM", distance = "1-Spearman",
                        reps = consensus_reps, pItem = 0.8, pFeature = 1,
                        k_range = "2-6", seed = 2026080151),
  sensitivity_method = list(algorithm = "hierarchical clustering",
                            distance = "Euclidean", linkage = "Ward.D2",
                            reps = consensus_reps, pItem = 0.8, pFeature = 1,
                            k_range = "2-6", seed = 2026080152),
  gate = list(minimum_cluster_fraction = 0.10, maximum_PAC = 0.25,
              minimum_mean_silhouette = 0.25, minimum_algorithm_ARI = 0.50,
              passed = gate_passed, selected_k = selected_k),
  outcome_blinding = paste(
    "No overall-survival, TIDE, Kim2018 response or other clinical endpoint",
    "was used in feature selection, clustering, gate evaluation or k selection."
  ),
  input_sha256 = as.list(setNames(vapply(inputs, digest, character(1),
                                         algo = "sha256", file = TRUE), basename(inputs)))
)
write_json(audit, file.path(audit_dir, "TCGA_frozen_state_consensus_clustering_v1.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat("\nConsensus clustering metrics:\n")
print(metrics_dt)
cat("\nGate passed:", gate_passed, " selected k:", selected_k, "\n")
if (gate_passed) print(cluster_summary)
