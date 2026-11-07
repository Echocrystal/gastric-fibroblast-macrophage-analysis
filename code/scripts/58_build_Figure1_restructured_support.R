#!/usr/bin/env Rscript

# Build the additional non-PLAT source data required by the restructured
# Figure 1. All calculations are R-only. CellChat pooled probabilities,
# donor-level ligand/receptor expression support and cross-cohort normalized
# receiver ranks are kept as separate estimands.

options(stringsAsFactors = FALSE, warn = 1)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
figure_lib <- file.path(project, "runtime", "R_figure")
cellchat_lib <- "I:/ECM-new/scRNA/20.rebuild_GSE163558/R_libs_4.4"
.libPaths(c(figure_lib, cellchat_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(hdf5r)
  library(Matrix)
})

source_dir <- file.path(project, "figures", "data", "Figure1_GSE206785_unbiased")
analysis_dir <- file.path(project, "analysis", "results", "figure1_gse206785_unbiased")
validation_dir <- file.path(project, "analysis", "results", "communication", "Figure7_redesign")
audit_dir <- file.path(project, "analysis", "audit")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

input_h5ad <- file.path(analysis_dir, "GSE206785_donor_balanced_all_lineage.h5ad")
network_file <- file.path(source_dir, "panel_e_cellchat_broad_network.tsv")
interaction_file <- file.path(source_dir, "panel_e_cellchat_interactions.tsv.gz")
discovery_receiver_file <- file.path(source_dir, "panel_e_cellchat_fibroblast_receiver_summary.tsv")
validation_cellchat_file <- file.path(validation_dir, "GSE183904_global_CellChat_scores.tsv")

required <- c(input_h5ad, network_file, interaction_file,
              discovery_receiver_file, validation_cellchat_file)
if (any(!file.exists(required))) {
  stop("Missing input: ", paste(required[!file.exists(required)], collapse = ", "))
}

read_h5ad_csr <- function(path) {
  f <- H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)

  if (!all(c("data", "indices", "indptr") %in% names(f[["X"]]))) {
    stop("The h5ad X slot is not CSR encoded")
  }
  data <- f[["X/data"]][]
  indices <- as.integer(f[["X/indices"]][])
  indptr <- as.integer(f[["X/indptr"]][])
  shape <- as.integer(h5attr(f[["X"]], "shape"))
  x_cell_gene <- new(
    "dgRMatrix", p = indptr, j = indices, x = as.numeric(data), Dim = shape
  )

  obs_index_name <- as.character(h5attr(f[["obs"]], "_index"))
  var_index_name <- as.character(h5attr(f[["var"]], "_index"))
  cell_ids <- as.character(f[[paste0("obs/", obs_index_name)]][])
  genes <- as.character(f[[paste0("var/", var_index_name)]][])

  read_obs_column <- function(name) {
    node <- f[[paste0("obs/", name)]]
    if (inherits(node, "H5Group")) {
      categories <- as.character(node[["categories"]][])
      codes <- as.integer(node[["codes"]][])
      out <- rep(NA_character_, length(codes))
      valid <- codes >= 0L
      out[valid] <- categories[codes[valid] + 1L]
      out
    } else {
      as.character(node[])
    }
  }

  meta <- data.table(
    cell_id = cell_ids,
    cell_type = read_obs_column("cell_type"),
    Patient = read_obs_column("Patient")
  )
  expr <- as(t(x_cell_gene), "CsparseMatrix")
  rownames(expr) <- genes
  colnames(expr) <- cell_ids
  list(expr = expr, meta = meta)
}

# -------------------------------------------------------------------------
# Panel c source: broad signaling roles from the same frozen CellChat matrix.
# Autocrine edges remain included, matching CellChat role summaries.
# -------------------------------------------------------------------------

network <- fread(network_file)
outgoing <- network[, .(
  outgoing_count = sum(interaction_count),
  outgoing_weight = sum(interaction_weight)
), by = .(cell_type = source)]
incoming <- network[, .(
  incoming_count = sum(interaction_count),
  incoming_weight = sum(interaction_weight)
), by = .(cell_type = target)]
roles <- merge(outgoing, incoming, by = "cell_type", all = TRUE)
roles[is.na(roles)] <- 0
roles[, outgoing_weight_rank := frank(-outgoing_weight, ties.method = "min")]
roles[, incoming_weight_rank := frank(-incoming_weight, ties.method = "min")]
setorder(roles, outgoing_weight_rank, incoming_weight_rank)
fwrite(roles, file.path(source_dir, "panel_c_cellchat_signaling_roles.tsv"), sep = "\t")

# -------------------------------------------------------------------------
# Panel f source: donor-level expression support for the frozen pooled
# Fibroblast -> immune CellChat interactions. This is not a donor CellChat
# probability. A pair is supported when all ligand and receptor components are
# detected in >=10% of the corresponding donor-level cell populations.
# -------------------------------------------------------------------------

message("Reading donor-balanced GSE206785 h5ad")
obj <- read_h5ad_csr(input_h5ad)
immune_targets <- c("B cells", "CD4+ T", "CD8+ T", "NK/ILC", "Mast", "Myeloid", "Plasma")
comm <- fread(interaction_file)
comm <- comm[
  source == "Fibroblast" & target %chin% immune_targets & prob > 0 & pval < 0.05
]
if (!nrow(comm)) stop("No frozen Fibroblast -> immune CellChat interactions")
comm[, pair_id := paste(target, interaction_name, sep = "::")]

split_complex <- function(x, available_genes) {
  alias_components <- list(
    "CD94:NKG2C" = c("KLRD1", "KLRC2")
  )
  if (x %in% names(alias_components)) {
    return(alias_components[[x]][alias_components[[x]] %chin% available_genes])
  }
  parts <- strsplit(as.character(x), "_", fixed = TRUE)[[1L]]
  parts[parts %chin% available_genes]
}

available_genes <- rownames(obj$expr)
ligand_components <- lapply(comm$ligand, split_complex, available_genes = available_genes)
receptor_components <- lapply(comm$receptor, split_complex, available_genes = available_genes)
if (any(lengths(ligand_components) == 0L) || any(lengths(receptor_components) == 0L)) {
  missing_ligands <- unique(comm$ligand[lengths(ligand_components) == 0L])
  missing_receptors <- unique(comm$receptor[lengths(receptor_components) == 0L])
  stop(
    "At least one frozen CellChat complex has no expression-matrix component. Ligands: ",
    paste(missing_ligands, collapse = ", "), "; receptors: ",
    paste(missing_receptors, collapse = ", ")
  )
}
comm[, ligand_components := vapply(ligand_components, paste, collapse = ";", character(1))]
comm[, receptor_components := vapply(receptor_components, paste, collapse = ";", character(1))]

genes_needed <- sort(unique(c(
  unlist(strsplit(comm$ligand_components, ";", fixed = TRUE)),
  unlist(strsplit(comm$receptor_components, ";", fixed = TRUE))
)))
expr <- obj$expr[genes_needed, , drop = FALSE]
meta <- copy(obj$meta)

cell_counts <- meta[, .N, by = .(Patient, cell_type)]
setnames(cell_counts, "N", "n_cells")
eligible <- merge(
  cell_counts[cell_type == "Fibroblast", .(Patient, n_fibroblast = n_cells)],
  cell_counts[cell_type %chin% immune_targets,
              .(Patient, target = cell_type, n_receiver = n_cells)],
  by = "Patient", all = FALSE
)
eligible <- eligible[n_fibroblast >= 20L & n_receiver >= 20L]
if (!nrow(eligible)) stop("No donor/receiver pair passes the >=20-cell gate")

groups_needed <- unique(rbindlist(list(
  eligible[, .(Patient, cell_type = "Fibroblast")],
  eligible[, .(Patient, cell_type = target)]
)))

metric_list <- vector("list", nrow(groups_needed))
for (i in seq_len(nrow(groups_needed))) {
  p <- groups_needed$Patient[i]
  g <- groups_needed$cell_type[i]
  idx <- which(meta$Patient == p & meta$cell_type == g)
  x <- expr[, idx, drop = FALSE]
  metric_list[[i]] <- data.table(
    Patient = p,
    cell_type = g,
    gene = rownames(expr),
    mean_expression = as.numeric(Matrix::rowMeans(x)),
    detection_fraction = as.numeric(Matrix::rowMeans(x > 0)),
    n_cells = length(idx)
  )
}
gene_metrics <- rbindlist(metric_list)

metric_lookup <- new.env(parent = emptyenv())
for (i in seq_len(nrow(gene_metrics))) {
  key <- paste(gene_metrics$Patient[i], gene_metrics$cell_type[i],
               gene_metrics$gene[i], sep = "||")
  metric_lookup[[key]] <- c(
    mean = gene_metrics$mean_expression[i],
    detection = gene_metrics$detection_fraction[i]
  )
}

lookup_components <- function(patient, cell_type, component_string) {
  components <- strsplit(component_string, ";", fixed = TRUE)[[1L]]
  vals <- lapply(components, function(gene) {
    metric_lookup[[paste(patient, cell_type, gene, sep = "||")]]
  })
  if (any(vapply(vals, is.null, logical(1)))) {
    return(c(min_mean = NA_real_, min_detection = NA_real_))
  }
  mat <- do.call(rbind, vals)
  c(min_mean = min(mat[, "mean"]), min_detection = min(mat[, "detection"]))
}

support_rows <- vector("list", nrow(eligible))
for (i in seq_len(nrow(eligible))) {
  p <- eligible$Patient[i]
  target_value <- eligible$target[i]
  pairs <- comm[target == target_value]
  z <- copy(pairs)
  ligand_values <- t(vapply(
    z$ligand_components,
    function(s) lookup_components(p, "Fibroblast", s),
    numeric(2)
  ))
  receptor_values <- t(vapply(
    z$receptor_components,
    function(s) lookup_components(p, target_value, s),
    numeric(2)
  ))
  z[, `:=`(
    Patient = p,
    n_fibroblast = eligible$n_fibroblast[i],
    n_receiver = eligible$n_receiver[i],
    ligand_min_mean = ligand_values[, "min_mean"],
    ligand_min_detection = ligand_values[, "min_detection"],
    receptor_min_mean = receptor_values[, "min_mean"],
    receptor_min_detection = receptor_values[, "min_detection"]
  )]
  z[, expression_gate_support :=
      ligand_min_detection >= 0.10 & receptor_min_detection >= 0.10]
  z[, donor_expression_potential := sqrt(pmax(ligand_min_mean, 0) *
                                         pmax(receptor_min_mean, 0))]
  support_rows[[i]] <- z
}
support_long <- rbindlist(support_rows, fill = TRUE)
setcolorder(support_long, c(
  "Patient", "source", "target", "ligand", "receptor", "pathway_name",
  "prob", "pval", "n_fibroblast", "n_receiver",
  "ligand_min_detection", "receptor_min_detection",
  "expression_gate_support", "donor_expression_potential"
))
fwrite(
  support_long,
  file.path(source_dir, "panel_f_donor_cellchat_pair_support.tsv.gz"),
  sep = "\t"
)

donor_summary <- support_long[, .(
  n_frozen_pairs = .N,
  supported_pairs = sum(expression_gate_support),
  pair_support_fraction = mean(expression_gate_support),
  pooled_weight_supported_fraction = sum(prob * expression_gate_support) / sum(prob),
  median_expression_potential = median(donor_expression_potential[expression_gate_support],
                                       na.rm = TRUE),
  n_fibroblast = unique(n_fibroblast),
  n_receiver = unique(n_receiver)
), by = .(Patient, target)]
donor_summary[!is.finite(median_expression_potential), median_expression_potential := NA_real_]
fwrite(
  donor_summary,
  file.path(source_dir, "panel_f_donor_receiver_support.tsv"),
  sep = "\t"
)

receiver_summary <- donor_summary[, .(
  eligible_donors = .N,
  median_pair_support_fraction = median(pair_support_fraction),
  q1_pair_support_fraction = quantile(pair_support_fraction, 0.25),
  q3_pair_support_fraction = quantile(pair_support_fraction, 0.75),
  median_weighted_support = median(pooled_weight_supported_fraction),
  q1_weighted_support = quantile(pooled_weight_supported_fraction, 0.25),
  q3_weighted_support = quantile(pooled_weight_supported_fraction, 0.75),
  minimum_weighted_support = min(pooled_weight_supported_fraction),
  maximum_weighted_support = max(pooled_weight_supported_fraction)
), by = target]
receiver_summary[, weighted_support_rank := frank(-median_weighted_support, ties.method = "min")]
setorder(receiver_summary, weighted_support_rank)
fwrite(
  receiver_summary,
  file.path(source_dir, "panel_f_donor_receiver_support_summary.tsv"),
  sep = "\t"
)

# -------------------------------------------------------------------------
# Panel g source: discovery and validation within-cohort receiver ranks.
# Raw weights are retained for audit, but the plotted scale is normalized to
# the maximum immune-receiver value within each cohort and metric.
# -------------------------------------------------------------------------

discovery <- fread(discovery_receiver_file)[target %chin% immune_targets]
discovery[, `:=`(
  cohort = "GSE206785",
  receiver_label = target
)]
discovery <- discovery[, .(
  cohort, receiver_label,
  interaction_count,
  interaction_weight
)]

validation_raw <- fread(validation_cellchat_file)
validation_immune <- c("B_cell", "T_cell", "NK", "Mast", "Myeloid", "Plasma")
validation <- validation_raw[
  source == "Fibroblast" & target %chin% validation_immune & cellchat_pvals < 0.05,
  .(
    interaction_count = .N,
    interaction_weight = sum(lr_probs)
  ),
  by = .(receiver_label = target)
]
validation[, cohort := "GSE183904"]
validation <- validation[, .(
  cohort, receiver_label,
  interaction_count,
  interaction_weight
)]

cross_cohort <- rbindlist(list(discovery, validation), use.names = TRUE)
cross_cohort[, `:=`(
  normalized_count = interaction_count / max(interaction_count),
  normalized_weight = interaction_weight / max(interaction_weight),
  count_rank = frank(-interaction_count, ties.method = "min"),
  weight_rank = frank(-interaction_weight, ties.method = "min"),
  receivers_displayed = .N
), by = cohort]
setorder(cross_cohort, cohort, weight_rank, count_rank)
fwrite(
  cross_cohort,
  file.path(source_dir, "panel_g_cross_cohort_receiver_comparison.tsv"),
  sep = "\t"
)

myeloid_rank <- cross_cohort[receiver_label == "Myeloid", .(
  cohort, receiver_label, interaction_count, interaction_weight,
  normalized_count, normalized_weight,
  count_rank, weight_rank, receivers_displayed
)]
fwrite(
  myeloid_rank,
  file.path(source_dir, "panel_g_cross_cohort_myeloid_rank.tsv"),
  sep = "\t"
)

manifest <- data.table(
  field = c(
    "discovery_dataset", "validation_dataset", "donor_gate",
    "pair_support_gate", "weighted_support_definition",
    "cross_cohort_definition", "claim_boundary"
  ),
  value = c(
    "GSE206785 donor-balanced tumour cells",
    "GSE183904 global broad-lineage CellChat",
    ">=20 Fibroblast and >=20 receiver cells in the same tumour patient",
    "all ligand and receptor complex components detected in >=10% of the corresponding donor populations",
    "fraction of the pooled CellChat probability sum carried by donor-supported frozen interactions",
    "count and weight normalized to the largest immune receiver separately within each cohort",
    "donor expression support and CellChat mass-action probabilities are proxies; neither establishes binding, spatial contact or causality"
  )
)
fwrite(
  manifest,
  file.path(audit_dir, "Figure1_restructured_support_manifest.tsv"),
  sep = "\t"
)

print(roles)
print(receiver_summary)
print(cross_cohort)
cat("Generated restructured Figure 1 source data successfully\n")
