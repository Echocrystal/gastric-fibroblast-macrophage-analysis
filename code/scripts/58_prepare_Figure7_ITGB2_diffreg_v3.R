# Prepare reference-style differential-regulation outputs for Figure 7 v3.
# The source RDS objects are never modified. ITGB2 itself is retained in the
# volcano table for provenance but excluded from downstream-gene hit counts.

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
base <- file.path(project, "analysis/results/virtual_knockout/GSE183904/C1QC_APOE_macrophage_receptor_screen_v2")
out <- file.path(project, "figures/source_data/Figure7_receiver_screen_v3")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

seeds <- c(2026072739, 2026072740, 2026072741)
target <- "ITGB2"

read_diffreg <- function(seed) {
  rds <- file.path(base, paste0("scTenifoldKnk_seed", seed, "__", target, ".rds"))
  x <- readRDS(rds)
  d <- x$diffRegulation
  d$seed <- seed
  d
}

diff_all <- do.call(rbind, lapply(seeds, read_diffreg))
write.table(diff_all, file.path(out, "Figure7B_ITGB2_diffreg_all_seeds.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

rep_seed <- seeds[1]
rep <- diff_all[diff_all$seed == rep_seed, ]
rep$target_self_hit <- rep$gene == target
rep$reference_hit <- abs(rep$Z) >= 2 & rep$p.adj < 0.01 & !rep$target_self_hit
write.table(rep, file.path(out, "Figure7B_ITGB2_diffreg_seed2026072739.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Gene-level cross-seed summary. Percentile is high for larger network distance.
diff_all$distance_percentile <- ave(diff_all$distance, diff_all$seed,
                                    FUN = function(z) rank(z, ties.method = "average") / length(z))
gene_ids <- setdiff(unique(diff_all$gene), target)
summary <- do.call(rbind, lapply(gene_ids, function(g) {
  q <- diff_all[diff_all$gene == g, ]
  data.frame(
    gene = g,
    median_Z = median(q$Z),
    min_Z = min(q$Z),
    max_Z = max(q$Z),
    median_FC = median(q$FC),
    median_log10_FC = median(log10(pmax(q$FC, 1e-12))),
    min_log10_FC = min(log10(pmax(q$FC, 1e-12))),
    max_log10_FC = max(log10(pmax(q$FC, 1e-12))),
    median_distance_percentile = median(q$distance_percentile),
    min_distance_percentile = min(q$distance_percentile),
    max_distance_percentile = max(q$distance_percentile),
    q05_seeds = sum(q$p.adj < 0.05),
    q01_seeds = sum(q$p.adj < 0.01),
    z2_seeds = sum(abs(q$Z) >= 2),
    stringsAsFactors = FALSE
  )
}))

# Stable ranking: recurrence first, then standardized effect and network distance.
summary <- summary[order(-summary$q05_seeds, -summary$median_Z,
                         -summary$median_distance_percentile, summary$gene), ]
write.table(summary, file.path(out, "Figure7C_ITGB2_recurrent_diffreg_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ITGB2-centred WT network for the representative seed. Include all strict
# differential-regulation hits plus the 12 strongest direct ITGB2 neighbours.
x <- readRDS(file.path(base, paste0("scTenifoldKnk_seed", rep_seed, "__", target, ".rds")))
M <- as.matrix(x$tensorNetworks$WT)
d <- x$diffRegulation
sig_genes <- d$gene[d$gene != target & abs(d$Z) >= 2 & d$p.adj < 0.01]
neighbours <- colnames(M)[order(abs(M[target, ]), decreasing = TRUE)[1:12]]
nodes <- unique(c(target, sig_genes, neighbours))

node_meta <- data.frame(
  gene = nodes,
  node_class = ifelse(nodes == target, "perturbed_target",
                      ifelse(nodes %in% sig_genes, "differential_regulation_hit", "direct_network_neighbour")),
  stringsAsFactors = FALSE
)
node_meta$Z <- d$Z[match(node_meta$gene, d$gene)]
node_meta$p_adj <- d$p.adj[match(node_meta$gene, d$gene)]
node_meta$FC <- d$FC[match(node_meta$gene, d$gene)]
node_meta$median_q05_seeds <- summary$q05_seeds[match(node_meta$gene, summary$gene)]
write.table(node_meta, file.path(out, "Figure7D_ITGB2_network_nodes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

edge_idx <- which(upper.tri(M[nodes, nodes, drop = FALSE]), arr.ind = TRUE)
edge_values <- M[nodes, nodes, drop = FALSE][edge_idx]
edges <- data.frame(
  from = nodes[edge_idx[, 1]],
  to = nodes[edge_idx[, 2]],
  weight = as.numeric(edge_values),
  abs_weight = abs(as.numeric(edge_values)),
  stringsAsFactors = FALSE
)
# Keep all direct target edges and strong links among differential-regulation hits.
edges$keep <- (edges$from == target | edges$to == target) |
  (edges$from %in% sig_genes & edges$to %in% sig_genes & edges$abs_weight >= 0.35)
edges <- edges[edges$keep, c("from", "to", "weight", "abs_weight")]
edges <- edges[order(-edges$abs_weight, edges$from, edges$to), ]
write.table(edges, file.path(out, "Figure7D_ITGB2_network_edges.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  paste0("representative_seed\t", rep_seed),
  paste0("target\t", target),
  paste0("n_genes\t", nrow(rep)),
  paste0("strict_hits_excluding_target\t", sum(rep$reference_hit)),
  paste0("strict_threshold\t|Z|>=2 and q<0.01"),
  paste0("network_nodes\t", nrow(node_meta)),
  paste0("network_edges\t", nrow(edges))
), file.path(out, "Figure7_v3_prepare_summary.tsv"))

message("Prepared Figure 7 v3 source tables in: ", out)
