#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else
  "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- if (length(args) >= 2L) normalizePath(args[[2]], winslash = "/", mustWork = TRUE) else
  "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
  library(SingleCellExperiment)
})

resource_dir <- file.path(project_dir, "resources", "ligand_receptor")
audit_dir <- file.path(project_dir, "analysis", "audit")
dir.create(resource_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

input_rds <- file.path(resource_dir, "lr_network_human_21122021.rds")
if (!file.exists(input_rds)) stop("Missing official NicheNet LR network: ", input_rds)

lr_raw <- as.data.table(readRDS(input_rds))
required <- c("from", "to", "database", "source")
if (!all(required %chin% names(lr_raw))) stop("Unexpected NicheNet LR schema")

lr_raw[, `:=`(
  ligand = toupper(trimws(as.character(from))),
  receptor = toupper(trimws(as.character(to))),
  database = as.character(database),
  source = as.character(source)
)]
lr_raw <- lr_raw[nzchar(ligand) & nzchar(receptor)]
lr <- lr_raw[, .(
  databases = paste(sort(unique(database)), collapse = ";"),
  sources = paste(sort(unique(source)), collapse = ";"),
  evidence_rows = .N
), by = .(ligand, receptor)]

# Presence in both single-cell datasets is recorded before analysis. It is not
# used to choose a clinical outcome or a winning axis.
g206_expr <- file.path(source_project, "data", "raw", "scrna", "GSE206785", "GSE206785_scgex.txt.gz")
con <- gzfile(g206_expr, open = "rt")
g206_header <- strsplit(readLines(con, n = 1L, warn = FALSE), ",", fixed = TRUE)[[1L]]
close(con)
if (!length(g206_header) || g206_header[[1L]] != "Cell") stop("Unexpected GSE206785 header")
g206_genes <- unique(toupper(g206_header[-1L]))

g183_sce <- file.path(source_project, "data", "intermediate", "scrna", "GSE183904",
                      "sample_sce", "sample1_secondary_qc_sce.rds")
if (!file.exists(g183_sce)) stop("Missing audited GSE183904 SCE checkpoint")
g183_genes <- unique(toupper(rownames(readRDS(g183_sce))))

genes <- sort(unique(c(lr$ligand, lr$receptor)))
gene_manifest <- data.table(
  gene_symbol = genes,
  role_ligand = genes %chin% lr$ligand,
  role_receptor = genes %chin% lr$receptor,
  available_GSE206785 = genes %chin% g206_genes,
  available_GSE183904 = genes %chin% g183_genes
)

benchmark_pairs <- data.table(
  ligand = c("CXCL12", "CXCL13", "CCL19", "CCL21", "TGFB1", "IL6", "CCL2", "CSF1", "SPP1"),
  receptor = c("CXCR4", "CXCR5", "CCR7", "CCR7", "TGFBR2", "IL6R", "CCR2", "CSF1R", "CD44"),
  benchmark_reason = c(
    "published gastric CAF/TLS benchmark",
    "canonical TLS/B-cell benchmark",
    "canonical TLS trafficking benchmark",
    "canonical TLS trafficking benchmark",
    "crowded CAF immunosuppression benchmark",
    "crowded CAF inflammatory benchmark",
    "crowded myeloid recruitment benchmark",
    "crowded macrophage-support benchmark",
    "crowded stromal-myeloid benchmark"
  )
)
lr <- merge(lr, benchmark_pairs, by = c("ligand", "receptor"), all.x = TRUE)
lr[, published_benchmark := !is.na(benchmark_reason)]
lr[, `:=`(
  ligand_available_GSE206785 = ligand %chin% g206_genes,
  receptor_available_GSE206785 = receptor %chin% g206_genes,
  ligand_available_GSE183904 = ligand %chin% g183_genes,
  receptor_available_GSE183904 = receptor %chin% g183_genes
)]
setorder(lr, ligand, receptor)

lr_tsv <- file.path(resource_dir, "lr_network_human_21122021_frozen.tsv")
gene_tsv <- file.path(resource_dir, "lr_gene_universe_frozen.tsv")
fwrite(lr, lr_tsv, sep = "\t", na = "")
fwrite(gene_manifest, gene_tsv, sep = "\t")

manifest <- list(
  frozen_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  source = "NicheNet human ligand-receptor network, Zenodo record 7074291",
  source_url = "https://zenodo.org/records/7074291/files/lr_network_human_21122021.rds",
  source_sha256 = digest(file = input_rds, algo = "sha256"),
  unique_pairs = nrow(lr),
  unique_ligands = uniqueN(lr$ligand),
  unique_receptors = uniqueN(lr$receptor),
  unique_genes = nrow(gene_manifest),
  genes_available_GSE206785 = sum(gene_manifest$available_GSE206785),
  genes_available_GSE183904 = sum(gene_manifest$available_GSE183904),
  pairs_complete_GSE206785 = sum(lr$ligand_available_GSE206785 & lr$receptor_available_GSE206785),
  pairs_complete_GSE183904 = sum(lr$ligand_available_GSE183904 & lr$receptor_available_GSE183904),
  output_lr_sha256 = digest(file = lr_tsv, algo = "sha256"),
  output_gene_sha256 = digest(file = gene_tsv, algo = "sha256"),
  selection_policy = paste(
    "All unique NicheNet human LR pairs were frozen before outcome analysis;",
    "benchmark flags do not remove pairs and are used only to prevent false novelty claims."
  )
)
write_json(manifest, file.path(audit_dir, "lr_universe_freeze_manifest.json"),
           pretty = TRUE, auto_unbox = TRUE)

cat(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), "\n")
