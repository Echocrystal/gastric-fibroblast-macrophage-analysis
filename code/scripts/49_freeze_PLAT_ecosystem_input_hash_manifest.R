options(stringsAsFactors = FALSE, warn = 1)

project <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"
source_project <- "I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
.libPaths(c(file.path(source_project, "runtime", "R_library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
})

manifest_file <- file.path(source_project, "data", "processed", "scrna", "GSE183904",
                           "GSE183904_sample_patient_clinical_manifest.tsv")
manifest <- fread(manifest_file)[include_primary_replication == TRUE &
                                   specimen_class == "Primary_Tumor"]
sce_files <- file.path(source_project, "data", "intermediate", "scrna", "GSE183904",
                       "sample_sce", paste0(manifest$sample_id, "_secondary_qc_sce.rds"))
names(sce_files) <- paste0("GSE183904_", manifest$sample_id)

other_files <- c(
  GSE183904_manifest = manifest_file,
  frozen_signatures = file.path(project, "resources", "gene_sets",
                                "GSE183904_PLAT_ecosystem_bulk_signatures.tsv"),
  TCGA_expression = "I:/ECM-new/ECM-public/data/processed/bulk/tcga_stad_full_symbol_log2cpm_412.rds",
  TCGA_clinical = "I:/ECM-new/ECM-public/results/known_subtypes_v2/TCGA_STAD_current412_with_official_TCGA2014_subtype.tsv.gz",
  GSE62254_expression = file.path(source_project, "data", "processed", "bulk",
                                  "GSE62254_projectable_expression.rds"),
  GSE84437_expression = file.path(source_project, "data", "processed", "bulk",
                                  "GSE84437_projectable_expression.rds"),
  GSE15459_expression = file.path(source_project, "data", "processed", "bulk",
                                  "GSE15459_projectable_expression.rds"),
  Kim2018_TPM = "I:/ECM-new/ECM-public/data/processed/ici/Kim2018/Kim2018_gene_symbol_TPM.tsv.gz",
  Kim2018_clinical = "I:/ECM-new/ECM-public/manifests/Kim2018_Nature_Table2_clinical_45.tsv"
)
all_files <- c(sce_files, other_files)
if (!all(file.exists(all_files))) stop("Input file absent")

hash_table <- data.table(
  input_id = names(all_files),
  path = normalizePath(all_files, winslash = "/", mustWork = TRUE),
  bytes = as.numeric(file.info(all_files)$size),
  sha256 = vapply(all_files, digest, character(1), algo = "sha256", file = TRUE)
)
hash_table[, input_class := fifelse(grepl("^GSE183904_sample", input_id),
                                    "single_cell_SCE", "cohort_or_manifest")]

audit_dir <- file.path(project, "analysis", "audit")
fwrite(hash_table, file.path(audit_dir, "PLAT_ecosystem_input_sha256_manifest.tsv"), sep = "\t")
write_json(
  list(generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
       n_single_cell_SCE = hash_table[input_class == "single_cell_SCE", .N],
       n_total_inputs = nrow(hash_table), inputs = hash_table),
  file.path(audit_dir, "PLAT_ecosystem_input_sha256_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA
)
cat("Frozen input hashes: ", nrow(hash_table), " files; ",
    hash_table[input_class == "single_cell_SCE", .N], " SCE objects.\n", sep = "")

