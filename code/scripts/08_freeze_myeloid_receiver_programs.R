#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else
  "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026"

suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
  library(jsonlite)
  library(digest)
})

resource_dir <- file.path(project_dir, "resources", "gene_sets")
audit_dir <- file.path(project_dir, "analysis", "audit")
dir.create(resource_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

# These pathways were fixed before inspecting their association with any candidate axis.
# They cover the plausible receiver phenotypes without assuming that transcript abundance
# is equivalent to phosphorylation or direct receptor activation.
primary_programs <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_HYPOXIA",
  "GOBP_CELL_MATRIX_ADHESION",
  "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_FOCAL_ADHESION_ASSEMBLY",
  "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
  "GOBP_MACROPHAGE_MIGRATION",
  "GOBP_PHAGOCYTOSIS",
  "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
  "REACTOME_INTEGRIN_SIGNALING",
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
  "REACTOME_FCGAMMA_RECEPTOR_FCGR_DEPENDENT_PHAGOCYTOSIS"
)

collections <- rbindlist(lapply(c("H", "C5", "C2"), function(collection_id) {
  as.data.table(msigdbr(species = "Homo sapiens", collection = collection_id))
}), fill = TRUE)

missing_programs <- setdiff(primary_programs, unique(collections$gs_name))
if (length(missing_programs)) {
  stop("Missing fixed MSigDB programs: ", paste(missing_programs, collapse = ", "))
}

primary <- unique(collections[gs_name %chin% primary_programs, .(
  program = gs_name,
  gene_symbol = toupper(gene_symbol),
  source = "MSigDB",
  evidence_tier = "primary_external",
  gs_collection,
  gs_subcollection,
  gs_description,
  gs_url
)])

# Identity programs are deliberately small and are used only to describe the receiver
# continuum. They are not eligible to become the selected downstream pathway.
identity_markers <- list(
  ID_MONO_FCN1_VCAN = c(
    "CD14", "VCAN", "FCN1", "S100A8", "S100A9", "S100A10",
    "CTSS", "LILRB1", "LILRB2", "CTSD", "TYMP", "MNDA"
  ),
  ID_MACRO_C1QC = c(
    "C1QA", "C1QB", "C1QC", "APOE", "APOC1", "CD163",
    "MERTK", "MARCO", "FOLR2", "VSIG4", "MSR1", "CTSD"
  ),
  ID_MACRO_SPP1_TREM2 = c(
    "SPP1", "TREM2", "GPNMB", "LPL", "LGALS3", "FABP5",
    "CTSB", "CTSD", "APOC1", "MMP9", "CHI3L1", "CD9"
  ),
  ID_DC_CD1C = c(
    "CD1C", "CLEC10A", "FCER1A", "CD1E", "HLA-DQA1",
    "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CST3"
  ),
  ID_INFLAMMATORY_MONO = c(
    "IL1B", "CXCL8", "CCL2", "CCL3", "CCL4", "TNF",
    "NFKBIA", "NFKBIZ", "PTGS2", "IL6", "S100A8", "S100A9"
  ),
  ID_TAM_IMMUNOREGULATORY = c(
    "IL10", "TGFB1", "CD274", "PDCD1LG2", "VSIR", "LGALS9",
    "LILRB2", "LILRB4", "SIGLEC10", "LAIR1", "HAVCR2", "IDO1"
  )
)

identity <- rbindlist(lapply(names(identity_markers), function(program_name) {
  data.table(
    program = program_name,
    gene_symbol = identity_markers[[program_name]],
    source = "a_priori_canonical_markers",
    evidence_tier = "exploratory_identity",
    gs_collection = NA_character_,
    gs_subcollection = NA_character_,
    gs_description = "Canonical marker program used only for receiver-state description",
    gs_url = NA_character_
  )
}))

programs <- unique(rbind(primary, identity, fill = TRUE))
setorder(programs, evidence_tier, program, gene_symbol)

output_path <- file.path(resource_dir, "myeloid_receiver_programs_frozen.tsv")
fwrite(programs, output_path, sep = "\t")

axis_genes <- data.table(
  axis = c(
    rep("ANGPTL2_ITGA5_ITGB1", 3L),
    rep("PLAU_PLAUR", 2L),
    rep("TNC_ITGA5_ITGB1", 3L),
    rep("HLAF_LILRB2", 2L)
  ),
  role = c(
    "ligand", "receptor_subunit", "receptor_subunit",
    "ligand", "receptor",
    "ligand", "receptor_subunit", "receptor_subunit",
    "ligand", "receptor"
  ),
  gene_symbol = c(
    "ANGPTL2", "ITGA5", "ITGB1",
    "PLAU", "PLAUR",
    "TNC", "ITGA5", "ITGB1",
    "HLA-F", "LILRB2"
  )
)
axis_path <- file.path(resource_dir, "candidate_axes_frozen.tsv")
fwrite(axis_genes, axis_path, sep = "\t")

audit <- list(
  generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
  msigdbr_package_version = as.character(packageVersion("msigdbr")),
  primary_programs = primary_programs,
  primary_program_count = length(primary_programs),
  exploratory_identity_program_count = length(identity_markers),
  unique_program_genes = uniqueN(programs$gene_symbol),
  circularity_guard = paste(
    "All candidate ligand and receptor genes are removed from every downstream",
    "program at scoring time; identity programs cannot be promoted as the selected pathway."
  ),
  transcriptomic_scope = paste(
    "Transcript programs are proxies for cell state and do not establish receptor binding,",
    "protein abundance, phosphorylation, or causal pathway activation."
  ),
  program_file_sha256 = digest(file = output_path, algo = "sha256"),
  axis_file_sha256 = digest(file = axis_path, algo = "sha256")
)
write_json(audit, file.path(audit_dir, "myeloid_receiver_program_freeze_manifest.json"),
           pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(audit, pretty = TRUE, auto_unbox = TRUE), "\n")

