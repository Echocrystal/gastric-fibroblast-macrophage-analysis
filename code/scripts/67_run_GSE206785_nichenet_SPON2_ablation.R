#!/usr/bin/env Rscript

## Outcome-blind SPON2-specific NicheNet ligand-ablation wrapper.
## The original four-candidate workflow is left unchanged; only the candidate,
## output directory, matched-ligand resource and audit filename are redirected.
orig <- "I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026/scripts/24_run_GSE206785_nichenet_ligand_ablation.R"
txt <- readLines(orig, warn = FALSE)
txt <- sub('candidates <- c("TGFBI", "ANGPTL2", "PLAU", "TNC")',
           'candidates <- c("SPON2")', txt, fixed = TRUE)
txt <- sub('result_dir <- file.path(project, "analysis", "results", "nichenet", "GSE206785")',
           'result_dir <- file.path(project, "analysis", "results", "nichenet", "GSE206785_SPON2")', txt, fixed = TRUE)
txt <- sub('GSE206785_nichenet_candidate_matched_ligands.tsv',
           'GSE206785_nichenet_SPON2_matched_ligands.tsv', txt, fixed = TRUE)
txt <- sub('GSE206785_NicheNet_ligand_ablation_manifest.json',
           'GSE206785_NicheNet_SPON2_ligand_ablation_manifest.json', txt, fixed = TRUE)
eval(parse(text = txt), envir = .GlobalEnv)
