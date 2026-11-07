# WJGO analysis code: local release staging

This directory is a **local staging package**, not a public repository or a
verified end-to-end rerun. It accompanies the manuscript on a PLAT-rich
fibroblast and C1QC+APOE+ macrophage context in gastric cancer. The current
submission manuscript and supplemental workbook were not edited to create it.

## Contents

- `code/`: selected copies of the analysis and figure-generation scripts.
  `code/CODE_FILE_MANIFEST.tsv` records each copied file and its SHA-256 hash.
- `gene_sets/`: the frozen 30-gene sender and 30-gene receiver signatures,
  their original source table, and fibroblast/myeloid state marker tables.
- `resources/`: 11 additional frozen gene/program and ligand-receptor TSVs
  referenced by selected scripts. `resources/RESOURCE_FILE_MANIFEST.tsv`
  records source paths, copy hashes, origins, and release recommendations;
  `resources/README.md` gives third-party attribution and reuse terms.
- `docs/`: data provenance, main and supplemental figure script maps, a
  random-seed registry, known documentation gaps, and publication instructions.

`manuscript_revisions_draft/` contains an internal supplemental-table revision
draft. It and the local staging helper `_stage_code.py` are excluded by
`.gitignore` and are not part of the proposed public code archive.

The two primary 30-gene signatures are already present in the submitted
`Supplemental_Tables_1-10.xlsx` sheet `S6_FigureS6a_TCGA_gene_annotati`.
The separately named TSV in `gene_sets/` makes them easier to retrieve. A
`signature_rank` is a selection rank, **not** a scoring weight: measured genes
are equally weighted when sample-level percentile ranks are averaged, then
scores are standardized within each cohort.

## What this package establishes

The files identify code, frozen gene lists, seeds, and public input sources
for the current analysis. Original code copies have deliberately not been
rewritten. Several scripts still contain Windows absolute paths and refer to
processed inputs in another local project. There is no single command that
recreates every submitted figure from a fresh download at this stage. See
`docs/RELEASE_READINESS.md` before describing this as a fully portable pipeline.

Some frozen program tables include MSigDB gene sets, and two ligand-receptor
tables derive from the CC BY 4.0 NicheNet-v2 Zenodo record. Their source and
attribution requirements are documented in `resources/README.md`. A repository
license chosen for the authors' code does not replace those source terms.

The submission's numerical figure data are supplied separately in
`Supplemental_Tables_1-10.xlsx`. Raw GEO/TCGA matrices, local Python/R
runtimes, previous submission packages, reviewer correspondence, and signed
documents are not included here.

## Before publishing the source archive

1. Confirm the selected script-to-output mapping in
   `docs/MAIN_FIGURE_CODE_MAP.tsv`, `docs/SUPPLEMENTAL_FIGURE_CODE_MAP.tsv`,
   and `code/CODE_FILE_MANIFEST.tsv`.
2. Review the copied files for credentials and personal data, choose a code
   license, and publish a tagged version. Put the actual URL and version in
   the manuscript only after publication.

For a portable rerun, replace local absolute paths with a project/data root,
document the public input preparation steps, record installable Python/R
environments, and rerun representative modules in a clean clone. Until then,
describe the repository as source-code availability, not a verified turnkey
pipeline.

Choose the statement in `docs/CODE_AVAILABILITY_TEXT.md` that matches the
actual release state at submission.
