# Gastric fibroblast–macrophage analysis source archive

This is a **source-code availability archive**, prepared on 3 November 2026.
It contains selected original analysis scripts, frozen gene lists, source maps,
and the completed spatial-endpoint and clinical-biomarker revision modules.
The archive has not been uploaded by the assistant. No public URL is claimed.

## Contents

- `code/`: 103 original staged scripts (47 Python, 56 R), with the original
  `CODE_FILE_MANIFEST.tsv` and script-to-figure maps in `docs/`.
- `gene_sets/`: the frozen primary 30-gene sender and 30-gene receiver lists,
  marker tables and signature-build helper.
- `resources/`: 11 frozen resource TSVs and source/reuse attribution.
- `revision_20261103/`: completed spatial composite versus ITGB2 comparison,
  clinical stratification, score–biomarker relationships and response-annotated
  candidate figures. Each module includes small frozen public-data-derived
  inputs and numeric outputs. See `docs/REVISION_MODULES.md`.
- `revision_signature_overlap_20261101/`: original extraction/analysis sources,
  60 primary genes, published reference lists, symbol-map snapshots and a
  portable recount of the seven final displayed comparisons.
- `docs/table11_draft/`: candidate Table 11 workbook and two TSV inventories,
  explicitly marked as a review draft with unresolved fields.
- `UPLOAD_FILE_MANIFEST.tsv`: byte sizes and SHA-256 for the upload files.

The original manuscript-source code includes historical Figure 4 / macrophage
analyses. **New C3/phagocytosis revision analyses are not included in this
archive version**; their ongoing diagnostic work must be reconciled before
claiming that this archive contains the final revision of every panel.

## Reproduction scope

The original 103 scripts are preserved without rewriting. Several use local
Windows absolute paths or processed matrices from a historical project.
There is no validated single-command clean-clone reproduction of all figures.
`docs/RELEASE_READINESS.md` records these limits. Syntax parsing and hashes do
not establish scientific rerun success.

The included revision modules use their colocated frozen snapshots by default.
They can be run module by module; spatial analysis reaggregates saved
section-level estimates and does **not** refit raw spot-level models. These
inputs support the reported revision calculations, not the entire upstream
data-processing pipeline. No local runtime, raw GEO/TCGA matrix, reviewer
letter, manuscript or signed form is included. The candidate Table 11 inventory
is the only newly added review workbook; it is explicitly nonfinal.

## Source attribution and licensing

Retain `resources/README.md` and its MSigDB and NicheNet-v2 attribution and
source terms. The Kim2018 clinical annotation inputs originate from the
published public Table 2 / frozen manuscript figure source data, as recorded
in their source manifests; this package adds fixed scores and configuration
labels. Primary-study source citations are in `docs/DATA_SOURCES.md`.

No repository-wide code license has been selected by the authors. This archive
does not add a license automatically. Source attribution does not replace
source-specific reuse terms.

## Upload

Upload the **contents of this directory**, so README, code, resources and docs
appear at the repository root. Do not upload the parent research directory.
After upload, verify the file tree and record the real repository URL and
commit hash. A version tag such as `v0.1-source-archive` may identify this source
availability snapshot; reserve a complete reproducibility claim for an
environment-locked clean-clone rerun of the finalized analyses.
