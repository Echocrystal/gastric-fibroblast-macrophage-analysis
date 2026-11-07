# Gastric fibroblast–macrophage analysis code

Source archive assembled on 7 November 2026 for the revised WJGO manuscript.
It contains the original staged analysis scripts, frozen signatures/resources,
and the later spatial, clinical-biomarker, signature-overlap and GSEA revisions.

## Contents

- `code/`: the 103 original staged R/Python scripts plus the preserved S10
  cosmetic assembly source. Original script-to-figure maps are in `docs/`.
- `gene_sets/`: the complete primary 30-gene sender and 30-gene receiver lists
  and state-marker source tables. Selection rank is not a scoring weight.
- `resources/`: frozen gene-program and ligand–receptor tables, with source
  attribution and third-party reuse terms in `resources/README.md`.
- `revision_signature_overlap_20261101/`: published reference lists, harmonised
  gene symbols and a portable recount of the seven displayed comparisons.
- `revision_20261103/`: completed revision analyses and their small frozen
  inputs; the directory name is historical and includes the 4 November updates.
  The current plotting/analysis entry points are listed in
  `docs/CURRENT_REVISION_MAP.md`.
- `UPLOAD_FILE_MANIFEST.tsv`: repository-relative file names, sizes and hashes.

The current S5 uses four GSEA panels. Earlier figure candidates remain as
historical sources; use the current revision map for the intended final plots.
The Overall spatial estimate pools GSE270678 main ST and GSE251950 participants.
ST_R is a sensitivity subset and is not added as a third independent cohort.

## Running and scope

See `docs/CURRENT_REVISION_MAP.md` for module commands, dependencies and output
scope. Selected revision modules run from the supplied small snapshots. The
archive also preserves historical scripts that require prepared matrices,
external resources or path configuration. No complete raw-data-to-all-figures
clean-clone run has been validated. Publishing this source archive establishes
code availability; it does not by itself establish full pipeline portability.

The new GSEA curves reconstruct running enrichment scores from the frozen
original ranking, gene sets and GSEA statistics; drawing curves does not rerun
the original permutation inference. The spatial pooling module starts from
frozen section/donor coefficients rather than raw Visium spots.

Raw GEO/TCGA matrices and local runtimes are obtained separately. Data sources
are documented in `docs/DATA_SOURCES.md`. The finalized Supplementary Table 11
is supplied with the journal supplementary materials; the older internal draft
is excluded here. No public repository URL or DOI is claimed before upload.

## Attribution and version

Retain the original source attributions and reuse terms. No repository-wide
license has been selected by the authors; third-party source terms still apply.
After uploading, record the repository URL and a release tag or commit in the
manuscript, Table 11 and reviewer response. An archival DOI can be obtained by
linking a GitHub release to Zenodo if needed.
