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
- `clinical_models/`: the complete frozen-input clinical reporting module,
  with existing-model terms, fit summaries and numerical QA; see
  `clinical_models/README.md`.
- `resources/ecm_direction/`, `resources/spatial_covariates/` and
  `resources/reporting_checks/`: additional fixed ECM definitions, frozen spot
  covariates/mappings and model-specific denominator/inclusion records prepared
  on 8 November 2026. Their status is distinguished from the public versions below.
- `code/reproduction/`: bounded reconstruction and verification adapters; see
  `docs/ECM_SPATIAL_RESOURCES.md` for commands, separately obtained raw inputs,
  completed checks and scope. The preserved original sources are in
  `code/frozen_sources/` and `code/provenance_ecm/`.
- `scripts/verify_existing_reporting_checks.py`: checks of frozen reporting
  resources without model refitting; see `resources/reporting_checks/README.md`.
- `docs/PUBLIC_BASE_AND_PREPARED_UPDATE.md`: immutable public base, prepared
  update status, merge instructions and byte-manifest corrections.
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
is excluded here. Public version identifiers and the status of the prepared
8 November additions are recorded below.

## Attribution and version

Retain the original source attributions and reuse terms. No repository-wide
license has been selected by the authors; third-party source terms still apply.

Public repository: https://github.com/Echocrystal/gastric-fibroblast-macrophage-analysis

- Versioned release `v1.0.0`:
  https://github.com/Echocrystal/gastric-fibroblast-macrophage-analysis/releases/tag/v1.0.0
  (commit `d3be5c7a2253b4802548613626a492b83a3f8e9e`).
- Table 11 resource additions:
  https://github.com/Echocrystal/gastric-fibroblast-macrophage-analysis/tree/d19a68ba00b1fcaeb6822af1268cb41ff8c2be42
- Complete clinical model reporting:
  https://github.com/Echocrystal/gastric-fibroblast-macrophage-analysis/tree/c0936462c87442279844f75a0a47f3e84baca8f8
  (public main HEAD verified on 8 November 2026).

This README correction and the additional ECM, spatial-covariate and
model-inclusion resources are part of the locally prepared 8 November update,
supplied with the revised submission as
`Supplementary_Reproducibility_Resources_20261108.zip`. They were not present in
the public versions listed above and have not yet been published as a new
repository commit or release. See `docs/PUBLIC_BASE_AND_PREPARED_UPDATE.md`.
No remote repository edits were performed during preparation of this update.
