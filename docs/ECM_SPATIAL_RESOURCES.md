# Frozen ECM direction and spatial covariates

Prepared locally on 8 November 2026 as additions to the study source archive.
This directory has not itself been published or assigned a new release tag.
Merge its relative paths into a working copy of
`Echocrystal/gastric-fibroblast-macrophage-analysis`; it adds resources and
adapters without replacing the historical statistical source files.

## Definitions and provenance

`resources/ecm_direction/tcga_ecm_centroid_gene_contrasts.tsv` preserves the
original 420-gene TCGA ECM centroid contrast from the preceding local ECM
project. No publication DOI is inferred for that model. Sorting the fixed
C2-minus-C1 contrast in descending order yields the first 50 C2-high genes H
and last 50 C1-high genes L. `frozen_ECM_highlow50.tsv` records their membership,
order and exact centroid values. The adjustment score is

`D = mean(gene Z for H) - mean(gene Z for L)`.

This is distinct from the current study's sender/receiver ecological score.
Patient-level and spatial representations use the same gene identities but
different data layers and standardisation domains:

- GSE183904: aggregate all secondary-QC cells per primary tumour. TMM-normalise
  the original complete 1,717-gene target-panel counts, calculate logCPM with
  prior count 1, then standardise each gene across the 26 primary tumours using
  sample SD. The provided primary all-cell score is the actual ecology
  covariate. It is not the fibroblast-intrinsic sensitivity score.
- Spatial: log1p(CP10K) using full-library totals; standardise each gene across
  retained measured spots within each section using sample SD. GSE270678
  retains requested genes (absent features enter as zero expression);
  GSE251950 intersects sets with measured features and requires at least two
  measured genes. Non-finite gene Z values, including zero variance, contribute
  zero to these precomputed ECM/broad-lineage means.
- The dynamic broad-myeloid mean in the current spatial model uses the eight
  genes in `spatial_model_broad_lineage_markers.tsv`, with at least six measured
  and informative genes. Its mean omits zero-variance genes, as in preserved
  script 35. This differs from zero filling in the precomputed upstream panels.

The coarse annotation source, actual model panels and frozen upstream score
gene map are supplied separately from fine-state annotation modules. Original
R/Python scoring sources are preserved byte-for-byte in `code/provenance_ecm/`.
They retain historical paths for provenance; use the new adapters for the
specified reproducible tasks.

## Small resources supplied

- The full 420-gene centroid contrast and two fixed 50-gene panels.
- The 26-patient sample map, original all-cell frozen scores and compact raw
  target-panel pseudobulk counts (44,642 rows; no individual-cell matrix).
- Actual spatial fibroblast, epithelial and myeloid marker definitions.
- Frozen measured/QC spot covariates, coordinates and identifiers: 63,060 main
  ST spots/33 sections; 86,606 ST_R spots/39 source sections; 30,856 GSE251950
  spots/10 source sections. Native source sections are preserved, including
  excluded samples. Original script 54 selects 33 main ST sections/19 donors,
  37 ST_R sections/16 overlapping donors and nine external primary tumours.
  The ST_R endpoint gate subsequently retains 15 donors; a manifest count of
  16 is therefore not a discrepancy in the 15-donor endpoint result.
- Native covariate numeric strings and original row order are retained.
  Library totals absent from the original ST_R table are added from the
  original H5 matrices. The map links every section to its original raw file.
- `GEO_raw_H5_file_index.tsv` supplies the SHA-256, byte size, GEO accession and
  expected relative location for all 82 source H5s. Original direct download
  URLs are retained where present in source records. GSE251950's recorded
  series archive is `GSE251950_RAW.tar`, containing the per-sample archives
  listed by GEO. No raw H5 matrix is included in this small resource archive.

`docs/ECM_SPATIAL_RESOURCE_MANIFEST.tsv` indexes these additions. Historical
analysis locks are supplied as provenance; their inclusion does not create a
new prespecification or change the post hoc status of the SPON2 analysis.

## Raw H5 layout

Obtain the original filtered-feature H5 files from GEO using the accessions,
names and hashes in `GEO_raw_H5_file_index.tsv`. Put them below a chosen raw
root, preserving this structure:

```text
RAW_ROOT/
  GSE270678/main_public_files/GSM..._filtered_feature_bc_matrix.h5
  GSE270678/ST_R_public_files/GSM..._filtered_feature_bc_matrix.h5
  GSE251950/samples/FOLDER/FOLDER_filtered_feature_bc_matrix.h5
```

The frozen spot selections and coordinates are already supplied, so these
reproduction commands do not require histology image files or a new tissue
selection. They do require the original H5 expression data to recalculate
sender/receiver/programme and ligand/receptor measurements.

## Commands

Run at the merged repository root. Replace `/path/to/RAW_ROOT` with the chosen
data directory. Python requires numpy, pandas, scipy and h5py. The R command
requires data.table and edgeR. A tested version record is supplied; it is not
an environment lock for all historical scripts in the repository.

```sh
Rscript code/reproduction/recompute_patient_ecm_direction.R . validation/patient_ecm

python code/reproduction/verify_frozen_spatial_covariates.py \
  --raw-root /path/to/RAW_ROOT --output-dir validation/spatial_covariates

python code/reproduction/run_frozen_spatial_models.py \
  --raw-root /path/to/RAW_ROOT --output-dir validation/spatial_inputs \
  --check-inputs

python code/reproduction/run_frozen_spatial_models.py \
  --raw-root /path/to/RAW_ROOT --output-dir validation/spatial_sections \
  --verify-first-sections \
  --expected-section-effects revision_20261103/01_spatial_endpoint/inputs/SPON2_ITGB2_section_effects.tsv
```

The adapter uses the preserved byte-identical scripts in `code/frozen_sources/`.
It changes module-level path bindings and injects the already loaded base
module, preserving statistical function bodies, fixed definitions and seeds.
`--repo-root` can override automatic discovery. `--definition-root` is optional
and only needed for an explicitly separate copy of the existing signature and
programme resources; a normal merged repository needs neither option.
`--check-inputs` verifies all required H5 sizes and SHA-256 hashes.

To run all original script-54 sections and its original aggregation/inference:

```sh
python code/reproduction/run_frozen_spatial_models.py \
  --raw-root /path/to/RAW_ROOT --output-dir reproduction/SPON2_context \
  --run-all
```

This last full run is exposed as an executable command, but was not executed
as part of this resource completion. Bounded validation details are below.
The new resource archive does not claim portability of every historical
single-cell, bulk and virtual-perturbation pipeline.

## Verification actually completed

- Recomputed all 26 patient-level ECM scores from the complete original target
  panel using TMM. Maximum absolute score difference: 3.11e-15.
- Checked identity `ECM = high50 - low50` in every supplied spatial row, unique
  section/barcode keys and exact H5-derived library totals for all 82 sections.
- Reconstructed ECM high/low and broad fibroblast/epithelial scores from H5s
  for D01T2, D01T2_2 and GSM7990473. Maximum absolute discrepancy: 5.56e-15.
- Recomputed all 19 existing section endpoints in each of these three cases
  against the previously frozen section-effects table: all denominators
  matched, maximum correlation discrepancy below 9.89e-17.
- Repeated those bounded model calculations in a relocated repository fixture
  using only signature/programme resources from the original public upload
  archive, these new resources and explicitly located raw H5s. No original
  local analysis directory or hidden I:-drive project root was used by the
  statistical functions.

Machine-readable bounded verification files are supplied in
`docs/validation_ecm_spatial/`. These checks validate the supplied resources and
selected existing section calculations. They do not turn same-cohort,
post-selection or post hoc analyses into independent validation.
