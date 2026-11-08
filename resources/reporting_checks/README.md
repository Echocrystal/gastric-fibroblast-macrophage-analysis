# Existing-model denominators and missing-data reporting

These small snapshots report the inputs and inclusion rules of analyses that
were already completed. No CAMERA, limma, Cox or response model was fitted for
this reporting update. The local extraction script reads the original frozen
inputs and stored fitted objects; the package verifier uses Python's standard
library and does not fit models.

## TCGA CAMERA

- `TCGA_CAMERA_frozen_covariates_412.tsv` preserves the 412 original expression
  sample IDs and their order, ecological score, raw and standardised covariates,
  and official four-class TCGA 2014 subtype coding used by original script 52.
- `TCGA_CAMERA_patient_inclusion.tsv` records every candidate sample in each
  model, its inclusion status and the exact missing model variables.
- `TCGA_CAMERA_model_denominators.tsv` records the existing models' n, design
  ranks, residual degrees of freedom and original gene-statistics n.
- The two `*_design_matrix.tsv` files contain the corresponding model matrices.
- `TCGA_CAMERA_model_variable_missingness.tsv`, `*_missingness_patterns.tsv`
  and `*_subtype_purity_counts.tsv` show the complete-case denominators and
  overlapping missingness. The 400-case purity model excludes 12 tumours.
  The 269-case full model excludes 143: 131 lack subtype only, five lack purity
  only, and seven lack both. Both lineage covariates and the exposure are present
  for all 412 tumours.
- `TCGA_CAMERA_original_input_hash_checks.tsv` records that the original
  expression, score and lineage input SHA-256 values exactly matched their
  archived analysis-audit values when these snapshots were exported. Raw
  expression is not duplicated here.

The primary and full designs have ranks 3 and 8, leaving 397 and 261 residual
degrees of freedom. The sample counts match both the original per-gene limma
statistics and the original analysis stdout. These are reporting checks, not a
new pathway analysis. Model covariates were not statistically imputed.

## Other TCGA models

`TCGA_survival_patient_inclusion.tsv` and `TCGA_survival_denominator_checks.tsv`
record the 21 existing original-projection and lineage-sensitivity survival
models. Their n and event counts were checked against the existing result
tables without fitting. Original script 45's saved native linked-data object
preserves character NA values that TSV export can represent as blanks; script
48's documented blank-to-NA cleaning was retained for its sensitivity models.

`TCGA_model_covariate_availability.tsv` reports individual field availability.
The `tcga_subtype` field used by the existing survival sensitivity (380 available;
STAD_CIN/STAD_EBV/STAD_GS/STAD_MSI/STAD_POLE) differs from the explicit official
`tcga2014_subtype` field used in CAMERA (274 available; CIN/GS/MSI/EBV).
Their denominators must not be substituted for one another. Field availability
also does not determine a multivariable model's joint complete-case n.

## Kim2018

`Kim2018_analysis_input_field_availability.tsv` describes the actual analysed
clinical input: MSI and EBV ISH are available for 45 patients and PD-L1 CPS for
41, with four CPS values missing. Treatment line, tumour anatomic site and
histology are absent fields in that input and were neither adjusted for nor
imputed. This does not claim that the original investigators never collected
these variables. An absent field is different from 45 documented item-level
missing entries.

`Kim2018_stored_fit_denominator_checks.tsv` was extracted from the 11 stored
fitted objects and checked against their existing summaries. All retain 45
patients, including 12 responders, with zero exclusions for their model
variables. CPS is not a covariate in those models, so missing CPS does not
exclude patients from them. `Kim2018_existing_model_patient_inclusion.tsv`
copies the existing inclusion records using the previously published
deidentified KIM001-KIM045 identifiers.

## Verification

From any working directory, run:

```text
python scripts/verify_existing_reporting_checks.py
```

The verifier checks unique IDs, complete-case inclusion and overlapping
missingness, independently calculates matrix ranks from the exported designs,
and checks consistency with existing denominator records. It prints JSON and
does not change any resource or statistical result.
