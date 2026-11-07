# Complete reporting of the existing Kim2018 response models

This package reproduces the **11 existing mean-bias-reduced logistic models** and exports their **31 coefficients**, including intercepts and adjustment covariates. It uses the frozen numerical model inputs and does not add model specifications, covariates, thresholds or new biological analyses.

## Run

Requires R with `data.table`, `brglm2` and `jsonlite` installed. The verified environment used R 4.4.0 and brglm2 1.0.1; see `outputs/sessionInfo.txt`.

From this directory:

```sh
Rscript export_existing_Kim2018_models.R
```

From another directory, pass the path to this package as the first argument:

```sh
Rscript path/to/clinical_models/export_existing_Kim2018_models.R path/to/clinical_models
```

All input/output paths are relative to the package directory. If packages are installed in a custom library, set `KIM2018_R_LIBRARIES` to its path (multiple paths may be separated by semicolons). No particular drive, user directory or external project checkout is required.

## Inputs

`inputs/Kim2018_model_inputs_deidentified.tsv` contains the same 45 frozen numerical model-input rows: response, MSI status and fixed scores. Original cohort identifiers have been replaced by arbitrary sequential identifiers `KIM001`–`KIM045`. The package contains no original patient identifiers, names, contact information or raw expression matrix. Numerical scores, response labels, MSI labels and row order are unchanged. `INPUT_FILE_MANIFEST.json` records input SHA256 checksums.

The three other small input tables are the original coefficient and Figure 6H numerical results, retained for exact reproduction checks. The script stops if the original 16 reported terms or the four Figure 6H rows do not reproduce within numerical tolerance.

## Outputs

- `outputs/Kim2018_existing_models_ALL_TERMS.tsv`: all 31 coefficients, SE, Wald z, exponentiated coefficient, 95% CI, P, original applicable BH q, sample/event counts and fit status.
- `outputs/Kim2018_existing_models_FIT_SUMMARY.tsv`: all 11 model formulas, complete-case sample counts, response events, missing-data exclusions and convergence diagnostics.
- `outputs/Kim2018_model_patient_inclusion.tsv`: model inclusion by the arbitrary row identifiers.
- `outputs/Figure6H_four_terms_complete.tsv`: the four figure estimates linked to their fitted models.
- `outputs/*REPRODUCTION_QA.tsv`: comparison with the original coefficient and figure tables.
- `outputs/Kim2018_existing_models_complete.json`: complete numerical results and software versions.
- `outputs/PUBLIC_PACKAGE_QA.json`: verification of the deidentified inputs and regenerated outputs.

The models use 45 patients, including 12 responders and 33 non-responders, with no missing-data exclusions for the variables in these models. All fits converged with full-rank designs and no fitting warnings. This package reproduces model reporting from frozen scores; it does not reproduce the earlier expression-to-score pipeline.

## Effect scales and multiplicity

For non-intercept terms, exponentiated coefficients are odds ratios; for the intercept, the value represents baseline odds. Continuous scores are standardized within the cohort. In the interaction model, each main effect is conditional on the other component having Z=0.

The four Figure 6H rows are selected terms from **three fitted models**. Original BH q values are preserved from the implemented families: four single-score exposure tests separately within univariable and MSI-adjusted analyses, and eight non-intercept terms across the three joint/interaction models. Blank q values indicate terms outside these original correction families; no new q value is imputed. The four displayed rows were not jointly corrected as a four-test family.

The average ecosystem score supplied to the models is `Z_cohort[(sender_Z + receiver_Z) / 2]`, including the final cohort standardisation. The script keeps the supplied frozen scores unchanged.

These single-arm response associations remain exploratory and do not estimate treatment benefit or establish incremental predictive value beyond established biomarkers.
