# Standalone spatial endpoint comparison for Reviewer 2

This review package recomputes aggregation and inference from frozen existing section-level coefficients. It does not re-fit spot-level models from raw spatial matrices and does not alter the manuscript or original tables.

## Execution

`python run_spatial_endpoint_comparison.py`

Dependencies: Python 3, numpy, pandas, matplotlib and threadpoolctl. Paths resolve relative to this package. Source inputs and original analysis scripts are snapshotted locally on first execution; subsequent runs use the frozen copies.

## Formula and analysis

SPON2 expression at each measured Visium spot was related to the exact ring-1 mean of either an ITGB2-receiver composite or ITGB2 alone under identical rank residualisation. The spot-level composite was [z(ITGB2) + z(receiver score)]/2; both terms were standardised within section before spatial averaging. Covariates were log1p total UMI, same-spot broad fibroblast and epithelial scores, ring-1 broad myeloid score, same-spot frozen ECM-direction score, fixed sender score and PLAT expression. Section correlations were Fisher-z averaged within donor and then equally across donors/patients. Intervals are percentile donor-bootstrap 95% CIs (9,999 draws); nominal P values are exhaustive two-sided donor sign-flip tests. Original BH q values retain their separate families: four candidate receptor-receiver composites for the composite, and three supporting endpoints for ITGB2 alone, separately within each dataset. ITGB2-alone q values do not derive from the four-axis candidate family. The additional six-test BH column is a new exploratory family for this standalone comparison and does not replace original inference. GSE270678 ST_R comprises repeated/extended sections from overlapping donors and is not independent validation. Section eligibility required SPON2 and ITGB2 detection in at least 1% of spots and at least 100 complete endpoint/covariate spots. These mixed-cell associations cannot establish direct binding, receptor activation or causality.

The original composite and candidate-calibration rows represent the same coefficient and P value. Their bootstrap intervals differ slightly because they were generated at different positions in the original sequential RNG stream. The comparison retains the original primary ring-1 distance-profile CI and attaches the identical-coefficient candidate-family q value. All 57 original meta-analysis groups are reconstructed in the original order with seed 2026082654 so that CI values are reproduced.

## Verification

Recomputed 817 donor-endpoint rows and 57 meta-analysis rows from 1,501 section-endpoint rows. All original donor estimates, pooled estimates, CIs, exact P values and original-family q values match within absolute tolerance 1e-12. See verification.json and input_provenance_sha256.tsv.

## Review results

| Dataset | Endpoint | n | Adjusted rho | 95% CI | Exact P | Original-family BH q |
|---|---|---:|---:|---|---:|---:|
| GSE270678_main_ST | ITGB2 + receiver composite | 19 | 0.022975 | 0.014140 to 0.032618 | 8.7738037e-05 | 0.00035095215 |
| GSE270678_main_ST | ITGB2 alone | 19 | 0.005142 | -0.002935 to 0.012882 | 0.23246384 | 0.23246384 |
| GSE270678_ST_R_sensitivity | ITGB2 + receiver composite | 15 | 0.019382 | 0.006484 to 0.033864 | 0.016418457 | 0.047119141 |
| GSE270678_ST_R_sensitivity | ITGB2 alone | 15 | 0.008866 | -0.001165 to 0.022091 | 0.15203857 | 0.15203857 |
| GSE251950_primary | ITGB2 + receiver composite | 9 | 0.045565 | 0.021696 to 0.072707 | 0.00390625 | 0.015625 |
| GSE251950_primary | ITGB2 alone | 9 | 0.024738 | 0.009378 to 0.041767 | 0.015625 | 0.0234375 |

The single-gene endpoint is not significant in the main cohort or ST_R, and is small and positive in GSE251950. The composite therefore has stronger statistical support than ITGB2 alone in the primary cohort, but this comparison does not quantify the percentage contribution of the receiver score or prove that it causes the composite result.

## Output files

- spatial_endpoint_comparison.tsv / .json: six-row comparison with separately labelled original and new exploratory q values.
- spatial_endpoint_donor_results.tsv: donor/patient values for the displayed two endpoints.
- spatial_endpoint_section_results.tsv: section values and original eligibility flags for those endpoints.
- all_recomputed_meta_results.tsv: all original endpoint groups used to reproduce RNG order and original q families.
- spatial_endpoint_comparison.png / .pdf: forest plot.
- verification.json: exact numerical checks; input_provenance_sha256.tsv: source snapshots and hashes.
