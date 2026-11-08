# SPON2-ITGB2 spatial candidate-context analysis lock - 2026-08-26

## Status

Locked before inspecting any endpoint from the new SPON2-ITGB2 spatial-neighbourhood analysis. This is a focused post hoc analysis motivated by the cross-cohort single-cell communication screen and the receiver-side ITGB2 perturbation result. No spatial endpoint, ring, covariate set, cohort or calibration family will be replaced after result inspection.

## Scientific question

Do Visium spots carrying a SPON2-rich PLAT-rich fibroblast context show proximity-conditioned co-expression with an ITGB2-positive C1QC+APOE+ macrophage receiver context, and is the receiver-side receptor signal locally associated with ECM-remodelling programmes?

## Cohort hierarchy

- Primary: GSE270678 main ST, 33 sections from 19 donors.
- Same-donor sensitivity: GSE270678 ST_R repeated or extended sections from 16 donors overlapping the primary cohort. ST_R is analysed separately and is not treated as independent replication.
- Independent secondary projection: nine primary gastric-tumour sections from nine patients in GSE251950; the matched peritoneal-metastasis section is excluded.
- Section effects are estimated first. Multiple sections are Fisher-z pooled within donor, and donors or patients are weighted equally.

## Fixed states and focal genes

- Fibroblast sender state: top 30 genes of `sender_CXCL14_POSTN_state_axis_excluded`.
- Macrophage receiver state: top 30 genes of `receiver_C1QC_APOE_specific_axis_excluded`.
- SPON2 is absent from the sender signature; ITGB2 is absent from the receiver signature.
- PLAT remains a single-gene fibroblast-state anchor and is not part of the sender score.
- Focal ligand or receptor genes are excluded from any functional programme used as an outcome.
- All signatures and programmes use the existing fixed gene lists; no gene is added or removed on the basis of spatial results.

## Expression and spatial graph

- Raw H5 counts are aligned to the existing QC-passing in-tissue spot table.
- Expression is transformed as `log1p(CP10K)` separately within each section.
- State and programme scores are means of within-section gene z scores.
- The graph uses exact 10x Visium hex-grid offsets. Ring 0 is the same spot; ring 1 contains the six immediate lattice neighbours; rings 2 and 3 are exact graph-distance shells rather than cumulative disks.
- Missing edge neighbours are not imputed; the mean of available spots within each shell is used.

## Coverage gate

- SPON2 and ITGB2 must be present in the feature matrix and detected in at least 1% of QC-passing spots in a section.
- At least 100 spots with complete endpoint and covariate data are required for a section-level coefficient.
- A donor is analyzable when at least one section passes the section gate.
- Formal primary-cohort inference requires at least 15 analyzable donors. Independent projection is interpreted only when at least seven of nine primary-tumour patients are analyzable.
- Coverage failures are reported without relaxing the gate.

## Primary endpoint

Within each section, calculate the residual Spearman correlation between same-spot SPON2 expression and the exact ring-1 mean of an ITGB2-receiver composite. The composite is the mean of within-section standardized ITGB2 expression and the standardized fixed C1QC+APOE+ receiver-state score.

Rank residualization includes:

1. log total UMI at the sender spot;
2. broad fibroblast score at the sender spot;
3. exact ring-1 broad myeloid score;
4. epithelial score at the sender spot;
5. fixed ECM-direction score at the sender spot;
6. fixed fibroblast sender-state score;
7. PLAT expression.

A positive coefficient indicates that higher SPON2 expression within comparable fibroblast/ECM contexts is associated with a stronger neighbouring ITGB2-positive macrophage receiver context.

## Prespecified supporting endpoints

1. Same-spot SPON2 versus the fixed sender-state score, adjusted for log UMI, broad fibroblast, broad myeloid, epithelial and ECM-direction scores and PLAT.
2. SPON2 versus exact ring-1 ITGB2 alone under the primary covariate set.
3. SPON2 versus exact ring-1 receiver-state score alone under the primary covariate set.
4. The primary ITGB2-receiver composite association at rings 0, 2 and 3 to characterize distance dependence.
5. Same-spot ITGB2 versus ECM degradation, myeloid-cell migration, phagocytosis and cell-matrix adhesion programme scores, adjusted for log UMI, broad myeloid, fibroblast, epithelial and ECM-direction scores and the fixed receiver-state score.
6. Exact ring-1 receiver-side ITGB2 versus the same four neighbouring functional programmes under ring-matched lineage covariates.

Functional programmes are scored only when at least 20 informative genes remain after excluding SPON2, ITGB2, PLAT and all sender/receiver signature genes.

## Candidate-axis calibration

The same ring-1 analysis is applied to four Figure 3 fibroblast-to-myeloid candidate contexts:

- SPON2-ITGB2;
- APP-CD74;
- NRG1-HLA-DPB1;
- CXCL14-CXCR4.

For each pair, the endpoint is same-spot ligand expression versus the ring-1 mean of a receptor-receiver composite constructed identically to the primary composite. The same covariate policy is used, with the focal ligand and receptor kept outside the fixed state scores. Benjamini-Hochberg correction is applied across these four candidate contexts. SPON2-ITGB2 is reported with its rank among the four candidates; a calibration result does not replace the prespecified primary endpoint.

## Inference and interpretation thresholds

- Section correlations are Fisher-z pooled within donor.
- Donor-level pooled correlations are summarized by equal-weight Fisher-z pooling, donor bootstrap 95% confidence intervals and exhaustive two-sided donor sign-flip tests.
- Primary support requires a positive pooled coefficient, a 95% confidence interval excluding zero, and a two-sided donor sign-flip P < 0.05 in GSE270678 main ST.
- Cross-dataset robustness is supported by a concordant positive effect in at least six of nine GSE251950 patients and a positive pooled estimate. Given the nine-patient sample size, the independent cohort is treated as directional replication rather than a second confirmatory test.
- ST_R is used only to assess same-donor section robustness.
- Directionally coherent findings that do not meet all formal thresholds are retained as exploratory calibration and placed in supplementary material.

## Evidence boundary

- Visium spots are mixed-cell capture units. The analysis tests proximity-conditioned transcriptomic co-expression, not cell-resolved contact.
- A positive result can strengthen the spatial plausibility of the SPON2-ITGB2 candidate context but cannot establish protein binding, receptor activation, signal mediation, recruitment or causality.
- The analysis remains explicitly post hoc because the candidate axis was selected from earlier single-cell and perturbation results.
