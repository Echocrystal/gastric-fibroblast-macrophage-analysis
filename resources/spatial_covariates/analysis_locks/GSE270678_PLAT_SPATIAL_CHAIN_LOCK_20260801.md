# GSE270678 PLAT spatial-chain analysis lock — 2026-08-01

## Status

Frozen before inspecting any GSE270678 PLAT–sender–receiver–ECM endpoint result.

## Cohort hierarchy

- Primary external spatial validation: 33 main-ST sections from 19 donors.
- Replicate/sensitivity cohort: 39 ST_R sections; overlapping donors are not treated as independent from main ST.
- Inference unit: donor. Section effects are estimated first and Fisher-z pooled within donor.

## Frozen states

- Sender: top 30 genes of `sender_CXCL14_POSTN_state_axis_excluded`.
- Receiver: top 30 genes of `receiver_C1QC_APOE_specific_axis_excluded`.
- PLAT is absent from the sender score.
- Candidate PLAT-axis genes are absent from the receiver score.
- State signatures were derived and frozen before this spatial validation.

## Frozen functional endpoints

- ECM degradation: `REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX` after excluding PLAT, all sender-signature genes, all receiver-signature genes and candidate receiver-axis genes.
- Cell–matrix adhesion: `GOBP_CELL_MATRIX_ADHESION` under the same exclusions.

## Primary three-link family

1. `L1_PLAT_sender_same`: same-spot PLAT versus PLAT-excluded sender state.
2. `L2_sender_receiver_r1`: sender state versus exact first-ring neighbouring receiver state.
3. `L3_receiver_ECMdeg_r0`: receiver state versus same-spot ECM-degradation programme.

The three endpoints form one family for exhaustive donor-level sign-flip maxT-FWER inference.

## Prespecified secondary endpoints

- Sender–receiver association at rings 0, 2 and 3.
- Receiver–ECM degradation association at rings 1–3.
- Receiver–cell–matrix adhesion association at rings 0–3.
- PLAT incremental association with ring-1 receiver, ECM degradation and cell–matrix adhesion after adjusting for the frozen sender/receiver states.

## Spatial graph

- Exact Visium hex-grid offsets.
- Ring 0: same spot.
- Ring 1: six immediate lattice neighbours.
- Rings 2 and 3: exact graph-distance shells, not cumulative disks.
- The mean of available spots in the specified shell is used; missing edge neighbours are not imputed.

## Within-section scoring

- Raw H5 counts aligned to the existing QC-passing spot table.
- `log1p(CP10K)` per gene.
- Programme scores are means of within-section gene z scores.
- At least 10 informative genes are required for each 30-gene state signature.
- At least 20 informative genes are required for each large functional programme.

## Covariate policy

All models use rank residualization and include log total UMI, frozen broad lineage scores and the frozen ECM direction as appropriate. Spatial receiver/programme outcomes use ring-matched broad-lineage lag covariates. No stage, Lauren or sampling-location covariate is fitted at the spot level.

## Evidence boundary

- Associations are spatial transcriptomic associations between mixed-cell capture spots.
- They do not establish PLAT secretion, receptor engagement, direct cell contact, protein activity, mediation or causality.
- ST_R is a repeated/extended section sensitivity analysis, not an independent external cohort.
