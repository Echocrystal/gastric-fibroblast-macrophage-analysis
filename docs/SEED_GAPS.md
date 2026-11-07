# Random seed audit for the WJGO code release

Source of truth: the current upload copies of `Manuscript.docx` and `Supplementary Information.docx` in `WJGOsubmit/00_Ready_to_Upload`, read on 2026-10-08. Script paths in `SEED_REGISTRY.tsv` are relative to the staged `code/` directory. The QC script is a byte-identical copy from a separate upstream project. This is a documentation audit, not a claim that all scripts are ready to run from a public repository.

## Confirmed coverage

- Supplementary Method 2 gives the `scDblFinder` seed rule `183904 + sample number`, matching the upstream QC script.
- Supplementary Methods 3, 4, 7 and 8 give the GSE206785 sampling/embedding, global atlas SVD/UMAP/Leiden and lineage PCA/Leiden/UMAP seeds. Their reported values match the current scripts.
- Supplementary Method 9 gives `20260804` for the Figure 3 LIANA permutations. Supplementary Method 22 gives the three scTenifoldKnk fitting seeds `2026072739`, `2026072740` and `2026072741`.
- Exact donor/patient sign-flip tests in the spatial code enumerate every sign pattern. They require no random seed. EdgeR fitting, fixed signature projection, rank tests, HC3, BH adjustment and Fisher-z pooling do not themselves require a random seed once inputs are fixed.

## Missing or incomplete manuscript documentation

1. Supplementary Method 5 and the corresponding main-text CellChat paragraph report 100 bootstrap resamples but omit seed `20260803` (`scripts/57_run_GSE206785_CellChat_broad_immune.R:35`).
2. Supplementary Method 9 says the 100-cell patient-state cap is deterministic. Its implementation samples with seeded `runif`, using `20260801` (`scripts/46_prepare_GSE183904_CCC_inputs.R:4,79-96`). The 9999 residual permutations for the four-state patient screen use `2026072737 + match(state_id, recurrent_states)` (`scripts/37_characterize_GSE183904_C1QC_APOE_receiver_state.R:223-235`); the frozen label order maps C1, C3, C4 and C5 to `2026072738` through `2026072741`.
3. Supplementary Method 10 reports a restricted-model wild bootstrap but omits its 9999 draws and the primary/sensitivity seeds `2026072742` and `2026072743` (`scripts/42_correct_GSE183904_state_abundance_inference.R:116-122,154-159`).
4. Supplementary Method 12 and the main text call the 10000 GSEA permutations fixed-seed without giving `20260801` (`figures/gen_fig4_patient_association.py:197-202,511-516`).
5. Supplementary Method 14 reports consensus-clustering stability without giving primary PAM seed `2026080151`, sensitivity Ward seed `2026080152`, or the default 1000 repetitions (`scripts/51_consensus_cluster_frozen_PLAT_ecosystem_states.R:45-60`). The script permits `CC_REPS` to override 1000; check the run manifest before claiming that exact repetition count for released outputs.
6. Supplementary Methods 20 and 21 report bootstrap confidence intervals but omit the 9999-draw seeds `2026080135` (GSE270678 main), `2026080637` (matched-neighborhood sensitivity), `2026080638` (GSE251950), and `2026082654` (SPON2-ITGB2 context). Their sign-flip P values are exhaustive, so no seed should be attached to the sign flips.
7. Supplementary Method 22 lists the virtual-knockout algorithm seeds but omits the earlier 30-cell-per-sample donor-cap sampling rule `2026072738 + i`, where `i` is the 1-based frozen manifest row (`scripts/38_prepare_GSE183904_C1QC_APOE_scTenifoldKnk.R:55-60`).
8. Supplementary Method 23 and Supplemental Figure 13 omit the 2000-permutation ITGB2 differential-regulation GSEA rule `2026081101 + algorithmic seed + zero-based programme index` (`scripts/69_run_ITGB2_diffreg_gsea.py:23-25,102-109`).
9. Supplementary Method 16 omits the conditional Monte Carlo chi-square procedure used for sparse configuration tables: 99999 replicates and a sequential seed beginning after `20260805`. The exact increment rule differs between `scripts/62_validate_clinicopathology_external_cohorts.R:257-266` and `scripts/63_analyze_cohort_specific_clinicopathology.R:158,232-233`; identify the output-producing script for each reported table before adding a seed sentence.
10. `scripts/54_prepare_Figure1_GSE206785_unbiased.py:445-453` also runs a global LIANA consensus analysis with 1000 permutations and seed `20260802`. Supplementary Method 3 records `20260802` for preparation, but the communication subsection does not identify this global LIANA run separately from the broad CellChat analysis.

Supplementary Method 24 currently says that random seeds are provided in the respective method subsections. That statement is broader than the documented coverage above. Figure-only random seeds, such as the Figure 4 descriptive bootstrap bands, can be included in a release README and need not crowd the inferential methods.

## Branch and release checks

- Scripts 59 and 60 use separate earlier TCGA-only Monte Carlo chi-square seed rules. The current multi-cohort scripts 62 and 63 must not be assigned the earlier `20260804` seed by convenience. Check source-data paths used by the final figures and tables before deciding whether either older branch belongs in the release.
- The QC source script is copied into staged `code/upstream_qc/`, but its ordered sample manifest is not included. The seed formula alone cannot recreate QC calls.
- Fixed seeds improve repeatability but are insufficient by themselves. Record exact input accessions/files, ordered cell/sample manifests, package versions, parameters and outputs for each released workflow. Some threaded numerical steps may differ slightly across software or hardware versions.
- Do not present the numeric seed for an exploratory or historical branch as if it were used for a different final analysis. The registry distinguishes `current_core`, `current_supporting`, `branch_or_legacy`, `figure_only` and `upstream_dependency` for this reason.
