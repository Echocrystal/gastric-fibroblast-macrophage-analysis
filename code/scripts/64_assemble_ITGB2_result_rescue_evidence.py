#!/usr/bin/env python
"""Assemble a traceable raw evidence table for the PLAT-rich/SPON2-ITGB2 story."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
AUDIT = PROJECT / "analysis" / "audit"
PLAT = PROJECT / "analysis" / "results" / "PLAT_axis" / "GSE183904"
VK = (
    PROJECT
    / "analysis"
    / "results"
    / "virtual_knockout"
    / "GSE183904"
    / "C1QC_APOE_macrophage_receptor_screen_v2"
)

OUTFILE = AUDIT / "ITGB2_axis_result_rescue_evidence_table_20260814.tsv"
REPORT = AUDIT / "ITGB2_axis_result_rescue_findings_20260814.md"


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


rows: list[dict[str, object]] = []


def add(
    layer: str,
    candidate: str,
    endpoint: str,
    cohort: str,
    unit: str,
    n: object,
    statistic: str,
    estimate: object,
    p_value: object = np.nan,
    fdr: object = np.nan,
    robustness: str = "",
    status: str = "",
    boundary: str = "",
    source: Path | str = "",
) -> None:
    rows.append(
        {
            "evidence_layer": layer,
            "candidate": candidate,
            "endpoint": endpoint,
            "cohort": cohort,
            "independent_unit": unit,
            "n": n,
            "statistic": statistic,
            "estimate": estimate,
            "p_value": p_value,
            "fdr": fdr,
            "robustness": robustness,
            "evidence_status": status,
            "interpretation_boundary": boundary,
            "source_file": str(source),
        }
    )


# Patient-level state association.
state_file = PLAT / "fine_state_axis" / "GSE183904_fine_state_PLAT_receiver_state_associations.tsv"
state = read(state_file)
r = state[
    (state["exposure"] == "sender_state_fraction_asin")
    & (state["endpoint"] == "receiver_state_fraction_asin")
].iloc[0]
add(
    "patient_level_state_association",
    "PLAT-rich fibroblast state -> C1QC+APOE+ macrophage state",
    "arcsine-sqrt transformed state fractions",
    "GSE183904",
    "patient",
    int(r["n"]),
    "covariate-adjusted partial Spearman rho",
    r["adjusted_rho"],
    r["adjusted_permutation_p"],
    r["adjusted_permutation_FDR_BH"],
    f"LOO rho {r['loo_adjusted_rho_min']:.3f}-{r['loo_adjusted_rho_max']:.3f}; positive {r['loo_positive_fraction']:.0%}",
    "permutation-significant primary patient-level association",
    "Observational abundance association; patients are the inferential units.",
    state_file,
)

specific_file = (
    PLAT
    / "receiver_state_characterization"
    / "GSE183904_sender_abundance_vs_recurrent_myeloid_states.tsv"
)
specific = read(specific_file)
r = specific[specific["receiver_state"] == "C1QC_APOE_macrophage"].iloc[0]
add(
    "patient_level_state_specificity",
    "PLAT-rich fibroblast state -> C1QC+APOE+ macrophage state",
    "comparison across four recurrent myeloid states",
    "GSE183904",
    "patient",
    int(r["n"]),
    "covariate-adjusted partial Spearman rho",
    r["adjusted_rho"],
    r["adjusted_permutation_p"],
    r["adjusted_permutation_FDR_BH"],
    f"LOO rho {r['loo_adjusted_rho_min']:.3f}-{r['loo_adjusted_rho_max']:.3f}; positive {r['loo_positive_fraction']:.0%}",
    "four-state FDR-significant",
    "Supports receiver-state specificity within the tested recurrent myeloid states.",
    specific_file,
)

corrected_file = (
    PLAT
    / "corrected_state_abundance_inference"
    / "GSE183904_corrected_state_abundance_inference.tsv"
)
corrected = read(corrected_file)
for r in corrected.itertuples(index=False):
    add(
        "patient_level_conservative_sensitivity",
        "PLAT-rich fibroblast state -> C1QC+APOE+ macrophage state",
        r.representation,
        "GSE183904",
        "patient",
        int(r.n),
        "covariate-adjusted partial Spearman rho",
        r.adjusted_rank_residual_rho,
        r.wild_bootstrap_t_p,
        np.nan,
        f"HC3 P={r.HC3_p:.4f}; LOO rho {r.loo_rho_min:.3f}-{r.loo_rho_max:.3f}; positive {r.loo_positive_fraction:.0%}",
        "directionally stable conservative sensitivity",
        "HC3 and wild-bootstrap tests are conservative sensitivity analyses at n=26.",
        corrected_file,
    )

# Cross-cohort continuous LR ranking.
rank_file = AUDIT / "posthoc_cross_cohort_lr_continuous_ranking_20260814.tsv"
ranked = read(rank_file)
for ligand, receptor in [("APP", "CD74"), ("RARRES2", "CCRL2"), ("SPON2", "ITGB2"), ("CXCL14", "CXCR4"), ("PLAT", "LRP1"), ("PLAT", "ITGAM")]:
    r = ranked[
        (ranked["ligand_complex"] == ligand)
        & (ranked["receptor_complex"] == receptor)
    ].iloc[0]
    add(
        "cross_cohort_ligand_receptor_ranking",
        f"{ligand}-{receptor}",
        "equal-weight continuous rank across four LIANA components",
        "GSE183904 + GSE206785",
        "shared ligand-receptor edge",
        len(ranked),
        "rank among all shared edges",
        int(r["continuous_rank_all_shared_edges"]),
        np.nan,
        np.nan,
        (
            f"ranks D-mag={r['discovery_magnitude_rank']:.4g}, D-spec={r['discovery_specificity_rank']:.4g}, "
            f"V-mag={r['validation_magnitude_rank']:.4g}, V-spec={r['validation_specificity_rank']:.4g}; "
            f"donors {r['discovery_supported_donors']:.0f}/{r['discovery_eligible_donors']:.0f} and "
            f"{r['validation_supported_donors']:.0f}/{r['validation_eligible_donors']:.0f}"
        ),
        "strictly replicated" if bool(r["strict_two_cohort_replication"]) else "post-hoc high-priority or sensitivity candidate",
        "Continuous ranking is post hoc and does not redefine the preregistered strict 0.05 replication gate.",
        rank_file,
    )

# Spatial matched-ligand evidence.
spatial_file = PROJECT / "tables" / "supplementary" / "Table_S7" / "current_source" / "FigureS7F_matched_ligands_source.tsv"
spatial = read(spatial_file)
for ligand in ["SPON2", "PLAT"]:
    for programme in [
        "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
        "GOBP_CELL_MATRIX_ADHESION",
    ]:
        r = spatial[(spatial["ligand"] == ligand) & (spatial["program"] == programme)].iloc[0]
        add(
            "spatial_association",
            ligand,
            programme,
            "spatial transcriptomics cohort",
            "patient section",
            int(r["k"]),
            "pooled Spearman rho",
            r["pooled_rho"],
            r["sign_test_one_sided_p"],
            np.nan,
            f"positive patients {int(r['positive_patients'])}/{int(r['k'])}; matched-ligand rank {int(r['rank_descending'])}",
            "positive spatial association",
            "Spatial co-localization/gradient association; no direct ligand-receptor binding inference.",
            spatial_file,
        )

# Multi-seed virtual perturbation.
seed_file = VK / "GSE183904_C1QC_APOE_scTenifoldKnk_program_enrichment_summary.tsv"
seed = read(seed_file)
seed = seed[seed["target"] == "ITGB2"]
for r in seed.itertuples(index=False):
    add(
        "virtual_perturbation_multiseed",
        "ITGB2",
        r.program,
        "GSE183904 C1QC+APOE+ macrophages",
        "random seed",
        int(r.seeds),
        "median enrichment AUC",
        r.median_enrichment_auc,
        r.maximum_wilcoxon_p,
        r.maximum_wilcoxon_FDR,
        f"AUC range {r.minimum_enrichment_auc:.3f}-{r.maximum_enrichment_auc:.3f}; significant seeds {int(r.seeds_FDR_below_0_05)}/{int(r.seeds)}",
        "stable" if r.seeds_FDR_below_0_05 >= 2 else "unsupported",
        "Unsigned differential-network displacement, not expression downregulation or experimental knockout.",
        seed_file,
    )

# Leave-one-donor-out virtual perturbation.
lodo_file = VK / "primary_receptor_leave_one_donor_out_seed2026072739" / "GSE183904_primary_receptor_LODO_program_summary.tsv"
lodo = read(lodo_file)
lodo = lodo[lodo["target"] == "ITGB2"]
for r in lodo.itertuples(index=False):
    add(
        "virtual_perturbation_leave_one_donor_out",
        "ITGB2",
        r.program,
        "GSE183904 C1QC+APOE+ macrophages",
        "leave-one-patient-out run",
        int(r.leave_one_donor_runs),
        "median enrichment AUC",
        r.median_auc,
        np.nan,
        np.nan,
        f"AUC range {r.minimum_auc:.3f}-{r.maximum_auc:.3f}; AUC>0.5 {r.fraction_auc_above_0_5:.0%}; FDR<0.05 {int(r.runs_FDR_below_0_05)}/{int(r.leave_one_donor_runs)}",
        "LODO stable" if r.fraction_auc_above_0_5 >= 0.8 else "LODO unstable",
        "Computational leave-one-patient-out robustness, not biological replication.",
        lodo_file,
    )

# Directional network connectivity audit.
direction_file = VK / "ITGB2_directional_network_audit" / "ITGB2_WT_KO_program_connectivity_summary.tsv"
direction = read(direction_file)
for r in direction.itertuples(index=False):
    add(
        "virtual_perturbation_directional_connectivity",
        "ITGB2",
        r.program,
        "GSE183904 C1QC+APOE+ macrophages",
        "random seed",
        int(r.seeds),
        "median proportional change in total absolute programme connectivity",
        r.median_percent_change_total_abs_strength,
        np.nan,
        np.nan,
        f"loss in {int(r.seeds_with_strength_loss)}/{int(r.seeds)} seeds; range {r.minimum_percent_change_total_abs_strength:.3f}% to {r.maximum_percent_change_total_abs_strength:.3f}%",
        "consistent inferred connectivity loss",
        "Network-connectivity loss, not transcriptional downregulation or causal inhibition.",
        direction_file,
    )

direction_control_file = VK / "ITGB2_directional_network_audit" / "Figure7G_ITGB2_directional_connectivity_summary.tsv"
direction_control = read(direction_control_file)
for r in direction_control.itertuples(index=False):
    add(
        "virtual_perturbation_directional_matched_controls_20",
        "ITGB2",
        r.program,
        "GSE183904 C1QC+APOE+ macrophages",
        "network-matched control gene",
        20,
        "target percent change in total absolute programme connectivity",
        r.target_percent_change,
        r.empirical_lower_tail_p,
        r.empirical_FDR_BH,
        f"controls with equal/greater loss {int(r.controls_at_least_as_negative)}/20; control median {r.control_median_percent_change:.3f}%",
        "directional loss not target-specific" if r.empirical_FDR_BH >= 0.05 else "directional loss target-specific",
        "Percent units; matched-control directionality test does not measure expression or ECM abundance.",
        direction_control_file,
    )

# Original and expanded matched-control calibration.
old_control_file = VK / "primary_receptor_network_matched_controls_seed2026072739" / "GSE183904_primary_receptors_vs_network_matched_controls.tsv"
old_control = read(old_control_file)
old_control = old_control[old_control["target"] == "ITGB2"]
for r in old_control.itertuples(index=False):
    add(
        "virtual_perturbation_matched_controls_20",
        "ITGB2",
        r.program,
        "GSE183904 C1QC+APOE+ macrophages",
        "network-matched control gene",
        20,
        "target enrichment AUC",
        r.target_auc,
        r.target_empirical_upper_tail_p,
        r.empirical_FDR_BH_within_target,
        f"controls >= target {int(r.controls_at_least_target)}/20; control median AUC {r.matched_control_median_auc:.3f}",
        "prespecified v2 calibration",
        "One-seed empirical calibration; minimum attainable P is 1/21.",
        old_control_file,
    )

expanded_file = VK / "ITGB2_expanded_60_network_matched_controls_seed2026072739_posthoc" / "GSE183904_ITGB2_vs_expanded_60_network_matched_controls.tsv"
expanded = read(expanded_file)
for r in expanded.itertuples(index=False):
    add(
        "virtual_perturbation_matched_controls_60",
        "ITGB2",
        r.program,
        "GSE183904 C1QC+APOE+ macrophages",
        "network-matched control gene",
        60,
        "target enrichment AUC",
        r.target_auc,
        r.target_empirical_upper_tail_p,
        r.empirical_FDR_BH_six_programmes,
        f"controls >= target {int(r.controls_at_least_target)}/60; control median AUC {r.matched_control_median_auc:.3f}; percentile {r.target_percentile_among_61:.1%}",
        "post-hoc exploratory expanded calibration",
        "Control set was frozen before new KO runs; expansion remains post hoc because v2 outcomes were already known.",
        expanded_file,
    )

evidence = pd.DataFrame(rows)
evidence.to_csv(OUTFILE, sep="\t", index=False)

report = """# ITGB2-axis result-rescue findings

1. The patient-level PLAT-rich fibroblast/C1QC+APOE+ macrophage association is significant by the existing adjusted permutation analysis (rho=0.543, P=0.0052, FDR=0.0104) and remains positive in every leave-one-patient-out analysis.
2. SPON2-ITGB2 ranks 5th of 292 shared cross-cohort edges and recurs in 22/22 discovery and 6/6 validation donors. It narrowly misses strict validation magnitude (0.0609) and is therefore a high-priority near-threshold candidate, not a strictly replicated edge.
3. SPON2 has positive spatial associations with ECM degradation (rho=0.0933; 9/9 patients positive; matched-ligand rank 6) and cell-matrix adhesion (rho=0.0685; 8/9 positive; rank 7), exceeding PLAT on both named endpoints.
4. ITGB2 virtual perturbation consistently prioritizes ECM degradation, myeloid migration and phagocytosis across three seeds. Leave-one-patient-out analyses retain AUC>0.5 in 26/26 runs for all three programmes, with FDR<0.05 in 25-26 of 26 runs. The directional audit shows only a small ECM-connectivity decrease (median -0.274% across seeds), which is not specific against 20 matched nodes (P=0.190; FDR=0.786); it therefore does not establish reduced ECM abundance or activity.
5. In the frozen post-hoc expansion to 60 matched controls, ITGB2 ranks first for migration (P=0.0164) and second for ECM degradation (P=0.0328). Six-programme BH correction yields FDR=0.0984 for both endpoints, supporting exploratory prioritization but not multiplicity-adjusted significance.
6. The defensible main axis is: PLAT-rich/CXCL14+POSTN+ fibroblast state -> SPON2-ITGB2 communication context -> ITGB2-sensitive macrophage ECM-remodelling network. PLAT remains a sender-state marker; direct PLAT-ITGB2 engagement is not supported.

Suggested next analyses/experiments:

1. Update Figure 4 to foreground the frozen patient-level permutation P and FDR, with HC3/wild bootstrap retained as sensitivity analyses.
2. Present the target-agnostic continuous LR ranking as exploratory support and retain the preregistered strict replication calls for APP-CD74 and RARRES2-CCRL2.
3. Present the 60-control ITGB2 calibration transparently with raw empirical P and six-programme FDR; do not add more controls or redefine endpoints after seeing these results.
4. For a future independent validation, prospectively declare ECM degradation as the primary perturbation endpoint and repeat ITGB2 calibration in another macrophage dataset or an independently frozen seed/control set.
5. Experimental follow-up should test SPON2-ITGB2 proximity/engagement and ITGB2-dependent macrophage migration and ECM-remodelling phenotypes; current computational data do not establish direct binding or causality.
"""
REPORT.write_text(report, encoding="utf-8")

print(f"Wrote {len(evidence)} evidence rows to {OUTFILE}")
print(f"Wrote findings to {REPORT}")
