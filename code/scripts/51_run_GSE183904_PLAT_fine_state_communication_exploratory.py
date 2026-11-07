#!/usr/bin/env python
"""Exploratory, patient-aware PLAT-directed fine-state communication audit.

This script does not recluster, relabel, modify signatures, or draw main figures.
It compares the frozen PLAT-rich fibroblast sender against every reviewed
GSE183904 myeloid state. Cluster 2 is retained only as a donor-skewed
sensitivity state. LIANA results are co-expression proxies, not evidence of
binding, directionality, or causality.
"""

from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path

import anndata as ad
import liana as li
import numpy as np
import pandas as pd


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
AUDIT = PROJECT / "analysis" / "audit" / "GSE183904_PLAT_fine_state_communication_exploratory_plan.json"
INPUT = PROJECT / "analysis" / "processed" / "GSE183904_fine_states"
REVIEW = PROJECT / "analysis" / "audit" / "GSE183904_fine_state_label_review.tsv"
OUT = PROJECT / "analysis" / "results" / "PLAT_axis" / "GSE183904" / "fine_state_communication_exploratory"
OUT.mkdir(parents=True, exist_ok=True)

FIB_FILE = INPUT / "GSE183904_Fibroblast_fine_state_cell_metadata.tsv.gz"
MYE_FILE = INPUT / "GSE183904_Myeloid_fine_state_cell_metadata.tsv.gz"
ECOLOGY_FILE = (
    PROJECT
    / "analysis"
    / "results"
    / "PLAT_axis"
    / "GSE183904"
    / "receiver_state_characterization"
    / "GSE183904_sender_abundance_vs_recurrent_myeloid_states.tsv"
)
PROGRAM_FILE = (
    PROJECT
    / "analysis"
    / "results"
    / "PLAT_axis"
    / "GSE183904"
    / "receiver_state_characterization"
    / "GSE183904_C1QC_APOE_program_contrasts.tsv"
)

GENES = ["PLAT", "ITGAM", "ITGB2", "LRP1"]
RESOURCE = pd.DataFrame(
    {"ligand": ["PLAT", "PLAT", "PLAT"], "receptor": ["ITGAM", "ITGB2", "LRP1"]}
)
SENDER_GROUP = "F1_PLATrich"
PRIMARY_RECEIVERS = ["C1", "C3", "C4", "C5"]
ALL_RECEIVERS = ["C1", "C2", "C3", "C4", "C5"]


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def bh_fdr(pvalues: pd.Series) -> np.ndarray:
    p = pvalues.to_numpy(dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(ranked, 1.0)
    return out


def sign_flip_p(diff: np.ndarray, n_perms: int, rng: np.random.Generator) -> float:
    """One-sided patient-blocked sign-flip P for a positive mean difference."""
    diff = np.asarray(diff, dtype=float)
    diff = diff[np.isfinite(diff)]
    if len(diff) == 0:
        return float("nan")
    observed = float(diff.mean())
    signs = rng.choice(np.array([-1.0, 1.0]), size=(n_perms, len(diff)), replace=True)
    permuted = (signs * diff[None, :]).mean(axis=1)
    return float((1 + np.sum(permuted >= observed)) / (n_perms + 1))


def bootstrap_ci(
    values: np.ndarray, n_boot: int, rng: np.random.Generator, statistic: str = "mean"
) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return float("nan"), float("nan")
    draws = rng.choice(values, size=(n_boot, len(values)), replace=True)
    estimates = draws.mean(axis=1) if statistic == "mean" else np.median(draws, axis=1)
    lo, hi = np.quantile(estimates, [0.025, 0.975])
    return float(lo), float(hi)


def read_inputs() -> tuple[pd.DataFrame, dict[str, object]]:
    plan = json.loads(AUDIT.read_text(encoding="utf-8"))
    fib = pd.read_csv(FIB_FILE, sep="\t")
    mye = pd.read_csv(MYE_FILE, sep="\t")
    review = pd.read_csv(REVIEW, sep="\t")

    sender_review = review[(review.lineage == "Fibroblast") & (review.cluster == 1)]
    if len(sender_review) != 1 or sender_review.iloc[0].refined_label != "CXCL14_POSTN_PLAT_rich_fibroblast":
        raise RuntimeError("Frozen fibroblast cluster-1 label review is missing or changed")
    myeloid_review = review[review.lineage == "Myeloid"].set_index("cluster")
    if set(myeloid_review.index.astype(int)) != {1, 2, 3, 4, 5}:
        raise RuntimeError("Expected reviewed myeloid clusters 1-5")

    sender = fib[fib.cluster == 1].copy()
    sender["ccc_group"] = SENDER_GROUP
    sender["refined_label"] = "CXCL14_POSTN_PLAT_rich_fibroblast"
    mye = mye[mye.cluster.isin([1, 2, 3, 4, 5])].copy()
    mye["ccc_group"] = "C" + mye.cluster.astype(str)
    mye["refined_label"] = mye.cluster.map(myeloid_review.refined_label)

    keep = [
        "cell_id",
        "Patient_ID",
        "lineage",
        "cluster",
        "ccc_group",
        "refined_label",
    ] + [f"expr_{g}" for g in GENES] + [f"detect_{g}" for g in GENES]
    data = pd.concat([sender[keep], mye[keep]], ignore_index=True)
    if data.cell_id.duplicated().any():
        raise RuntimeError("Duplicated cell identifiers after sender/receiver assembly")
    if data[[f"expr_{g}" for g in GENES]].isna().any().any():
        raise RuntimeError("Missing PLAT-axis expression values")
    return data, plan


def to_anndata(data: pd.DataFrame) -> ad.AnnData:
    obs = data[["cell_id", "Patient_ID", "ccc_group", "refined_label"]].copy()
    obs = obs.set_index("cell_id", drop=False)
    obs["ccc_group"] = obs.ccc_group.astype("category")
    x = data[[f"expr_{g}" for g in GENES]].to_numpy(dtype=np.float32)
    return ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=pd.Index(GENES, name="gene")))


def run_liana(
    data: pd.DataFrame,
    receivers: list[str],
    expr_prop: float,
    n_perms: int,
    min_cells: int,
    seed: int,
    n_jobs: int,
) -> pd.DataFrame:
    groups = [SENDER_GROUP] + receivers
    selected = data[data.ccc_group.isin(groups)].copy()
    present = selected.groupby("ccc_group").size()
    receivers = [x for x in receivers if present.get(x, 0) >= min_cells]
    if present.get(SENDER_GROUP, 0) < min_cells or not receivers:
        return pd.DataFrame()
    selected = selected[selected.ccc_group.isin([SENDER_GROUP] + receivers)].copy()
    adata = to_anndata(selected)
    pairs = pd.DataFrame([(SENDER_GROUP, x) for x in receivers], columns=["source", "target"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = li.mt.rank_aggregate(
            adata,
            groupby="ccc_group",
            groupby_pairs=pairs,
            resource=RESOURCE,
            expr_prop=expr_prop,
            min_cells=min_cells,
            use_raw=False,
            n_perms=n_perms,
            seed=seed,
            n_jobs=n_jobs,
            return_all_lrs=True,
            inplace=False,
            verbose=False,
        )
    if result is None or result.empty:
        return pd.DataFrame()
    result = result.copy()
    result["combined_worst_rank"] = result[["magnitude_rank", "specificity_rank"]].max(axis=1)
    result["combined_geometric_rank"] = np.sqrt(
        result.magnitude_rank.clip(lower=np.finfo(float).tiny)
        * result.specificity_rank.clip(lower=np.finfo(float).tiny)
    )
    return result


def pooled_liana(data: pd.DataFrame, plan: dict[str, object]) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    thresholds = [
        float(plan["thresholds"]["primary_expression_proportion"]),
        float(plan["thresholds"]["conventional_expression_proportion_sensitivity"]),
    ]
    for threshold in thresholds:
        for analysis_set, receivers in [
            ("primary_recurrent", PRIMARY_RECEIVERS),
            ("all_with_C2_sensitivity", ALL_RECEIVERS),
        ]:
            log(f"pooled LIANA: {analysis_set}, expr_prop={threshold:.2f}")
            result = run_liana(
                data,
                receivers,
                threshold,
                int(plan["thresholds"]["permutations"]),
                int(plan["thresholds"]["minimum_cells_per_patient_state"]),
                int(plan["thresholds"]["random_seed"]),
                4,
            )
            result.insert(0, "analysis_set", analysis_set)
            result.insert(1, "expr_prop", threshold)
            rows.append(result)
    answer = pd.concat(rows, ignore_index=True)
    answer.to_csv(OUT / "pooled_PLAT_directed_LIANA_results.tsv", sep="\t", index=False)
    return answer


def per_patient_liana(data: pd.DataFrame, plan: dict[str, object]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[pd.DataFrame] = []
    errors: list[dict[str, object]] = []
    thresholds = [
        float(plan["thresholds"]["primary_expression_proportion"]),
        float(plan["thresholds"]["conventional_expression_proportion_sensitivity"]),
    ]
    min_cells = int(plan["thresholds"]["minimum_cells_per_patient_state"])
    n_perms = int(plan["thresholds"]["permutations"])
    seed = int(plan["thresholds"]["random_seed"])
    counts = data.groupby(["Patient_ID", "ccc_group"], observed=True).size().unstack(fill_value=0)
    eligible = counts.index[counts.get(SENDER_GROUP, pd.Series(0, index=counts.index)) >= min_cells]
    log(f"patient-level LIANA: {len(eligible)} sender-eligible patients")
    for patient_index, patient in enumerate(sorted(eligible)):
        patient_data = data[data.Patient_ID == patient].copy()
        available = [x for x in PRIMARY_RECEIVERS if counts.loc[patient].get(x, 0) >= min_cells]
        if not available:
            continue
        for threshold in thresholds:
            try:
                result = run_liana(
                    patient_data,
                    available,
                    threshold,
                    n_perms,
                    min_cells,
                    seed + patient_index,
                    1,
                )
                if result.empty:
                    errors.append(
                        {"Patient_ID": patient, "expr_prop": threshold, "error": "empty_result"}
                    )
                    continue
                result.insert(0, "Patient_ID", patient)
                result.insert(1, "expr_prop", threshold)
                rows.append(result)
            except Exception as exc:  # output is audited; failures are not silently discarded
                errors.append(
                    {"Patient_ID": patient, "expr_prop": threshold, "error": repr(exc)}
                )
    answer = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    error_table = pd.DataFrame(errors, columns=["Patient_ID", "expr_prop", "error"])
    answer.to_csv(OUT / "per_patient_liana_results.tsv.gz", sep="\t", index=False, compression="gzip")
    error_table.to_csv(OUT / "per_patient_liana_errors.tsv", sep="\t", index=False)
    return answer, error_table


def expression_potentials(data: pd.DataFrame, plan: dict[str, object]) -> pd.DataFrame:
    min_cells = int(plan["thresholds"]["minimum_cells_per_patient_state"])
    records: list[dict[str, object]] = []
    sender = data[data.ccc_group == SENDER_GROUP]
    receivers = data[data.ccc_group.isin(ALL_RECEIVERS)]
    for patient in sorted(data.Patient_ID.unique()):
        s = sender[sender.Patient_ID == patient]
        if len(s) < min_cells:
            continue
        s_mean = float(s.expr_PLAT.mean())
        s_det = float(s.detect_PLAT.mean())
        for receiver in ALL_RECEIVERS:
            r = receivers[(receivers.Patient_ID == patient) & (receivers.ccc_group == receiver)]
            if len(r) < min_cells:
                continue
            receptor_values = {
                "alphaMbeta2": {
                    "mean": float(min(r.expr_ITGAM.mean(), r.expr_ITGB2.mean())),
                    "detection": float(min(r.detect_ITGAM.mean(), r.detect_ITGB2.mean())),
                },
                "LRP1": {
                    "mean": float(r.expr_LRP1.mean()),
                    "detection": float(r.detect_LRP1.mean()),
                },
            }
            for branch, values in receptor_values.items():
                records.append(
                    {
                        "Patient_ID": patient,
                        "receiver": receiver,
                        "receiver_label": r.refined_label.iloc[0],
                        "primary_receiver": receiver in PRIMARY_RECEIVERS,
                        "branch": branch,
                        "n_sender_cells": len(s),
                        "n_receiver_cells": len(r),
                        "sender_PLAT_mean": s_mean,
                        "sender_PLAT_detection": s_det,
                        "receiver_mean": values["mean"],
                        "receiver_detection": values["detection"],
                        "expression_potential": float(np.sqrt(max(s_mean, 0) * max(values["mean"], 0))),
                        "detection_potential": float(s_det * values["detection"]),
                        "gate_5pct": bool(s_det >= 0.05 and values["detection"] >= 0.05),
                        "gate_10pct": bool(s_det >= 0.10 and values["detection"] >= 0.10),
                    }
                )
    answer = pd.DataFrame(records)
    answer.to_csv(
        OUT / "per_patient_PLAT_expression_potentials.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    return answer


def summarize_potentials(potentials: pd.DataFrame, plan: dict[str, object]) -> pd.DataFrame:
    rng = np.random.default_rng(int(plan["thresholds"]["random_seed"]))
    n_perms = int(plan["thresholds"]["permutations"])
    n_boot = int(plan["thresholds"]["bootstrap_replicates"])
    recurrent = potentials[potentials.primary_receiver].copy()
    recurrent["within_patient_rank"] = recurrent.groupby(
        ["Patient_ID", "branch"], observed=True
    ).expression_potential.rank(method="average", ascending=False)
    recurrent["is_top_within_patient"] = recurrent.within_patient_rank == 1
    other_mean = recurrent.groupby(["Patient_ID", "branch"], observed=True).expression_potential.transform(
        lambda x: (x.sum() - x) / (len(x) - 1)
    )
    recurrent["difference_vs_other_state_mean"] = recurrent.expression_potential - other_mean

    rows: list[dict[str, object]] = []
    for (receiver, branch), block in recurrent.groupby(["receiver", "branch"], observed=True):
        diff = block.difference_vs_other_state_mean.to_numpy(dtype=float)
        ci_lo, ci_hi = bootstrap_ci(diff, n_boot, rng, statistic="mean")
        rows.append(
            {
                "receiver": receiver,
                "receiver_label": block.receiver_label.iloc[0],
                "primary_receiver": True,
                "branch": branch,
                "eligible_patients": block.Patient_ID.nunique(),
                "median_expression_potential": block.expression_potential.median(),
                "median_detection_potential": block.detection_potential.median(),
                "gate_5pct_fraction": block.gate_5pct.mean(),
                "gate_10pct_fraction": block.gate_10pct.mean(),
                "mean_within_patient_rank": block.within_patient_rank.mean(),
                "top_within_patient_fraction": block.is_top_within_patient.mean(),
                "mean_difference_vs_other_states": diff.mean(),
                "mean_difference_bootstrap_ci_low": ci_lo,
                "mean_difference_bootstrap_ci_high": ci_hi,
                "direction_consistency": float(np.mean(diff > 0)),
                "one_sided_signflip_p": sign_flip_p(diff, n_perms, rng),
            }
        )

    # C2 is explicitly descriptive and never enters the recurrent-state null.
    c2 = potentials[potentials.receiver == "C2"]
    for branch, block in c2.groupby("branch", observed=True):
        rows.append(
            {
                "receiver": "C2",
                "receiver_label": block.receiver_label.iloc[0],
                "primary_receiver": False,
                "branch": branch,
                "eligible_patients": block.Patient_ID.nunique(),
                "median_expression_potential": block.expression_potential.median(),
                "median_detection_potential": block.detection_potential.median(),
                "gate_5pct_fraction": block.gate_5pct.mean(),
                "gate_10pct_fraction": block.gate_10pct.mean(),
                "mean_within_patient_rank": np.nan,
                "top_within_patient_fraction": np.nan,
                "mean_difference_vs_other_states": np.nan,
                "mean_difference_bootstrap_ci_low": np.nan,
                "mean_difference_bootstrap_ci_high": np.nan,
                "direction_consistency": np.nan,
                "one_sided_signflip_p": np.nan,
            }
        )
    answer = pd.DataFrame(rows)
    answer["signflip_FDR_BH"] = np.nan
    primary_index = answer.index[answer.primary_receiver]
    answer.loc[primary_index, "signflip_FDR_BH"] = bh_fdr(
        answer.loc[primary_index, "one_sided_signflip_p"]
    )
    answer.to_csv(OUT / "PLAT_edge_receiver_summary.tsv", sep="\t", index=False)
    return answer


def summarize_patient_liana(results: pd.DataFrame) -> pd.DataFrame:
    if results.empty:
        answer = pd.DataFrame()
    else:
        z = results.copy()
        z["dual_rank_lt_0_10"] = (z.magnitude_rank < 0.10) & (z.specificity_rank < 0.10)
        answer = (
            z.groupby(["expr_prop", "target", "receptor_complex"], observed=True)
            .agg(
                eligible_patients=("Patient_ID", "nunique"),
                median_magnitude_rank=("magnitude_rank", "median"),
                median_specificity_rank=("specificity_rank", "median"),
                median_combined_worst_rank=("combined_worst_rank", "median"),
                dual_rank_lt_0_10_fraction=("dual_rank_lt_0_10", "mean"),
            )
            .reset_index()
        )
    answer.to_csv(OUT / "per_patient_LIANA_receiver_summary.tsv", sep="\t", index=False)
    return answer


def build_threshold_sensitivity(
    pooled: pd.DataFrame, potential_summary: pd.DataFrame
) -> pd.DataFrame:
    primary = pooled[pooled.analysis_set == "primary_recurrent"].copy()
    primary = primary[
        [
            "expr_prop",
            "target",
            "receptor_complex",
            "magnitude_rank",
            "specificity_rank",
            "combined_worst_rank",
        ]
    ].rename(columns={"target": "receiver"})
    gates = potential_summary[
        ["receiver", "branch", "gate_5pct_fraction", "gate_10pct_fraction"]
    ].copy()
    rows = []
    for threshold in [0.05, 0.10]:
        x = primary[np.isclose(primary.expr_prop, threshold)].copy()
        x["expression_gate_support_fraction"] = np.nan
        for idx, row in x.iterrows():
            branch = "LRP1" if row.receptor_complex == "LRP1" else "alphaMbeta2"
            match = gates[(gates.receiver == row.receiver) & (gates.branch == branch)]
            if not match.empty:
                col = "gate_5pct_fraction" if threshold == 0.05 else "gate_10pct_fraction"
                x.loc[idx, "expression_gate_support_fraction"] = match.iloc[0][col]
        rows.append(x)
    answer = pd.concat(rows, ignore_index=True)
    answer.to_csv(OUT / "threshold_sensitivity_5pct_vs_10pct.tsv", sep="\t", index=False)
    return answer


def build_receiver_ranking(
    potential_summary: pd.DataFrame,
    pooled: pd.DataFrame,
    plan: dict[str, object],
) -> pd.DataFrame:
    primary = potential_summary[potential_summary.primary_receiver].copy()
    rank = (
        primary.groupby("receiver", observed=True)
        .agg(
            receiver_label=("receiver_label", "first"),
            mean_rank_across_branches=("mean_within_patient_rank", "mean"),
            mean_top_fraction_across_branches=("top_within_patient_fraction", "mean"),
            minimum_direction_consistency=("direction_consistency", "min"),
            significant_branches_FDR_lt_0_10=("signflip_FDR_BH", lambda x: int((x < 0.10).sum())),
            branches_gate5_consistency_ge_0_60=("gate_5pct_fraction", lambda x: int((x >= 0.60).sum())),
        )
        .reset_index()
    )

    pooled5 = pooled[
        (pooled.analysis_set == "primary_recurrent") & np.isclose(pooled.expr_prop, 0.05)
    ].copy()
    pooled_score = (
        pooled5.groupby("target", observed=True)
        .agg(
            pooled_best_magnitude_rank=("magnitude_rank", "min"),
            pooled_best_specificity_rank=("specificity_rank", "min"),
            pooled_best_combined_worst_rank=("combined_worst_rank", "min"),
        )
        .reset_index()
        .rename(columns={"target": "receiver"})
    )
    rank = rank.merge(pooled_score, on="receiver", how="left")

    if ECOLOGY_FILE.exists():
        ecology = pd.read_csv(ECOLOGY_FILE, sep="\t")
        label_to_receiver = dict(zip(rank.receiver_label, rank.receiver))
        ecology["receiver"] = ecology.receiver_state.map(label_to_receiver)
        ecology = ecology[
            ["receiver", "adjusted_rho", "adjusted_permutation_p", "adjusted_permutation_FDR_BH"]
        ].rename(
            columns={
                "adjusted_rho": "ecology_adjusted_rho",
                "adjusted_permutation_p": "ecology_permutation_p",
                "adjusted_permutation_FDR_BH": "ecology_permutation_FDR_BH",
            }
        )
        rank = rank.merge(ecology, on="receiver", how="left")

    rank = rank.sort_values(
        ["mean_rank_across_branches", "pooled_best_combined_worst_rank"],
        ascending=[True, True],
    ).reset_index(drop=True)
    rank.insert(0, "exploratory_rank", np.arange(1, len(rank) + 1))
    rank["passes_relaxed_direction_gate"] = (
        rank.minimum_direction_consistency
        >= float(plan["thresholds"]["minimum_direction_consistency"])
    )
    rank.to_csv(OUT / "receiver_state_exploratory_ranking.tsv", sep="\t", index=False)
    return rank


def write_summary(
    data: pd.DataFrame,
    plan: dict[str, object],
    pooled: pd.DataFrame,
    per_patient: pd.DataFrame,
    errors: pd.DataFrame,
    potential_summary: pd.DataFrame,
    ranking: pd.DataFrame,
) -> None:
    recurrent = potential_summary[potential_summary.primary_receiver]
    branch_leaders = {}
    for branch, block in recurrent.groupby("branch", observed=True):
        leader = block.sort_values(
            ["mean_within_patient_rank", "median_expression_potential"],
            ascending=[True, False],
        ).iloc[0]
        branch_leaders[branch] = {
            "receiver": leader.receiver,
            "receiver_label": leader.receiver_label,
            "mean_within_patient_rank": float(leader.mean_within_patient_rank),
            "direction_consistency": float(leader.direction_consistency),
            "signflip_FDR_BH": float(leader.signflip_FDR_BH),
        }
    unique_leaders = {x["receiver"] for x in branch_leaders.values()}
    if len(unique_leaders) == 1:
        decision = "single_receiver_candidate"
    else:
        decision = "distributed_C1_C3_receiver_landscape" if unique_leaders <= {"C1", "C3"} else "distributed_receiver_landscape"

    programme_note = None
    if PROGRAM_FILE.exists():
        programs = pd.read_csv(PROGRAM_FILE, sep="\t")
        strong = programs[
            (programs.positive_patients == programs.n_patients)
            & (programs.paired_wilcoxon_FDR_BH < 0.10)
        ]
        programme_note = {
            "C3_all_patient_positive_programmes": strong.program.tolist(),
            "count": int(len(strong)),
            "interpretation": "Orthogonal C3 receiver-state characterization; not used to define clusters or tune communication thresholds.",
        }

    summary = {
        "analysis_id": plan["analysis_id"],
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "software": {
            "python": sys.version,
            "liana": li.__version__,
            "anndata": ad.__version__,
            "pandas": pd.__version__,
            "numpy": np.__version__,
        },
        "input_cells": int(len(data)),
        "input_patients": int(data.Patient_ID.nunique()),
        "sender_cells": int((data.ccc_group == SENDER_GROUP).sum()),
        "pooled_liana_rows": int(len(pooled)),
        "per_patient_liana_rows": int(len(per_patient)),
        "per_patient_liana_error_rows": int(len(errors)),
        "branch_leaders": branch_leaders,
        "decision": decision,
        "top_integrated_receiver": ranking.iloc[0].receiver if not ranking.empty else None,
        "programme_note": programme_note,
        "claim_boundary": "PLAT-directed results are OmniPath-supported ligand-receptor co-expression potentials. They do not establish binding, spatial contact, directionality, or causality.",
        "figure_implication": "Show all recurrent states before highlighting any receiver; retain C2 only as a visibly donor-skewed sensitivity state.",
    }
    (OUT / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )


def main() -> None:
    t0 = time.time()
    data, plan = read_inputs()
    log(f"assembled {len(data)} cells from {data.Patient_ID.nunique()} patients")
    pooled = pooled_liana(data, plan)
    per_patient, errors = per_patient_liana(data, plan)
    summarize_patient_liana(per_patient)
    potentials = expression_potentials(data, plan)
    potential_summary = summarize_potentials(potentials, plan)
    build_threshold_sensitivity(pooled, potential_summary)
    ranking = build_receiver_ranking(potential_summary, pooled, plan)
    write_summary(data, plan, pooled, per_patient, errors, potential_summary, ranking)
    log(f"complete in {time.time() - t0:.1f} seconds; outputs: {OUT}")


if __name__ == "__main__":
    main()
