#!/usr/bin/env python
"""All-state PLAT-directed communication competition in GSE206785."""

from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import anndata as ad
import liana as li
import numpy as np
import pandas as pd


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
PLAN_FILE = PROJECT / "analysis" / "audit" / "GSE206785_PLAT_all_myeloid_communication_exploratory_plan.json"
H5AD = (
    PROJECT
    / "analysis"
    / "processed"
    / "cell_communication"
    / "GSE206785_tumor_fibroblast_myeloid_liana_input.h5ad"
)
MAPPING = (
    PROJECT
    / "analysis"
    / "results"
    / "PLAT_axis"
    / "GSE183904"
    / "fine_state_communication_exploratory"
    / "cross_cohort_state_mapping"
    / "cross_cohort_state_similarity.tsv"
)
OUT = (
    PROJECT
    / "analysis"
    / "results"
    / "PLAT_axis"
    / "GSE206785"
    / "all_myeloid_communication_exploratory"
)
OUT.mkdir(parents=True, exist_ok=True)

SENDER = "F13_Activated_CTHRC1"
GENES = ["PLAT", "ITGAM", "ITGB2", "LRP1"]
RESOURCE = pd.DataFrame(
    {"ligand": ["PLAT", "PLAT", "PLAT"], "receptor": ["ITGAM", "ITGB2", "LRP1"]}
)


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def run_liana(
    adata: ad.AnnData,
    receivers: list[str],
    threshold: float,
    plan: dict[str, object],
) -> pd.DataFrame:
    selected = adata[adata.obs.Annotation.astype(str).isin([SENDER] + receivers)].copy()
    selected.obs["Annotation"] = selected.obs.Annotation.astype("category")
    pairs = pd.DataFrame([(SENDER, x) for x in receivers], columns=["source", "target"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = li.mt.rank_aggregate(
            selected,
            groupby="Annotation",
            groupby_pairs=pairs,
            resource=RESOURCE,
            expr_prop=threshold,
            min_cells=int(plan["thresholds"]["minimum_cells_per_patient_state"]),
            use_raw=False,
            n_perms=int(plan["thresholds"]["permutations"]),
            seed=int(plan["thresholds"]["random_seed"]),
            n_jobs=4,
            return_all_lrs=True,
            inplace=False,
            verbose=False,
        )
    result = result.copy()
    result.insert(0, "expr_prop", threshold)
    result["combined_worst_rank"] = result[["magnitude_rank", "specificity_rank"]].max(axis=1)
    result["combined_geometric_rank"] = np.sqrt(
        result.magnitude_rank.clip(lower=np.finfo(float).tiny)
        * result.specificity_rank.clip(lower=np.finfo(float).tiny)
    )
    return result


def patient_potentials(
    adata: ad.AnnData,
    receivers: list[str],
    plan: dict[str, object],
) -> pd.DataFrame:
    min_cells = int(plan["thresholds"]["minimum_cells_per_patient_state"])
    obs = adata.obs.copy()
    x = adata.X
    gene_index = {gene: adata.var_names.get_loc(gene) for gene in GENES}
    annotation = obs.Annotation.astype(str).to_numpy()
    patient_values = obs.Patient.astype(str).to_numpy()
    records = []
    for patient in sorted(np.unique(patient_values)):
        sender_ix = np.flatnonzero((patient_values == patient) & (annotation == SENDER))
        if len(sender_ix) < min_cells:
            continue
        s_expr = np.asarray(x[sender_ix, gene_index["PLAT"]].todense()).ravel()
        s_mean = float(s_expr.mean())
        s_det = float(np.mean(s_expr > 0))
        for receiver in receivers:
            receiver_ix = np.flatnonzero((patient_values == patient) & (annotation == receiver))
            if len(receiver_ix) < min_cells:
                continue
            receiver_stats = {}
            for gene in ["ITGAM", "ITGB2", "LRP1"]:
                values = np.asarray(x[receiver_ix, gene_index[gene]].todense()).ravel()
                receiver_stats[gene] = {
                    "mean": float(values.mean()),
                    "detection": float(np.mean(values > 0)),
                }
            branches = {
                "alphaMbeta2": {
                    "mean": min(receiver_stats["ITGAM"]["mean"], receiver_stats["ITGB2"]["mean"]),
                    "detection": min(
                        receiver_stats["ITGAM"]["detection"],
                        receiver_stats["ITGB2"]["detection"],
                    ),
                },
                "LRP1": receiver_stats["LRP1"],
            }
            for branch, stats in branches.items():
                records.append(
                    {
                        "Patient": patient,
                        "receiver": receiver,
                        "branch": branch,
                        "n_sender_cells": len(sender_ix),
                        "n_receiver_cells": len(receiver_ix),
                        "sender_PLAT_mean": s_mean,
                        "sender_PLAT_detection": s_det,
                        "receiver_mean": stats["mean"],
                        "receiver_detection": stats["detection"],
                        "expression_potential": float(np.sqrt(max(s_mean, 0) * max(stats["mean"], 0))),
                        "detection_potential": float(s_det * stats["detection"]),
                        "gate_5pct": bool(s_det >= 0.05 and stats["detection"] >= 0.05),
                        "gate_10pct": bool(s_det >= 0.10 and stats["detection"] >= 0.10),
                    }
                )
    answer = pd.DataFrame(records)
    answer["within_patient_rank"] = answer.groupby(
        ["Patient", "branch"], observed=True
    ).expression_potential.rank(method="average", ascending=False)
    n_available = answer.groupby(["Patient", "branch"], observed=True).receiver.transform("size")
    answer["within_patient_rank_fraction"] = answer.within_patient_rank / n_available
    answer["is_top_within_patient"] = answer.within_patient_rank == 1
    answer.to_csv(
        OUT / "GSE206785_per_patient_PLAT_expression_potentials.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    return answer


def main() -> None:
    plan = json.loads(PLAN_FILE.read_text(encoding="utf-8"))
    adata = ad.read_h5ad(H5AD)
    if not set(GENES).issubset(adata.var_names):
        raise RuntimeError("PLAT resource genes missing from GSE206785 communication matrix")
    obs = adata.obs
    myeloid = obs[obs.Type.astype(str) == "Myeloid"]
    counts = myeloid.Annotation.astype(str).value_counts()
    patients = myeloid.groupby("Annotation", observed=True).Patient.nunique()
    receivers = sorted(
        state
        for state in counts.index
        if counts[state] >= int(plan["thresholds"]["minimum_total_state_cells"])
        and patients.get(state, 0) >= int(plan["thresholds"]["minimum_state_patients"])
    )
    log(f"eligible receiver states: {len(receivers)}")

    pooled_parts = []
    for threshold in [
        float(plan["thresholds"]["primary_expression_proportion"]),
        float(plan["thresholds"]["conventional_expression_proportion_sensitivity"]),
    ]:
        log(f"pooled LIANA expr_prop={threshold:.2f}")
        pooled_parts.append(run_liana(adata, receivers, threshold, plan))
    pooled = pd.concat(pooled_parts, ignore_index=True)
    pooled.to_csv(OUT / "GSE206785_all_myeloid_PLAT_directed_LIANA.tsv", sep="\t", index=False)

    potentials = patient_potentials(adata, receivers, plan)
    summary = (
        potentials.groupby(["receiver", "branch"], observed=True)
        .agg(
            eligible_patients=("Patient", "nunique"),
            median_expression_potential=("expression_potential", "median"),
            median_detection_potential=("detection_potential", "median"),
            gate_5pct_fraction=("gate_5pct", "mean"),
            gate_10pct_fraction=("gate_10pct", "mean"),
            mean_within_patient_rank_fraction=("within_patient_rank_fraction", "mean"),
            top_within_patient_fraction=("is_top_within_patient", "mean"),
        )
        .reset_index()
    )
    state_sizes = (
        myeloid.groupby("Annotation", observed=True)
        .agg(
            total_cells=("Annotation", "size"),
            total_patients=("Patient", "nunique"),
            maximum_patient_fraction=("Patient", lambda x: x.value_counts().max() / len(x)),
        )
        .reset_index()
        .rename(columns={"Annotation": "receiver"})
    )
    summary = summary.merge(state_sizes, on="receiver", how="left")

    pooled5 = pooled[np.isclose(pooled.expr_prop, 0.05)].copy()
    pooled_state = (
        pooled5.groupby("target", observed=True)
        .agg(
            pooled_best_magnitude_rank=("magnitude_rank", "min"),
            pooled_best_specificity_rank=("specificity_rank", "min"),
            pooled_best_combined_worst_rank=("combined_worst_rank", "min"),
        )
        .reset_index()
        .rename(columns={"target": "receiver"})
    )
    summary = summary.merge(pooled_state, on="receiver", how="left")

    mapping = pd.read_csv(MAPPING, sep="\t")
    best_mapping = (
        mapping.sort_values(["query_state", "spearman_rho"], ascending=[True, False])
        .groupby("query_state", as_index=False)
        .head(1)[
            ["query_state", "reference_state", "spearman_rho", "bootstrap_top_probability"]
        ]
        .rename(
            columns={
                "query_state": "receiver",
                "reference_state": "best_GSE183904_match",
                "spearman_rho": "mapping_spearman_rho",
                "bootstrap_top_probability": "mapping_bootstrap_top_probability",
            }
        )
    )
    summary = summary.merge(best_mapping, on="receiver", how="left")
    summary.to_csv(OUT / "GSE206785_PLAT_receiver_state_summary.tsv", sep="\t", index=False)

    leaders = {}
    for branch, block in summary.groupby("branch", observed=True):
        eligible = block[block.eligible_patients >= 5].copy()
        leader = eligible.sort_values(
            ["mean_within_patient_rank_fraction", "pooled_best_combined_worst_rank"],
            ascending=[True, True],
        ).iloc[0]
        leaders[branch] = {
            "receiver": leader.receiver,
            "eligible_patients": int(leader.eligible_patients),
            "best_GSE183904_match": leader.best_GSE183904_match,
            "mean_within_patient_rank_fraction": float(leader.mean_within_patient_rank_fraction),
            "gate_5pct_fraction": float(leader.gate_5pct_fraction),
            "gate_10pct_fraction": float(leader.gate_10pct_fraction),
        }
    m02 = summary[summary.receiver == "M02_Mono_CD14"]
    result_summary = {
        "analysis_id": plan["analysis_id"],
        "eligible_receiver_states": receivers,
        "sender_cells": int((obs.Annotation.astype(str) == SENDER).sum()),
        "sender_patients": int(obs.loc[obs.Annotation.astype(str) == SENDER, "Patient"].nunique()),
        "branch_leaders": leaders,
        "M02_branch_rows": m02.to_dict(orient="records"),
        "claim_boundary": plan["claim_boundary"],
    }
    (OUT / "analysis_summary.json").write_text(
        json.dumps(result_summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    log(f"complete: {OUT}")


if __name__ == "__main__":
    main()
