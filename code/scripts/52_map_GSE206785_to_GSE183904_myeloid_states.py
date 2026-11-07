#!/usr/bin/env python
"""Axis-excluded, donor-equal cross-cohort myeloid state mapping."""

from __future__ import annotations

import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import rankdata, spearmanr


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
PLAN_FILE = PROJECT / "analysis" / "audit" / "GSE206785_to_GSE183904_myeloid_state_mapping_plan.json"
INPUT183 = PROJECT / "analysis" / "processed" / "communication" / "GSE183904"
INPUT206 = (
    PROJECT
    / "analysis"
    / "processed"
    / "cell_communication"
    / "GSE206785_tumor_fibroblast_myeloid_liana_input.h5ad"
)
MARKERS = (
    PROJECT
    / "analysis"
    / "results"
    / "PLAT_axis"
    / "GSE183904"
    / "fine_states"
    / "GSE183904_Myeloid_fine_state_top_markers.tsv"
)
OUT = (
    PROJECT
    / "analysis"
    / "results"
    / "PLAT_axis"
    / "GSE183904"
    / "fine_state_communication_exploratory"
    / "cross_cohort_state_mapping"
)
OUT.mkdir(parents=True, exist_ok=True)

REFERENCE_MAP = {
    "M1_Unresolved": "C1_FCN1_AREG_inflammatory_monocyte",
    "M3_C1QC_APOE_target": "C3_C1QC_APOE_macrophage",
    "M4_NeutrophilLike": "C4_IL1B_CXCL8_inflammatory_myeloid",
    "M5_C1QC_boundary": "C5_HLAII_C1QC_macrophage",
}
CLUSTER_TO_REFERENCE = {
    1: "C1_FCN1_AREG_inflammatory_monocyte",
    3: "C3_C1QC_APOE_macrophage",
    4: "C4_IL1B_CXCL8_inflammatory_myeloid",
    5: "C5_HLAII_C1QC_macrophage",
}
AXIS_GENES = {"PLAT", "ANXA2", "ITGAM", "ITGB2", "LRP1", "ILK"}
STRESS_GENES = {"FOS", "JUN", "JUNB", "EGR1", "HSPA1A", "HSPA1B", "DUSP1"}
MIN_CELLS_PER_DONOR_STATE = 5


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def allowed_gene(gene: str) -> bool:
    gene = str(gene)
    return not (
        gene in AXIS_GENES
        or gene in STRESS_GENES
        or gene.startswith("MT-")
        or gene.startswith("RPL")
        or gene.startswith("RPS")
    )


def normalize_log_cp10k(x: sparse.spmatrix) -> sparse.csr_matrix:
    x = x.tocsr().astype(np.float32)
    lib = np.asarray(x.sum(axis=1)).ravel().astype(float)
    if np.any(lib <= 0):
        raise RuntimeError("Zero-library cells in GSE183904 mapping input")
    answer = x.multiply((1e4 / lib)[:, None]).tocsr()
    answer.data = np.log1p(answer.data)
    return answer


def donor_equal_centroids(
    x: sparse.spmatrix | np.ndarray,
    obs: pd.DataFrame,
    groups: list[str],
    group_col: str,
    donor_col: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    coverage = []
    group_values = obs[group_col].astype(str).to_numpy()
    donor_values = obs[donor_col].astype(str).to_numpy()
    for group in groups:
        group_mask = group_values == group
        donors = sorted(np.unique(donor_values[group_mask]))
        donor_profiles = []
        donor_cells = []
        for donor in donors:
            ix = np.flatnonzero(group_mask & (donor_values == donor))
            if len(ix) < MIN_CELLS_PER_DONOR_STATE:
                continue
            profile = np.asarray(x[ix].mean(axis=0)).ravel()
            donor_profiles.append(profile)
            donor_cells.append(len(ix))
        if not donor_profiles:
            continue
        rows.append(np.mean(np.vstack(donor_profiles), axis=0))
        coverage.append(
            {
                "state": group,
                "cells": int(group_mask.sum()),
                "patients_total": int(len(donors)),
                "patients_ge5": int(len(donor_profiles)),
                "minimum_cells_among_eligible_patients": int(min(donor_cells)),
                "maximum_cells_among_eligible_patients": int(max(donor_cells)),
            }
        )
    centroid = pd.DataFrame(np.vstack(rows), index=[x["state"] for x in coverage])
    return centroid, pd.DataFrame(coverage)


def zscore_genes(frame: pd.DataFrame) -> pd.DataFrame:
    mean = frame.mean(axis=0)
    sd = frame.std(axis=0, ddof=0)
    return (frame - mean) / sd.replace(0, np.nan)


def main() -> None:
    plan = json.loads(PLAN_FILE.read_text(encoding="utf-8"))
    rng = np.random.default_rng(int(plan["expression_mapping"]["random_seed"]))
    n_boot = int(plan["expression_mapping"]["bootstrap_replicates"])

    log("read GSE206785 communication-feature matrix")
    query = ad.read_h5ad(INPUT206)
    query_obs = query.obs.copy()
    query_mask = query_obs.Type.astype(str).eq("Myeloid").to_numpy()
    query_obs = query_obs.loc[query_mask].copy()
    query_x = query.X[query_mask]
    query_x = query_x.tocsr() if sparse.issparse(query_x) else sparse.csr_matrix(query_x)

    query_state_counts = query_obs.Annotation.astype(str).value_counts()
    query_patient_counts = (
        query_obs.groupby("Annotation", observed=True).Patient.nunique().astype(int)
    )
    query_groups = sorted(
        state
        for state in query_state_counts.index
        if query_state_counts[state] >= 50 and query_patient_counts.get(state, 0) >= 3
    )
    if "M02_Mono_CD14" not in query_groups:
        raise RuntimeError("Prespecified M02_Mono_CD14 query failed coverage gate")

    log("read GSE183904 full count bundle")
    stem = INPUT183 / "GSE183904_fine_full"
    ref_x_full = mmread(str(stem) + "_raw_counts_genes_by_cells.mtx").T.tocsr()
    ref_genes = pd.read_csv(str(stem) + "_genes.tsv", sep="\t").iloc[:, 0].astype(str)
    ref_obs = pd.read_csv(str(stem) + "_cells.tsv", sep="\t")
    if ref_x_full.shape != (len(ref_obs), len(ref_genes)):
        raise RuntimeError("GSE183904 full bundle shape mismatch")
    ref_mask = ref_obs.ccc_group.isin(REFERENCE_MAP).to_numpy()
    ref_obs = ref_obs.loc[ref_mask].copy()
    ref_obs["reference_state"] = ref_obs.ccc_group.map(REFERENCE_MAP)
    ref_x_full = ref_x_full[ref_mask]

    shared = sorted(
        gene
        for gene in set(query.var_names.astype(str)).intersection(ref_genes)
        if allowed_gene(gene)
    )
    if len(shared) < 500:
        raise RuntimeError(f"Too few axis-excluded shared genes: {len(shared)}")
    query_gene_index = pd.Index(query.var_names.astype(str)).get_indexer(shared)
    ref_gene_index = pd.Index(ref_genes).get_indexer(shared)
    query_x = query_x[:, query_gene_index]
    ref_x = normalize_log_cp10k(ref_x_full)[:, ref_gene_index]

    log(f"build donor-equal centroids over {len(shared)} shared genes")
    q_centroid, q_coverage = donor_equal_centroids(
        query_x, query_obs, query_groups, "Annotation", "Patient"
    )
    r_centroid, r_coverage = donor_equal_centroids(
        ref_x,
        ref_obs,
        list(REFERENCE_MAP.values()),
        "reference_state",
        "Patient_ID",
    )
    q_centroid.columns = shared
    r_centroid.columns = shared
    q_coverage.insert(0, "cohort", "GSE206785")
    r_coverage.insert(0, "cohort", "GSE183904")
    pd.concat([q_coverage, r_coverage], ignore_index=True).to_csv(
        OUT / "state_mapping_coverage.tsv", sep="\t", index=False
    )

    q_sd = q_centroid.std(axis=0, ddof=0)
    r_sd = r_centroid.std(axis=0, ddof=0)
    informative = q_sd.index[(q_sd > 0.05) & (r_sd > 0.05)].tolist()
    if len(informative) < 200:
        raise RuntimeError(f"Too few informative genes after variance gate: {len(informative)}")
    qz = zscore_genes(q_centroid[informative])
    rz = zscore_genes(r_centroid[informative])

    records = []
    rank_cache_q = {state: rankdata(qz.loc[state].to_numpy()) for state in qz.index}
    rank_cache_r = {state: rankdata(rz.loc[state].to_numpy()) for state in rz.index}
    bootstrap_indices = rng.integers(0, len(informative), size=(n_boot, len(informative)))
    top_counts = {q: {r: 0 for r in rz.index} for q in qz.index}
    bootstrap_store: dict[tuple[str, str], np.ndarray] = {}
    for query_state in qz.index:
        pair_boot = {}
        qrank = rank_cache_q[query_state]
        for ref_state in rz.index:
            rrank = rank_cache_r[ref_state]
            observed = float(spearmanr(qz.loc[query_state], rz.loc[ref_state]).statistic)
            boot = np.empty(n_boot, dtype=float)
            for i, ix in enumerate(bootstrap_indices):
                a = qrank[ix]
                b = rrank[ix]
                sa = a.std(ddof=0)
                sb = b.std(ddof=0)
                boot[i] = np.corrcoef(a, b)[0, 1] if sa > 0 and sb > 0 else np.nan
            pair_boot[ref_state] = boot
            bootstrap_store[(query_state, ref_state)] = boot
            records.append(
                {
                    "query_state": query_state,
                    "reference_state": ref_state,
                    "n_informative_genes": len(informative),
                    "spearman_rho": observed,
                    "gene_bootstrap_ci_low": float(np.nanquantile(boot, 0.025)),
                    "gene_bootstrap_ci_high": float(np.nanquantile(boot, 0.975)),
                }
            )
        matrix = np.vstack([pair_boot[r] for r in rz.index]).T
        winners = np.nanargmax(matrix, axis=1)
        for winner in winners:
            top_counts[query_state][rz.index[winner]] += 1
    mapping = pd.DataFrame(records)
    mapping["bootstrap_top_probability"] = [
        top_counts[row.query_state][row.reference_state] / n_boot
        for row in mapping.itertuples(index=False)
    ]
    mapping.to_csv(OUT / "cross_cohort_state_similarity.tsv", sep="\t", index=False)

    marker_table = pd.read_csv(MARKERS, sep="\t")
    marker_records = []
    qz_all = zscore_genes(q_centroid[shared])
    for cluster, reference_state in CLUSTER_TO_REFERENCE.items():
        marker_block = marker_table[(marker_table.cluster == cluster) & (marker_table.mean_log_difference > 0)].copy()
        marker_block = marker_block[marker_block.gene_symbol.map(allowed_gene)]
        marker_block = marker_block[marker_block.gene_symbol.isin(shared)].head(
            int(plan["marker_mapping"]["markers_per_state"])
        )
        if len(marker_block) < int(plan["marker_mapping"]["minimum_available_markers"]):
            continue
        genes = marker_block.gene_symbol.tolist()
        weights = marker_block.mean_log_difference.to_numpy(dtype=float)
        weights = weights / weights.sum()
        for query_state in q_centroid.index:
            score = float(np.nansum(qz_all.loc[query_state, genes].to_numpy(dtype=float) * weights))
            marker_records.append(
                {
                    "query_state": query_state,
                    "reference_state": reference_state,
                    "available_markers": len(genes),
                    "marker_score": score,
                    "markers": ";".join(genes),
                }
            )
    marker_scores = pd.DataFrame(marker_records)
    marker_scores["marker_rank_within_query"] = marker_scores.groupby(
        "query_state", observed=True
    ).marker_score.rank(method="min", ascending=False)
    marker_scores.to_csv(OUT / "cross_cohort_marker_projection.tsv", sep="\t", index=False)

    m02 = mapping[mapping.query_state == "M02_Mono_CD14"].sort_values("spearman_rho", ascending=False)
    m02_marker = marker_scores[marker_scores.query_state == "M02_Mono_CD14"].sort_values(
        "marker_score", ascending=False
    )
    summary = {
        "analysis_id": plan["analysis_id"],
        "shared_axis_excluded_genes": len(shared),
        "informative_genes": len(informative),
        "query_states": len(q_centroid),
        "reference_states": len(r_centroid),
        "M02_best_centroid_match": m02.iloc[0].reference_state,
        "M02_best_centroid_spearman_rho": float(m02.iloc[0].spearman_rho),
        "M02_best_centroid_bootstrap_top_probability": float(m02.iloc[0].bootstrap_top_probability),
        "M02_best_marker_projection": m02_marker.iloc[0].reference_state,
        "M02_best_marker_score": float(m02_marker.iloc[0].marker_score),
        "claim_boundary": plan["claim_boundary"],
        "feature_space_limitation": "GSE206785 input contains 2,029 communication-feature genes rather than an unbiased whole-transcriptome matrix.",
    }
    (OUT / "M02_mapping_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    pd.DataFrame({"gene": informative}).to_csv(
        OUT / "axis_excluded_informative_gene_universe.tsv", sep="\t", index=False
    )
    log(f"complete: {summary}")


if __name__ == "__main__":
    main()
