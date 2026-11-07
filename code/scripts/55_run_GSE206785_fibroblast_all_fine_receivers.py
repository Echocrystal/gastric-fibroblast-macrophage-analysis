#!/usr/bin/env python3
"""Resolve fibroblast outgoing communication across all GSE206785 fine states.

This is the unbiased receiver-resolution check required before focusing on the
myeloid compartment.  It evaluates every eligible fibroblast author state
against every eligible non-fibroblast author state in the same donor-balanced
tumour-cell object and then aggregates targets back to their broad lineage.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import liana as li
import numpy as np
import pandas as pd


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
OUT = PROJECT / "analysis/results/figure1_gse206785_unbiased"
FIGDATA = PROJECT / "figures/data/Figure1_GSE206785_unbiased"
AUDIT = PROJECT / "analysis/audit"
INPUT = OUT / "GSE206785_donor_balanced_all_lineage.h5ad"
SEED = 20260802
MIN_CELLS = 10
N_PERMS = 1000
EXPR_PROP = 0.10


def main() -> None:
    a = ad.read_h5ad(INPUT)
    counts = a.obs["Annotation"].value_counts()
    eligible = sorted(counts[counts.ge(MIN_CELLS)].index.astype(str))
    fibro_states = sorted(
        set(a.obs.loc[a.obs["Type"].eq("Fibroblast"), "Annotation"].astype(str)).intersection(eligible)
    )
    target_states = sorted(set(eligible).difference(fibro_states))
    pairs = pd.MultiIndex.from_product(
        [fibro_states, target_states], names=["source", "target"]
    ).to_frame(index=False)

    a.obs["fine_group"] = a.obs["Annotation"].astype("category")
    result = li.mt.rank_aggregate(
        a,
        groupby="fine_group",
        groupby_pairs=pairs,
        resource_name="consensus",
        expr_prop=EXPR_PROP,
        min_cells=MIN_CELLS,
        use_raw=False,
        n_perms=N_PERMS,
        seed=SEED,
        n_jobs=4,
        inplace=False,
        verbose=True,
    )
    result.to_csv(OUT / "GSE206785_fibroblast_to_all_fine_LIANA.tsv.gz", sep="\t", index=False)

    annotation_type = (
        a.obs[["Annotation", "cell_type"]]
        .drop_duplicates()
        .set_index("Annotation")["cell_type"]
        .to_dict()
    )
    z = result.copy()
    z["target_lineage"] = z["target"].map(annotation_type)
    z["robust"] = z["magnitude_rank"].le(0.05) & z["specificity_rank"].le(0.05)
    z["mag_strength"] = -np.log10(np.clip(z["magnitude_rank"], 1e-12, 1.0))
    z["spec_strength"] = -np.log10(np.clip(z["specificity_rank"], 1e-12, 1.0))
    z["joint_strength"] = np.minimum(z["mag_strength"], z["spec_strength"])
    z["robust_weight"] = np.where(z["robust"], z["joint_strength"], 0.0)

    pair_summary = (
        z.groupby(["source", "target", "target_lineage"], observed=True)
        .agg(
            n_tested=("ligand_complex", "size"),
            n_robust=("robust", "sum"),
            robust_weight=("robust_weight", "sum"),
            best_joint_strength=("joint_strength", "max"),
            top10_joint_strength=("joint_strength", lambda x: float(np.mean(np.sort(np.asarray(x))[-10:]))),
        )
        .reset_index()
    )
    pair_summary.to_csv(FIGDATA / "panel_h_fibroblast_all_fine_receiver_pairs.tsv", sep="\t", index=False)

    lineage_summary = (
        pair_summary.groupby("target_lineage", observed=True)
        .agg(
            tested_state_pairs=("target", "size"),
            state_pairs_with_robust=("n_robust", lambda x: int(np.sum(np.asarray(x) > 0))),
            robust_lr_total=("n_robust", "sum"),
            robust_weight_total=("robust_weight", "sum"),
            median_top10_joint=("top10_joint_strength", "median"),
            maximum_joint=("best_joint_strength", "max"),
        )
        .reset_index()
        .sort_values(["robust_weight_total", "median_top10_joint"], ascending=False)
        .reset_index(drop=True)
    )
    lineage_summary["receiver_rank"] = np.arange(1, len(lineage_summary) + 1)
    lineage_summary.to_csv(FIGDATA / "panel_h_fine_resolved_receiver_lineage_summary.tsv", sep="\t", index=False)

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": "GSE206785",
        "analysis": "donor-balanced fibroblast-state to all fine receiver-state LIANA consensus",
        "fibroblast_states": fibro_states,
        "target_states": len(target_states),
        "tested_state_pairs": len(pairs),
        "n_perms": N_PERMS,
        "expr_prop": EXPR_PROP,
        "robust_rule": "magnitude_rank <= 0.05 AND specificity_rank <= 0.05",
        "claim_boundary": "fine-state predicted communication; no binding or causal claim",
    }
    (AUDIT / "Figure1_GSE206785_fine_receiver_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(lineage_summary.to_string(index=False))
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
