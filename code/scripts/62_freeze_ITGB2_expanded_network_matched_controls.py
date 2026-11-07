#!/usr/bin/env python
"""Freeze a post-hoc expanded ITGB2 matched-control set before new KO runs.

The selection uses only pre-knockout expression, detection and common WT-network
features. No virtual-knockout outcome file is read. The first 20 selected genes
must reproduce the previously frozen v2 ITGB2 controls exactly.
"""

from __future__ import annotations

from datetime import datetime
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
AUDIT = PROJECT / "analysis" / "audit"
RESOURCE = PROJECT / "resources" / "gene_sets"
RESULT = (
    PROJECT
    / "analysis"
    / "results"
    / "virtual_knockout"
    / "GSE183904"
    / "C1QC_APOE_macrophage_receptor_screen_v2"
)

FEATURE_FILE = AUDIT / "GSE183904_primary_receptors_v2_scTenifoldKnk_control_matching_universe.tsv"
OLD_CONTROL_FILE = RESOURCE / "GSE183904_primary_receptors_v2_scTenifoldKnk_network_matched_controls.tsv"
OUTFILE = RESOURCE / "GSE183904_ITGB2_expanded_60_network_matched_controls_posthoc.tsv"
AUDIT_FILE = AUDIT / "GSE183904_ITGB2_expanded_60_network_matched_controls_posthoc.json"

TARGET = "ITGB2"
N_CONTROLS = 60
MATCH_VARIABLES = [
    "mean_raw_count",
    "detection_fraction",
    "mean_total_abs_strength",
    "mean_nonzero_degree",
]
AXIS_GENES = {
    "CD74",
    "CCRL2",
    "CXCR4",
    "ITGB2",
    "LRP1",
    "ITGAM",
    "CMKLR1",
    "ANXA2",
    "TREM2",
    "PLAT",
    "ILK",
}


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def transform(name: str, values: np.ndarray) -> np.ndarray:
    values = values.astype(float)
    if name == "detection_fraction":
        values = np.clip(values, 1e-4, 1 - 1e-4)
        return np.log(values / (1 - values))
    return np.log1p(values)


def main() -> None:
    features = pd.read_csv(FEATURE_FILE, sep="\t")
    target = features[features["gene_symbol"] == TARGET].copy()
    if len(target) != 1:
        raise RuntimeError("Expected one ITGB2 feature row")
    target_membership = int(target.iloc[0]["program_membership_count"])

    technical = features["technical"].astype(str).str.lower().eq("true")
    pool = features[
        ~features["gene_symbol"].isin(AXIS_GENES)
        & ~technical
        & (features["patient_count"] >= 3)
        & (features["detection_fraction"] >= 0.05)
        & ((features["program_membership_count"] - target_membership).abs() <= 1)
    ].copy()
    if len(pool) < N_CONTROLS:
        raise RuntimeError(f"Only {len(pool)} eligible controls for requested {N_CONTROLS}")

    squared_distance = np.zeros(len(pool), dtype=float)
    for name in MATCH_VARIABLES:
        combined = np.concatenate(
            [target[name].to_numpy(dtype=float), pool[name].to_numpy(dtype=float)]
        )
        transformed = transform(name, combined)
        center = float(np.median(transformed))
        scale = float(np.median(np.abs(transformed - center)))
        if not np.isfinite(scale) or scale <= 0:
            scale = float(np.std(transformed, ddof=1))
        if not np.isfinite(scale) or scale <= 0:
            raise RuntimeError(f"Degenerate match feature: {name}")
        target_z = (transformed[0] - center) / scale
        pool_z = (transformed[1:] - center) / scale
        pool[f"z_{name}"] = pool_z
        pool[f"target_z_{name}"] = target_z
        squared_distance += (pool_z - target_z) ** 2

    pool["feature_distance"] = np.sqrt(squared_distance)
    pool["membership_distance"] = (
        pool["program_membership_count"] - target_membership
    ).abs()
    pool["total_match_distance"] = (
        pool["feature_distance"] + 2 * pool["membership_distance"]
    )
    pool = pool.sort_values(
        ["total_match_distance", "gene_symbol"], kind="mergesort"
    ).reset_index(drop=True)
    pool.insert(0, "expanded_control_rank", np.arange(1, len(pool) + 1))

    old = pd.read_csv(OLD_CONTROL_FILE, sep="\t")
    old_itgb2 = old[old["target"] == TARGET].sort_values("control_rank")
    old_genes = old_itgb2["control_gene"].tolist()
    reproduced_genes = pool.head(len(old_genes))["gene_symbol"].tolist()
    if old_genes != reproduced_genes:
        raise RuntimeError(
            "Expanded matching did not reproduce the old first 20 controls exactly"
        )

    selected = pool.head(N_CONTROLS).copy()
    completed_dir = RESULT / "primary_receptor_network_matched_controls_seed2026072739"
    selected["seed2026072739_fit_already_exists"] = selected["gene_symbol"].map(
        lambda gene: (completed_dir / f"scTenifoldKnk_control__{gene}.rds").exists()
    )
    selected["analysis_status"] = "post_hoc_exploratory_frozen_before_expanded_KO_outcomes"
    keep = [
        "expanded_control_rank",
        "gene_symbol",
        "mean_raw_count",
        "detection_fraction",
        "patient_count",
        "mean_total_abs_strength",
        "mean_nonzero_degree",
        "strength_seed_CV",
        "program_membership_count",
        "feature_distance",
        "membership_distance",
        "total_match_distance",
        "seed2026072739_fit_already_exists",
        "analysis_status",
    ]
    selected[keep].to_csv(OUTFILE, sep="\t", index=False)

    exact_membership = int((selected["membership_distance"] == 0).sum())
    audit = {
        "generated_at": datetime.now().astimezone().isoformat(),
        "analysis_status": "post_hoc_exploratory_frozen_before_expanded_KO_outcomes",
        "target": TARGET,
        "controls_requested": N_CONTROLS,
        "eligible_pool_size": int(len(pool)),
        "selection_reads_KO_outcomes": False,
        "match_variables": MATCH_VARIABLES,
        "program_membership_rule": "absolute difference from ITGB2 <= 1",
        "exact_program_membership_controls_selected": exact_membership,
        "first_20_reproduce_previous_v2_manifest": True,
        "seed2026072739_fits_already_available": int(
            selected["seed2026072739_fit_already_exists"].sum()
        ),
        "seed2026072739_fits_still_required": int(
            (~selected["seed2026072739_fit_already_exists"]).sum()
        ),
        "empirical_upper_tail_resolution_with_60_controls": {
            "zero_controls_at_least_target": 1 / 61,
            "one_control_at_least_target": 2 / 61,
            "two_controls_at_least_target": 3 / 61,
        },
        "input_sha256": {
            "feature_file": file_sha256(FEATURE_FILE),
            "old_control_file": file_sha256(OLD_CONTROL_FILE),
        },
        "output_file": str(OUTFILE),
        "claim_boundary": (
            "Expanded controls improve empirical-null resolution only. The analysis is "
            "post hoc, does not guarantee significance and does not establish causality."
        ),
    }
    AUDIT_FILE.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    print(f"Eligible controls: {len(pool)}")
    print(f"Frozen controls: {len(selected)}")
    print(f"Existing seed-1 fits: {audit['seed2026072739_fits_already_available']}")
    print(f"Missing seed-1 fits: {audit['seed2026072739_fits_still_required']}")
    print(
        selected[
            [
                "expanded_control_rank",
                "gene_symbol",
                "total_match_distance",
                "program_membership_count",
                "seed2026072739_fit_already_exists",
            ]
        ]
        .head(25)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
