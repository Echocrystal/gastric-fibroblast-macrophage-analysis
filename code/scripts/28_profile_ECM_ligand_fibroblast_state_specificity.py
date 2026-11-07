#!/usr/bin/env python3
"""Profile scanned ECM ligands across GSE206785 tumor fibroblast states.

This is a state-localization gate, not a differential-expression discovery
test.  The 72-ligand universe was frozen by script 25 before this analysis.
Both cell-weighted and patient-equal summaries are reported.  The paired
contrast is restricted to patients containing both F13 CTHRC1 CAFs and other
tumor fibroblasts, avoiding a cell-level pseudo-replication P value.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
SOURCE = Path(
    r"I:\ECM-new\ECM_DEG_WGCNA_CoxLASSO_public_2026\ecm_immune_ecotype_public_2026"
)
EXPR = SOURCE / "data/raw/scrna/GSE206785/GSE206785_scgex.txt.gz"
META = SOURCE / "data/raw/scrna/GSE206785/GSE206785_metadata.txt.gz"
RANKING = (
    PROJECT
    / "analysis/results/axis_scan/GSE206785/GSE206785_ECM_ligand_axis_ranking.tsv"
)
OUT = PROJECT / "analysis/results/axis_scan/GSE206785/state_specificity"
AUDIT = PROJECT / "analysis/audit"
OUT.mkdir(parents=True, exist_ok=True)
AUDIT.mkdir(parents=True, exist_ok=True)

F13 = "F13_Activated_CTHRC1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


metadata = pd.read_csv(META, sep="\t", dtype=str)
ranking = pd.read_csv(RANKING, sep="\t")
ligands = ranking["ligand"].astype(str).str.upper().drop_duplicates().tolist()
if len(ligands) != 72:
    raise RuntimeError(f"Expected frozen 72-ligand universe, observed {len(ligands)}")

# pandas still scans the compressed CSV, but materializes only 73 columns.
expression = pd.read_csv(EXPR, usecols=["Cell", *ligands])
if len(expression) != len(metadata):
    raise RuntimeError(f"Expression/metadata row mismatch: {len(expression)} vs {len(metadata)}")

suffix = "-" + metadata["Sample"].astype(str)
if not np.all([cell.endswith(sample) for cell, sample in zip(expression["Cell"], suffix)]):
    raise RuntimeError("Expression rows are not aligned to metadata sample suffixes")

is_fib = metadata["Tissue"].eq("Tumor") & metadata["Type"].eq("Fibroblast")
meta = metadata.loc[is_fib, ["Patient", "Annotation"]].reset_index(drop=True)
values = expression.loc[is_fib, ligands].reset_index(drop=True)
raw = np.rint(np.expm1(values.to_numpy(dtype=np.float64))).astype(np.int32)
if np.any(raw < 0):
    raise RuntimeError("Negative reconstructed counts")

raw_df = pd.DataFrame(raw, columns=ligands)
long = pd.concat([meta, raw_df], axis=1).melt(
    id_vars=["Patient", "Annotation"], var_name="ligand", value_name="raw_count"
)
long["detected"] = long["raw_count"].gt(0)

state = (
    long.groupby(["Annotation", "ligand"], observed=True)
    .agg(
        cells=("raw_count", "size"),
        patients=("Patient", "nunique"),
        cell_weighted_mean_raw=("raw_count", "mean"),
        cell_weighted_detection=("detected", "mean"),
        count_sum=("raw_count", "sum"),
    )
    .reset_index()
)

patient_state = (
    long.groupby(["Patient", "Annotation", "ligand"], observed=True)
    .agg(
        patient_state_mean_raw=("raw_count", "mean"),
        patient_state_detection=("detected", "mean"),
        cells=("raw_count", "size"),
    )
    .reset_index()
)
patient_equal = (
    patient_state.groupby(["Annotation", "ligand"], observed=True)
    .agg(
        patient_equal_mean_raw=("patient_state_mean_raw", "mean"),
        patient_equal_detection=("patient_state_detection", "mean"),
    )
    .reset_index()
)
state = state.merge(patient_equal, on=["Annotation", "ligand"], how="left")
state["fraction_of_gene_counts_across_tumor_fibroblasts"] = state["count_sum"] / state.groupby(
    "ligand", observed=True
)["count_sum"].transform("sum").replace(0, np.nan)

contrast_rows: list[dict[str, object]] = []
for ligand in ligands:
    subset = long.loc[long["ligand"].eq(ligand)].copy()
    subset["compartment"] = np.where(subset["Annotation"].eq(F13), "F13", "other")
    patient_compartment = (
        subset.groupby(["Patient", "compartment"], observed=True)
        .agg(mean_raw=("raw_count", "mean"), detection=("detected", "mean"), cells=("raw_count", "size"))
        .reset_index()
    )
    wide_mean = patient_compartment.pivot(index="Patient", columns="compartment", values="mean_raw")
    wide_det = patient_compartment.pivot(index="Patient", columns="compartment", values="detection")
    paired = wide_mean.dropna(subset=["F13", "other"])
    paired_det = wide_det.loc[paired.index]
    mean_diff = paired["F13"] - paired["other"]
    det_diff = paired_det["F13"] - paired_det["other"]
    try:
        p_mean = float(wilcoxon(mean_diff, alternative="two-sided").pvalue) if np.any(mean_diff) else 1.0
    except ValueError:
        p_mean = np.nan
    try:
        p_det = float(wilcoxon(det_diff, alternative="two-sided").pvalue) if np.any(det_diff) else 1.0
    except ValueError:
        p_det = np.nan

    f13_cells = subset.loc[subset["compartment"].eq("F13"), "raw_count"]
    other_cells = subset.loc[subset["compartment"].eq("other"), "raw_count"]
    contrast_rows.append(
        {
            "ligand": ligand,
            "paired_patients": int(len(paired)),
            "F13_cell_weighted_mean_raw": float(f13_cells.mean()),
            "other_fibroblast_cell_weighted_mean_raw": float(other_cells.mean()),
            "F13_vs_other_log2_mean_ratio_pseudocount_0.1": float(
                np.log2((f13_cells.mean() + 0.1) / (other_cells.mean() + 0.1))
            ),
            "paired_patient_mean_raw_difference": float(mean_diff.mean()),
            "paired_patient_median_raw_difference": float(mean_diff.median()),
            "paired_patient_mean_expression_wilcoxon_p": p_mean,
            "paired_patient_mean_detection_difference": float(det_diff.mean()),
            "paired_patient_detection_wilcoxon_p": p_det,
            "paired_patients_F13_mean_higher": int((mean_diff > 0).sum()),
            "paired_patients_F13_mean_lower": int((mean_diff < 0).sum()),
        }
    )

contrast = pd.DataFrame(contrast_rows)
contrast["paired_patient_mean_expression_fdr_BH"] = np.minimum(
    1.0,
    contrast["paired_patient_mean_expression_wilcoxon_p"].rank(method="min")
    * len(contrast)
    / contrast["paired_patient_mean_expression_wilcoxon_p"].notna().sum(),
)
# Recalculate BH monotonically in sorted P-value order.
for col, out_col in [
    ("paired_patient_mean_expression_wilcoxon_p", "paired_patient_mean_expression_fdr_BH"),
    ("paired_patient_detection_wilcoxon_p", "paired_patient_detection_fdr_BH"),
]:
    p = contrast[col].to_numpy(dtype=float)
    valid = np.isfinite(p)
    order = np.argsort(p[valid])
    adjusted = np.full(len(p), np.nan)
    pv = p[valid][order]
    qv = np.minimum.accumulate((pv * valid.sum() / np.arange(1, valid.sum() + 1))[::-1])[::-1]
    valid_indices = np.flatnonzero(valid)[order]
    adjusted[valid_indices] = np.minimum(qv, 1.0)
    contrast[out_col] = adjusted

contrast = contrast.merge(
    ranking[["ligand", "scan_rank", "evidence_score", "minimum_rho", "minimum_normalized_AP"]],
    on="ligand",
    how="left",
)
contrast = contrast.sort_values(
    ["F13_vs_other_log2_mean_ratio_pseudocount_0.1", "scan_rank"], ascending=[False, True]
)

state.to_csv(OUT / "GSE206785_ECM_ligand_fibroblast_state_localization.tsv", sep="\t", index=False)
patient_state.to_csv(
    OUT / "GSE206785_ECM_ligand_patient_state_pseudobulk.tsv.gz",
    sep="\t",
    index=False,
    compression={"method": "gzip", "mtime": 0},
)
contrast.to_csv(OUT / "GSE206785_ECM_ligand_F13_state_specificity.tsv", sep="\t", index=False)

audit = {
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "dataset": "GSE206785",
    "frozen_ligand_count": len(ligands),
    "tumor_fibroblast_cells": int(is_fib.sum()),
    "tumor_fibroblast_patients": int(meta["Patient"].nunique()),
    "F13_cells": int(meta["Annotation"].eq(F13).sum()),
    "F13_patients": int(meta.loc[meta["Annotation"].eq(F13), "Patient"].nunique()),
    "source_value_semantics": "author-processed log1p integer counts",
    "reconstruction": "round(expm1(value))",
    "input_sha256": {
        "expression": sha256_file(EXPR),
        "metadata": sha256_file(META),
        "frozen_ranking": sha256_file(RANKING),
    },
    "interpretation": "State localization and patient-paired contrast; not a causal ligand effect.",
}
(AUDIT / "GSE206785_ECM_ligand_fibroblast_state_specificity.json").write_text(
    json.dumps(audit, indent=2), encoding="utf-8"
)

focus = ["PLAT", "HSPG2", "TIMP1", "DCN", "THBS1", "MFGE8", "LAMC3", "VCAN"]
print(contrast.loc[contrast["ligand"].isin(focus)].to_string(index=False))
print("\nTop 20 F13-localized ligands from frozen universe:")
print(contrast.head(20).to_string(index=False))
