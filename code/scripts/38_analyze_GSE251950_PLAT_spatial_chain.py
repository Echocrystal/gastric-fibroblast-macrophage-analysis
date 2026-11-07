"""Independent GSE251950 projection of the frozen PLAT spatial chain.

Coverage gates and endpoints are frozen in
analysis/design/PLAT_SPATIAL_ENRICHMENT_LOCK_20260806.md.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
PUBLIC = Path(r"I:\ECM-new\ECM-public")
BASE_SCRIPT = PROJECT / "scripts" / "35_analyze_GSE270678_PLAT_spatial_chain.py"
LOCK_FILE = PROJECT / "analysis" / "design" / "PLAT_SPATIAL_ENRICHMENT_LOCK_20260806.md"
RAW = PUBLIC / "data" / "raw" / "spatial" / "GSE251950" / "samples"
SPOT_FILE = PUBLIC / "data" / "processed" / "spatial" / "GSE251950" / "GSE251950_spot_ecology_scores.tsv.gz"
MANIFEST_FILE = PUBLIC / "data" / "processed" / "spatial" / "GSE251950" / "GSE251950_spatial_sample_manifest.tsv"
OUT = PROJECT / "analysis" / "results" / "spatial" / "GSE251950_PLAT_frozen_chain_v1"
OUT.mkdir(parents=True, exist_ok=True)

SPEC = importlib.util.spec_from_file_location("plat_chain_base", BASE_SCRIPT)
M = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(M)

SEED = 2026080638
RNG = np.random.default_rng(SEED)
MIN_ANALYZABLE_PATIENTS = 8


def h5_totals(h5_path: Path, barcodes: list[str]) -> np.ndarray:
    with h5py.File(h5_path, "r") as fh:
        m = fh["matrix"]
        h5_barcodes = np.asarray(m["barcodes"]).astype("U")
        matrix = csc_matrix(
            (m["data"][...], m["indices"][...], m["indptr"][...]),
            shape=tuple(m["shape"][...]),
        )
    lookup = {b: i for i, b in enumerate(h5_barcodes)}
    idx = np.array([lookup.get(b, -1) for b in barcodes], dtype=int)
    if np.any(idx < 0):
        raise RuntimeError(f"QC barcode missing from {h5_path.name}")
    return np.asarray(matrix[:, idx].sum(axis=0)).ravel().astype(float)


def coverage_count(norm: pd.DataFrame, genes: list[str]) -> tuple[int, int]:
    available = [g for g in genes if g in norm.index]
    if not available:
        return 0, 0
    x = norm.loc[available].to_numpy(float)
    _, informative = M.zscore_genes(x)
    return len(available), int(informative.sum())


def bootstrap_ci(values: np.ndarray, draws: int = 9999) -> tuple[float, float]:
    r = np.asarray(values, float)
    r = r[np.isfinite(r)]
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    sampled = RNG.choice(z, size=(draws, len(z)), replace=True)
    return tuple(np.quantile(np.tanh(sampled.mean(axis=1)), [0.025, 0.975]).tolist())


def main():
    sender_genes, receiver_genes, deg_genes, adh_genes = M.load_signatures()
    needed = sorted(set(sender_genes) | set(receiver_genes) | set(deg_genes) | set(adh_genes) | set(M.MYELOID_GENES) | {"PLAT"})
    manifest = pd.read_csv(MANIFEST_FILE, sep="\t")
    manifest = manifest.loc[manifest.primary_analysis.astype(bool)].copy()
    if len(manifest) != 9 or manifest.patient_id.nunique() != 9:
        raise RuntimeError(f"Expected nine primary patients, observed {len(manifest)}/{manifest.patient_id.nunique()}")
    spots = pd.read_csv(SPOT_FILE, sep="\t")
    spots = spots.loc[spots.patient_id.isin(manifest.patient_id)].copy()

    slide_rows, coverage_rows = [], []
    cached: dict[str, tuple[pd.DataFrame, pd.DataFrame]] = {}
    for i, row in enumerate(manifest.itertuples(index=False), 1):
        print(f"[coverage {i}/9] {row.sample_title}", flush=True)
        base = spots.loc[spots.GSM.eq(row.GSM)].sort_values("barcode").reset_index(drop=True)
        base = base.rename(columns={"total_UMI": "total_counts"})
        folder = RAW / row.folder_id
        h5_path = folder / f"{row.folder_id}_filtered_feature_bc_matrix.h5"
        totals = h5_totals(h5_path, base.barcode.tolist())
        base["total_counts"] = totals
        norm = M.read_selected_expression(h5_path, base.barcode.tolist(), needed, totals)
        s_avail, s_info = coverage_count(norm, sender_genes)
        r_avail, r_info = coverage_count(norm, receiver_genes)
        d_avail, d_info = coverage_count(norm, deg_genes)
        a_avail, a_info = coverage_count(norm, adh_genes)
        plat_present = "PLAT" in norm.index
        gate = bool(plat_present and s_info >= 10 and r_info >= 10 and d_info >= 20)
        coverage_rows.append({
            "GSM": row.GSM, "sample_title": row.sample_title, "patient_id": row.patient_id,
            "n_spots": len(base), "PLAT_present": plat_present,
            "PLAT_detected_fraction": float(np.mean(norm.loc["PLAT"].to_numpy(float) > 0)) if plat_present else np.nan,
            "sender_available": s_avail, "sender_informative": s_info,
            "receiver_available": r_avail, "receiver_informative": r_info,
            "ECMdeg_available": d_avail, "ECMdeg_informative": d_info,
            "adhesion_available": a_avail, "adhesion_informative": a_info,
            "coverage_gate": gate,
        })
        cached[row.GSM] = (base, norm)

    coverage = pd.DataFrame(coverage_rows)
    analyzable = coverage.loc[coverage.coverage_gate, "GSM"].tolist()
    coverage.to_csv(OUT / "GSE251950_PLAT_frozen_chain_gene_coverage.tsv", sep="\t", index=False)
    if len(analyzable) < MIN_ANALYZABLE_PATIENTS:
        audit = {
            "dataset": "GSE251950", "primary_patients": 9, "analyzable_patients": len(analyzable),
            "coverage_gate_passed": False, "minimum_required": MIN_ANALYZABLE_PATIENTS,
            "warning": "Formal chain not run because the frozen coverage gate failed.",
        }
        (OUT / "GSE251950_PLAT_frozen_chain_audit.json").write_text(
            json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        print("\nCOVERAGE GATE FAILED")
        print(coverage.to_string(index=False))
        return

    for i, row in enumerate(manifest.loc[manifest.GSM.isin(analyzable)].itertuples(index=False), 1):
        print(f"[analysis {i}/{len(analyzable)}] {row.sample_title}", flush=True)
        base, norm = cached[row.GSM]
        effects, _ = M.endpoint_rows_for_slide(base, norm, sender_genes, receiver_genes, deg_genes, adh_genes)
        for effect in effects:
            slide_rows.append({
                "GSM": row.GSM, "sample_title": row.sample_title, "patient_id": row.patient_id, **effect,
            })

    slide = pd.DataFrame(slide_rows)
    donor = slide.rename(columns={"patient_id": "donor_id", "rho": "rho"})[
        ["donor_id", "endpoint", "ring", "rho", "GSM", "n_spots"]
    ].copy()
    donor["n_sections"] = 1
    donor["median_section_spots"] = donor.n_spots

    meta_rows = []
    for (endpoint, ring), d in donor.groupby(["endpoint", "ring"], sort=False):
        x = d.rho.to_numpy(float)
        lo, hi = bootstrap_ci(x)
        meta_rows.append({
            "endpoint": endpoint, "ring": ring, "n_patients": len(x),
            "pooled_rho": M.fisher_pool(x), "ci_low": lo, "ci_high": hi,
            "positive_patients": int(np.sum(x > 0)), "negative_patients": int(np.sum(x < 0)),
            "exact_signflip_p": M.exact_signflip_p(x), "primary_family": endpoint in M.PRIMARY,
        })
    meta = pd.DataFrame(meta_rows)
    wide = donor.pivot(index="donor_id", columns="endpoint", values="rho").reset_index()
    adjusted = M.maxT_primary(wide)
    meta["primary_maxT_FWER_p"] = np.nan
    for endpoint in M.PRIMARY:
        meta.loc[meta.endpoint.eq(endpoint), "primary_maxT_FWER_p"] = adjusted[endpoint]

    slide.to_csv(OUT / "GSE251950_PLAT_frozen_chain_slide_effects.tsv", sep="\t", index=False)
    donor.to_csv(OUT / "GSE251950_PLAT_frozen_chain_patient_effects.tsv", sep="\t", index=False)
    meta.to_csv(OUT / "GSE251950_PLAT_frozen_chain_patient_meta.tsv", sep="\t", index=False)

    audit = {
        "seed": SEED, "dataset": "GSE251950 primary gastric tumours", "primary_patients": 9,
        "analyzable_patients": len(analyzable), "coverage_gate_passed": True,
        "excluded_samples": coverage.loc[~coverage.coverage_gate, "sample_title"].tolist(),
        "primary_endpoints": M.PRIMARY,
        "primary_maxT_statistic": "joint patient-sign-flip maximum absolute studentized Fisher-z mean",
        "maxT_complete_patients": int(adjusted["n_complete"]),
        "lock_file": str(LOCK_FILE),
        "warning": "Mixed-cell Visium associations; no direct contact, recruitment, mediation or causality claim.",
    }
    (OUT / "GSE251950_PLAT_frozen_chain_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("\nCOVERAGE")
    print(coverage.to_string(index=False))
    print("\nPRIMARY RESULTS")
    print(meta.loc[meta.endpoint.isin(M.PRIMARY)].sort_values("endpoint").to_string(index=False))
    print("\nRING PROFILES")
    print(meta.loc[meta.endpoint.str.startswith(("L2_", "L3_"))].sort_values("endpoint").to_string(index=False))
    print("\nINCREMENTAL")
    print(meta.loc[meta.endpoint.str.startswith("PLAT_incremental")].sort_values("endpoint").to_string(index=False))


if __name__ == "__main__":
    main()
