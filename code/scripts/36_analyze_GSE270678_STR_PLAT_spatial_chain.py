"""ST_R repeated/extended-section sensitivity for the frozen PLAT spatial chain."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.stats import spearmanr


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
PUBLIC = Path(r"I:\ECM-new\ECM-public")
BASE_SCRIPT = PROJECT / "scripts" / "35_analyze_GSE270678_PLAT_spatial_chain.py"
SPEC = importlib.util.spec_from_file_location("plat_chain_main", BASE_SCRIPT)
M = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(M)

RAW = PUBLIC / "data" / "raw" / "spatial" / "GSE270678" / "ST_R_public_files"
SPOT_FILE = (
    PUBLIC / "data" / "processed" / "spatial" / "GSE270678" / "v4"
    / "GSE270678_ST_R_project_compartment_spot_scores_v4.tsv.gz"
)
MANIFEST_FILE = PUBLIC / "manifests" / "v4" / "GSE270678_spatial_sample_manifest.tsv"
MAIN_DONOR_FILE = (
    PROJECT / "analysis" / "results" / "spatial" / "GSE270678_PLAT_ecosystem_v1"
    / "GSE270678_PLAT_spatial_chain_donor_effects.tsv"
)
OUT = PROJECT / "analysis" / "results" / "spatial" / "GSE270678_PLAT_ecosystem_v1" / "ST_R_sensitivity"
OUT.mkdir(parents=True, exist_ok=True)


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
        raise RuntimeError(f"ST_R QC barcode missing from {h5_path.name}")
    return np.asarray(matrix[:, idx].sum(axis=0)).ravel().astype(float)


def main():
    sender_genes, receiver_genes, deg_genes, adh_genes = M.load_signatures()
    needed = sorted(
        set(sender_genes) | set(receiver_genes) | set(deg_genes) | set(adh_genes)
        | set(M.MYELOID_GENES) | {"PLAT"}
    )
    main_donor = pd.read_csv(MAIN_DONOR_FILE, sep="\t")
    overlapping = sorted(main_donor.donor_id.unique().tolist())

    spots = pd.read_csv(SPOT_FILE, sep="\t")
    spots = spots.loc[spots.in_tissue.astype(bool) & spots.donor_id.isin(overlapping)].copy()
    manifest_all = pd.read_csv(MANIFEST_FILE, sep="\t")
    str_all = manifest_all.loc[
        (manifest_all.platform_role == "ST_R")
        & manifest_all.is_replicate_or_extension.astype(bool)
    ].drop_duplicates("sample_id")
    excluded_extension = str_all.loc[
        ~str_all.donor_id.isin(overlapping), ["sample_id", "donor_id"]
    ].copy()
    manifest = str_all.loc[str_all.donor_id.isin(overlapping)].copy()
    if len(manifest) != 37 or manifest.donor_id.nunique() != 16:
        raise RuntimeError(
            f"Expected 37 overlapping-donor ST_R sections / 16 donors, got {len(manifest)} / {manifest.donor_id.nunique()}"
        )

    slide_rows = []
    coverage_rows = []
    for i, row in enumerate(manifest.itertuples(index=False), 1):
        print(f"[{i:02d}/{len(manifest)}] {row.sample_id}", flush=True)
        base = spots.loc[spots.sample_id.eq(row.sample_id)].sort_values("barcode").reset_index(drop=True)
        h5_path = RAW / row.filtered_h5_file
        totals = h5_totals(h5_path, base.barcode.tolist())
        base["total_counts"] = totals
        norm = M.read_selected_expression(h5_path, base.barcode.tolist(), needed, totals)
        effects, coverage = M.endpoint_rows_for_slide(
            base, norm, sender_genes, receiver_genes, deg_genes, adh_genes
        )
        for effect in effects:
            slide_rows.append({
                "sample_id": row.sample_id,
                "main_slice_id": row.main_slice_id,
                "donor_id": row.donor_id,
                "stage": row.stage_geo,
                "lauren": row.lauren_geo,
                "sampling_location": row.sampling_location_geo,
                **effect,
            })
        coverage_rows.append({
            "sample_id": row.sample_id, "main_slice_id": row.main_slice_id,
            "donor_id": row.donor_id, "QC_spots": len(base), **coverage,
        })

    slide = pd.DataFrame(slide_rows)
    coverage = pd.DataFrame(coverage_rows)
    donor = (
        slide.groupby(["donor_id", "endpoint", "ring"], as_index=False)
        .agg(
            rho=("rho", lambda x: M.fisher_pool(x.to_numpy(float))),
            n_sections=("sample_id", "nunique"),
            median_section_spots=("n_spots", "median"),
        )
    )
    meta_rows = []
    for (endpoint, ring), d in donor.groupby(["endpoint", "ring"], sort=False):
        values = d.rho.to_numpy(float)
        lo, hi = M.bootstrap_ci(values)
        meta_rows.append({
            "endpoint": endpoint, "ring": ring,
            "n_donors": int(np.isfinite(values).sum()),
            "pooled_rho": M.fisher_pool(values), "ci_low": lo, "ci_high": hi,
            "positive_donors": int(np.sum(values > 0)),
            "negative_donors": int(np.sum(values < 0)),
            "exact_signflip_p": M.exact_signflip_p(values),
            "primary_family": endpoint in M.PRIMARY,
        })
    meta = pd.DataFrame(meta_rows)
    wide = donor.pivot(index="donor_id", columns="endpoint", values="rho").reset_index()
    adjusted = M.maxT_primary(wide)
    meta["primary_maxT_FWER_p"] = np.nan
    for endpoint in M.PRIMARY:
        meta.loc[meta.endpoint.eq(endpoint), "primary_maxT_FWER_p"] = adjusted[endpoint]

    comparison_rows = []
    merged = main_donor.merge(
        donor, on=["donor_id", "endpoint", "ring"], suffixes=("_main", "_STR")
    )
    for (endpoint, ring), d in merged.groupby(["endpoint", "ring"], sort=False):
        x = d.rho_main.to_numpy(float)
        y = d.rho_STR.to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(y)
        corr = spearmanr(x[ok], y[ok]) if ok.sum() >= 3 else (np.nan, np.nan)
        comparison_rows.append({
            "endpoint": endpoint, "ring": ring, "n_overlap_donors": int(ok.sum()),
            "same_sign_donors": int(np.sum(np.sign(x[ok]) == np.sign(y[ok]))),
            "main_STR_spearman_rho": float(corr.statistic) if ok.sum() >= 3 else np.nan,
            "main_STR_spearman_p_descriptive": float(corr.pvalue) if ok.sum() >= 3 else np.nan,
            "main_pooled_rho_overlap": M.fisher_pool(x[ok]),
            "STR_pooled_rho_overlap": M.fisher_pool(y[ok]),
        })
    comparison = pd.DataFrame(comparison_rows)

    slide.to_csv(OUT / "GSE270678_STR_PLAT_spatial_chain_slide_effects.tsv", sep="\t", index=False)
    donor.to_csv(OUT / "GSE270678_STR_PLAT_spatial_chain_donor_effects.tsv", sep="\t", index=False)
    meta.to_csv(OUT / "GSE270678_STR_PLAT_spatial_chain_patient_meta.tsv", sep="\t", index=False)
    coverage.to_csv(OUT / "GSE270678_STR_PLAT_spatial_chain_gene_coverage.tsv", sep="\t", index=False)
    comparison.to_csv(OUT / "GSE270678_main_vs_STR_donor_concordance.tsv", sep="\t", index=False)

    audit = {
        "dataset": "GSE270678 ST_R repeated/extended sections",
        "sections": int(slide.sample_id.nunique()),
        "overlapping_donors": int(slide.donor_id.nunique()),
        "excluded_extension_only_donors": sorted(excluded_extension.donor_id.unique().tolist()),
        "excluded_extension_only_sections": sorted(excluded_extension.sample_id.unique().tolist()),
        "primary_endpoints": M.PRIMARY,
        "primary_maxT_statistic": "joint donor-sign-flip maximum absolute studentized Fisher-z mean",
        "warning": "ST_R shares donors with main ST and is a repeated-section sensitivity analysis, not an independent cohort.",
    }
    (OUT / "GSE270678_STR_PLAT_spatial_chain_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print("\nST_R PRIMARY RESULTS")
    print(meta.loc[meta.endpoint.isin(M.PRIMARY)].sort_values("endpoint").to_string(index=False))
    print("\nMAIN vs ST_R PRIMARY CONCORDANCE")
    print(comparison.loc[comparison.endpoint.isin(M.PRIMARY)].sort_values("endpoint").to_string(index=False))
    print("\nST_R RING PROFILES")
    print(meta.loc[meta.endpoint.str.startswith(("L2_", "L3_"))].sort_values("endpoint").to_string(index=False))
    print("\nST_R INCREMENTAL")
    print(meta.loc[meta.endpoint.str.startswith("PLAT_incremental")].sort_values("endpoint").to_string(index=False))


if __name__ == "__main__":
    main()
