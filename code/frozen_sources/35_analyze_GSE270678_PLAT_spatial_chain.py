"""Frozen external spatial-chain analysis for the PLAT gastric-cancer ecosystem.

Primary cohort: GSE270678 main ST, 33 sections / 19 donors.
The script performs no plotting and writes traceable section-, donor- and
meta-level TSV outputs. See the matching 2026-08-01 analysis lock.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.stats import rankdata


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
PUBLIC = Path(r"I:\ECM-new\ECM-public")
RAW = PUBLIC / "data" / "raw" / "spatial" / "GSE270678" / "main_public_files"
SPOT_FILE = (
    PUBLIC
    / "data"
    / "processed"
    / "spatial"
    / "GSE270678"
    / "v4"
    / "GSE270678_main_ST_project_compartment_spot_scores_v4.tsv.gz"
)
MANIFEST_FILE = PUBLIC / "manifests" / "v4" / "GSE270678_spatial_sample_manifest.tsv"
SIGNATURE_FILE = PROJECT / "resources" / "gene_sets" / "GSE183904_PLAT_ecosystem_bulk_signatures.tsv"
PROGRAM_FILE = PROJECT / "resources" / "gene_sets" / "myeloid_receiver_programs_frozen.tsv"
LOCK_FILE = PROJECT / "analysis" / "design" / "GSE270678_PLAT_SPATIAL_CHAIN_LOCK_20260801.md"
OUT = PROJECT / "analysis" / "results" / "spatial" / "GSE270678_PLAT_ecosystem_v1"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 2026080135
RNG = np.random.default_rng(SEED)

SENDER_ID = "sender_CXCL14_POSTN_state_axis_excluded"
RECEIVER_ID = "receiver_C1QC_APOE_specific_axis_excluded"
MYELOID_GENES = ["LST1", "FCER1G", "TYROBP", "AIF1", "CTSS", "CD68", "CSF1R", "LILRB1"]
AXIS_EXCLUDE = {"PLAT", "ANXA2", "ITGAM", "ITGB2", "LRP1"}
DEG_PROGRAM = "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX"
ADH_PROGRAM = "GOBP_CELL_MATRIX_ADHESION"
PRIMARY = ["L1_PLAT_sender_same", "L2_sender_receiver_r1", "L3_receiver_ECMdeg_r0"]

HEX_STEPS = ((-1, -1), (-1, 1), (0, -2), (0, 2), (1, -1), (1, 1))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def exact_ring_offsets(max_ring: int = 3) -> dict[int, list[tuple[int, int]]]:
    distance = {(0, 0): 0}
    frontier = {(0, 0)}
    for ring in range(1, max_ring + 1):
        new_frontier = set()
        for r, c in frontier:
            for dr, dc in HEX_STEPS:
                node = (r + dr, c + dc)
                if node not in distance:
                    distance[node] = ring
                    new_frontier.add(node)
        frontier = new_frontier
    return {r: sorted([x for x, d in distance.items() if d == r]) for r in range(max_ring + 1)}


RING_OFFSETS = exact_ring_offsets(3)


def spatial_shell(values: np.ndarray, rows: np.ndarray, cols: np.ndarray, ring: int) -> np.ndarray:
    if ring == 0:
        return np.asarray(values, dtype=float).copy()
    keys = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(rows, cols))}
    out = np.full(len(values), np.nan, dtype=float)
    offsets = RING_OFFSETS[ring]
    for i, (r, c) in enumerate(zip(rows, cols)):
        idx = [keys.get((int(r + dr), int(c + dc))) for dr, dc in offsets]
        idx = [j for j in idx if j is not None and np.isfinite(values[j])]
        if idx:
            out[i] = float(np.mean(values[idx]))
    return out


def zscore_genes(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    means = np.nanmean(matrix, axis=1, keepdims=True)
    sds = np.nanstd(matrix, axis=1, ddof=1, keepdims=True)
    informative = np.isfinite(sds[:, 0]) & (sds[:, 0] > 0)
    z = np.zeros_like(matrix, dtype=float)
    z[informative] = (matrix[informative] - means[informative]) / sds[informative]
    return z, informative


def module_score(norm: pd.DataFrame, genes: list[str], minimum: int) -> tuple[np.ndarray, int, int]:
    available = [g for g in genes if g in norm.index]
    if len(available) < minimum:
        raise RuntimeError(f"Only {len(available)} genes available; required {minimum}")
    x = norm.loc[available].to_numpy(float)
    z, informative = zscore_genes(x)
    if int(informative.sum()) < minimum:
        raise RuntimeError(f"Only {int(informative.sum())} informative genes; required {minimum}")
    return np.nanmean(z[informative], axis=0), len(available), int(informative.sum())


def residual_rank_correlation(x: np.ndarray, y: np.ndarray, covariates: list[np.ndarray]) -> tuple[int, float]:
    arrays = [np.asarray(x, float), np.asarray(y, float)] + [np.asarray(c, float) for c in covariates]
    ok = np.logical_and.reduce([np.isfinite(a) for a in arrays])
    n = int(ok.sum())
    if n < 100:
        return n, np.nan
    xr = rankdata(arrays[0][ok], method="average")
    yr = rankdata(arrays[1][ok], method="average")
    ranked_cov = []
    for c in arrays[2:]:
        cr = rankdata(c[ok], method="average")
        if np.nanstd(cr) > 0:
            ranked_cov.append(cr)
    design = np.ones((n, 1 + len(ranked_cov)), dtype=float)
    if ranked_cov:
        design[:, 1:] = np.column_stack(ranked_cov)
    rx = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
    ry = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
    if np.std(rx) == 0 or np.std(ry) == 0:
        return n, np.nan
    return n, float(np.corrcoef(rx, ry)[0, 1])


def fisher_pool(values: np.ndarray) -> float:
    r = np.asarray(values, dtype=float)
    r = r[np.isfinite(r)]
    if not len(r):
        return np.nan
    return float(np.tanh(np.mean(np.arctanh(np.clip(r, -0.999999, 0.999999)))))


def bootstrap_ci(values: np.ndarray, draws: int = 9999) -> tuple[float, float]:
    r = np.asarray(values, dtype=float)
    r = r[np.isfinite(r)]
    if len(r) < 2:
        return np.nan, np.nan
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    sampled = RNG.choice(z, size=(draws, len(z)), replace=True)
    pooled = np.tanh(sampled.mean(axis=1))
    return tuple(np.quantile(pooled, [0.025, 0.975]).tolist())


def exact_signflip_p(values: np.ndarray) -> float:
    r = np.asarray(values, dtype=float)
    r = r[np.isfinite(r)]
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    observed = abs(float(z.mean()))
    n = len(z)
    count = 0
    total = 1 << n
    for mask in range(total):
        signs = np.fromiter((1.0 if (mask >> i) & 1 else -1.0 for i in range(n)), float, count=n)
        if abs(float(np.mean(z * signs))) >= observed - 1e-15:
            count += 1
    return count / total


def maxT_primary(donor_wide: pd.DataFrame) -> dict[str, float]:
    complete = donor_wide.dropna(subset=PRIMARY).copy()
    z = np.arctanh(np.clip(complete[PRIMARY].to_numpy(float), -0.999999, 0.999999))
    n = z.shape[0]
    # The three links have genuinely different between-donor variances, so a
    # raw maximum of absolute means is not exchangeable across endpoints.
    # Studentize each endpoint before taking maxT; signs remain jointly flipped
    # by donor to preserve cross-link dependence.
    observed = np.abs(z.mean(axis=0)) / np.maximum(z.std(axis=0, ddof=1) / np.sqrt(n), 1e-12)
    exceed = np.zeros(len(PRIMARY), dtype=np.int64)
    total = 1 << n
    for mask in range(total):
        signs = np.fromiter((1.0 if (mask >> i) & 1 else -1.0 for i in range(n)), float, count=n)[:, None]
        signed = z * signs
        null_t = np.abs(signed.mean(axis=0)) / np.maximum(signed.std(axis=0, ddof=1) / np.sqrt(n), 1e-12)
        null_max = float(np.max(null_t))
        exceed += null_max >= observed - 1e-15
    return {endpoint: exceed[i] / total for i, endpoint in enumerate(PRIMARY)} | {"n_complete": n}


def load_signatures() -> tuple[list[str], list[str], list[str], list[str]]:
    sig = pd.read_csv(SIGNATURE_FILE, sep="\t")
    sender = (
        sig.loc[(sig.signature_id == SENDER_ID) & (sig.signature_rank <= 30)]
        .sort_values("signature_rank")
        .gene_symbol.tolist()
    )
    receiver = (
        sig.loc[(sig.signature_id == RECEIVER_ID) & (sig.signature_rank <= 30)]
        .sort_values("signature_rank")
        .gene_symbol.tolist()
    )
    programs = pd.read_csv(PROGRAM_FILE, sep="\t")
    excluded = set(sender) | set(receiver) | AXIS_EXCLUDE
    deg = programs.loc[programs.program.eq(DEG_PROGRAM), "gene_symbol"].drop_duplicates().tolist()
    adh = programs.loc[programs.program.eq(ADH_PROGRAM), "gene_symbol"].drop_duplicates().tolist()
    deg = [g for g in deg if g not in excluded]
    adh = [g for g in adh if g not in excluded]
    return sender, receiver, deg, adh


def read_selected_expression(h5_path: Path, barcodes: list[str], genes: list[str], totals: np.ndarray) -> pd.DataFrame:
    with h5py.File(h5_path, "r") as fh:
        m = fh["matrix"]
        names = np.asarray(m["features"]["name"]).astype("U")
        h5_barcodes = np.asarray(m["barcodes"]).astype("U")
        matrix = csc_matrix((m["data"][...], m["indices"][...], m["indptr"][...]), shape=tuple(m["shape"][...]))
    barcode_index = {b: i for i, b in enumerate(h5_barcodes)}
    col_idx = np.array([barcode_index.get(b, -1) for b in barcodes], dtype=int)
    if np.any(col_idx < 0):
        missing = np.asarray(barcodes)[col_idx < 0][:5]
        raise RuntimeError(f"QC barcodes missing from H5 {h5_path.name}: {missing}")
    gene_rows: dict[str, list[int]] = {}
    for i, gene in enumerate(names):
        if gene in genes:
            gene_rows.setdefault(gene, []).append(i)
    out = np.zeros((len(gene_rows), len(barcodes)), dtype=float)
    labels = []
    for k, (gene, rows) in enumerate(gene_rows.items()):
        out[k] = np.asarray(matrix[rows][:, col_idx].sum(axis=0)).ravel()
        labels.append(gene)
    norm = np.log1p(out * (1e4 / np.maximum(np.asarray(totals, float), 1.0))[None, :])
    return pd.DataFrame(norm, index=labels, columns=barcodes)


def endpoint_rows_for_slide(base: pd.DataFrame, norm: pd.DataFrame, sender_genes, receiver_genes, deg_genes, adh_genes):
    sender, s_avail, s_info = module_score(norm, sender_genes, 10)
    receiver, r_avail, r_info = module_score(norm, receiver_genes, 10)
    ecm_deg, d_avail, d_info = module_score(norm, deg_genes, 20)
    adhesion, a_avail, a_info = module_score(norm, adh_genes, 20)
    myeloid, m_avail, m_info = module_score(norm, MYELOID_GENES, 6)
    plat = norm.loc["PLAT"].to_numpy(float)
    log_umi = np.log1p(base.total_counts.to_numpy(float))
    fibro = base.MARKER_Fibroblast.to_numpy(float)
    epithelial = base.MARKER_Epithelial.to_numpy(float)
    ecm_direction = base.frozen_ECM_direction.to_numpy(float)
    rows = base.array_row.to_numpy(int)
    cols = base.array_col.to_numpy(int)

    shell = {}
    for ring in range(4):
        shell[("receiver", ring)] = spatial_shell(receiver, rows, cols, ring)
        shell[("myeloid", ring)] = spatial_shell(myeloid, rows, cols, ring)
        shell[("fibro", ring)] = spatial_shell(fibro, rows, cols, ring)
        shell[("deg", ring)] = spatial_shell(ecm_deg, rows, cols, ring)
        shell[("adh", ring)] = spatial_shell(adhesion, rows, cols, ring)

    endpoints: list[tuple[str, np.ndarray, np.ndarray, list[np.ndarray], int]] = []
    endpoints.append((
        "L1_PLAT_sender_same", plat, sender,
        [log_umi, fibro, myeloid, epithelial, ecm_direction], 0,
    ))
    for ring in range(4):
        endpoints.append((
            f"L2_sender_receiver_r{ring}", sender, shell[("receiver", ring)],
            [log_umi, fibro, shell[("myeloid", ring)], epithelial, ecm_direction, plat], ring,
        ))
        endpoints.append((
            f"L3_receiver_ECMdeg_r{ring}", receiver, shell[("deg", ring)],
            [log_umi, myeloid, shell[("fibro", ring)], epithelial, ecm_direction, sender, plat], ring,
        ))
        endpoints.append((
            f"L3_receiver_adhesion_r{ring}", receiver, shell[("adh", ring)],
            [log_umi, myeloid, shell[("fibro", ring)], epithelial, ecm_direction, sender, plat], ring,
        ))
    endpoints.extend([
        (
            "PLAT_incremental_receiver_r1", plat, shell[("receiver", 1)],
            [sender, log_umi, fibro, shell[("myeloid", 1)], epithelial, ecm_direction], 1,
        ),
        (
            "PLAT_incremental_ECMdeg_r1", plat, shell[("deg", 1)],
            [sender, receiver, log_umi, fibro, shell[("myeloid", 1)], epithelial, ecm_direction], 1,
        ),
        (
            "PLAT_incremental_adhesion_r1", plat, shell[("adh", 1)],
            [sender, receiver, log_umi, fibro, shell[("myeloid", 1)], epithelial, ecm_direction], 1,
        ),
    ])

    effect_rows = []
    for endpoint, x, y, cov, ring in endpoints:
        n, rho = residual_rank_correlation(x, y, cov)
        effect_rows.append({"endpoint": endpoint, "ring": ring, "n_spots": n, "rho": rho})
    coverage = {
        "sender_available": s_avail, "sender_informative": s_info,
        "receiver_available": r_avail, "receiver_informative": r_info,
        "ECMdeg_available": d_avail, "ECMdeg_informative": d_info,
        "adhesion_available": a_avail, "adhesion_informative": a_info,
        "myeloid_available": m_avail, "myeloid_informative": m_info,
        "PLAT_detected_fraction": float(np.mean(plat > 0)),
    }
    return effect_rows, coverage


def main():
    sender_genes, receiver_genes, deg_genes, adh_genes = load_signatures()
    needed = sorted(set(sender_genes) | set(receiver_genes) | set(deg_genes) | set(adh_genes) | set(MYELOID_GENES) | {"PLAT"})

    spots = pd.read_csv(SPOT_FILE, sep="\t")
    spots = spots.loc[spots.in_tissue.astype(bool)].copy()
    manifest = pd.read_csv(MANIFEST_FILE, sep="\t")
    manifest = manifest.loc[(manifest.platform_role == "ST") & (~manifest.is_replicate_or_extension.astype(bool))].copy()
    manifest = manifest.drop_duplicates("sample_id")
    if len(manifest) != 33 or spots.donor_id.nunique() != 19:
        raise RuntimeError(f"Expected 33 sections / 19 donors, got {len(manifest)} / {spots.donor_id.nunique()}")

    slide_rows = []
    coverage_rows = []
    for i, row in enumerate(manifest.itertuples(index=False), 1):
        sample = row.sample_id
        print(f"[{i:02d}/33] {sample}", flush=True)
        base = spots.loc[spots.sample_id.eq(sample)].sort_values("barcode").reset_index(drop=True)
        h5_path = RAW / row.filtered_h5_file
        if not h5_path.exists():
            raise FileNotFoundError(h5_path)
        norm = read_selected_expression(
            h5_path, base.barcode.tolist(), needed, base.total_counts.to_numpy(float)
        )
        effects, coverage = endpoint_rows_for_slide(
            base, norm, sender_genes, receiver_genes, deg_genes, adh_genes
        )
        for effect in effects:
            slide_rows.append({
                "sample_id": sample, "donor_id": row.donor_id,
                "stage": row.stage_geo, "lauren": row.lauren_geo,
                "sampling_location": row.sampling_location_geo,
                **effect,
            })
        coverage_rows.append({
            "sample_id": sample, "donor_id": row.donor_id,
            "QC_spots": len(base), **coverage,
        })

    slide = pd.DataFrame(slide_rows)
    coverage = pd.DataFrame(coverage_rows)
    donor = (
        slide.groupby(["donor_id", "endpoint", "ring"], as_index=False)
        .agg(
            rho=("rho", lambda x: fisher_pool(x.to_numpy(float))),
            n_sections=("sample_id", "nunique"),
            median_section_spots=("n_spots", "median"),
        )
    )

    meta_rows = []
    for (endpoint, ring), d in donor.groupby(["endpoint", "ring"], sort=False):
        values = d.rho.to_numpy(float)
        lo, hi = bootstrap_ci(values)
        meta_rows.append({
            "endpoint": endpoint, "ring": ring, "n_donors": int(np.isfinite(values).sum()),
            "pooled_rho": fisher_pool(values), "ci_low": lo, "ci_high": hi,
            "positive_donors": int(np.sum(values > 0)), "negative_donors": int(np.sum(values < 0)),
            "exact_signflip_p": exact_signflip_p(values),
            "primary_family": endpoint in PRIMARY,
        })
    meta = pd.DataFrame(meta_rows)

    wide = donor.pivot(index="donor_id", columns="endpoint", values="rho").reset_index()
    maxT = maxT_primary(wide)
    meta["primary_maxT_FWER_p"] = np.nan
    for endpoint in PRIMARY:
        meta.loc[meta.endpoint.eq(endpoint), "primary_maxT_FWER_p"] = maxT[endpoint]

    slide.to_csv(OUT / "GSE270678_PLAT_spatial_chain_slide_effects.tsv", sep="\t", index=False)
    donor.to_csv(OUT / "GSE270678_PLAT_spatial_chain_donor_effects.tsv", sep="\t", index=False)
    meta.to_csv(OUT / "GSE270678_PLAT_spatial_chain_patient_meta.tsv", sep="\t", index=False)
    coverage.to_csv(OUT / "GSE270678_PLAT_spatial_chain_gene_coverage.tsv", sep="\t", index=False)

    audit = {
        "seed": SEED,
        "dataset": "GSE270678 main ST",
        "sections": int(slide.sample_id.nunique()),
        "donors": int(slide.donor_id.nunique()),
        "primary_endpoints": PRIMARY,
        "primary_maxT_statistic": "joint donor-sign-flip maximum absolute studentized Fisher-z mean",
        "maxT_complete_donors": int(maxT["n_complete"]),
        "ring_offsets": {str(k): v for k, v in RING_OFFSETS.items()},
        "sender_genes": sender_genes,
        "receiver_genes": receiver_genes,
        "ECMdeg_genes_after_exclusion": len(deg_genes),
        "adhesion_genes_after_exclusion": len(adh_genes),
        "input_sha256": {
            SPOT_FILE.name: sha256(SPOT_FILE),
            MANIFEST_FILE.name: sha256(MANIFEST_FILE),
            SIGNATURE_FILE.name: sha256(SIGNATURE_FILE),
            PROGRAM_FILE.name: sha256(PROGRAM_FILE),
            LOCK_FILE.name: sha256(LOCK_FILE),
        },
        "warning": "Mixed-cell Visium spot associations; no cell contact, receptor engagement, mediation or causality claim.",
    }
    with (OUT / "GSE270678_PLAT_spatial_chain_audit.json").open("w", encoding="utf-8") as fh:
        json.dump(audit, fh, indent=2, ensure_ascii=False)

    print("\nPRIMARY RESULTS")
    print(meta.loc[meta.endpoint.isin(PRIMARY)].sort_values("endpoint").to_string(index=False))
    print("\nRING PROFILES")
    print(meta.loc[meta.endpoint.str.startswith(("L2_", "L3_"))].sort_values(["endpoint"]).to_string(index=False))
    print("\nINCREMENTAL")
    print(meta.loc[meta.endpoint.str.startswith("PLAT_incremental")].sort_values("endpoint").to_string(index=False))


def refresh_cached_inference():
    """Recompute the locked donor-level joint inference without rereading H5 files."""
    donor_file = OUT / "GSE270678_PLAT_spatial_chain_donor_effects.tsv"
    meta_file = OUT / "GSE270678_PLAT_spatial_chain_patient_meta.tsv"
    audit_file = OUT / "GSE270678_PLAT_spatial_chain_audit.json"
    donor = pd.read_csv(donor_file, sep="\t")
    meta = pd.read_csv(meta_file, sep="\t")
    wide = donor.pivot(index="donor_id", columns="endpoint", values="rho").reset_index()
    adjusted = maxT_primary(wide)
    meta["primary_maxT_FWER_p"] = np.nan
    for endpoint in PRIMARY:
        meta.loc[meta.endpoint.eq(endpoint), "primary_maxT_FWER_p"] = adjusted[endpoint]
    meta.to_csv(meta_file, sep="\t", index=False)
    audit = json.loads(audit_file.read_text(encoding="utf-8"))
    audit["primary_maxT_statistic"] = "joint donor-sign-flip maximum absolute studentized Fisher-z mean"
    audit["primary_maxT_complete_donors"] = int(adjusted["n_complete"])
    audit_file.write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    print(meta.loc[meta.endpoint.isin(PRIMARY)].sort_values("endpoint").to_string(index=False))


if __name__ == "__main__":
    if "--inference-only" in sys.argv:
        refresh_cached_inference()
    else:
        main()
