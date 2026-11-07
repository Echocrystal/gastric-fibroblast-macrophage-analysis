"""Frozen PLAT-positive matched-neighbourhood and conditional-attenuation analyses.

See analysis/design/PLAT_SPATIAL_ENRICHMENT_LOCK_20260806.md.
Main ST and same-donor ST_R sensitivity are analysed separately.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.spatial import cKDTree


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
PUBLIC = Path(r"I:\ECM-new\ECM-public")
BASE_SCRIPT = PROJECT / "scripts" / "35_analyze_GSE270678_PLAT_spatial_chain.py"
LOCK_FILE = PROJECT / "analysis" / "design" / "PLAT_SPATIAL_ENRICHMENT_LOCK_20260806.md"
MANIFEST_FILE = PUBLIC / "manifests" / "v4" / "GSE270678_spatial_sample_manifest.tsv"
OUT = PROJECT / "analysis" / "results" / "spatial" / "GSE270678_PLAT_ecosystem_v1" / "matched_neighborhood_v1"
OUT.mkdir(parents=True, exist_ok=True)

SPEC = importlib.util.spec_from_file_location("plat_chain_base", BASE_SCRIPT)
M = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(M)

SEED = 2026080637
RNG = np.random.default_rng(SEED)
CALIPER = 1.5
MIN_PAIRS = 15
MAX_POST_SMD = 0.25
MATCH_COVARIATES = ["sender", "log_umi", "fibroblast", "myeloid", "epithelial", "ecm_direction"]
MATCH_ENDPOINTS = [
    "matched_receiver_r0", "matched_receiver_r1", "matched_receiver_r2", "matched_receiver_r3",
    "matched_ECMdeg_r1", "matched_adhesion_r1",
]
PRIMARY_MATCHED = "matched_receiver_r1"
LADDER_MODELS = ["M0_base", "M1_plus_sender", "M2_plus_receiver_r1", "M3_plus_sender_receiver_r1"]

DATASETS = {
    "Main_ST": {
        "raw": PUBLIC / "data" / "raw" / "spatial" / "GSE270678" / "main_public_files",
        "spot": PUBLIC / "data" / "processed" / "spatial" / "GSE270678" / "v4"
        / "GSE270678_main_ST_project_compartment_spot_scores_v4.tsv.gz",
        "role": "ST", "replicate": False, "expected_sections": 33, "expected_donors": 19,
        "formal_min_donors": 15,
    },
    "ST_R_sensitivity": {
        "raw": PUBLIC / "data" / "raw" / "spatial" / "GSE270678" / "ST_R_public_files",
        "spot": PUBLIC / "data" / "processed" / "spatial" / "GSE270678" / "v4"
        / "GSE270678_ST_R_project_compartment_spot_scores_v4.tsv.gz",
        "role": "ST_R", "replicate": True, "expected_sections": 37, "expected_donors": 16,
        "formal_min_donors": 12,
    },
}


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


def standardize(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    sd = np.nanstd(x, ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        return np.zeros_like(x)
    return (x - np.nanmean(x)) / sd


def exact_signflip_mean_p(values: np.ndarray) -> float:
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    observed = abs(float(np.mean(x)))
    total = 1 << len(x)
    exceed = 0
    for mask in range(total):
        signs = np.fromiter((1.0 if (mask >> i) & 1 else -1.0 for i in range(len(x))), float, count=len(x))
        if abs(float(np.mean(x * signs))) >= observed - 1e-15:
            exceed += 1
    return exceed / total


def bootstrap_mean_ci(values: np.ndarray, draws: int = 9999) -> tuple[float, float]:
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    sampled = RNG.choice(x, size=(draws, len(x)), replace=True)
    return tuple(np.quantile(sampled.mean(axis=1), [0.025, 0.975]).tolist())


def standardized_mean_difference(x: np.ndarray, group: np.ndarray) -> float:
    a = np.asarray(x, float)[group]
    b = np.asarray(x, float)[~group]
    pooled = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2.0)
    if not np.isfinite(pooled) or pooled <= 0:
        return 0.0
    return float((np.mean(a) - np.mean(b)) / pooled)


def deterministic_match(zcov: np.ndarray, treated: np.ndarray, strata: np.ndarray):
    """Global-shortest greedy 1:1 matching without replacement within strata."""
    selected_pairs: list[tuple[int, int, float]] = []
    for stratum in sorted(np.unique(strata).tolist()):
        ti = np.flatnonzero(treated & (strata == stratum))
        ci = np.flatnonzero((~treated) & (strata == stratum))
        if not len(ti) or not len(ci):
            continue
        tree = cKDTree(zcov[ci])
        k = min(20, len(ci))
        dist, nn = tree.query(zcov[ti], k=k)
        if k == 1:
            dist = dist[:, None]
            nn = nn[:, None]
        candidates = []
        for i in range(len(ti)):
            for j in range(k):
                d = float(dist[i, j])
                if np.isfinite(d) and d <= CALIPER:
                    candidates.append((d, int(ti[i]), int(ci[int(nn[i, j])])))
        candidates.sort(key=lambda x: (x[0], x[1], x[2]))
        used_t, used_c = set(), set()
        for d, t, c in candidates:
            if t not in used_t and c not in used_c:
                selected_pairs.append((t, c, d))
                used_t.add(t)
                used_c.add(c)
    return selected_pairs


def compute_scores(base: pd.DataFrame, norm: pd.DataFrame, sender_genes, receiver_genes, deg_genes, adh_genes):
    sender, _, _ = M.module_score(norm, sender_genes, 10)
    receiver, _, _ = M.module_score(norm, receiver_genes, 10)
    deg, _, _ = M.module_score(norm, deg_genes, 20)
    adhesion, _, _ = M.module_score(norm, adh_genes, 20)
    myeloid, _, _ = M.module_score(norm, M.MYELOID_GENES, 6)
    plat = norm.loc["PLAT"].to_numpy(float)
    rows = base.array_row.to_numpy(int)
    cols = base.array_col.to_numpy(int)

    receiver_z = standardize(receiver)
    deg_z = standardize(deg)
    adhesion_z = standardize(adhesion)
    shell_receiver = {r: M.spatial_shell(receiver_z, rows, cols, r) for r in range(4)}
    shell_deg = {r: M.spatial_shell(deg_z, rows, cols, r) for r in range(4)}
    shell_adh = {r: M.spatial_shell(adhesion_z, rows, cols, r) for r in range(4)}
    shell_myeloid = {r: M.spatial_shell(myeloid, rows, cols, r) for r in range(4)}

    scores = {
        "PLAT": plat,
        "sender": sender,
        "receiver": receiver,
        "receiver_r1_raw": M.spatial_shell(receiver, rows, cols, 1),
        "log_umi": np.log1p(base.total_counts.to_numpy(float)),
        "fibroblast": base.MARKER_Fibroblast.to_numpy(float),
        "myeloid": myeloid,
        "myeloid_r1": shell_myeloid[1],
        "epithelial": base.MARKER_Epithelial.to_numpy(float),
        "ecm_direction": base.frozen_ECM_direction.to_numpy(float),
        "stromal_dominance": base.stromal_dominance.astype(bool).to_numpy(),
        "deg_r1_raw": M.spatial_shell(deg, rows, cols, 1),
        "adh_r1_raw": M.spatial_shell(adhesion, rows, cols, 1),
    }
    outcomes = {f"matched_receiver_r{r}": shell_receiver[r] for r in range(4)}
    outcomes["matched_ECMdeg_r1"] = shell_deg[1]
    outcomes["matched_adhesion_r1"] = shell_adh[1]
    return scores, outcomes


def match_slide(base: pd.DataFrame, scores: dict[str, np.ndarray], outcomes: dict[str, np.ndarray]):
    eligible = (
        np.isfinite(scores["receiver_r1_raw"])
        & (scores["sender"] >= np.nanmedian(scores["sender"]))
        & (scores["fibroblast"] >= np.nanmedian(scores["fibroblast"]))
    )
    for cov in MATCH_COVARIATES:
        eligible &= np.isfinite(scores[cov])
    idx = np.flatnonzero(eligible)
    if not len(idx):
        return [], {}, {"n_eligible": 0, "n_treated": 0, "n_controls": 0, "n_pairs": 0, "max_abs_post_smd": np.nan}

    cov = np.column_stack([scores[x][idx] for x in MATCH_COVARIATES])
    zcov = np.column_stack([standardize(cov[:, j]) for j in range(cov.shape[1])])
    treated = scores["PLAT"][idx] > 0
    strata = scores["stromal_dominance"][idx].astype(int)

    pre_smd = {x: standardized_mean_difference(cov[:, j], treated) for j, x in enumerate(MATCH_COVARIATES)}
    local_pairs = deterministic_match(zcov, treated, strata)
    if not local_pairs:
        return [], {}, {
            "n_eligible": len(idx), "n_treated": int(treated.sum()), "n_controls": int((~treated).sum()),
            "n_pairs": 0, "max_abs_post_smd": np.nan, **{f"pre_smd_{k}": v for k, v in pre_smd.items()},
        }
    pairs = [(int(idx[t]), int(idx[c]), d) for t, c, d in local_pairs]
    t_global = np.array([x[0] for x in pairs], int)
    c_global = np.array([x[1] for x in pairs], int)
    matched_cov = np.vstack([
        np.column_stack([scores[x][t_global] for x in MATCH_COVARIATES]),
        np.column_stack([scores[x][c_global] for x in MATCH_COVARIATES]),
    ])
    matched_group = np.r_[np.ones(len(pairs), dtype=bool), np.zeros(len(pairs), dtype=bool)]
    post_smd = {
        x: standardized_mean_difference(matched_cov[:, j], matched_group)
        for j, x in enumerate(MATCH_COVARIATES)
    }
    effects = {
        endpoint: float(np.nanmean(np.asarray(values)[t_global] - np.asarray(values)[c_global]))
        for endpoint, values in outcomes.items()
    }
    diag = {
        "n_eligible": len(idx), "n_treated": int(treated.sum()), "n_controls": int((~treated).sum()),
        "n_pairs": len(pairs), "median_match_distance": float(np.median([x[2] for x in pairs])),
        "max_match_distance": float(np.max([x[2] for x in pairs])),
        "max_abs_post_smd": float(max(abs(v) for v in post_smd.values())),
        **{f"pre_smd_{k}": v for k, v in pre_smd.items()},
        **{f"post_smd_{k}": v for k, v in post_smd.items()},
    }
    pair_rows = []
    for pair_id, (t, c, distance) in enumerate(pairs, 1):
        row = {
            "pair_id": pair_id, "treated_barcode": base.barcode.iloc[t], "control_barcode": base.barcode.iloc[c],
            "match_distance": distance, "treated_PLAT": scores["PLAT"][t], "control_PLAT": scores["PLAT"][c],
        }
        for endpoint, values in outcomes.items():
            row[f"treated_{endpoint}"] = values[t]
            row[f"control_{endpoint}"] = values[c]
            row[f"difference_{endpoint}"] = values[t] - values[c]
        pair_rows.append(row)
    return pair_rows, effects, diag


def ladder_slide(scores: dict[str, np.ndarray]):
    base = [scores["log_umi"], scores["fibroblast"], scores["myeloid_r1"], scores["epithelial"], scores["ecm_direction"]]
    models = {
        "M0_base": base,
        "M1_plus_sender": base + [scores["sender"]],
        "M2_plus_receiver_r1": base + [scores["receiver_r1_raw"]],
        "M3_plus_sender_receiver_r1": base + [scores["sender"], scores["receiver_r1_raw"]],
    }
    outcomes = {"ECMdeg_r1": scores["deg_r1_raw"], "adhesion_r1": scores["adh_r1_raw"]}
    rows = []
    for endpoint, outcome in outcomes.items():
        for model, covariates in models.items():
            n, rho = M.residual_rank_correlation(scores["PLAT"], outcome, covariates)
            rows.append({"endpoint": endpoint, "model": model, "n_spots": n, "rho": rho})
    return rows


def aggregate_matched(slide: pd.DataFrame, dataset: str, formal_min: int):
    qualified = slide.loc[slide.quality_gate].copy()
    donor = qualified.groupby(["donor_id", "endpoint"], as_index=False).agg(
        effect=("effect", "mean"), n_sections=("sample_id", "nunique"), median_pairs=("n_pairs", "median")
    )
    meta_rows = []
    for endpoint, d in donor.groupby("endpoint", sort=False):
        x = d.effect.to_numpy(float)
        lo, hi = bootstrap_mean_ci(x)
        meta_rows.append({
            "dataset": dataset, "endpoint": endpoint, "n_donors": len(x),
            "mean_paired_difference": float(np.mean(x)), "ci_low": lo, "ci_high": hi,
            "positive_donors": int(np.sum(x > 0)), "negative_donors": int(np.sum(x < 0)),
            "exact_signflip_p": exact_signflip_mean_p(x),
            "primary_endpoint": endpoint == PRIMARY_MATCHED,
            "formal_donor_gate": len(x) >= formal_min,
        })
    return donor, pd.DataFrame(meta_rows)


def aggregate_ladder(slide: pd.DataFrame, dataset: str):
    donor = slide.groupby(["donor_id", "endpoint", "model"], as_index=False).agg(
        rho=("rho", lambda x: M.fisher_pool(x.to_numpy(float))), n_sections=("sample_id", "nunique")
    )
    meta_rows = []
    for (endpoint, model), d in donor.groupby(["endpoint", "model"], sort=False):
        x = d.rho.to_numpy(float)
        lo, hi = M.bootstrap_ci(x)
        meta_rows.append({
            "dataset": dataset, "endpoint": endpoint, "model": model, "n_donors": len(x),
            "pooled_rho": M.fisher_pool(x), "ci_low": lo, "ci_high": hi,
            "positive_donors": int(np.sum(x > 0)), "negative_donors": int(np.sum(x < 0)),
            "exact_signflip_p": M.exact_signflip_p(x),
        })
    meta = pd.DataFrame(meta_rows)
    for endpoint in meta.endpoint.unique():
        m0 = float(meta.loc[(meta.endpoint == endpoint) & (meta.model == "M0_base"), "pooled_rho"].iloc[0])
        if m0 != 0:
            mask = meta.endpoint == endpoint
            meta.loc[mask, "attenuation_vs_M0_fraction"] = 1.0 - meta.loc[mask, "pooled_rho"] / m0
    return donor, meta


def main():
    sender_genes, receiver_genes, deg_genes, adh_genes = M.load_signatures()
    needed = sorted(set(sender_genes) | set(receiver_genes) | set(deg_genes) | set(adh_genes) | set(M.MYELOID_GENES) | {"PLAT"})
    manifest_all = pd.read_csv(MANIFEST_FILE, sep="\t")
    main_donors = set(manifest_all.loc[(manifest_all.platform_role == "ST") & (~manifest_all.is_replicate_or_extension.astype(bool)), "donor_id"])

    all_match_slide, all_pairs, all_ladder_slide = [], [], []
    dataset_audit = {}
    for dataset, cfg in DATASETS.items():
        spots = pd.read_csv(cfg["spot"], sep="\t")
        spots = spots.loc[spots.in_tissue.astype(bool)].copy()
        manifest = manifest_all.loc[
            (manifest_all.platform_role == cfg["role"])
            & (manifest_all.is_replicate_or_extension.astype(bool) == cfg["replicate"])
        ].drop_duplicates("sample_id")
        if dataset == "ST_R_sensitivity":
            manifest = manifest.loc[manifest.donor_id.isin(main_donors)].copy()
            spots = spots.loc[spots.donor_id.isin(main_donors)].copy()
        if len(manifest) != cfg["expected_sections"] or manifest.donor_id.nunique() != cfg["expected_donors"]:
            raise RuntimeError(f"{dataset}: unexpected sections/donors {len(manifest)}/{manifest.donor_id.nunique()}")

        for i, row in enumerate(manifest.itertuples(index=False), 1):
            print(f"[{dataset} {i:02d}/{len(manifest)}] {row.sample_id}", flush=True)
            base = spots.loc[spots.sample_id.eq(row.sample_id)].sort_values("barcode").reset_index(drop=True)
            h5_path = cfg["raw"] / row.filtered_h5_file
            totals = h5_totals(h5_path, base.barcode.tolist())
            base["total_counts"] = totals
            norm = M.read_selected_expression(h5_path, base.barcode.tolist(), needed, totals)
            scores, outcomes = compute_scores(base, norm, sender_genes, receiver_genes, deg_genes, adh_genes)
            pair_rows, effects, diag = match_slide(base, scores, outcomes)
            quality = (
                diag.get("n_pairs", 0) >= MIN_PAIRS
                and np.isfinite(diag.get("max_abs_post_smd", np.nan))
                and diag["max_abs_post_smd"] <= MAX_POST_SMD
            )
            for endpoint in MATCH_ENDPOINTS:
                all_match_slide.append({
                    "dataset": dataset, "sample_id": row.sample_id, "donor_id": row.donor_id,
                    "endpoint": endpoint, "effect": effects.get(endpoint, np.nan),
                    "quality_gate": quality, **diag,
                })
            for pair in pair_rows:
                all_pairs.append({"dataset": dataset, "sample_id": row.sample_id, "donor_id": row.donor_id, **pair})
            for result in ladder_slide(scores):
                all_ladder_slide.append({"dataset": dataset, "sample_id": row.sample_id, "donor_id": row.donor_id, **result})

        dataset_audit[dataset] = {
            "sections": int(len(manifest)), "donors": int(manifest.donor_id.nunique()),
            "formal_min_donors": cfg["formal_min_donors"],
        }

    matched_slide = pd.DataFrame(all_match_slide)
    pairs = pd.DataFrame(all_pairs)
    ladder_slide_df = pd.DataFrame(all_ladder_slide)
    matched_donors, matched_meta, ladder_donors, ladder_meta = [], [], [], []
    for dataset, cfg in DATASETS.items():
        d, m = aggregate_matched(matched_slide.loc[matched_slide.dataset == dataset], dataset, cfg["formal_min_donors"])
        d["dataset"] = dataset
        matched_donors.append(d)
        matched_meta.append(m)
        ld, lm = aggregate_ladder(ladder_slide_df.loc[ladder_slide_df.dataset == dataset], dataset)
        ld["dataset"] = dataset
        ladder_donors.append(ld)
        ladder_meta.append(lm)

    matched_donor = pd.concat(matched_donors, ignore_index=True)
    matched_meta = pd.concat(matched_meta, ignore_index=True)
    ladder_donor = pd.concat(ladder_donors, ignore_index=True)
    ladder_meta = pd.concat(ladder_meta, ignore_index=True)

    matched_slide.to_csv(OUT / "GSE270678_PLAT_matched_neighborhood_slide_effects.tsv", sep="\t", index=False)
    pairs.to_csv(OUT / "GSE270678_PLAT_matched_neighborhood_pairs.tsv.gz", sep="\t", index=False, compression="gzip")
    matched_donor.to_csv(OUT / "GSE270678_PLAT_matched_neighborhood_donor_effects.tsv", sep="\t", index=False)
    matched_meta.to_csv(OUT / "GSE270678_PLAT_matched_neighborhood_meta.tsv", sep="\t", index=False)
    ladder_slide_df.to_csv(OUT / "GSE270678_PLAT_conditional_ladder_slide_effects.tsv", sep="\t", index=False)
    ladder_donor.to_csv(OUT / "GSE270678_PLAT_conditional_ladder_donor_effects.tsv", sep="\t", index=False)
    ladder_meta.to_csv(OUT / "GSE270678_PLAT_conditional_ladder_meta.tsv", sep="\t", index=False)

    audit = {
        "seed": SEED, "lock_file": str(LOCK_FILE), "caliper": CALIPER,
        "minimum_pairs_per_section": MIN_PAIRS, "maximum_post_match_abs_SMD": MAX_POST_SMD,
        "matching_covariates": MATCH_COVARIATES, "primary_matched_endpoint": PRIMARY_MATCHED,
        "datasets": dataset_audit,
        "warning": "Mixed-cell Visium matched association; no recruitment, mediation or causality claim.",
    }
    (OUT / "GSE270678_PLAT_matched_neighborhood_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("\nMATCHING QUALITY")
    q = matched_slide.drop_duplicates(["dataset", "sample_id"])
    print(q.groupby("dataset").agg(
        sections=("sample_id", "nunique"), qualified=("quality_gate", "sum"),
        median_pairs=("n_pairs", "median"), median_max_post_smd=("max_abs_post_smd", "median")
    ).to_string())
    print("\nMATCHED META")
    print(matched_meta.to_string(index=False))
    print("\nCONDITIONAL LADDER")
    print(ladder_meta.to_string(index=False))


if __name__ == "__main__":
    main()
