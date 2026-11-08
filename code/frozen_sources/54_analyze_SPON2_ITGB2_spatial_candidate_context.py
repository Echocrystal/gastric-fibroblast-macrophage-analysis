"""Locked post hoc spatial analysis of the SPON2-ITGB2 candidate context."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
PUBLIC = Path(r"I:\ECM-new\ECM-public")
BASE_SCRIPT = PROJECT / "scripts" / "35_analyze_GSE270678_PLAT_spatial_chain.py"
LOCK_FILE = PROJECT / "analysis" / "design" / "SPON2_ITGB2_SPATIAL_CANDIDATE_CONTEXT_LOCK_20260826.md"
SIGNATURE_FILE = PROJECT / "resources" / "gene_sets" / "GSE183904_PLAT_ecosystem_bulk_signatures.tsv"
PROGRAM_FILE = PROJECT / "resources" / "gene_sets" / "myeloid_receiver_programs_frozen.tsv"
MANIFEST_270678 = PUBLIC / "manifests" / "v4" / "GSE270678_spatial_sample_manifest.tsv"
MANIFEST_251950 = PUBLIC / "data" / "processed" / "spatial" / "GSE251950" / "GSE251950_spatial_sample_manifest.tsv"
OUT = PROJECT / "analysis" / "results" / "spatial" / "SPON2_ITGB2_candidate_context_v1"
OUT.mkdir(parents=True, exist_ok=True)

SPEC = importlib.util.spec_from_file_location("spatial_chain_base", BASE_SCRIPT)
M = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(M)

SEED = 2026082654
RNG = np.random.default_rng(SEED)
MIN_DETECTION = 0.01
MIN_COMPLETE_SPOTS = 100
MIN_MAIN_DONORS = 15
MIN_EXTERNAL_PATIENTS = 7

SENDER_ID = "sender_CXCL14_POSTN_state_axis_excluded"
RECEIVER_ID = "receiver_C1QC_APOE_specific_axis_excluded"
AXES = [
    ("SPON2_ITGB2", "SPON2", "ITGB2"),
    ("APP_CD74", "APP", "CD74"),
    ("NRG1_HLA_DPB1", "NRG1", "HLA-DPB1"),
    ("CXCL14_CXCR4", "CXCL14", "CXCR4"),
]
PROGRAMS = {
    "ECM_degradation": "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
    "myeloid_migration": "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
    "phagocytosis": "GOBP_PHAGOCYTOSIS",
    "cell_matrix_adhesion": "GOBP_CELL_MATRIX_ADHESION",
}
FOCAL_GENES = {g for _, ligand, receptor in AXES for g in (ligand, receptor)} | {"PLAT"}

DATASETS = {
    "GSE270678_main_ST": {
        "raw": PUBLIC / "data" / "raw" / "spatial" / "GSE270678" / "main_public_files",
        "spot": PUBLIC / "data" / "processed" / "spatial" / "GSE270678" / "v4"
        / "GSE270678_main_ST_project_compartment_spot_scores_v4.tsv.gz",
        "kind": "270678", "role": "ST", "replicate": False,
        "expected_sections": 33, "expected_donors": 19,
    },
    "GSE270678_ST_R_sensitivity": {
        "raw": PUBLIC / "data" / "raw" / "spatial" / "GSE270678" / "ST_R_public_files",
        "spot": PUBLIC / "data" / "processed" / "spatial" / "GSE270678" / "v4"
        / "GSE270678_ST_R_project_compartment_spot_scores_v4.tsv.gz",
        "kind": "270678", "role": "ST_R", "replicate": True,
        "expected_sections": 37, "expected_donors": 16,
    },
    "GSE251950_primary": {
        "raw": PUBLIC / "data" / "raw" / "spatial" / "GSE251950" / "samples",
        "spot": PUBLIC / "data" / "processed" / "spatial" / "GSE251950" / "GSE251950_spot_ecology_scores.tsv.gz",
        "kind": "251950", "expected_sections": 9, "expected_donors": 9,
    },
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def standardize(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    sd = np.nanstd(x, ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        return np.full_like(x, np.nan)
    return (x - np.nanmean(x)) / sd


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


def load_fixed_definitions():
    sig = pd.read_csv(SIGNATURE_FILE, sep="\t")
    sender = sig.loc[(sig.signature_id == SENDER_ID) & (sig.signature_rank <= 30)].sort_values("signature_rank").gene_symbol.tolist()
    receiver = sig.loc[(sig.signature_id == RECEIVER_ID) & (sig.signature_rank <= 30)].sort_values("signature_rank").gene_symbol.tolist()
    if FOCAL_GENES.intersection(sender) or FOCAL_GENES.intersection(receiver):
        raise RuntimeError("A focal ligand/receptor is present in a fixed state signature")
    p = pd.read_csv(PROGRAM_FILE, sep="\t")
    excluded = set(sender) | set(receiver) | FOCAL_GENES
    programmes = {}
    for short, full in PROGRAMS.items():
        genes = p.loc[p.program.eq(full), "gene_symbol"].drop_duplicates().tolist()
        programmes[short] = [g for g in genes if g not in excluded]
    return sender, receiver, programmes


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    restored = np.empty_like(q)
    restored[order] = np.minimum(q, 1.0)
    out[ok] = restored
    return out


def exact_signflip_p(values: np.ndarray) -> float:
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return np.nan
    z = np.arctanh(np.clip(x, -0.999999, 0.999999))
    observed = abs(float(z.mean()))
    total = 1 << len(z)
    exceed = 0
    chunk = 65536
    bit_positions = np.arange(len(z), dtype=np.uint64)
    for start in range(0, total, chunk):
        masks = np.arange(start, min(start + chunk, total), dtype=np.uint64)[:, None]
        signs = (((masks >> bit_positions) & 1).astype(float) * 2.0) - 1.0
        null = np.abs((signs @ z) / len(z))
        exceed += int(np.sum(null >= observed - 1e-15))
    return exceed / total


def bootstrap_ci(values: np.ndarray, draws: int = 9999) -> tuple[float, float]:
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return np.nan, np.nan
    z = np.arctanh(np.clip(x, -0.999999, 0.999999))
    sampled = RNG.choice(z, size=(draws, len(z)), replace=True)
    return tuple(np.quantile(np.tanh(sampled.mean(axis=1)), [0.025, 0.975]).tolist())


def dataset_manifest(name: str, cfg: dict) -> pd.DataFrame:
    if cfg["kind"] == "270678":
        m = pd.read_csv(MANIFEST_270678, sep="\t")
        m = m.loc[
            (m.platform_role == cfg["role"])
            & (m.is_replicate_or_extension.astype(bool) == cfg["replicate"])
        ].drop_duplicates("sample_id").copy()
        if cfg["replicate"]:
            primary = pd.read_csv(MANIFEST_270678, sep="\t")
            main_donors = set(primary.loc[(primary.platform_role == "ST") & (~primary.is_replicate_or_extension.astype(bool)), "donor_id"])
            m = m.loc[m.donor_id.isin(main_donors)].copy()
        m["section_id"] = m.sample_id
        return m
    m = pd.read_csv(MANIFEST_251950, sep="\t")
    m = m.loc[m.primary_analysis.astype(bool)].copy()
    m["section_id"] = m.GSM
    m["donor_id"] = m.patient_id
    return m


def section_inputs(cfg: dict, row, spots: pd.DataFrame):
    if cfg["kind"] == "270678":
        base = spots.loc[spots.sample_id.eq(row.sample_id) & spots.in_tissue.astype(bool)].sort_values("barcode").reset_index(drop=True)
        h5_path = cfg["raw"] / row.filtered_h5_file
    else:
        base = spots.loc[spots.GSM.eq(row.GSM)].sort_values("barcode").reset_index(drop=True)
        base = base.rename(columns={"total_UMI": "total_counts"})
        folder = cfg["raw"] / row.folder_id
        h5_path = folder / f"{row.folder_id}_filtered_feature_bc_matrix.h5"
    totals = h5_totals(h5_path, base.barcode.tolist())
    base["total_counts"] = totals
    return base, h5_path


def module(norm: pd.DataFrame, genes: list[str], minimum: int) -> tuple[np.ndarray, int, int]:
    return M.module_score(norm, genes, minimum)


def analyze_section(base: pd.DataFrame, norm: pd.DataFrame, sender_genes, receiver_genes, programmes):
    sender, s_avail, s_info = module(norm, sender_genes, 10)
    receiver, r_avail, r_info = module(norm, receiver_genes, 10)
    myeloid, m_avail, m_info = module(norm, M.MYELOID_GENES, 6)
    programme_scores = {}
    programme_coverage = {}
    for short, genes in programmes.items():
        score, available, informative = module(norm, genes, 20)
        programme_scores[short] = score
        programme_coverage[short] = (available, informative)

    rows = base.array_row.to_numpy(int)
    cols = base.array_col.to_numpy(int)
    log_umi = np.log1p(base.total_counts.to_numpy(float))
    fibro = base.MARKER_Fibroblast.to_numpy(float)
    epithelial = base.MARKER_Epithelial.to_numpy(float)
    ecm_direction = base.frozen_ECM_direction.to_numpy(float)
    plat = norm.loc["PLAT"].to_numpy(float)
    genes = {g: norm.loc[g].to_numpy(float) for g in FOCAL_GENES if g in norm.index}
    shell_myeloid = {ring: M.spatial_shell(myeloid, rows, cols, ring) for ring in range(4)}
    receiver_z = standardize(receiver)

    effects = []
    def add(endpoint, category, x, y, covariates, ring, axis=""):
        n, rho = M.residual_rank_correlation(x, y, covariates)
        effects.append({"endpoint": endpoint, "category": category, "axis": axis, "ring": ring, "n_spots": n, "rho": rho})

    primary_cov = [log_umi, fibro, shell_myeloid[1], epithelial, ecm_direction, sender, plat]
    itgb2_z = standardize(genes["ITGB2"])
    composite = np.nanmean(np.vstack([itgb2_z, receiver_z]), axis=0)
    for ring in range(4):
        ring_cov = [log_umi, fibro, shell_myeloid[ring], epithelial, ecm_direction, sender, plat]
        add(
            f"SPON2_ITGB2_receiver_composite_r{ring}", "distance_profile",
            genes["SPON2"], M.spatial_shell(composite, rows, cols, ring), ring_cov, ring,
            "SPON2_ITGB2",
        )
    add("SPON2_sender_same", "supporting", genes["SPON2"], sender,
        [log_umi, fibro, myeloid, epithelial, ecm_direction, plat], 0)
    add("SPON2_ITGB2_r1", "supporting", genes["SPON2"], M.spatial_shell(genes["ITGB2"], rows, cols, 1), primary_cov, 1)
    add("SPON2_receiver_r1", "supporting", genes["SPON2"], M.spatial_shell(receiver, rows, cols, 1), primary_cov, 1)

    for axis, ligand, receptor in AXES:
        receptor_composite = np.nanmean(np.vstack([standardize(genes[receptor]), receiver_z]), axis=0)
        add(
            f"candidate_{axis}_r1", "candidate_calibration",
            genes[ligand], M.spatial_shell(receptor_composite, rows, cols, 1),
            primary_cov, 1, axis,
        )

    receiver_cov_r0 = [log_umi, myeloid, fibro, epithelial, ecm_direction, receiver, sender, plat, genes["SPON2"]]
    receiver_cov_r1 = [log_umi, shell_myeloid[1], fibro, epithelial, ecm_direction, receiver, sender, plat, genes["SPON2"]]
    for short, score in programme_scores.items():
        add(f"ITGB2_{short}_r0", "receiver_function", genes["ITGB2"], score, receiver_cov_r0, 0)
        add(
            f"ITGB2_{short}_r1", "receiver_function",
            M.spatial_shell(genes["ITGB2"], rows, cols, 1), M.spatial_shell(score, rows, cols, 1),
            receiver_cov_r1, 1,
        )

    detection = {g: float(np.mean(genes[g] > 0)) if g in genes else 0.0 for g in sorted(FOCAL_GENES)}
    gate = all(g in genes for g in ("SPON2", "ITGB2")) and detection["SPON2"] >= MIN_DETECTION and detection["ITGB2"] >= MIN_DETECTION
    coverage = {
        "section_gate": bool(gate), "sender_available": s_avail, "sender_informative": s_info,
        "receiver_available": r_avail, "receiver_informative": r_info,
        "myeloid_available": m_avail, "myeloid_informative": m_info,
        **{f"{k}_detected_fraction": v for k, v in detection.items()},
    }
    for short, (available, informative) in programme_coverage.items():
        coverage[f"{short}_available"] = available
        coverage[f"{short}_informative"] = informative
    return effects, coverage


def aggregate(section: pd.DataFrame):
    donor = (
        section.loc[section.section_gate & (section.n_spots >= MIN_COMPLETE_SPOTS)]
        .groupby(["dataset", "donor_id", "endpoint", "category", "axis", "ring"], as_index=False)
        .agg(rho=("rho", lambda x: M.fisher_pool(x.to_numpy(float))), n_sections=("section_id", "nunique"), median_spots=("n_spots", "median"))
    )
    meta_rows = []
    for keys, d in donor.groupby(["dataset", "endpoint", "category", "axis", "ring"], sort=False, dropna=False):
        dataset, endpoint, category, axis, ring = keys
        x = d.rho.to_numpy(float)
        lo, hi = bootstrap_ci(x)
        meta_rows.append({
            "dataset": dataset, "endpoint": endpoint, "category": category, "axis": axis, "ring": ring,
            "n_donors": int(np.isfinite(x).sum()), "pooled_rho": M.fisher_pool(x), "ci_low": lo, "ci_high": hi,
            "positive_donors": int(np.sum(x > 0)), "negative_donors": int(np.sum(x < 0)),
            "exact_signflip_p": exact_signflip_p(x),
        })
    meta = pd.DataFrame(meta_rows)
    meta["candidate_BH_FDR"] = np.nan
    meta["supporting_BH_FDR"] = np.nan
    for dataset, d in meta.loc[meta.category.eq("candidate_calibration")].groupby("dataset"):
        idx = d.index.to_numpy()
        meta.loc[idx, "candidate_BH_FDR"] = bh_adjust(d.exact_signflip_p.to_numpy(float))
        ranks = d.pooled_rho.rank(method="min", ascending=False)
        meta.loc[idx, "candidate_effect_rank"] = ranks.to_numpy(float)
    for (dataset, category), d in meta.loc[
        meta.category.isin(["supporting", "receiver_function"])
    ].groupby(["dataset", "category"]):
        idx = d.index.to_numpy()
        meta.loc[idx, "supporting_BH_FDR"] = bh_adjust(d.exact_signflip_p.to_numpy(float))
    return donor, meta


def write_report(meta: pd.DataFrame, coverage: pd.DataFrame):
    lines = ["# SPON2-ITGB2 spatial candidate-context results", "", "## Coverage", ""]
    for dataset, d in coverage.groupby("dataset"):
        lines.append(f"- {dataset}: {d.section_id.nunique()} sections, {d.loc[d.section_gate, 'donor_id'].nunique()} donors/patients passed the SPON2/ITGB2 detection gate.")
    lines.extend(["", "## Primary ring-1 endpoint", ""])
    p = meta.loc[meta.endpoint.eq("SPON2_ITGB2_receiver_composite_r1")]
    for row in p.itertuples(index=False):
        lines.append(
            f"- {row.dataset}: n={row.n_donors}, pooled rho={row.pooled_rho:.6f}, 95% CI {row.ci_low:.6f} to {row.ci_high:.6f}, "
            f"positive={row.positive_donors}/{row.n_donors}, sign-flip P={row.exact_signflip_p:.6g}."
        )
    lines.extend(["", "## Candidate calibration", ""])
    c = meta.loc[meta.category.eq("candidate_calibration")].sort_values(["dataset", "candidate_effect_rank"])
    for row in c.itertuples(index=False):
        lines.append(
            f"- {row.dataset}, {row.axis}: rho={row.pooled_rho:.6f}, rank={row.candidate_effect_rank:.0f}/4, "
            f"P={row.exact_signflip_p:.6g}, FDR={row.candidate_BH_FDR:.6g}."
        )
    lines.extend([
        "", "## Evidence boundary", "",
        "These are mixed-cell Visium proximity-conditioned transcriptomic associations. They do not establish direct ligand-receptor binding, receptor activation, mediation or causality.",
    ])
    (OUT / "SPON2_ITGB2_spatial_result_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(section: pd.DataFrame, coverage: pd.DataFrame):
    donor, meta = aggregate(section)
    section.to_csv(OUT / "SPON2_ITGB2_section_effects.tsv", sep="\t", index=False)
    donor.to_csv(OUT / "SPON2_ITGB2_donor_effects.tsv", sep="\t", index=False)
    meta.to_csv(OUT / "SPON2_ITGB2_meta_summary.tsv", sep="\t", index=False)
    coverage.to_csv(OUT / "SPON2_ITGB2_gene_endpoint_coverage.tsv", sep="\t", index=False)
    meta.loc[meta.category.eq("distance_profile")].to_csv(
        OUT / "SPON2_ITGB2_ring_distance_profiles.tsv", sep="\t", index=False
    )
    meta.loc[meta.category.eq("candidate_calibration")].to_csv(
        OUT / "SPON2_ITGB2_candidate_axis_calibration.tsv", sep="\t", index=False
    )
    write_report(meta, coverage)
    return donor, meta


def main():
    sender_genes, receiver_genes, programmes = load_fixed_definitions()
    needed = sorted(set(sender_genes) | set(receiver_genes) | set(M.MYELOID_GENES) | FOCAL_GENES | {g for genes in programmes.values() for g in genes})
    all_effects, all_coverage = [], []
    input_files = [Path(__file__), BASE_SCRIPT, LOCK_FILE, SIGNATURE_FILE, PROGRAM_FILE, MANIFEST_270678, MANIFEST_251950]

    for dataset, cfg in DATASETS.items():
        manifest = dataset_manifest(dataset, cfg)
        spots = pd.read_csv(cfg["spot"], sep="\t")
        input_files.append(cfg["spot"])
        if len(manifest) != cfg["expected_sections"] or manifest.donor_id.nunique() != cfg["expected_donors"]:
            raise RuntimeError(f"{dataset}: expected {cfg['expected_sections']}/{cfg['expected_donors']} sections/donors, observed {len(manifest)}/{manifest.donor_id.nunique()}")
        for i, row in enumerate(manifest.itertuples(index=False), 1):
            print(f"[{dataset} {i:02d}/{len(manifest)}] {row.section_id}", flush=True)
            base, h5_path = section_inputs(cfg, row, spots)
            totals = base.total_counts.to_numpy(float)
            norm = M.read_selected_expression(h5_path, base.barcode.tolist(), needed, totals)
            missing = sorted(FOCAL_GENES.difference(norm.index))
            if missing:
                raise RuntimeError(f"{dataset}/{row.section_id}: missing focal genes {missing}")
            effects, coverage = analyze_section(base, norm, sender_genes, receiver_genes, programmes)
            for effect in effects:
                all_effects.append({"dataset": dataset, "section_id": row.section_id, "donor_id": row.donor_id, **effect})
            all_coverage.append({"dataset": dataset, "section_id": row.section_id, "donor_id": row.donor_id, "n_spots": len(base), **coverage})

    section = pd.DataFrame(all_effects)
    coverage = pd.DataFrame(all_coverage)
    section = section.merge(coverage[["dataset", "section_id", "section_gate"]], on=["dataset", "section_id"], how="left")
    donor, meta = write_outputs(section, coverage)

    main_donors = int(coverage.loc[(coverage.dataset == "GSE270678_main_ST") & coverage.section_gate, "donor_id"].nunique())
    external = int(coverage.loc[(coverage.dataset == "GSE251950_primary") & coverage.section_gate, "donor_id"].nunique())
    audit = {
        "seed": SEED, "analysis_status": "locked_post_hoc", "output_directory": str(OUT),
        "primary_endpoint": "SPON2_ITGB2_receiver_composite_r1",
        "coverage_thresholds": {"minimum_gene_detection_fraction": MIN_DETECTION, "minimum_complete_spots": MIN_COMPLETE_SPOTS},
        "main_formal_gate": {"required": MIN_MAIN_DONORS, "observed": main_donors, "passed": main_donors >= MIN_MAIN_DONORS},
        "external_directional_gate": {"required": MIN_EXTERNAL_PATIENTS, "observed": external, "passed": external >= MIN_EXTERNAL_PATIENTS},
        "candidate_axes": [axis for axis, _, _ in AXES],
        "spatial_graph": {"type": "exact_Visium_hex_shells", "rings": [0, 1, 2, 3], "offsets": {str(k): v for k, v in M.RING_OFFSETS.items()}},
        "input_sha256": {str(path): sha256(path) for path in input_files},
        "warning": "Mixed-cell Visium proximity-conditioned co-expression; no direct binding, receptor activation, mediation or causality claim.",
    }
    (OUT / "SPON2_ITGB2_analysis_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")

    print("\nPRIMARY ENDPOINT")
    print(meta.loc[meta.endpoint.eq("SPON2_ITGB2_receiver_composite_r1")].to_string(index=False))
    print("\nCANDIDATE CALIBRATION")
    print(meta.loc[meta.category.eq("candidate_calibration")].sort_values(["dataset", "candidate_effect_rank"]).to_string(index=False))
    print("\nRECEIVER FUNCTION")
    print(meta.loc[meta.category.eq("receiver_function")].sort_values(["dataset", "endpoint"]).to_string(index=False))


def refresh_cached_outputs():
    section = pd.read_csv(OUT / "SPON2_ITGB2_section_effects.tsv", sep="\t", keep_default_na=False)
    coverage = pd.read_csv(OUT / "SPON2_ITGB2_gene_endpoint_coverage.tsv", sep="\t")
    _, meta = write_outputs(section, coverage)
    audit_file = OUT / "SPON2_ITGB2_analysis_audit.json"
    audit = json.loads(audit_file.read_text(encoding="utf-8"))
    audit["input_sha256"][str(Path(__file__))] = sha256(Path(__file__))
    audit_file.write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    print(meta.loc[meta.category.isin(["candidate_calibration", "receiver_function"])].to_string(index=False))


if __name__ == "__main__":
    if "--aggregate-only" in sys.argv:
        refresh_cached_outputs()
    else:
        main()
