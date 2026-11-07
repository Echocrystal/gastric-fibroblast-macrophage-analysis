#!/usr/bin/env python3
"""Prepare an unbiased, GSE206785-centred source-data package for Figure 1.

This script deliberately separates three estimands:

1. a sender-side ECM-ligand universe defined across all eligible tumour
   fibroblasts, without using a receiver or outcome;
2. PLAT localisation across lineages and recurrent fibroblast states; and
3. predicted cell-cell communication inferred after donor-stratified sampling
   of all major tumour lineages.

The GEO matrix contains author-processed log1p(integer counts), cells in rows.
The script streams the dense text file and never treats cells as biological
replicates.  Plotting is performed later, exclusively in R.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import liana as li
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
SOURCE = Path(r"I:\ECM-new\ECM-public")
RAW = SOURCE / "data/raw/scrna/GSE206785"
OUT = PROJECT / "analysis/results/figure1_gse206785_unbiased"
FIGDATA = PROJECT / "figures/data/Figure1_GSE206785_unbiased"
AUDIT = PROJECT / "analysis/audit"
for directory in (OUT, FIGDATA, AUDIT):
    directory.mkdir(parents=True, exist_ok=True)

SEED = 20260802
RNG = np.random.default_rng(SEED)
COMM_CAP_PER_PATIENT_TYPE = 50
FIBRO_CAP_PER_PATIENT = 150
EXPR_PROP = 0.10
MIN_PATIENTS = 5
MIN_CELLS = 10
N_PERMS = 1000

META_PATH = RAW / "GSE206785_metadata.txt.gz"
EXPR_PATH = RAW / "GSE206785_scgex.txt.gz"
LR_MANIFEST = PROJECT / "resources/ligand_receptor/lr_gene_universe_frozen.tsv"
OLD_MEMBERSHIP = PROJECT / "figures/data/Figure1B_four_set_gene_membership.tsv.gz"
FINE_EXPR = PROJECT / "analysis/processed/GSE206785_tumor_normal_fine_state_lr_expression.tsv.gz"
STATE_SUMMARY = (
    PROJECT
    / "analysis/results/cell_communication/Figure1_GSE206785/"
      "GSE206785_LIANA_consensus_state_balanced_source_target_summary.tsv"
)
MATRISOME = SOURCE / "gene_sets/naba_matrisome_locked_membership.tsv"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def percentile(series: pd.Series) -> pd.Series:
    return series.rank(method="average", pct=True)


def capped_indices(frame: pd.DataFrame, group_cols: list[str], cap: int) -> np.ndarray:
    selected: list[int] = []
    for _, group in frame.groupby(group_cols, sort=True, observed=True):
        idx = group.index.to_numpy(dtype=int)
        if len(idx) > cap:
            idx = np.sort(RNG.choice(idx, size=cap, replace=False))
        selected.extend(idx.tolist())
    return np.asarray(sorted(selected), dtype=int)


def clean_gene_for_hvg(gene: str) -> bool:
    return not bool(
        re.match(r"^(MT-|RPL|RPS|HBA|HBB|HBD|HBG|MTRNR|RP11-|RP[0-9]+-)", gene)
    )


def complex_detected(
    matrix_bool: sparse.csr_matrix,
    rows: np.ndarray,
    complex_name: str,
    gene_to_col: dict[str, int],
) -> float:
    parts = [p for p in str(complex_name).split("_") if p]
    cols = [gene_to_col.get(p) for p in parts]
    if not parts or any(c is None for c in cols) or len(rows) == 0:
        return math.nan
    values = matrix_bool[rows][:, cols].toarray()
    return float(np.mean(np.all(values, axis=1)))


def aggregate_liana(result: pd.DataFrame) -> pd.DataFrame:
    robust = result[
        result["magnitude_rank"].le(0.05)
        & result["specificity_rank"].le(0.05)
    ].copy()
    robust["joint_rank"] = robust[["magnitude_rank", "specificity_rank"]].max(axis=1)
    robust["joint_support"] = -np.log10(np.clip(robust["joint_rank"], 1e-12, 1.0))
    tested = (
        result.groupby(["source", "target"], observed=True)
        .agg(
            n_tested=("ligand_complex", "size"),
            best_magnitude_rank=("magnitude_rank", "min"),
            best_specificity_rank=("specificity_rank", "min"),
        )
        .reset_index()
    )
    support = (
        robust.groupby(["source", "target"], observed=True)
        .agg(n_robust=("ligand_complex", "size"), evidence_weight=("joint_support", "sum"))
        .reset_index()
    )
    out = tested.merge(support, on=["source", "target"], how="left")
    out[["n_robust", "evidence_weight"]] = out[["n_robust", "evidence_weight"]].fillna(0)
    out["n_robust"] = out["n_robust"].astype(int)
    return out


def main() -> None:
    metadata = pd.read_csv(META_PATH, sep="\t", dtype=str)
    expected = ["Sample", "Patient", "Tissue", "Platform", "Subtype", "Type", "Annotation"]
    if metadata.columns.tolist() != expected or metadata.isna().any().any():
        raise RuntimeError("Unexpected or incomplete GSE206785 metadata")

    low_quality = metadata["Annotation"].str.contains("lowQual|UnIdent", case=False, regex=True)
    valid = metadata["Tissue"].eq("Tumor") & ~low_quality & ~metadata["Type"].eq("Glial")
    eligible = metadata.loc[valid].copy()
    comm_indices = capped_indices(
        eligible, ["Patient", "Type"], COMM_CAP_PER_PATIENT_TYPE
    )
    fibro_indices = capped_indices(
        eligible.loc[eligible["Type"].eq("Fibroblast")], ["Patient"], FIBRO_CAP_PER_PATIENT
    )
    selected_union = set(comm_indices.tolist()) | set(fibro_indices.tolist())
    comm_position = {idx: pos for pos, idx in enumerate(comm_indices.tolist())}

    with gzip.open(EXPR_PATH, "rt", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split(",")
    if header[0] != "Cell" or len(header[1:]) != len(set(header[1:])):
        raise RuntimeError("Unexpected or duplicated GSE206785 expression header")
    all_genes = np.asarray(header[1:], dtype=object)
    gene_to_all_col = {g: i for i, g in enumerate(all_genes.tolist())}

    lr_manifest = pd.read_csv(LR_MANIFEST, sep="\t")
    lr_genes = sorted(
        set(
            lr_manifest.loc[lr_manifest["available_GSE206785"].astype(bool), "gene_symbol"]
            .astype(str)
        )
    )

    # State-variable LR genes give the display embedding biological breadth
    # without reading every gene for every cell. They do not enter candidate ranking.
    fine_small = pd.read_csv(
        FINE_EXPR,
        sep="\t",
        usecols=["Tissue", "Type", "Annotation", "gene_symbol", "mean_log_expression"],
    )
    fine_small = fine_small[
        fine_small["Tissue"].eq("Tumor")
        & ~fine_small["Annotation"].str.contains("lowQual|UnIdent", case=False, regex=True)
    ]
    state_means = (
        fine_small.groupby(["Annotation", "gene_symbol"], observed=True)["mean_log_expression"]
        .mean()
        .reset_index()
    )
    gene_var = state_means.groupby("gene_symbol", observed=True)["mean_log_expression"].var()
    variable_genes = [
        g for g in gene_var.sort_values(ascending=False).index if clean_gene_for_hvg(str(g))
    ][:900]
    canonical = [
        "EPCAM", "KRT8", "KRT18", "KRT19", "COL1A1", "COL1A2", "DCN", "LUM",
        "CTHRC1", "PLAT", "PECAM1", "VWF", "RGS5", "MCAM", "PTPRC", "LST1",
        "TYROBP", "FCER1G", "CD3D", "CD3E", "CD4", "CD8A", "NKG7", "GNLY",
        "MS4A1", "CD79A", "MZB1", "JCHAIN", "TPSAB1", "KIT",
    ]
    analysis_genes = sorted(
        set(lr_genes).union(variable_genes).union(canonical).intersection(gene_to_all_col)
    )
    analysis_cols = np.asarray([gene_to_all_col[g] for g in analysis_genes], dtype=int)
    X = np.zeros((len(comm_indices), len(analysis_genes)), dtype=np.float32)
    cell_ids = np.empty(len(comm_indices), dtype=object)

    fibro_patients = sorted(metadata.loc[fibro_indices, "Patient"].unique().tolist())
    fibro_patient_to_row = {p: i for i, p in enumerate(fibro_patients)}
    fibro_sum = np.zeros((len(fibro_patients), len(all_genes)), dtype=np.float64)
    fibro_detected = np.zeros((len(fibro_patients), len(all_genes)), dtype=np.int32)
    fibro_n = np.zeros(len(fibro_patients), dtype=np.int32)
    fibro_index_set = set(fibro_indices.tolist())

    row_cursor = 0
    alignment_mismatches = 0
    with gzip.open(EXPR_PATH, "rt", encoding="utf-8", newline="") as handle:
        handle.readline()
        for row_idx, line in enumerate(handle):
            if row_idx not in selected_union:
                continue
            split_at = line.find(",")
            if split_at < 0:
                raise RuntimeError(f"Malformed expression row {row_idx}")
            cell_id = line[:split_at]
            observed_sample = cell_id.rsplit("-", 1)[-1]
            if observed_sample != metadata.at[row_idx, "Sample"]:
                alignment_mismatches += 1
            values = np.fromstring(line[split_at + 1 :], sep=",", dtype=np.float32)
            if len(values) != len(all_genes) or not np.isfinite(values).all():
                raise RuntimeError(f"Invalid expression row {row_idx}: {len(values)} values")
            if row_idx in comm_position:
                pos = comm_position[row_idx]
                X[pos] = np.rint(np.expm1(values[analysis_cols])).astype(np.float32)
                cell_ids[pos] = cell_id
            if row_idx in fibro_index_set:
                p = fibro_patient_to_row[metadata.at[row_idx, "Patient"]]
                fibro_sum[p] += values
                fibro_detected[p] += values > 0
                fibro_n[p] += 1
            row_cursor += 1
    if row_cursor != len(selected_union) or alignment_mismatches:
        raise RuntimeError(
            f"Selected-row audit failed rows={row_cursor}/{len(selected_union)}, "
            f"mismatches={alignment_mismatches}"
        )

    fibro_mean_by_patient = fibro_sum / fibro_n[:, None]
    fibro_detection_by_patient = fibro_detected / fibro_n[:, None]
    fibro_metrics = pd.DataFrame(
        {
            "gene_symbol": all_genes,
            "patient_equal_mean_log_expression": fibro_mean_by_patient.mean(axis=0),
            "patient_equal_detection": fibro_detection_by_patient.mean(axis=0),
            "patients_detected": (fibro_detected > 0).sum(axis=0),
            "patients_evaluated": len(fibro_patients),
        }
    )
    fibro_metrics["expression_gate"] = (
        fibro_metrics["patient_equal_detection"].ge(EXPR_PROP)
        & fibro_metrics["patients_detected"].ge(MIN_PATIENTS)
    )
    fibro_metrics.to_csv(FIGDATA / "panel_b_tumour_fibroblast_expression_metrics.tsv.gz", sep="\t", index=False)

    membership_old = pd.read_csv(OLD_MEMBERSHIP, sep="\t")
    prior_sets = {
        name: set(membership_old.loc[membership_old[name].astype(bool), "gene_symbol"].astype(str))
        for name in ["Matrisome", "LR prior", "NicheNet"]
    }
    broad_fibro_set = set(
        fibro_metrics.loc[fibro_metrics["expression_gate"], "gene_symbol"].astype(str)
    )
    sets = {**prior_sets, "Tumour fibroblast expression": broad_fibro_set}
    union = sorted(set().union(*sets.values()))
    membership = pd.DataFrame({"gene_symbol": union})
    for name, genes in sets.items():
        membership[name] = membership["gene_symbol"].isin(genes)
    membership["membership_code"] = membership[list(sets)].astype(int).astype(str).agg("".join, axis=1)
    membership["exclusive_region"] = membership.apply(
        lambda row: " & ".join([name for name in sets if bool(row[name])]), axis=1
    )
    membership.to_csv(FIGDATA / "panel_b_four_set_membership.tsv.gz", sep="\t", index=False)
    all_codes = [f"{i:04b}" for i in range(1, 16)]
    region_counts = membership["membership_code"].value_counts().reindex(all_codes, fill_value=0)
    region_table = pd.DataFrame(
        {
            "membership_code": all_codes,
            "exclusive_gene_count": region_counts.to_numpy(dtype=int),
            "exclusive_region": [
                " & ".join([name for bit, name in zip(code, sets) if bit == "1"])
                for code in all_codes
            ],
        }
    )
    region_table.to_csv(FIGDATA / "panel_b_four_set_exclusive_regions.tsv", sep="\t", index=False)
    candidates = sorted(set.intersection(*sets.values()))

    # Receiver-unbiased ranking across all tumour fibroblasts.
    fine_candidates = pd.read_csv(FINE_EXPR, sep="\t")
    fine_candidates = fine_candidates[
        fine_candidates["Tissue"].eq("Tumor")
        & ~fine_candidates["Annotation"].str.contains("lowQual|UnIdent", case=False, regex=True)
        & fine_candidates["gene_symbol"].isin(candidates)
        & fine_candidates["n_cells"].ge(5)
    ].copy()
    fine_candidates["weighted_mean"] = (
        fine_candidates["mean_log_expression"] * fine_candidates["n_cells"]
    )
    by_patient_type = (
        fine_candidates.groupby(["Patient", "Type", "gene_symbol"], observed=True)
        .agg(weighted_sum=("weighted_mean", "sum"), cells=("n_cells", "sum"))
        .reset_index()
    )
    by_patient_type["mean_log_expression"] = by_patient_type["weighted_sum"] / by_patient_type["cells"]
    type_means = (
        by_patient_type.groupby(["Type", "gene_symbol"], observed=True)
        .agg(patient_equal_mean=("mean_log_expression", "mean"), patients=("Patient", "nunique"))
        .reset_index()
    )
    common_types = set(type_means.loc[type_means["patients"].ge(5), "Type"])
    fib_type = type_means[type_means["Type"].eq("Fibroblast")][
        ["gene_symbol", "patient_equal_mean"]
    ].rename(columns={"patient_equal_mean": "fibroblast_type_mean"})
    other_type = (
        type_means[type_means["Type"].isin(common_types - {"Fibroblast"})]
        .groupby("gene_symbol", observed=True)["patient_equal_mean"]
        .max()
        .rename("best_other_type_mean")
        .reset_index()
    )
    ranking = pd.DataFrame({"ligand": candidates}).merge(
        fibro_metrics.rename(columns={"gene_symbol": "ligand"}), on="ligand", how="left"
    ).merge(
        fib_type.rename(columns={"gene_symbol": "ligand"}), on="ligand", how="left"
    ).merge(
        other_type.rename(columns={"gene_symbol": "ligand"}), on="ligand", how="left"
    ).fillna(0)
    ranking["fibroblast_lineage_specificity_log2_ratio"] = np.log2(
        (ranking["fibroblast_type_mean"] + 0.05) / (ranking["best_other_type_mean"] + 0.05)
    )
    rank_metrics = [
        "patient_equal_mean_log_expression",
        "patient_equal_detection",
        "patients_detected",
        "fibroblast_lineage_specificity_log2_ratio",
    ]
    for metric in rank_metrics:
        ranking[f"percentile_{metric}"] = percentile(ranking[metric])
    ranking["evidence_score"] = ranking[[f"percentile_{m}" for m in rank_metrics]].mean(axis=1)

    mat = pd.read_csv(MATRISOME, sep="\t")
    class_map: dict[str, str] = {}
    for set_name, label in [
        ("NABA_COLLAGENS", "Collagen"),
        ("NABA_ECM_GLYCOPROTEINS", "ECM glycoprotein"),
        ("NABA_ECM_REGULATORS", "ECM regulator"),
        ("NABA_PROTEOGLYCANS", "Proteoglycan"),
    ]:
        for gene in mat.loc[mat["gs_name"].eq(set_name), "gene_symbol"].astype(str):
            class_map.setdefault(gene, label)
    ranking["matrisome_class"] = ranking["ligand"].map(class_map).fillna("Other matrisome")
    ranking = ranking.sort_values(
        ["evidence_score", "patient_equal_detection", "ligand"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    ranking["scan_rank"] = np.arange(1, len(ranking) + 1)
    ranking.to_csv(FIGDATA / "panel_c_receiver_unbiased_ligand_ranking.tsv", sep="\t", index=False)

    # PLAT localisation from patient-level, cell-count-weighted state summaries.
    plat = pd.read_csv(FINE_EXPR, sep="\t")
    plat = plat[
        plat["Tissue"].eq("Tumor")
        & plat["gene_symbol"].eq("PLAT")
        & plat["n_cells"].ge(5)
    ].copy()
    plat["detected_cells"] = plat["n_cells"] * plat["detected_fraction"]
    plat["weighted_mean"] = plat["n_cells"] * plat["mean_log_expression"]
    ptype = (
        plat.groupby(["Patient", "Type"], observed=True)
        .agg(cells=("n_cells", "sum"), detected=("detected_cells", "sum"), weighted=("weighted_mean", "sum"))
        .reset_index()
    )
    ptype["detection"] = ptype["detected"] / ptype["cells"]
    ptype["mean_log_expression"] = ptype["weighted"] / ptype["cells"]
    lineage = (
        ptype.groupby("Type", observed=True)
        .agg(
            mean_log_expression=("mean_log_expression", "mean"),
            detection_fraction=("detection", "mean"),
            patients=("Patient", "nunique"),
        )
        .reset_index()
        .sort_values("mean_log_expression", ascending=False)
    )
    lineage.to_csv(FIGDATA / "panel_e_PLAT_lineage_localisation.tsv", sep="\t", index=False)

    pfib = plat[plat["Type"].eq("Fibroblast")].copy()
    fib_state = (
        pfib.groupby("Annotation", observed=True)
        .agg(
            cells=("n_cells", "sum"),
            patients=("Patient", "nunique"),
            count_proxy=("reconstructed_count_sum", "sum"),
            cell_weighted_mean=("weighted_mean", "sum"),
            detected_cells=("detected_cells", "sum"),
        )
        .reset_index()
    )
    fib_state["cell_weighted_mean"] /= fib_state["cells"]
    fib_state["detection_fraction"] = fib_state["detected_cells"] / fib_state["cells"]
    fib_state["fraction_of_PLAT_count_proxy"] = fib_state["count_proxy"] / fib_state["count_proxy"].sum()
    fib_state.to_csv(FIGDATA / "panel_e_PLAT_fibroblast_state_localisation.tsv", sep="\t", index=False)

    # Build the donor-balanced all-lineage AnnData object.
    obs = metadata.loc[comm_indices].copy()
    obs.insert(0, "cell_id", cell_ids)
    obs.index = pd.Index(cell_ids.astype(str), name="cell_id_index")
    type_map = {"B": "B cells", "Innate lymphoid": "NK/ILC"}
    obs["cell_type"] = obs["Type"].replace(type_map)
    counts = sparse.csr_matrix(X)
    adata = ad.AnnData(counts, obs=obs, var=pd.DataFrame(index=pd.Index(analysis_genes, name="gene_symbol")))
    adata.layers["counts"] = counts.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=min(1200, adata.n_vars))
    sc.pp.pca(adata, mask_var="highly_variable", n_comps=30, random_state=SEED)
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30, random_state=SEED)
    sc.tl.umap(adata, min_dist=0.35, random_state=SEED)
    adata.write_h5ad(OUT / "GSE206785_donor_balanced_all_lineage.h5ad", compression="gzip")
    plat_col = analysis_genes.index("PLAT")
    plat_values = np.asarray(adata.X[:, plat_col].toarray()).ravel()
    umap = obs[["cell_id", "Patient", "cell_type", "Annotation"]].copy()
    umap["UMAP1"] = adata.obsm["X_umap"][:, 0]
    umap["UMAP2"] = adata.obsm["X_umap"][:, 1]
    umap["PLAT_log1p_CP10K"] = plat_values
    umap.to_csv(FIGDATA / "panel_d_global_umap.tsv.gz", sep="\t", index=False)

    fib_mask = obs["Type"].eq("Fibroblast").to_numpy()
    fib = adata[fib_mask].copy()
    sc.pp.highly_variable_genes(fib, flavor="seurat", n_top_genes=min(900, fib.n_vars))
    sc.pp.pca(fib, mask_var="highly_variable", n_comps=25, random_state=SEED)
    sc.pp.neighbors(fib, n_neighbors=15, n_pcs=25, random_state=SEED)
    sc.tl.umap(fib, min_dist=0.30, random_state=SEED)
    fib_plat = np.asarray(fib[:, "PLAT"].X.toarray()).ravel()
    fib_umap = fib.obs[["cell_id", "Patient", "Annotation"]].copy()
    fib_umap["UMAP1"] = fib.obsm["X_umap"][:, 0]
    fib_umap["UMAP2"] = fib.obsm["X_umap"][:, 1]
    fib_umap["PLAT_log1p_CP10K"] = fib_plat
    fib_umap.to_csv(FIGDATA / "panel_e_fibroblast_umap.tsv.gz", sep="\t", index=False)

    # LIANA consensus: magnitude and specificity remain separate columns.
    result = li.mt.rank_aggregate(
        adata,
        groupby="cell_type",
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
    result.to_csv(OUT / "GSE206785_global_LIANA_consensus.tsv.gz", sep="\t", index=False)
    global_summary = aggregate_liana(result)

    robust = result[
        result["magnitude_rank"].le(0.05)
        & result["specificity_rank"].le(0.05)
    ].copy()
    robust["joint_rank"] = robust[["magnitude_rank", "specificity_rank"]].max(axis=1)
    robust["joint_support"] = -np.log10(np.clip(robust["joint_rank"], 1e-12, 1.0))

    # Donor co-detection support for fibroblast outgoing robust pairs.
    bool_matrix = (adata.layers["counts"] > 0).tocsr()
    gene_to_col = {g: i for i, g in enumerate(analysis_genes)}
    obs_reset = adata.obs.reset_index(drop=True)
    pair_support_rows: list[dict[str, object]] = []
    fib_robust = robust[robust["source"].eq("Fibroblast")]
    for row in fib_robust.itertuples(index=False):
        eligible_donors = 0
        supported_donors = 0
        for patient in sorted(obs_reset["Patient"].unique()):
            source_rows = np.flatnonzero(
                obs_reset["Patient"].eq(patient).to_numpy()
                & obs_reset["cell_type"].eq("Fibroblast").to_numpy()
            )
            target_rows = np.flatnonzero(
                obs_reset["Patient"].eq(patient).to_numpy()
                & obs_reset["cell_type"].eq(row.target).to_numpy()
            )
            if len(source_rows) < MIN_CELLS or len(target_rows) < MIN_CELLS:
                continue
            eligible_donors += 1
            ldet = complex_detected(bool_matrix, source_rows, row.ligand_complex, gene_to_col)
            rdet = complex_detected(bool_matrix, target_rows, row.receptor_complex, gene_to_col)
            if np.isfinite(ldet) and np.isfinite(rdet) and ldet >= EXPR_PROP and rdet >= EXPR_PROP:
                supported_donors += 1
        pair_support_rows.append(
            {
                "source": row.source,
                "target": row.target,
                "ligand_complex": row.ligand_complex,
                "receptor_complex": row.receptor_complex,
                "magnitude_rank": row.magnitude_rank,
                "specificity_rank": row.specificity_rank,
                "joint_support": row.joint_support,
                "eligible_donors": eligible_donors,
                "supported_donors": supported_donors,
                "donor_support_fraction": supported_donors / eligible_donors if eligible_donors else math.nan,
            }
        )
    pair_support = pd.DataFrame(pair_support_rows)
    pair_support.to_csv(FIGDATA / "panel_g_fibroblast_outgoing_pair_support.tsv", sep="\t", index=False)
    if len(pair_support):
        donor_target = (
            pair_support.groupby("target", observed=True)
            .agg(
                median_donor_support=("donor_support_fraction", "median"),
                recurrent_pairs=("donor_support_fraction", lambda x: int(np.sum(np.asarray(x) >= 0.50))),
            )
            .reset_index()
        )
        global_summary = global_summary.merge(
            donor_target, on="target", how="left", suffixes=("", "_fibroblast")
        )
    global_summary.to_csv(FIGDATA / "panel_f_global_communication_matrix.tsv", sep="\t", index=False)
    fibro_out = global_summary[global_summary["source"].eq("Fibroblast")].copy()
    fibro_out = fibro_out.sort_values(["evidence_weight", "n_robust"], ascending=False).reset_index(drop=True)
    fibro_out["receiver_rank"] = np.arange(1, len(fibro_out) + 1)
    fibro_out.to_csv(FIGDATA / "panel_g_fibroblast_receiver_ranking.tsv", sep="\t", index=False)

    # Reuse the already audited state-balanced GSE206785 fine-state result.
    state = pd.read_csv(STATE_SUMMARY, sep="\t")
    state.to_csv(FIGDATA / "panel_h_fibroblast_myeloid_state_landscape.tsv", sep="\t", index=False)

    plat_rank = ranking.loc[ranking["ligand"].eq("PLAT")]
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": "GSE206785",
        "core_claim_boundary": (
            "predicted fibroblast-myeloid communication plus independent PLAT fibroblast localisation; "
            "no PLAT-mediated signaling or causal claim"
        ),
        "seed": SEED,
        "eligible_tumour_cells": int(valid.sum()),
        "communication_cells": int(len(comm_indices)),
        "communication_cap_per_patient_type": COMM_CAP_PER_PATIENT_TYPE,
        "fibroblast_expression_cells": int(len(fibro_indices)),
        "fibroblast_cap_per_patient": FIBRO_CAP_PER_PATIENT,
        "major_lineages": sorted(obs["cell_type"].unique().tolist()),
        "analysis_genes": len(analysis_genes),
        "candidate_set_sizes": {name: len(value) for name, value in sets.items()},
        "four_way_candidates": len(candidates),
        "plat_in_candidate_universe": "PLAT" in candidates,
        "plat_rank": int(plat_rank["scan_rank"].iloc[0]) if len(plat_rank) else None,
        "plat_score": float(plat_rank["evidence_score"].iloc[0]) if len(plat_rank) else None,
        "liana": {
            "version": li.__version__,
            "expr_prop": EXPR_PROP,
            "min_cells": MIN_CELLS,
            "n_perms": N_PERMS,
            "robust_rule": "magnitude_rank <= 0.05 AND specificity_rank <= 0.05",
        },
        "input_sha256": {
            "expression": sha256_file(EXPR_PATH),
            "metadata": sha256_file(META_PATH),
            "lr_manifest": sha256_file(LR_MANIFEST),
            "state_summary": sha256_file(STATE_SUMMARY),
        },
        "notes": [
            "All plotted cells are tumour cells with author lowQual/UnIdent labels excluded.",
            "Glial cells are excluded from the communication atlas because of sparse coverage.",
            "Candidate ranking is sender-side and does not use myeloid identity, M02, survival, or treatment response.",
            "The fine-state landscape reuses the audited GSE206785 state-balanced LIANA result.",
        ],
    }
    (AUDIT / "Figure1_GSE206785_unbiased_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    print("\nFibroblast outgoing receiver ranking:\n", fibro_out.to_string(index=False))


if __name__ == "__main__":
    main()
