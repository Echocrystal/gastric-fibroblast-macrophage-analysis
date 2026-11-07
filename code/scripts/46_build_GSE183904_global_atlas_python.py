from __future__ import annotations

import hashlib
import json
import os
import re
import sys
from pathlib import Path

PROJECT = Path(r"I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026")
PY_RUNTIME = PROJECT / "runtime" / "python_plotting"
sys.path.insert(0, str(PY_RUNTIME))

import h5py
import igraph as ig
import leidenalg
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from sklearn.utils.sparsefuncs import inplace_column_scale
import umap


SEED = 20260730
N_HVG = 3000
N_PCS = 50
N_NEIGHBORS = 30
LEIDEN_RESOLUTION = 0.75

RAW_DIR = PROJECT / "analysis" / "processed" / "GSE183904_global_atlas" / "raw_h5"
PROC_DIR = PROJECT / "analysis" / "processed" / "GSE183904_global_atlas"
RESULT_DIR = PROJECT / "analysis" / "results" / "PLAT_axis" / "GSE183904" / "global_atlas"
AUDIT_DIR = PROJECT / "analysis" / "audit"
for directory in (PROC_DIR, RESULT_DIR, AUDIT_DIR):
    directory.mkdir(parents=True, exist_ok=True)

GENES_FILE = PROC_DIR / "GSE183904_global_atlas_genes.tsv"
CELLS_FILE = PROC_DIR / "GSE183904_global_atlas_cells.tsv.gz"
EXPORT_SUMMARY_FILE = PROC_DIR / "GSE183904_global_atlas_export_summary.tsv"

AXIS_GENES = {"PLAT", "ANXA2", "ITGAM", "ITGB2", "LRP1", "ILK"}
CELL_CYCLE_GENES = {"MKI67", "TOP2A", "TYMS", "UBE2C", "CENPF"}

MARKER_SETS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Endothelial": ["PECAM1", "VWF", "EMCN", "KDR"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"],
    "Mural": ["RGS5", "CSPG4", "MCAM", "ACTA2", "TAGLN"],
    "T cell": ["CD3D", "CD3E", "TRAC", "LTB"],
    "NK": ["NKG7", "GNLY", "KLRD1"],
    "B cell": ["MS4A1", "CD79A", "CD37", "CD74"],
    "Plasma": ["MZB1", "JCHAIN", "SDC1", "IGKC"],
    "Myeloid": ["LST1", "TYROBP", "FCER1G", "CTSS"],
    "Mast": ["TPSAB1", "TPSB2", "KIT", "MS4A2"],
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_csc(h5_file: Path) -> tuple[sp.csc_matrix, list[str]]:
    with h5py.File(h5_file, "r") as handle:
        shape = tuple(int(x) for x in handle["shape_genes_by_cells"][:])
        indptr = handle["indptr"][:].astype(np.int64, copy=False)
        indices = handle["indices"][:].astype(np.int32, copy=False)
        data = handle["data"][:].astype(np.float32, copy=False)
        cell_ids = [x.decode("utf-8") if isinstance(x, bytes) else str(x)
                    for x in handle["cell_id"][:]]
    return sp.csc_matrix((data, indices, indptr), shape=shape), cell_ids


def log_cp10k(counts: sp.csc_matrix) -> sp.csc_matrix:
    lib = np.asarray(counts.sum(axis=0)).ravel().astype(np.float64)
    lib[~np.isfinite(lib) | (lib <= 0)] = 1.0
    scale = (1.0e4 / lib).astype(np.float32)
    out = counts.copy()
    inplace_column_scale(out, scale)
    np.log1p(out.data, out=out.data)
    return out


def eligible_gene(gene: str) -> bool:
    if gene in AXIS_GENES or gene in CELL_CYCLE_GENES:
        return False
    return re.match(r"^(MT-|RPL|RPS|HB)", gene) is None


def main() -> None:
    np.random.seed(SEED)
    genes = pd.read_csv(GENES_FILE, sep="\t")["gene_symbol"].astype(str).tolist()
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    cells = pd.read_csv(CELLS_FILE, sep="\t")
    sample_summary = pd.read_csv(EXPORT_SUMMARY_FILE, sep="\t")
    if len(cells) != 109788 or cells["Patient_ID"].nunique() != 26:
        raise RuntimeError("Unexpected GSE183904 primary-tumour cell universe")

    print("Pass 1/3: patient-aware gene variability", flush=True)
    variance_sum = np.zeros(len(genes), dtype=np.float64)
    variance_n = np.zeros(len(genes), dtype=np.int32)
    detected_patients = np.zeros(len(genes), dtype=np.int32)
    pooled_detected_cells = np.zeros(len(genes), dtype=np.int64)

    for row in sample_summary.itertuples(index=False):
        h5_file = RAW_DIR / f"{row.sample_id}_counts_csc.h5"
        mat, ids = load_csc(h5_file)
        expected = cells.loc[cells["sample_id"].eq(row.sample_id), "cell_id"].tolist()
        if ids != expected:
            raise RuntimeError(f"Cell-order mismatch for {row.sample_id}")
        x = log_cp10k(mat)
        n = x.shape[1]
        s1 = np.asarray(x.sum(axis=1)).ravel().astype(np.float64)
        x2 = x.copy()
        np.square(x2.data, out=x2.data)
        s2 = np.asarray(x2.sum(axis=1)).ravel().astype(np.float64)
        var = np.maximum((s2 - (s1 * s1) / n) / max(n - 1, 1), 0.0)
        ok = np.isfinite(var)
        variance_sum[ok] += var[ok]
        variance_n[ok] += 1
        nnz = np.diff(x.tocsr().indptr)
        detected_patients += (nnz > 0).astype(np.int32)
        pooled_detected_cells += nnz.astype(np.int64)
        print(f"  {row.sample_id}: {n:,} cells", flush=True)
        del mat, x, x2

    mean_within_patient_variance = variance_sum / np.maximum(variance_n, 1)
    pooled_detection = pooled_detected_cells / len(cells)
    eligible = np.array([
        eligible_gene(g) and detected_patients[i] >= 5 and pooled_detection[i] >= 0.01
        for i, g in enumerate(genes)
    ])
    eligible_idx = np.flatnonzero(eligible)
    order = eligible_idx[np.argsort(-mean_within_patient_variance[eligible_idx], kind="stable")]
    hvg_idx = order[:N_HVG]
    hvg_genes = [genes[i] for i in hvg_idx]
    if len(hvg_genes) != N_HVG:
        raise RuntimeError("Insufficient eligible genes for global atlas")
    if AXIS_GENES.intersection(hvg_genes):
        raise RuntimeError("Axis gene leaked into global clustering features")

    hvg_table = pd.DataFrame({
        "gene_symbol": genes,
        "mean_within_patient_variance": mean_within_patient_variance,
        "contributing_patients": variance_n,
        "detected_patients": detected_patients,
        "pooled_detection_fraction": pooled_detection,
        "eligible": eligible,
        "selected_HVG": np.isin(np.arange(len(genes)), hvg_idx),
    })
    hvg_table.to_csv(RESULT_DIR / "GSE183904_global_atlas_patient_aware_HVGs.tsv",
                     sep="\t", index=False)

    print("Pass 2/3: normalized HVG matrix and Python embedding", flush=True)
    sample_matrices: list[sp.csr_matrix] = []
    offsets: dict[str, tuple[int, int]] = {}
    offset = 0
    for row in sample_summary.itertuples(index=False):
        mat, _ = load_csc(RAW_DIR / f"{row.sample_id}_counts_csc.h5")
        sub = mat.tocsr()[hvg_idx, :].tocsc()
        x = log_cp10k(sub).T.tocsr()
        x.data = x.data.astype(np.float32, copy=False)
        sample_matrices.append(x)
        offsets[row.sample_id] = (offset, offset + x.shape[0])
        offset += x.shape[0]
        print(f"  {row.sample_id}: matrix {x.shape[0]:,} x {x.shape[1]:,}", flush=True)
        del mat, sub

    X = sp.vstack(sample_matrices, format="csr", dtype=np.float32)
    del sample_matrices
    if X.shape != (len(cells), N_HVG):
        raise RuntimeError(f"Unexpected HVG matrix shape: {X.shape}")

    mean = np.asarray(X.mean(axis=0)).ravel()
    mean_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    sd = np.sqrt(np.maximum(mean_sq - mean * mean, 1e-6)).astype(np.float32)
    inplace_column_scale(X, (1.0 / sd).astype(np.float32))

    svd = TruncatedSVD(n_components=N_PCS, algorithm="randomized", n_iter=7,
                       random_state=SEED)
    pcs = svd.fit_transform(X).astype(np.float32)
    del X

    patient = cells["Patient_ID"].astype(str).to_numpy()
    for patient_id in pd.unique(patient):
        idx = np.flatnonzero(patient == patient_id)
        pcs[idx] -= pcs[idx].mean(axis=0, keepdims=True)
    pc_sd = pcs.std(axis=0, ddof=1)
    pc_sd[~np.isfinite(pc_sd) | (pc_sd <= 0)] = 1.0
    pcs /= pc_sd

    reducer = umap.UMAP(
        n_neighbors=N_NEIGHBORS,
        n_components=2,
        min_dist=0.25,
        metric="euclidean",
        random_state=SEED,
        transform_seed=SEED,
        low_memory=True,
        n_jobs=1,
        verbose=True,
    )
    embedding = reducer.fit_transform(pcs).astype(np.float32)
    graph = reducer.graph_.maximum(reducer.graph_.T).tocsr()
    graph.setdiag(0)
    graph.eliminate_zeros()

    print("Leiden clustering from the Python UMAP neighbour graph", flush=True)
    ig_graph = ig.Graph.Weighted_Adjacency(
        graph, mode="undirected", attr="weight", loops=False
    )
    partition = leidenalg.find_partition(
        ig_graph,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=LEIDEN_RESOLUTION,
        seed=SEED,
        n_iterations=20,
    )
    cluster = np.asarray(partition.membership, dtype=np.int32) + 1
    del graph, ig_graph, partition

    cells["global_cluster"] = cluster
    cells["UMAP1"] = embedding[:, 0]
    cells["UMAP2"] = embedding[:, 1]
    cells.to_csv(PROC_DIR / "GSE183904_global_atlas_cell_metadata.tsv.gz",
                 sep="\t", index=False, compression="gzip")
    np.savez_compressed(
        PROC_DIR / "GSE183904_global_atlas_PCs_UMAP.npz",
        pcs=pcs, umap=embedding, cluster=cluster,
        explained_variance_ratio=svd.explained_variance_ratio_,
    )

    cross = pd.crosstab(cells["global_cluster"], cells["coarse_label"])
    cross_fraction = cross.div(cross.sum(axis=1), axis=0)
    cluster_summary = pd.DataFrame({
        "global_cluster": cross.index,
        "n_cells": cross.sum(axis=1).to_numpy(),
        "n_patients": cells.groupby("global_cluster")["Patient_ID"].nunique().reindex(cross.index).to_numpy(),
        "majority_coarse_label": cross.idxmax(axis=1).to_numpy(),
        "majority_fraction": cross_fraction.max(axis=1).to_numpy(),
    })
    cluster_summary.to_csv(RESULT_DIR / "GSE183904_global_atlas_cluster_summary.tsv",
                           sep="\t", index=False)
    cross.reset_index().to_csv(
        RESULT_DIR / "GSE183904_global_atlas_cluster_by_coarse_label_counts.tsv",
        sep="\t", index=False
    )

    print("Pass 3/3: canonical-marker summaries by global cluster", flush=True)
    marker_genes = list(dict.fromkeys(g for values in MARKER_SETS.values() for g in values))
    missing_markers = [g for g in marker_genes if g not in gene_to_idx]
    if missing_markers:
        raise RuntimeError(f"Missing canonical markers: {missing_markers}")
    marker_idx = np.array([gene_to_idx[g] for g in marker_genes], dtype=np.int64)
    clusters = np.sort(np.unique(cluster))
    cluster_to_col = {int(c): i for i, c in enumerate(clusters)}
    sum_expr = np.zeros((len(marker_genes), len(clusters)), dtype=np.float64)
    detected = np.zeros_like(sum_expr)
    n_cells_cluster = np.zeros(len(clusters), dtype=np.int64)

    for row in sample_summary.itertuples(index=False):
        start, end = offsets[row.sample_id]
        sample_cluster = cluster[start:end]
        mat, _ = load_csc(RAW_DIR / f"{row.sample_id}_counts_csc.h5")
        x = log_cp10k(mat.tocsr()[marker_idx, :].tocsc()).tocsr()
        for c in np.unique(sample_cluster):
            cell_idx = np.flatnonzero(sample_cluster == c)
            j = cluster_to_col[int(c)]
            block = x[:, cell_idx]
            sum_expr[:, j] += np.asarray(block.sum(axis=1)).ravel()
            detected[:, j] += np.diff(block.indptr)
            n_cells_cluster[j] += len(cell_idx)
        print(f"  markers: {row.sample_id}", flush=True)
        del mat, x

    records = []
    marker_to_lineage = {
        gene: lineage for lineage, members in MARKER_SETS.items() for gene in members
    }
    for i, gene in enumerate(marker_genes):
        for j, c in enumerate(clusters):
            records.append({
                "global_cluster": int(c),
                "lineage": marker_to_lineage[gene],
                "gene_symbol": gene,
                "mean_log1p_CP10K": sum_expr[i, j] / n_cells_cluster[j],
                "detected_fraction": detected[i, j] / n_cells_cluster[j],
                "n_cells": int(n_cells_cluster[j]),
            })
    marker_table = pd.DataFrame.from_records(records)
    marker_table.to_csv(RESULT_DIR / "GSE183904_global_atlas_marker_summary.tsv",
                        sep="\t", index=False)

    audit = {
        "generated_at": pd.Timestamp.now(tz="Asia/Shanghai").isoformat(),
        "dataset": "GSE183904",
        "primary_tumours": int(cells["sample_id"].nunique()),
        "patients": int(cells["Patient_ID"].nunique()),
        "cells": int(len(cells)),
        "genes": int(len(genes)),
        "global_clusters": int(len(clusters)),
        "backend": "Python 3.12",
        "normalization": "log1p(CP10K)",
        "feature_selection": {
            "strategy": "mean within-patient variance",
            "hvg_count": N_HVG,
            "minimum_detection_fraction": 0.01,
            "minimum_detected_patients": 5,
            "excluded_axis_genes": sorted(AXIS_GENES),
            "excluded_cell_cycle_genes": sorted(CELL_CYCLE_GENES),
            "excluded_prefixes": ["MT-", "RPL", "RPS", "HB"],
        },
        "dimension_reduction": {
            "scaled_sparse_truncated_svd_components": N_PCS,
            "patient_centered_PC_space": True,
            "umap_neighbors": N_NEIGHBORS,
            "umap_min_dist": 0.25,
            "seed": SEED,
        },
        "clustering": {
            "method": "Leiden RBConfiguration on UMAP fuzzy neighbour graph",
            "resolution": LEIDEN_RESOLUTION,
            "iterations": 20,
        },
        "annotation_note": (
            "Broad labels were prespecified marker-score annotations and were not used "
            "to construct the unsupervised graph, clusters, or UMAP."
        ),
        "inference_note": "Cell-level clusters and UMAP are descriptive; inferential analyses use patients.",
        "input_sha256": {
            "genes": sha256(GENES_FILE),
            "cells": sha256(CELLS_FILE),
            "export_summary": sha256(EXPORT_SUMMARY_FILE),
        },
    }
    with (AUDIT_DIR / "GSE183904_global_atlas_python.json").open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, ensure_ascii=False)

    print(
        f"Global atlas complete: {len(cells):,} cells, {len(clusters)} clusters; "
        f"PLAT excluded from {N_HVG} clustering features.",
        flush=True,
    )


if __name__ == "__main__":
    main()
