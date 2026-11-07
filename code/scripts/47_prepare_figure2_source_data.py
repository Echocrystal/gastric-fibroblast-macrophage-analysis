from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

PROJECT = Path(r"I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026")
sys.path.insert(0, str(PROJECT / "runtime" / "python_plotting"))

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.utils.sparsefuncs import inplace_column_scale


PROC = PROJECT / "analysis" / "processed"
RES = PROJECT / "analysis" / "results"
RAW_H5 = PROC / "GSE183904_global_atlas" / "raw_h5"
OUT = PROJECT / "figures" / "source_data" / "Figure2"
OUT.mkdir(parents=True, exist_ok=True)


FIB_LABELS = {
    1: "CXCL14+POSTN+PLAT-rich fibroblast",
    2: "Mural/pericyte contamination",
    3: "Immune contamination",
    4: "Activated matrix fibroblast (donor-skewed)",
    5: "Epithelial/mesothelial contamination",
    6: "PI16+ adventitial fibroblast (donor-restricted)",
}

FIB_MARKERS = [
    "CXCL14", "POSTN", "PLAT", "CFD", "MFAP4", "MMP2", "TCF21", "PDGFRA",
    "CTHRC1", "FAP", "COL1A1", "RGS5", "ACTA2", "SRGN", "CCL5", "KRT19", "PI16",
]

GLOBAL_MARKER_SETS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT19"],
    "Endothelial": ["PECAM1", "VWF", "KDR"],
    "Fibroblast": ["COL1A1", "DCN", "LUM"],
    "Mural": ["RGS5", "ACTA2", "TAGLN"],
    "T_cell": ["CD3D", "CD3E", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1"],
    "B_cell": ["MS4A1", "CD79A", "CD74"],
    "Plasma": ["MZB1", "JCHAIN", "SDC1"],
    "Myeloid": ["LST1", "TYROBP", "FCER1G"],
    "Mast": ["TPSAB1", "KIT", "MS4A2"],
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_csc(path: Path) -> tuple[sp.csc_matrix, list[str]]:
    with h5py.File(path, "r") as h:
        shape = tuple(int(x) for x in h["shape_genes_by_cells"][:])
        mat = sp.csc_matrix(
            (
                h["data"][:].astype(np.float32, copy=False),
                h["indices"][:].astype(np.int32, copy=False),
                h["indptr"][:].astype(np.int64, copy=False),
            ),
            shape=shape,
        )
        cells = [x.decode("utf-8") if isinstance(x, bytes) else str(x)
                 for x in h["cell_id"][:]]
    return mat, cells


def log_cp10k(counts: sp.csc_matrix, library_size: np.ndarray | None = None) -> sp.csc_matrix:
    if library_size is None:
        library_size = np.asarray(counts.sum(axis=0)).ravel()
    library_size = np.asarray(library_size, dtype=float)
    library_size[~np.isfinite(library_size) | (library_size <= 0)] = 1.0
    out = counts.copy()
    inplace_column_scale(out, (1.0e4 / library_size).astype(np.float32))
    np.log1p(out.data, out=out.data)
    return out


def zscore(values: pd.Series) -> pd.Series:
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        return pd.Series(np.zeros(len(values)), index=values.index)
    return (values - values.mean()) / sd


def main() -> None:
    written: list[Path] = []

    global_meta_file = PROC / "GSE183904_global_atlas" / "GSE183904_global_atlas_cell_metadata.tsv.gz"
    global_meta = pd.read_csv(global_meta_file, sep="\t")
    ab_file = OUT / "Figure2_panels_ab_global_atlas_cells.tsv.gz"
    global_meta[[
        "cell_id", "sample_id", "Patient_ID", "coarse_label", "global_cluster", "UMAP1", "UMAP2"
    ]].to_csv(ab_file, sep="\t", index=False, compression="gzip")
    written.append(ab_file)

    genes_file = PROC / "GSE183904_global_atlas" / "GSE183904_global_atlas_genes.tsv"
    genes = pd.read_csv(genes_file, sep="\t")["gene_symbol"].astype(str).tolist()
    gene_to_idx = {g: i for i, g in enumerate(genes)}

    # Conventional broad-lineage marker dot plot. Marker genes are ordered by
    # their expected lineage, while expression is measured in every broad cell type.
    global_marker_genes = list(dict.fromkeys(
        gene for members in GLOBAL_MARKER_SETS.values() for gene in members
    ))
    missing_global = [g for g in global_marker_genes if g not in gene_to_idx]
    if missing_global:
        raise RuntimeError(f"Missing global canonical markers: {missing_global}")
    global_marker_idx = np.array([gene_to_idx[g] for g in global_marker_genes], dtype=int)
    broad_labels = list(GLOBAL_MARKER_SETS)
    broad_to_col = {label: i for i, label in enumerate(broad_labels)}
    global_sums = np.zeros((len(global_marker_genes), len(broad_labels)), dtype=float)
    global_detections = np.zeros_like(global_sums)
    global_cells = np.zeros(len(broad_labels), dtype=int)

    for sample_id, atlas_sample in global_meta.groupby("sample_id", sort=False):
        mat, cell_ids = load_csc(RAW_H5 / f"{sample_id}_counts_csc.h5")
        if cell_ids != atlas_sample["cell_id"].tolist():
            raise RuntimeError(f"Global atlas cell order mismatch for {sample_id}")
        library_size = np.asarray(mat.sum(axis=0)).ravel()
        marker_mat = log_cp10k(mat.tocsr()[global_marker_idx, :].tocsc(), library_size).tocsr()
        labels = atlas_sample["coarse_label"].astype(str).to_numpy()
        for label in broad_labels:
            mask = np.flatnonzero(labels == label)
            if not len(mask):
                continue
            block = marker_mat[:, mask]
            j = broad_to_col[label]
            global_sums[:, j] += np.asarray(block.sum(axis=1)).ravel()
            global_detections[:, j] += np.diff(block.indptr)
            global_cells[j] += len(mask)

    global_records = []
    marker_group = {
        gene: lineage for lineage, members in GLOBAL_MARKER_SETS.items() for gene in members
    }
    for i, gene in enumerate(global_marker_genes):
        means = global_sums[i, :] / global_cells
        zs = (means - means.mean()) / max(means.std(ddof=0), 1e-8)
        for j, cell_type in enumerate(broad_labels):
            global_records.append({
                "cell_type": cell_type,
                "gene_symbol": gene,
                "marker_group": marker_group[gene],
                "mean_log1p_CP10K": means[j],
                "scaled_mean": zs[j],
                "detected_fraction": global_detections[i, j] / global_cells[j],
                "n_cells": global_cells[j],
            })
    broad_markers = pd.DataFrame.from_records(global_records)
    c_file = OUT / "Figure2_panel_c_broad_lineage_canonical_markers.tsv"
    broad_markers.to_csv(c_file, sep="\t", index=False)
    written.append(c_file)

    fib_file = (
        PROC / "GSE183904_fine_states" /
        "GSE183904_Fibroblast_fine_state_cell_metadata.tsv.gz"
    )
    cxcl_file = (
        PROC / "GSE183904_fine_states" /
        "GSE183904_CXCL14_APOE_cell_features.tsv.gz"
    )
    fib = pd.read_csv(fib_file, sep="\t")
    cxcl = pd.read_csv(cxcl_file, sep="\t", usecols=["cell_id", "expr_CXCL14", "detect_CXCL14"])
    fib = fib.merge(cxcl, on="cell_id", how="left", validate="one_to_one")
    if fib["expr_CXCL14"].isna().any():
        raise RuntimeError("Missing CXCL14 values after fibroblast merge")
    fib["refined_label"] = fib["cluster"].map(FIB_LABELS)
    dg_file = OUT / "Figure2_panels_dg_fibroblast_cells.tsv.gz"
    fib[[
        "cell_id", "sample_id", "Patient_ID", "cluster", "refined_label",
        "recurrent_cluster", "UMAP1", "UMAP2", "expr_PLAT", "expr_CXCL14", "expr_POSTN",
        "detect_PLAT", "detect_CXCL14", "detect_POSTN",
    ]].to_csv(dg_file, sep="\t", index=False, compression="gzip")
    written.append(dg_file)

    missing = [g for g in FIB_MARKERS if g not in gene_to_idx]
    if missing:
        raise RuntimeError(f"Missing fibroblast markers: {missing}")
    marker_idx = np.array([gene_to_idx[g] for g in FIB_MARKERS], dtype=int)
    sums = np.zeros((len(FIB_MARKERS), 6), dtype=float)
    detections = np.zeros_like(sums)
    state_cells = np.zeros(6, dtype=int)

    for sample_id, fib_sample in fib.groupby("sample_id", sort=False):
        mat, cell_ids = load_csc(RAW_H5 / f"{sample_id}_counts_csc.h5")
        lookup = {cell_id: i for i, cell_id in enumerate(cell_ids)}
        positions = np.array([lookup[x] for x in fib_sample["cell_id"]], dtype=int)
        library_size = np.asarray(mat.sum(axis=0)).ravel()
        marker_mat = log_cp10k(mat.tocsr()[marker_idx, :].tocsc(), library_size)
        marker_mat = marker_mat[:, positions].tocsr()
        clusters = fib_sample["cluster"].to_numpy(dtype=int)
        for cluster in range(1, 7):
            mask = np.flatnonzero(clusters == cluster)
            if not len(mask):
                continue
            block = marker_mat[:, mask]
            sums[:, cluster - 1] += np.asarray(block.sum(axis=1)).ravel()
            detections[:, cluster - 1] += np.diff(block.indptr)
            state_cells[cluster - 1] += len(mask)

    records = []
    for i, gene in enumerate(FIB_MARKERS):
        raw_means = sums[i, :] / state_cells
        zs = (raw_means - raw_means.mean()) / max(raw_means.std(ddof=0), 1e-8)
        for cluster in range(1, 7):
            records.append({
                "cluster": cluster,
                "refined_label": FIB_LABELS[cluster],
                "gene_symbol": gene,
                "mean_log1p_CP10K": raw_means[cluster - 1],
                "scaled_mean": zs[cluster - 1],
                "detected_fraction": detections[i, cluster - 1] / state_cells[cluster - 1],
                "n_cells": state_cells[cluster - 1],
            })
    h = pd.DataFrame.from_records(records)
    h_file = OUT / "Figure2_panel_h_fibroblast_state_markers.tsv"
    h.to_csv(h_file, sep="\t", index=False)
    written.append(h_file)

    fib_summary_file = (
        RES / "PLAT_axis" / "GSE183904" / "fine_states" /
        "GSE183904_Fibroblast_fine_state_cluster_summary.tsv"
    )
    recurrence = pd.read_csv(fib_summary_file, sep="\t")
    recurrence["refined_label"] = recurrence["cluster"].map(FIB_LABELS)
    i_file = OUT / "Figure2_panel_i_fibroblast_state_recurrence.tsv"
    recurrence[[
        "cluster", "refined_label", "n_cells", "n_patients",
        "patients_with_at_least_5_cells", "max_patient_fraction", "recurrent_cluster",
        "PLAT_mean_log1p_CP10K", "PLAT_detected_fraction",
    ]].to_csv(i_file, sep="\t", index=False)
    written.append(i_file)

    patient_state = fib.groupby(["Patient_ID", "cluster"], as_index=False).size().rename(
        columns={"size": "n_cells"}
    )
    patient_universe = pd.MultiIndex.from_product(
        [sorted(fib["Patient_ID"].unique()), range(1, 7)],
        names=["Patient_ID", "cluster"],
    ).to_frame(index=False)
    patient_state = patient_universe.merge(
        patient_state, on=["Patient_ID", "cluster"], how="left"
    )
    patient_state["n_cells"] = patient_state["n_cells"].fillna(0).astype(int)
    patient_totals = fib.groupby("Patient_ID").size().rename("patient_fibroblast_cells")
    patient_state = patient_state.merge(patient_totals, on="Patient_ID", how="left")
    patient_state["within_patient_fraction"] = (
        patient_state["n_cells"] / patient_state["patient_fibroblast_cells"]
    )
    patient_state["refined_label"] = patient_state["cluster"].map(FIB_LABELS)
    i_patient_file = OUT / "Figure2_panel_i_fibroblast_state_by_patient.tsv"
    patient_state.to_csv(i_patient_file, sep="\t", index=False)
    written.append(i_patient_file)

    specificity_dir = RES / "axis_scan" / "GSE206785" / "state_specificity"
    selected_states = [
        "F08_Activated_VSTM2A", "F13_Activated_CTHRC1", "F14_Activated_ADAM28",
        "F09_Quiescent_IGF1", "F06_Quiescent_",
    ]
    pb_file = specificity_dir / "GSE206785_ECM_ligand_patient_state_pseudobulk.tsv.gz"
    pb_all = pd.read_csv(pb_file, sep="\t")
    j = pb_all[
        pb_all["ligand"].eq("PLAT") & pb_all["Annotation"].isin(selected_states)
    ].copy()
    j["Annotation"] = pd.Categorical(j["Annotation"], selected_states, ordered=True)
    j = j.sort_values(["Annotation", "Patient"])
    j_file = OUT / "Figure2_panel_j_GSE206785_PLAT_state_localization.tsv"
    j.to_csv(j_file, sep="\t", index=False)
    written.append(j_file)

    pb = pb_all.copy()
    pb = pb[pb["ligand"].eq("PLAT")].copy()
    f13 = pb[pb["Annotation"].eq("F13_Activated_CTHRC1")][
        ["Patient", "patient_state_mean_raw", "cells"]
    ].rename(columns={"patient_state_mean_raw": "F13_mean_raw", "cells": "F13_cells"})
    other = pb[~pb["Annotation"].eq("F13_Activated_CTHRC1")].copy()
    other["weighted_expression"] = other["patient_state_mean_raw"] * other["cells"]
    other = other.groupby("Patient", as_index=False).agg(
        weighted_expression=("weighted_expression", "sum"),
        other_cells=("cells", "sum"),
    )
    other["other_fibroblast_mean_raw"] = other["weighted_expression"] / other["other_cells"]
    paired = f13.merge(other, on="Patient", how="inner", validate="one_to_one")
    paired["F13_minus_other"] = paired["F13_mean_raw"] - paired["other_fibroblast_mean_raw"]
    if len(paired) != 17:
        raise RuntimeError(f"Expected 17 paired GSE206785 patients, observed {len(paired)}")
    specificity = pd.read_csv(
        specificity_dir / "GSE206785_ECM_ligand_F13_state_specificity.tsv", sep="\t"
    )
    plat_specificity = specificity.loc[specificity["ligand"].eq("PLAT")]
    if len(plat_specificity) != 1:
        raise RuntimeError("Expected exactly one PLAT row in GSE206785 specificity results")
    stats = plat_specificity.iloc[0]
    paired["wilcoxon_p"] = float(stats["paired_patient_mean_expression_wilcoxon_p"])
    paired["fdr_BH"] = float(stats["paired_patient_mean_expression_fdr_BH"])
    paired["median_difference"] = float(stats["paired_patient_median_raw_difference"])
    paired["positive_pairs"] = int(stats["paired_patients_F13_mean_higher"])
    paired["paired_patients"] = int(stats["paired_patients"])
    k_file = OUT / "Figure2_panel_k_GSE206785_patient_paired_PLAT.tsv"
    paired.to_csv(k_file, sep="\t", index=False)
    written.append(k_file)

    manifest = {
        "generated_at": pd.Timestamp.now(tz="Asia/Shanghai").isoformat(),
        "figure": "Figure 2",
        "files": [
            {"name": path.name, "sha256": sha256(path), "bytes": path.stat().st_size}
            for path in written
        ],
        "coordinate_identity": (
            "Panels a and b are generated from the same rows and UMAP1/UMAP2 columns "
            "in Figure2_panels_ab_global_atlas_cells.tsv.gz."
        ),
        "clustering_exclusion": "PLAT and the prespecified PLAT-axis genes were excluded from HVGs.",
    }
    manifest_file = OUT / "Figure2_source_data_manifest.json"
    with manifest_file.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
    print(f"Prepared {len(written)} Figure 2 source-data tables in {OUT}")


if __name__ == "__main__":
    main()
