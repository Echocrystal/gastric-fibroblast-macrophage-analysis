#!/usr/bin/env python3
"""Generate donor-aware cell-communication source data for gastric-cancer Fig. 1.

The analysis is deliberately consensus-first.  It reconstructs exact integer
counts from the author-provided log1p matrix, normalizes with full-library sizes,
and evaluates every expression-qualified fibroblast-to-myeloid state pair.  The
pooled LIANA results are descriptive communication proxies; patient-level
ligand/receptor co-detection and a state-balanced downsampling run provide the
biological-replicate and abundance-sensitivity layers used in the figure.

PLAT is never forced into the result.  Direct PLAT edges are reported only when
they occur in the queried resource; resource membership is audited separately.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import liana as li
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
SOURCE = Path(
    r"I:\ECM-new\ECM_DEG_WGCNA_CoxLASSO_public_2026\ecm_immune_ecotype_public_2026"
)
EXPR = SOURCE / "data/raw/scrna/GSE206785/GSE206785_scgex.txt.gz"
META = SOURCE / "data/raw/scrna/GSE206785/GSE206785_metadata.txt.gz"
FROZEN_LR = PROJECT / "resources/ligand_receptor/lr_network_human_21122021_frozen.tsv"

OUT = PROJECT / "analysis/results/cell_communication/Figure1_GSE206785"
PROCESSED = PROJECT / "analysis/processed/cell_communication"
AUDIT = PROJECT / "analysis/audit"
for directory in (OUT, PROCESSED, AUDIT):
    directory.mkdir(parents=True, exist_ok=True)

SEED = 20260801
EXPR_PROP = 0.10
MIN_STATE_CELLS = 20
MIN_STATE_PATIENTS = 3
MIN_DONOR_STATE_CELLS = 5
N_PERMS = int(os.environ.get("CCC_N_PERMS", "1000"))
DOWNSAMPLE_N_PERMS = int(os.environ.get("CCC_DOWNSAMPLE_N_PERMS", "250"))
DOWNSAMPLE_CELLS_PER_STATE = int(os.environ.get("CCC_DOWNSAMPLE_CELLS", "150"))
N_JOBS = max(1, min(int(os.environ.get("CCC_N_JOBS", "4")), os.cpu_count() or 1))

F13 = "F13_Activated_CTHRC1"
M02 = "M02_Mono_CD14"
H5AD = PROCESSED / "GSE206785_tumor_fibroblast_myeloid_liana_input.h5ad"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def split_complex(symbol: str) -> list[str]:
    """Split LIANA complex labels into their gene subunits."""
    return [part for part in str(symbol).split("_") if part]


def select_resource_genes(resource: pd.DataFrame) -> set[str]:
    genes: set[str] = set()
    for column in ("ligand", "receptor"):
        for symbol in resource[column].astype(str):
            genes.update(split_complex(symbol))
    return genes


def build_input() -> tuple[ad.AnnData, pd.DataFrame, dict[str, pd.DataFrame]]:
    resources = {
        name: li.resource.select_resource(name).drop_duplicates().reset_index(drop=True)
        for name in ("consensus", "cellphonedb", "cellchatdb")
    }
    frozen = pd.read_csv(FROZEN_LR, sep="\t", usecols=["ligand", "receptor"])
    resources["frozen_omnipath"] = frozen.drop_duplicates().reset_index(drop=True)

    metadata = pd.read_csv(META, sep="\t", dtype=str)
    expected = ["Sample", "Patient", "Tissue", "Platform", "Subtype", "Type", "Annotation"]
    if metadata.columns.tolist() != expected:
        raise RuntimeError(f"Unexpected metadata columns: {metadata.columns.tolist()}")

    selected = metadata[
        metadata["Tissue"].eq("Tumor")
        & metadata["Type"].isin(["Fibroblast", "Myeloid"])
        & ~metadata["Annotation"].str.contains("lowQual", case=False, na=False)
    ].copy()
    state_summary = (
        selected.groupby(["Type", "Annotation"], observed=True)
        .agg(n_cells=("Annotation", "size"), n_patients=("Patient", "nunique"))
        .reset_index()
    )
    kept_states = state_summary.loc[
        state_summary["n_cells"].ge(MIN_STATE_CELLS)
        & state_summary["n_patients"].ge(MIN_STATE_PATIENTS),
        "Annotation",
    ].tolist()
    if F13 not in kept_states or M02 not in kept_states:
        raise RuntimeError("Frozen F13 and M02 states did not pass the prespecified coverage gate")
    selected = selected[selected["Annotation"].isin(kept_states)].copy()
    selected_rows = set(selected.index.astype(int))

    gene_universe = {"PLAT"}
    for resource in resources.values():
        gene_universe.update(select_resource_genes(resource))

    if H5AD.exists():
        adata = ad.read_h5ad(H5AD)
        if set(adata.obs_names) != set(selected.index.astype(str)):
            raise RuntimeError("Cached AnnData cell set does not match the frozen metadata gate")
        return adata, state_summary, resources

    rows: list[sparse.csr_matrix] = []
    obs_records: list[dict[str, object]] = []
    library_sizes: list[float] = []
    selected_gene_names: list[str] | None = None
    selected_gene_indices: np.ndarray | None = None
    row_count = 0
    residual_max = 0.0

    with gzip.open(EXPR, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split(",")
        if not header or header[0] != "Cell":
            raise RuntimeError("Expression header does not begin with Cell")
        all_genes = np.asarray(header[1:], dtype=object)
        if len(set(all_genes.tolist())) != len(all_genes):
            raise RuntimeError("Expression matrix contains duplicated gene symbols")
        selected_gene_indices = np.asarray(
            [i for i, gene in enumerate(all_genes) if gene in gene_universe], dtype=int
        )
        selected_gene_names = all_genes[selected_gene_indices].astype(str).tolist()
        if "PLAT" not in selected_gene_names:
            raise RuntimeError("PLAT is absent from the extracted gene universe")

        for row_idx, line in enumerate(handle):
            row_count += 1
            if row_idx not in selected_rows:
                continue
            line = line.rstrip("\r\n")
            first_comma = line.find(",")
            if first_comma <= 0:
                raise RuntimeError(f"Malformed expression row {row_idx}")
            cell_id = line[:first_comma]
            values = np.fromstring(line[first_comma + 1 :], dtype=np.float64, sep=",")
            if values.size != len(all_genes):
                raise RuntimeError(f"Gene count mismatch at expression row {row_idx}")
            counts = np.rint(np.expm1(values))
            residual = np.abs(values - np.log1p(counts))
            if np.any(counts > 0):
                residual_max = max(residual_max, float(residual[counts > 0].max()))
            libsize = float(counts.sum())
            if not math.isfinite(libsize) or libsize <= 0:
                raise RuntimeError(f"Invalid library size at expression row {row_idx}")
            selected_counts = counts[selected_gene_indices].astype(np.float32, copy=False)
            rows.append(sparse.csr_matrix(selected_counts.reshape(1, -1)))
            record = metadata.loc[row_idx].to_dict()
            record["cell_id"] = cell_id
            record["expression_row_index_zero_based"] = row_idx
            record["metadata_row_index"] = str(row_idx)
            obs_records.append(record)
            library_sizes.append(libsize)

    if row_count != len(metadata):
        raise RuntimeError(f"Expression/metadata row mismatch: {row_count} vs {len(metadata)}")
    if selected_gene_names is None or selected_gene_indices is None:
        raise RuntimeError("Expression matrix header was not parsed")
    if len(rows) != len(selected):
        raise RuntimeError(f"Extracted {len(rows)} cells but expected {len(selected)}")

    counts = sparse.vstack(rows, format="csr", dtype=np.float32)
    scale = np.asarray(1e4 / np.asarray(library_sizes, dtype=np.float64), dtype=np.float32)
    normalized = counts.multiply(scale[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data)

    obs = pd.DataFrame(obs_records).set_index("metadata_row_index", drop=True)
    var = pd.DataFrame(index=pd.Index(selected_gene_names, name="gene_symbol"))
    adata = ad.AnnData(X=normalized, obs=obs, var=var)
    adata.layers["counts_lr_subset"] = counts
    adata.obs["full_library_size"] = np.asarray(library_sizes, dtype=np.float64)
    adata.uns["extraction"] = {
        "source_semantics": "author-provided log1p integer counts",
        "count_reconstruction": "round(expm1(value))",
        "normalization": "log1p(1e4 * selected_gene_count / full_library_count)",
        "max_nonzero_reconstruction_residual": residual_max,
        "min_state_cells": MIN_STATE_CELLS,
        "min_state_patients": MIN_STATE_PATIENTS,
    }
    adata.write_h5ad(H5AD, compression="gzip")
    return adata, state_summary, resources


def source_target_pairs(adata: ad.AnnData) -> tuple[pd.DataFrame, list[str], list[str]]:
    fib = sorted(adata.obs.loc[adata.obs["Type"].eq("Fibroblast"), "Annotation"].unique())
    mye = sorted(adata.obs.loc[adata.obs["Type"].eq("Myeloid"), "Annotation"].unique())
    pairs = pd.DataFrame([(s, t) for s in fib for t in mye], columns=["source", "target"])
    return pairs, fib, mye


def run_liana(
    adata: ad.AnnData,
    groupby_pairs: pd.DataFrame,
    resource_name: str,
    resource: pd.DataFrame | None,
    n_perms: int,
    key: str,
) -> pd.DataFrame:
    kwargs: dict[str, object] = {
        "adata": adata,
        "groupby": "Annotation",
        "expr_prop": EXPR_PROP,
        "min_cells": MIN_STATE_CELLS,
        "groupby_pairs": groupby_pairs,
        "use_raw": False,
        "n_perms": n_perms,
        "seed": SEED,
        "n_jobs": N_JOBS,
        "return_all_lrs": False,
        "inplace": False,
        "verbose": True,
    }
    if resource is None:
        kwargs["resource_name"] = resource_name
    else:
        kwargs["resource"] = resource[["ligand", "receptor"]].drop_duplicates()
    result = li.mt.rank_aggregate(**kwargs)
    if result is None or result.empty:
        raise RuntimeError(f"LIANA returned no results for {key}")
    result = result.copy()
    result.insert(0, "resource", key)
    result["combined_worst_rank"] = result[["magnitude_rank", "specificity_rank"]].max(axis=1)
    result["combined_geometric_rank"] = np.sqrt(
        result["magnitude_rank"].clip(lower=np.finfo(float).tiny)
        * result["specificity_rank"].clip(lower=np.finfo(float).tiny)
    )
    result.to_csv(OUT / f"GSE206785_LIANA_{key}_full.tsv.gz", sep="\t", index=False)
    return result


def summarize_source_target(result: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for (source, target), block in result.groupby(["source", "target"], observed=True):
        block = block.sort_values(["combined_worst_rank", "combined_geometric_rank"])
        top = block.head(10)
        magnitude_strength = -np.log10(top["magnitude_rank"].clip(lower=1e-12))
        specificity_strength = -np.log10(top["specificity_rank"].clip(lower=1e-12))
        records.append(
            {
                "source": source,
                "target": target,
                "n_expression_qualified_interactions": len(block),
                "n_dual_rank_le_0_05": int(
                    ((block["magnitude_rank"] <= 0.05) & (block["specificity_rank"] <= 0.05)).sum()
                ),
                "best_combined_worst_rank": float(block["combined_worst_rank"].min()),
                "top10_mean_magnitude_strength": float(magnitude_strength.mean()),
                "top10_mean_specificity_strength": float(specificity_strength.mean()),
                "top10_joint_strength": float(np.sqrt(magnitude_strength.mean() * specificity_strength.mean())),
            }
        )
    return pd.DataFrame(records)


def component_score(values: pd.Series, complex_name: str) -> float:
    subunits = split_complex(complex_name)
    if any(gene not in values.index for gene in subunits):
        return float("nan")
    return float(min(values[gene] for gene in subunits))


def donor_support(adata: ad.AnnData, interactions: pd.DataFrame) -> pd.DataFrame:
    relevant_genes: set[str] = set()
    for col in ("ligand_complex", "receptor_complex"):
        for symbol in interactions[col].astype(str):
            relevant_genes.update(split_complex(symbol))
    relevant_genes = relevant_genes.intersection(adata.var_names)
    gene_list = sorted(relevant_genes)
    gene_index = adata.var_names.get_indexer(gene_list)

    records: list[dict[str, object]] = []
    patients = sorted(adata.obs["Patient"].unique())
    for patient in patients:
        source_mask = (
            adata.obs["Patient"].eq(patient).to_numpy()
            & adata.obs["Annotation"].eq(F13).to_numpy()
        )
        target_mask = (
            adata.obs["Patient"].eq(patient).to_numpy()
            & adata.obs["Annotation"].eq(M02).to_numpy()
        )
        n_source = int(source_mask.sum())
        n_target = int(target_mask.sum())
        if n_source < MIN_DONOR_STATE_CELLS or n_target < MIN_DONOR_STATE_CELLS:
            continue
        source_x = adata[source_mask, gene_index].X
        target_x = adata[target_mask, gene_index].X
        source_mean = pd.Series(np.asarray(source_x.mean(axis=0)).ravel(), index=gene_list)
        target_mean = pd.Series(np.asarray(target_x.mean(axis=0)).ravel(), index=gene_list)
        source_detect = pd.Series(np.asarray((source_x > 0).mean(axis=0)).ravel(), index=gene_list)
        target_detect = pd.Series(np.asarray((target_x > 0).mean(axis=0)).ravel(), index=gene_list)

        for row in interactions.itertuples(index=False):
            ligand_mean = component_score(source_mean, row.ligand_complex)
            receptor_mean = component_score(target_mean, row.receptor_complex)
            ligand_detect = component_score(source_detect, row.ligand_complex)
            receptor_detect = component_score(target_detect, row.receptor_complex)
            magnitude = math.sqrt(max(ligand_mean, 0) * max(receptor_mean, 0))
            records.append(
                {
                    "Patient": patient,
                    "source": F13,
                    "target": M02,
                    "ligand_complex": row.ligand_complex,
                    "receptor_complex": row.receptor_complex,
                    "n_source_cells": n_source,
                    "n_target_cells": n_target,
                    "ligand_mean_log_normalized": ligand_mean,
                    "receptor_mean_log_normalized": receptor_mean,
                    "ligand_detected_fraction": ligand_detect,
                    "receptor_detected_fraction": receptor_detect,
                    "coexpression_magnitude": magnitude,
                    "passes_10pct_codetection": bool(
                        ligand_detect >= EXPR_PROP and receptor_detect >= EXPR_PROP
                    ),
                }
            )
    donor = pd.DataFrame(records)
    donor.to_csv(OUT / "GSE206785_F13_M02_interaction_donor_support_long.tsv.gz", sep="\t", index=False)
    summary = (
        donor.groupby(["ligand_complex", "receptor_complex"], observed=True)
        .agg(
            informative_patients=("Patient", "nunique"),
            codetection_patients=("passes_10pct_codetection", "sum"),
            codetection_fraction=("passes_10pct_codetection", "mean"),
            median_coexpression_magnitude=("coexpression_magnitude", "median"),
            min_source_cells=("n_source_cells", "min"),
            min_target_cells=("n_target_cells", "min"),
        )
        .reset_index()
    )
    return summary


def balanced_downsample(adata: ad.AnnData) -> ad.AnnData:
    rng = np.random.default_rng(SEED)
    selected_indices: list[int] = []
    annotation = adata.obs["Annotation"].to_numpy()
    for state in sorted(np.unique(annotation)):
        idx = np.flatnonzero(annotation == state)
        take = min(len(idx), DOWNSAMPLE_CELLS_PER_STATE)
        selected_indices.extend(rng.choice(idx, size=take, replace=False).tolist())
    return adata[np.asarray(sorted(selected_indices), dtype=int)].copy()


def main() -> None:
    adata, state_summary, resources = build_input()
    state_summary.to_csv(OUT / "GSE206785_Fig1_CCC_state_coverage.tsv", sep="\t", index=False)
    pairs, fibro_states, myeloid_states = source_target_pairs(adata)
    pairs.to_csv(OUT / "GSE206785_Fig1_CCC_tested_state_pairs.tsv", sep="\t", index=False)

    results: dict[str, pd.DataFrame] = {}
    results["consensus"] = run_liana(
        adata, pairs, "consensus", None, N_PERMS, "consensus"
    )
    results["cellphonedb"] = run_liana(
        adata, pairs, "cellphonedb", None, N_PERMS, "cellphonedb"
    )
    results["cellchatdb"] = run_liana(
        adata, pairs, "cellchatdb", None, N_PERMS, "cellchatdb"
    )
    results["frozen_omnipath"] = run_liana(
        adata, pairs, "consensus", resources["frozen_omnipath"], N_PERMS, "frozen_omnipath"
    )

    source_target = summarize_source_target(results["consensus"])
    source_target.to_csv(
        OUT / "GSE206785_LIANA_consensus_source_target_summary.tsv", sep="\t", index=False
    )

    f13_m02 = results["consensus"].query("source == @F13 and target == @M02").copy()
    f13_m02 = f13_m02.sort_values(["combined_worst_rank", "combined_geometric_rank"])
    if f13_m02.empty:
        raise RuntimeError("Consensus resource yielded no F13-to-M02 interactions")

    membership_records: list[dict[str, object]] = []
    all_pairs: set[tuple[str, str]] = set()
    for name, resource in resources.items():
        membership = set(map(tuple, resource[["ligand", "receptor"]].astype(str).to_numpy()))
        all_pairs.update(membership)
        for ligand, receptor in membership:
            if ligand == "PLAT" or receptor == "PLAT":
                membership_records.append(
                    {"resource": name, "ligand": ligand, "receptor": receptor}
                )
    pd.DataFrame(membership_records).to_csv(
        OUT / "GSE206785_PLAT_resource_membership.tsv", sep="\t", index=False
    )

    resource_pair_sets = {
        name: set(map(tuple, frame[["ligand", "receptor"]].astype(str).to_numpy()))
        for name, frame in resources.items()
    }
    f13_m02["support_consensus_resource"] = [
        (l, r) in resource_pair_sets["consensus"]
        for l, r in zip(f13_m02["ligand_complex"], f13_m02["receptor_complex"])
    ]
    f13_m02["support_cellphonedb_resource"] = [
        (l, r) in resource_pair_sets["cellphonedb"]
        for l, r in zip(f13_m02["ligand_complex"], f13_m02["receptor_complex"])
    ]
    f13_m02["support_cellchatdb_resource"] = [
        (l, r) in resource_pair_sets["cellchatdb"]
        for l, r in zip(f13_m02["ligand_complex"], f13_m02["receptor_complex"])
    ]
    f13_m02["support_frozen_omnipath_resource"] = [
        (l, r) in resource_pair_sets["frozen_omnipath"]
        for l, r in zip(f13_m02["ligand_complex"], f13_m02["receptor_complex"])
    ]
    independent_cols = ["support_cellphonedb_resource", "support_cellchatdb_resource"]
    f13_m02["independent_curated_resource_count"] = f13_m02[independent_cols].sum(axis=1)

    # Evaluate every expression-qualified F13->M02 interaction at donor level.
    # This deliberately retains low-ranked PLAT rows so their expression
    # availability is visible even when their pooled communication rank is weak.
    donor = donor_support(adata, f13_m02)
    f13_m02 = f13_m02.merge(
        donor, on=["ligand_complex", "receptor_complex"], how="left", validate="one_to_one"
    )
    f13_m02.to_csv(
        OUT / "GSE206785_LIANA_consensus_F13_to_M02_ranked.tsv", sep="\t", index=False
    )

    downsampled = balanced_downsample(adata)
    down_pairs, _, _ = source_target_pairs(downsampled)
    down = run_liana(
        downsampled,
        down_pairs,
        "consensus",
        None,
        DOWNSAMPLE_N_PERMS,
        "consensus_state_balanced",
    )
    down_summary = summarize_source_target(down)
    down_summary.to_csv(
        OUT / "GSE206785_LIANA_consensus_state_balanced_source_target_summary.tsv",
        sep="\t",
        index=False,
    )
    merged = source_target.merge(
        down_summary,
        on=["source", "target"],
        suffixes=("_full", "_balanced"),
        validate="one_to_one",
    )
    rho = spearmanr(
        merged["top10_joint_strength_full"], merged["top10_joint_strength_balanced"]
    )

    plat_rows = f13_m02[
        f13_m02["ligand_complex"].astype(str).str.contains(r"(?:^|_)PLAT(?:_|$)", regex=True)
        | f13_m02["receptor_complex"].astype(str).str.contains(r"(?:^|_)PLAT(?:_|$)", regex=True)
    ].copy()
    plat_rows.to_csv(
        OUT / "GSE206785_LIANA_consensus_F13_to_M02_PLAT_rows.tsv", sep="\t", index=False
    )

    audit = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": "GSE206785",
        "purpose": "Figure 1 consensus fibroblast-to-myeloid communication analysis",
        "core_claim_boundary": (
            "co-expression proxy and predicted receiver programme; no binding or causal claim"
        ),
        "versions": {
            package: version(package)
            for package in ("numpy", "pandas", "scipy", "anndata", "scanpy", "liana")
        },
        "parameters": {
            "expr_prop": EXPR_PROP,
            "min_state_cells": MIN_STATE_CELLS,
            "min_state_patients": MIN_STATE_PATIENTS,
            "n_perms": N_PERMS,
            "downsample_n_perms": DOWNSAMPLE_N_PERMS,
            "downsample_cells_per_state": DOWNSAMPLE_CELLS_PER_STATE,
            "seed": SEED,
            "n_jobs": N_JOBS,
        },
        "input": {
            "expression_sha256": sha256_file(EXPR),
            "metadata_sha256": sha256_file(META),
            "frozen_lr_sha256": sha256_file(FROZEN_LR),
            "cells": int(adata.n_obs),
            "genes": int(adata.n_vars),
            "fibroblast_states": fibro_states,
            "myeloid_states": myeloid_states,
        },
        "sensitivity": {
            "source_target_joint_strength_full_vs_balanced_spearman_rho": float(rho.statistic),
            "source_target_joint_strength_full_vs_balanced_spearman_p": float(rho.pvalue),
        },
        "plat_consensus_rows": int(len(plat_rows)),
        "plat_has_independent_cellphonedb_or_cellchatdb_resource_support": bool(
            (plat_rows.get("independent_curated_resource_count", pd.Series(dtype=float)) > 0).any()
        ),
    }
    (AUDIT / "Figure1_GSE206785_cell_communication_manifest.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
