#!/usr/bin/env python
"""Run all-five-myeloid-state communication analyses for Figure 3.

The analysis is intentionally state-first: the sender is the F1 PLAT-rich
fibroblast state, but the ligand-receptor repertoire is not restricted to PLAT
as a ligand.  C1-C5 are all retained and evaluated with the same thresholds.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import anndata as ad
import liana as li
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
INPUT = PROJECT / "analysis" / "processed" / "communication" / "GSE183904"
OUT = PROJECT / "analysis" / "results" / "communication" / "Figure3_all5_raw"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 20260804
N_PERMS = 1000
EXPR_PROP = 0.10
MIN_CELLS = 10
N_JOBS = 4
SOURCE = "F1_PLATrich"
TARGETS = [
    "M1_Unresolved",
    "M2_SPP1_APOE_sensitivity",
    "M3_C1QC_APOE_target",
    "M4_NeutrophilLike",
    "M5_C1QC_boundary",
]


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def read_bundle(prefix: Path) -> ad.AnnData:
    log(f"read {prefix.name}")
    matrix = mmread(str(prefix) + "_raw_counts_genes_by_cells.mtx").T.tocsr()
    genes = pd.read_csv(str(prefix) + "_genes.tsv", sep="\t").iloc[:, 0].astype(str)
    obs = pd.read_csv(str(prefix) + "_cells.tsv", sep="\t").set_index("cell_id", drop=False)
    if matrix.shape != (len(obs), len(genes)):
        raise ValueError(f"Shape mismatch: {matrix.shape} vs {(len(obs), len(genes))}")
    if obs.index.duplicated().any() or genes.duplicated().any():
        raise ValueError("Duplicated cell identifiers or gene symbols")
    data = ad.AnnData(matrix, obs=obs, var=pd.DataFrame(index=genes.to_numpy()))
    data.obs["ccc_group"] = data.obs["ccc_group"].astype("category")
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)
    return data


def run_liana(data: ad.AnnData, outfile: Path) -> pd.DataFrame:
    if outfile.exists() and outfile.stat().st_size > 0:
        log(f"reuse {outfile.name}")
        return pd.read_csv(outfile, sep="\t")
    log(f"run LIANA rank_aggregate on {data.n_obs} cells and {data.obs.ccc_group.nunique()} groups")
    result = li.mt.rank_aggregate(
        data,
        groupby="ccc_group",
        resource_name="consensus",
        expr_prop=EXPR_PROP,
        min_cells=MIN_CELLS,
        use_raw=False,
        n_perms=N_PERMS,
        seed=SEED,
        n_jobs=N_JOBS,
        inplace=False,
        verbose=True,
    )
    result.to_csv(outfile, sep="\t", index=False)
    return result


def run_cellchat(data: ad.AnnData, outfile: Path) -> pd.DataFrame:
    if outfile.exists() and outfile.stat().st_size > 0:
        log(f"reuse {outfile.name}")
        return pd.read_csv(outfile, sep="\t")
    log(f"run LIANA CellChat method on {data.n_obs} cells")
    result = li.mt.cellchat(
        data,
        groupby="ccc_group",
        resource_name="consensus",
        expr_prop=EXPR_PROP,
        min_cells=MIN_CELLS,
        use_raw=False,
        n_perms=N_PERMS,
        seed=SEED,
        n_jobs=N_JOBS,
        inplace=False,
        verbose=True,
    )
    result.to_csv(outfile, sep="\t", index=False)
    return result


def select_display_pairs(rank: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Choose a target-balanced display set without privileging C3.

    For each C1-C5 target, take the six strongest consensus-magnitude pairs
    among interactions with specificity_rank <= 0.10. The union is ordered by
    its best magnitude rank across the five targets. The two best observed
    PLAT-ligand identities are appended so the marker gene remains visible
    without restricting the repertoire to PLAT-only edges.
    """
    focus = rank[(rank.source == SOURCE) & rank.target.isin(TARGETS)].copy()
    focus["pair"] = focus.ligand_complex.astype(str) + "–" + focus.receptor_complex.astype(str)
    selected = []
    for target in TARGETS:
        state = focus[focus.target == target].sort_values(
            ["magnitude_rank", "specificity_rank"], kind="stable"
        )
        eligible = state[state.specificity_rank <= 0.10]
        if len(eligible) < 6:
            eligible = state
        selected.extend(eligible.head(6).pair.tolist())

    plat = focus[focus.ligand_complex.astype(str) == "PLAT"].sort_values(
        ["magnitude_rank", "specificity_rank"], kind="stable"
    )
    selected.extend(plat.drop_duplicates("pair").head(2).pair.tolist())
    selected = list(dict.fromkeys(selected))

    order = (
        focus[focus.pair.isin(selected)]
        .groupby("pair", observed=True)
        .agg(
            best_magnitude_rank=("magnitude_rank", "min"),
            best_specificity_rank=("specificity_rank", "min"),
            n_states_observed=("target", "nunique"),
        )
        .reset_index()
        .sort_values(["best_magnitude_rank", "best_specificity_rank", "pair"], kind="stable")
    )
    order["display_order"] = np.arange(1, len(order) + 1)
    matrix = focus[focus.pair.isin(order.pair)].merge(
        order[["pair", "display_order"]], on="pair", how="left"
    )
    matrix["target_order"] = matrix.target.map({name: i + 1 for i, name in enumerate(TARGETS)})
    matrix = matrix.sort_values(["display_order", "target_order"])
    matrix.to_csv(OUT / "panel_c_LIANA_all5_display_matrix.tsv", sep="\t", index=False)
    order.to_csv(OUT / "panel_c_LIANA_display_pair_order.tsv", sep="\t", index=False)
    return matrix, order


def summarize_cellchat(cellchat: pd.DataFrame) -> pd.DataFrame:
    focus = cellchat[(cellchat.source == SOURCE) & cellchat.target.isin(TARGETS)].copy()
    focus["supported"] = focus.cellchat_pvals <= 0.05
    focus["supported_probability"] = np.where(focus.supported, focus.lr_probs, 0.0)
    summary = (
        focus.groupby("target", observed=True)
        .agg(
            incoming_probability=("supported_probability", "sum"),
            supported_lr_pairs=("supported", "sum"),
            tested_lr_pairs=("ligand_complex", "size"),
            best_probability=("lr_probs", "max"),
        )
        .reindex(TARGETS)
        .reset_index()
    )
    summary.to_csv(OUT / "panel_d_CellChat_all5_incoming.tsv", sep="\t", index=False)
    focus.to_csv(OUT / "panel_d_CellChat_all5_focused_raw.tsv", sep="\t", index=False)
    return summary


def component_indices(data: ad.AnnData, complex_name: str) -> list[int] | None:
    genes = str(complex_name).split("_")
    if not set(genes).issubset(data.var_names):
        return None
    return [data.var_names.get_loc(gene) for gene in genes]


def donor_expression_potential(full: ad.AnnData, display: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    identities = display[["pair", "ligand_complex", "receptor_complex"]].drop_duplicates("pair")
    rows = []
    donor_values = full.obs.Patient_ID.astype(str).to_numpy()
    group_values = full.obs.ccc_group.astype(str).to_numpy()
    for donor in sorted(np.unique(donor_values)):
        source_ix = np.flatnonzero((donor_values == donor) & (group_values == SOURCE))
        if len(source_ix) < MIN_CELLS:
            continue
        for target in TARGETS:
            target_ix = np.flatnonzero((donor_values == donor) & (group_values == target))
            if len(target_ix) < MIN_CELLS:
                continue
            for row in identities.itertuples(index=False):
                ligand_ix = component_indices(full, row.ligand_complex)
                receptor_ix = component_indices(full, row.receptor_complex)
                if ligand_ix is None or receptor_ix is None:
                    continue
                ligand_detect = np.asarray((full.X[source_ix][:, ligand_ix] > 0).mean(axis=0)).ravel()
                receptor_detect = np.asarray((full.X[target_ix][:, receptor_ix] > 0).mean(axis=0)).ravel()
                ligand_mean = np.asarray(full.X[source_ix][:, ligand_ix].mean(axis=0)).ravel()
                receptor_mean = np.asarray(full.X[target_ix][:, receptor_ix].mean(axis=0)).ravel()
                supported = bool((ligand_detect >= EXPR_PROP).all() and (receptor_detect >= EXPR_PROP).all())
                potential = float(np.sqrt(max(ligand_mean.min(), 0.0) * max(receptor_mean.min(), 0.0)))
                rows.append(
                    {
                        "Patient_ID": donor,
                        "target": target,
                        "pair": row.pair,
                        "n_source_cells": len(source_ix),
                        "n_target_cells": len(target_ix),
                        "ligand_min_detection": float(ligand_detect.min()),
                        "receptor_min_detection": float(receptor_detect.min()),
                        "expression_gate_support": supported,
                        "lr_expression_potential": potential,
                    }
                )
    pair_level = pd.DataFrame(rows)
    pair_level.to_csv(OUT / "panel_e_all5_donor_pair_potentials.tsv", sep="\t", index=False)
    programme = (
        pair_level[pair_level.expression_gate_support]
        .groupby(["Patient_ID", "target"], observed=True)
        .agg(
            lr_programme_potential=("lr_expression_potential", "median"),
            n_supported_display_pairs=("pair", "nunique"),
            n_source_cells=("n_source_cells", "first"),
            n_target_cells=("n_target_cells", "first"),
        )
        .reset_index()
    )
    programme.to_csv(OUT / "panel_e_all5_donor_programme_potential.tsv", sep="\t", index=False)
    return pair_level, programme


def main() -> None:
    versions = {
        "python": sys.version,
        "liana": li.__version__,
        "anndata": ad.__version__,
        "scanpy": sc.__version__,
        "seed": SEED,
        "n_perms": N_PERMS,
        "expr_prop": EXPR_PROP,
        "min_cells": MIN_CELLS,
        "selection_rule": "top six LIANA magnitude pairs per target among specificity_rank <= 0.10, union; append top two PLAT pair identities",
    }
    (OUT / "software_and_parameters.json").write_text(json.dumps(versions, indent=2), encoding="utf-8")

    balanced = read_bundle(INPUT / "GSE183904_fine_all5_donor_balanced")
    rank = run_liana(balanced, OUT / "GSE183904_fine_all5_LIANA_rank_aggregate.tsv")
    cellchat = run_cellchat(balanced, OUT / "GSE183904_fine_all5_CellChat_scores.tsv")
    display, order = select_display_pairs(rank)
    summarize_cellchat(cellchat)
    log(f"selected {len(order)} ligand-receptor identities for the all-five-state display")
    del balanced

    full = read_bundle(INPUT / "GSE183904_fine_all5_full")
    _, programme = donor_expression_potential(full, display)
    log(f"wrote {len(programme)} eligible donor-by-state programme estimates")
    log("Figure 3 all-five-state communication analysis complete")


if __name__ == "__main__":
    main()
