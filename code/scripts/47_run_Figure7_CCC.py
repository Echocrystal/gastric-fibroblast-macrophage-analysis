#!/usr/bin/env python
"""Run the preregistered Figure 7 communication analyses.

The script writes method-level source tables only. Figure drawing is handled by
a separate matplotlib script after numerical and claim-boundary review.
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
from scipy.sparse import hstack


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
INPUT183 = PROJECT / "analysis" / "processed" / "communication" / "GSE183904"
INPUT206 = PROJECT / "analysis" / "processed" / "virtual_knockout" / "GSE206785"
OUT = PROJECT / "analysis" / "results" / "communication" / "Figure7_redesign"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 20260801
N_PERMS = 1000
EXPR_PROP = 0.10
MIN_CELLS = 10
N_JOBS = 4


def log(msg: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{stamp}] {msg}", flush=True)


def read_bundle(prefix: Path, group_col: str) -> ad.AnnData:
    log(f"read matrix {prefix.name}")
    x = mmread(str(prefix) + "_raw_counts_genes_by_cells.mtx").T.tocsr()
    genes = pd.read_csv(str(prefix) + "_genes.tsv", sep="\t").iloc[:, 0].astype(str)
    obs = pd.read_csv(str(prefix) + "_cells.tsv", sep="\t")
    if "cell_id" not in obs:
        raise ValueError(f"Missing cell_id in {prefix}_cells.tsv")
    obs = obs.set_index("cell_id", drop=False)
    if x.shape != (len(obs), len(genes)):
        raise ValueError(f"Shape mismatch for {prefix}: {x.shape}, {len(obs)}, {len(genes)}")
    if obs.index.duplicated().any():
        raise ValueError(f"Duplicated cell identifiers in {prefix}")
    if genes.duplicated().any():
        raise ValueError(f"Duplicated gene symbols in {prefix}")
    a = ad.AnnData(x, obs=obs, var=pd.DataFrame(index=genes.to_numpy()))
    a.obs[group_col] = a.obs[group_col].astype("category")
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def read_gse206785() -> ad.AnnData:
    pieces = []
    obs_parts = []
    genes_ref = None
    specs = [
        ("sender_CTHRC1_CAF", "F13_CTHRC1_sender"),
        ("receiver_CD14_monocyte", "M02_CD14_receiver"),
    ]
    for stem, group in specs:
        prefix = INPUT206 / f"GSE206785_{stem}_full_raw_counts"
        x = mmread(str(prefix) + ".mtx.gz").tocsc()
        genes = pd.read_csv(str(prefix) + "_genes.tsv", sep="\t").iloc[:, 0].astype(str)
        cells = pd.read_csv(str(prefix) + "_cells.tsv", sep="\t").iloc[:, 0].astype(str)
        meta = pd.read_csv(str(prefix) + "_metadata.tsv", sep="\t")
        meta = meta.set_index("cell_id").loc[cells].copy()
        meta["ccc_group"] = group
        if genes_ref is None:
            genes_ref = genes
        elif not genes_ref.equals(genes):
            raise ValueError("GSE206785 sender/receiver gene-order mismatch")
        if x.shape != (len(genes), len(cells)):
            raise ValueError(f"GSE206785 shape mismatch for {stem}")
        pieces.append(x)
        obs_parts.append(meta)
    x = hstack(pieces, format="csr").T
    obs = pd.concat(obs_parts, axis=0)
    if obs.index.duplicated().any():
        raise ValueError("Duplicated GSE206785 cell identifiers")
    a = ad.AnnData(x, obs=obs, var=pd.DataFrame(index=genes_ref.to_numpy()))
    a.obs["ccc_group"] = a.obs["ccc_group"].astype("category")
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def run_method(
    a: ad.AnnData,
    groupby: str,
    method_name: str,
    resource_name: str,
    outfile: Path,
) -> pd.DataFrame:
    if outfile.exists() and outfile.stat().st_size > 0:
        log(f"reuse {outfile.name}")
        return pd.read_csv(outfile, sep="\t")
    method = {
        "rank_aggregate": li.mt.rank_aggregate,
        "cellphonedb": li.mt.cellphonedb,
        "cellchat": li.mt.cellchat,
    }[method_name]
    log(f"run {method_name}, resource={resource_name}, cells={a.n_obs}, groups={a.obs[groupby].nunique()}")
    t0 = time.time()
    result = method(
        a,
        groupby=groupby,
        resource_name=resource_name,
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
    log(f"wrote {outfile.name}: {len(result)} rows, {time.time() - t0:.1f} s")
    return result


def aggregate_global_matrix(rank: pd.DataFrame, outfile: Path) -> None:
    z = rank.copy()
    z["robust"] = (z["magnitude_rank"] <= 0.05) & (z["specificity_rank"] <= 0.05)
    z["evidence_weight"] = np.where(
        z["robust"], -np.log10(z["magnitude_rank"].clip(lower=1e-12)), 0.0
    )
    ans = (
        z.groupby(["source", "target"], observed=True)
        .agg(
            n_tested=("ligand_complex", "size"),
            n_robust=("robust", "sum"),
            evidence_weight=("evidence_weight", "sum"),
            best_magnitude_rank=("magnitude_rank", "min"),
            best_specificity_rank=("specificity_rank", "min"),
        )
        .reset_index()
    )
    ans.to_csv(outfile, sep="\t", index=False)


def aggregate_cell_roles(cellchat: pd.DataFrame, outfile: Path) -> None:
    z = cellchat.copy()
    z["supported"] = z["cellchat_pvals"] <= 0.05
    z["supported_prob"] = np.where(z["supported"], z["lr_probs"], 0.0)
    outgoing = (
        z.groupby("source", observed=True)
        .agg(outgoing_probability=("supported_prob", "sum"), outgoing_edges=("supported", "sum"))
        .reset_index()
        .rename(columns={"source": "cell_group"})
    )
    incoming = (
        z.groupby("target", observed=True)
        .agg(incoming_probability=("supported_prob", "sum"), incoming_edges=("supported", "sum"))
        .reset_index()
        .rename(columns={"target": "cell_group"})
    )
    outgoing.merge(incoming, on="cell_group", how="outer").fillna(0).to_csv(
        outfile, sep="\t", index=False
    )


def donor_pair_support(
    a: ad.AnnData,
    pooled: pd.DataFrame,
    donor_col: str,
    group_col: str,
    source_group: str,
    target_group: str,
    outfile: Path,
    top_n: int = 40,
) -> pd.DataFrame:
    focus_all = pooled[(pooled.source == source_group) & (pooled.target == target_group)].copy()
    focus_top = focus_all.sort_values(["magnitude_rank", "specificity_rank"]).head(top_n)
    # PLAT rows are always audited even when they do not enter the top-N pooled
    # edges; this prevents an apparently absent donor result from being confused
    # with a failed expression gate.
    mandatory_ligands = {"PLAT", "APP", "NRG1", "CXCL14", "SPON2", "RARRES2"}
    focus_mandatory = focus_all[focus_all["ligand_complex"].astype(str).isin(mandatory_ligands)]
    focus = pd.concat([focus_top, focus_mandatory], ignore_index=True).drop_duplicates(
        ["source", "target", "ligand_complex", "receptor_complex"]
    )
    rows = []
    donors = sorted(a.obs[donor_col].astype(str).unique())
    for donor in donors:
        obs = a.obs
        source_ix = np.flatnonzero(
            (obs[donor_col].astype(str).to_numpy() == donor)
            & (obs[group_col].astype(str).to_numpy() == source_group)
        )
        target_ix = np.flatnonzero(
            (obs[donor_col].astype(str).to_numpy() == donor)
            & (obs[group_col].astype(str).to_numpy() == target_group)
        )
        if len(source_ix) < MIN_CELLS or len(target_ix) < MIN_CELLS:
            continue
        for r in focus.itertuples(index=False):
            ligands = str(r.ligand_complex).split("_")
            receptors = str(r.receptor_complex).split("_")
            if not set(ligands + receptors).issubset(a.var_names):
                continue
            lidx = [a.var_names.get_loc(g) for g in ligands]
            ridx = [a.var_names.get_loc(g) for g in receptors]
            ldetect = np.asarray((a.X[source_ix][:, lidx] > 0).mean(axis=0)).ravel()
            rdetect = np.asarray((a.X[target_ix][:, ridx] > 0).mean(axis=0)).ravel()
            lmean = np.asarray(a.X[source_ix][:, lidx].mean(axis=0)).ravel()
            rmean = np.asarray(a.X[target_ix][:, ridx].mean(axis=0)).ravel()
            support = bool((ldetect >= EXPR_PROP).all() and (rdetect >= EXPR_PROP).all())
            magnitude = float(np.sqrt(max(lmean.min(), 0) * max(rmean.min(), 0)))
            rows.append(
                {
                    "donor": donor,
                    "source": source_group,
                    "target": target_group,
                    "ligand_complex": r.ligand_complex,
                    "receptor_complex": r.receptor_complex,
                    "n_source_cells": len(source_ix),
                    "n_target_cells": len(target_ix),
                    "ligand_min_detection": float(ldetect.min()),
                    "receptor_min_detection": float(rdetect.min()),
                    "expression_gate_support": support,
                    "donor_magnitude": magnitude,
                }
            )
    ans = pd.DataFrame(rows)
    ans.to_csv(outfile, sep="\t", index=False)
    summary = (
        ans.groupby(["source", "target", "ligand_complex", "receptor_complex"], observed=True)
        .agg(
            eligible_donors=("donor", "nunique"),
            supported_donors=("expression_gate_support", "sum"),
            median_donor_magnitude=("donor_magnitude", "median"),
            min_donor_magnitude=("donor_magnitude", "min"),
        )
        .reset_index()
    )
    summary["support_fraction"] = summary.supported_donors / summary.eligible_donors
    summary.to_csv(outfile.with_name(outfile.stem + "_summary.tsv"), sep="\t", index=False)
    return summary


def focused_consensus_table(
    pooled_rank: pd.DataFrame,
    cpdb: pd.DataFrame,
    cellchat: pd.DataFrame,
    resource_runs: dict[str, pd.DataFrame],
    donor_summary: pd.DataFrame,
    outfile: Path,
) -> None:
    keys = ["source", "target", "ligand_complex", "receptor_complex"]
    z = pooled_rank.copy()
    z = z[z.source.str.startswith("F") & z.target.str.startswith("M")]
    z = z.merge(
        cpdb[keys + ["lr_means", "cellphone_pvals"]], on=keys, how="left"
    ).merge(
        cellchat[keys + ["lr_probs", "cellchat_pvals"]], on=keys, how="left"
    )
    for resource, table in resource_runs.items():
        present = table[keys].drop_duplicates().assign(**{f"present_{resource}": True})
        z = z.merge(present, on=keys, how="left")
        z[f"present_{resource}"] = z[f"present_{resource}"].fillna(False)
    z = z.merge(donor_summary, on=keys, how="left")
    z.to_csv(outfile, sep="\t", index=False)


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
    }
    (OUT / "software_and_parameters.json").write_text(json.dumps(versions, indent=2), encoding="utf-8")

    global_a = read_bundle(INPUT183 / "GSE183904_global_donor_balanced", "coarse_label")
    global_rank = run_method(
        global_a, "coarse_label", "rank_aggregate", "consensus", OUT / "GSE183904_global_LIANA_rank_aggregate.tsv"
    )
    global_cc = run_method(
        global_a, "coarse_label", "cellchat", "consensus", OUT / "GSE183904_global_CellChat_scores.tsv"
    )
    aggregate_global_matrix(global_rank, OUT / "panel_a_global_communication_matrix.tsv")
    aggregate_cell_roles(global_cc, OUT / "panel_b_global_cell_roles.tsv")
    del global_a

    fine_a = read_bundle(INPUT183 / "GSE183904_fine_donor_balanced", "ccc_group")
    fine_rank = run_method(
        fine_a, "ccc_group", "rank_aggregate", "consensus", OUT / "GSE183904_fine_LIANA_rank_aggregate.tsv"
    )
    fine_cpdb = run_method(
        fine_a, "ccc_group", "cellphonedb", "consensus", OUT / "GSE183904_fine_CellPhoneDB_scores.tsv"
    )
    fine_cc = run_method(
        fine_a, "ccc_group", "cellchat", "consensus", OUT / "GSE183904_fine_CellChat_scores.tsv"
    )
    resource_runs = {}
    for resource in ["cellphonedb", "cellchatdb"]:
        resource_runs[resource] = run_method(
            fine_a,
            "ccc_group",
            "rank_aggregate",
            resource,
            OUT / f"GSE183904_fine_LIANA_rank_aggregate_{resource}.tsv",
        )

    full_a = read_bundle(INPUT183 / "GSE183904_fine_full", "ccc_group")
    donor_summary = donor_pair_support(
        full_a,
        fine_rank,
        "Patient_ID",
        "ccc_group",
        "F1_PLATrich",
        "M3_C1QC_APOE_target",
        OUT / "GSE183904_F1_to_M3_donor_pair_support.tsv",
    )
    del full_a

    focused_consensus_table(
        fine_rank,
        fine_cpdb,
        fine_cc,
        resource_runs,
        donor_summary,
        OUT / "panel_c_d_focused_fibroblast_myeloid_consensus.tsv",
    )
    del fine_a

    val_a = read_gse206785()
    val_rank = run_method(
        val_a, "ccc_group", "rank_aggregate", "consensus", OUT / "GSE206785_focused_LIANA_rank_aggregate.tsv"
    )
    val_cpdb = run_method(
        val_a, "ccc_group", "cellphonedb", "consensus", OUT / "GSE206785_focused_CellPhoneDB_scores.tsv"
    )
    val_cc = run_method(
        val_a, "ccc_group", "cellchat", "consensus", OUT / "GSE206785_focused_CellChat_scores.tsv"
    )
    val_donor = donor_pair_support(
        val_a,
        val_rank,
        "Patient",
        "ccc_group",
        "F13_CTHRC1_sender",
        "M02_CD14_receiver",
        OUT / "GSE206785_F13_to_M02_donor_pair_support.tsv",
    )

    # Cross-cohort joins use ligand-receptor identities; state labels remain cohort-specific.
    key_lr = ["ligand_complex", "receptor_complex"]
    disc = donor_summary.copy().rename(columns={c: f"GSE183904_{c}" for c in donor_summary.columns if c not in key_lr})
    val = val_donor.copy().rename(columns={c: f"GSE206785_{c}" for c in val_donor.columns if c not in key_lr})
    cross = disc.merge(val, on=key_lr, how="outer")
    cross.to_csv(OUT / "panel_f_cross_cohort_donor_recurrence.tsv", sep="\t", index=False)

    # Preserve validation method tables for later panel-level consensus assembly.
    log(f"GSE206785 validation rows: LIANA={len(val_rank)}, CPDB={len(val_cpdb)}, CellChat={len(val_cc)}")
    log("Figure 7 communication inference complete")


if __name__ == "__main__":
    main()
