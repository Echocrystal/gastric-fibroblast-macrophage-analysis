#!/usr/bin/env python
"""Prepare receiver-unbiased Figure 1 candidate and communication sources.

The former Figure 1 universe used an M02-specific receptor-availability gate.
This revision removes that gate. Candidate construction and ranking are
therefore sender-side only; the receiver lineage is nominated later by the
independent, coarse-lineage communication analysis.
"""

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
DATA = PROJECT / "figures" / "data"
DATA.mkdir(parents=True, exist_ok=True)

MEMBERSHIP = DATA / "Figure1B_four_set_gene_membership.tsv.gz"
EXPRESSION = (
    PROJECT
    / "analysis/processed/GSE206785_tumor_normal_fine_state_lr_expression.tsv.gz"
)
MATRISOME = Path(r"I:\ECM-new\ECM-public\gene_sets\naba_matrisome_locked_membership.tsv")
GLOBAL_MATRIX = (
    PROJECT
    / "analysis/results/communication/Figure7_redesign/"
      "panel_a_global_communication_matrix.tsv"
)

F13 = "F13_Activated_CTHRC1"
MIN_CELLS = 5
MIN_STATE_PATIENTS = 6
DETECTION_GATE = 0.10


def percentile(series: pd.Series) -> pd.Series:
    """Average-tie percentile in (0, 1], deterministic across platforms."""
    return series.rank(method="average", pct=True)


membership = pd.read_csv(MEMBERSHIP, sep="\t")
set_cols = ["Matrisome", "LR prior", "NicheNet", "F13 expression"]
candidate_genes = sorted(
    membership.loc[membership[set_cols].all(axis=1), "gene_symbol"].astype(str)
)
if len(candidate_genes) != 88:
    raise ValueError(f"Expected 88 four-way candidates, found {len(candidate_genes)}")

expr = pd.read_csv(
    EXPRESSION,
    sep="\t",
    usecols=[
        "Patient", "Tissue", "Type", "Annotation", "gene_symbol", "n_cells",
        "mean_log_expression", "detected_fraction",
    ],
)
fib = expr[
    expr["Tissue"].eq("Tumor")
    & expr["Type"].eq("Fibroblast")
    & expr["gene_symbol"].isin(candidate_genes)
    & expr["n_cells"].ge(MIN_CELLS)
].copy()

coverage = (
    fib[["Patient", "Annotation"]]
    .drop_duplicates()
    .groupby("Annotation", observed=True)["Patient"]
    .nunique()
)
eligible_states = coverage[coverage.ge(MIN_STATE_PATIENTS)].index
fib = fib[fib["Annotation"].isin(eligible_states)].copy()

state = (
    fib.groupby(["Annotation", "gene_symbol"], observed=True)
    .agg(
        patient_equal_mean=("mean_log_expression", "mean"),
        patient_equal_detection=("detected_fraction", "mean"),
        patients=("Patient", "nunique"),
        patient_recurrence=(
            "detected_fraction", lambda x: float(np.mean(np.asarray(x) >= DETECTION_GATE))
        ),
    )
    .reset_index()
)

sender = state[state["Annotation"].eq(F13)].copy()
other = (
    state[~state["Annotation"].eq(F13)]
    .groupby("gene_symbol", observed=True)
    .agg(
        best_other_mean=("patient_equal_mean", "max"),
        best_other_detection=("patient_equal_detection", "max"),
    )
    .reset_index()
)

ranking = (
    pd.DataFrame({"ligand": candidate_genes})
    .merge(sender.drop(columns="Annotation"), left_on="ligand", right_on="gene_symbol", how="left")
    .drop(columns="gene_symbol")
    .merge(other, left_on="ligand", right_on="gene_symbol", how="left")
    .drop(columns="gene_symbol")
    .fillna(0)
)
ranking["sender_specificity_log2_ratio"] = np.log2(
    (ranking["patient_equal_mean"] + 0.05) / (ranking["best_other_mean"] + 0.05)
)

# The score is explicitly receiver-unbiased. It rewards expression level,
# detection, patient recurrence and specificity of the published CTHRC1+
# activated-fibroblast state relative to other recurrent fibroblast states.
metrics = [
    "patient_equal_mean",
    "patient_equal_detection",
    "patient_recurrence",
    "sender_specificity_log2_ratio",
]
for metric in metrics:
    ranking[f"percentile_{metric}"] = percentile(ranking[metric])
ranking["evidence_score"] = ranking[[f"percentile_{x}" for x in metrics]].mean(axis=1)

mat = pd.read_csv(MATRISOME, sep="\t")
priority_sets = [
    ("NABA_COLLAGENS", "Collagen"),
    ("NABA_ECM_GLYCOPROTEINS", "ECM glycoprotein"),
    ("NABA_ECM_REGULATORS", "ECM regulator"),
    ("NABA_PROTEOGLYCANS", "Proteoglycan"),
]
class_map = {}
for set_name, label in priority_sets:
    genes = set(mat.loc[mat["gs_name"].eq(set_name), "gene_symbol"].astype(str))
    for gene in candidate_genes:
        if gene in genes and gene not in class_map:
            class_map[gene] = label
ranking["matrisome_class"] = ranking["ligand"].map(class_map).fillna("Other matrisome")

ranking = ranking.sort_values(
    ["evidence_score", "patient_equal_detection", "ligand"],
    ascending=[False, False, True],
).reset_index(drop=True)
ranking["scan_rank"] = np.arange(1, len(ranking) + 1)

ranking.to_csv(DATA / "Figure1C_ranked_88_sender_ligands_source.tsv", sep="\t", index=False)
state.to_csv(DATA / "Figure1C_sender_state_expression_source.tsv", sep="\t", index=False)

universe = pd.DataFrame(
    [
        ["NABA Matrisome", 513, "unique genes"],
        ["Ligand-receptor prior", 1226, "unique ligands"],
        ["NicheNet ligand-target prior", 1226, "unique ligands"],
        ["F13 activated-fibroblast expression", 6591, "expressed genes"],
        ["Four-way sender-side intersection", 88, "ECM-associated ligands"],
    ],
    columns=["step", "n", "unit"],
)
universe.to_csv(DATA / "Figure1B_candidate_universe_source.tsv", sep="\t", index=False)

global_matrix = pd.read_csv(GLOBAL_MATRIX, sep="\t")
fib_out = (
    global_matrix[global_matrix["source"].eq("Fibroblast")]
    .sort_values(["evidence_weight", "n_robust"], ascending=[False, False])
    .reset_index(drop=True)
)
if set(fib_out["target"]) != {
    "Epithelial", "Endothelial", "Fibroblast", "Mural", "T_cell", "NK",
    "B_cell", "Plasma", "Myeloid", "Mast",
}:
    raise ValueError("Unexpected coarse receiver-lineage set")
fib_out["receiver_rank"] = np.arange(1, len(fib_out) + 1)
fib_out.to_csv(DATA / "Figure1D_fibroblast_all_lineage_communication_source.tsv", sep="\t", index=False)

plat = ranking.loc[ranking["ligand"].eq("PLAT")].iloc[0]
print(
    f"Prepared {len(ranking)} receiver-unbiased candidates; "
    f"PLAT rank {int(plat.scan_rank)}/{len(ranking)}, score={plat.evidence_score:.6f}."
)
print(fib_out[["receiver_rank", "target", "n_robust", "evidence_weight"]].to_string(index=False))
