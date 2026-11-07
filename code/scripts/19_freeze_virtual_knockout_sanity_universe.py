#!/usr/bin/env python3
"""Freeze a compact, donor-capped GRN universe for virtual-KO calibration.

The full 25,287-gene matrices are retained as immutable inputs. This script
creates a deterministic sanity-analysis layer only: at most 30 unique cells
per donor, 250 HVGs per target state, four pre-existing ECM/myeloid program
families, and forced candidate/benchmark genes. No outcome, survival, or ICI
information enters the selection.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
INPUT = PROJECT / "analysis/processed/virtual_knockout/GSE206785"
OUT = INPUT / "sanity_v1"
RESOURCE = PROJECT / "resources/gene_sets"
AUDIT = PROJECT / "analysis/audit"
OUT.mkdir(parents=True, exist_ok=True)
AUDIT.mkdir(parents=True, exist_ok=True)

SEED = 20260726
DONOR_CAP = 30
N_HVG_PER_STATE = 250
MIN_CELL_FRACTION = 0.05
MIN_PATIENTS = 3

STATES = {
    "sender_CTHRC1_CAF": "F13_Activated_CTHRC1",
    "receiver_CD14_monocyte": "M02_Mono_CD14",
}

PRIMARY_PROGRAMS = [
    "GOBP_CELL_MATRIX_ADHESION",
    "GOBP_FOCAL_ADHESION_ASSEMBLY",
    "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY",
    "REACTOME_INTEGRIN_SIGNALING",
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
]

AXIS_GENES = [
    "TGFBI",
    "ANGPTL2",
    "PLAU",
    "TNC",
    "ITGA5",
    "ITGB1",
    "PLAUR",
]

# These genes are not used to select an axis; they ensure that established
# focal-adhesion machinery and candidate downstream TFs remain observable.
MECHANISM_GUARD_GENES = [
    "CTHRC1",
    "PTK2",
    "SRC",
    "TLN1",
    "VCL",
    "ILK",
    "FERMT3",
    "RHOA",
    "RAC1",
    "ROCK1",
    "ROCK2",
    "RELA",
    "STAT3",
    "JUN",
    "FOS",
    "HIF1A",
]

CIRCULARITY_EXCLUSION = set(
    AXIS_GENES + ["CTHRC1", "HLA-F", "LILRB2"]
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_state(label: str) -> tuple[sparse.csc_matrix, pd.Series, pd.DataFrame]:
    prefix = INPUT / f"GSE206785_{label}_full_raw_counts"
    matrix = mmread(str(prefix) + ".mtx.gz").tocsc().astype(np.int32)
    genes = pd.read_csv(str(prefix) + "_genes.tsv", sep="\t")["gene_symbol"]
    metadata = pd.read_csv(str(prefix) + "_metadata.tsv", sep="\t", dtype=str)
    if matrix.shape != (len(genes), len(metadata)):
        raise RuntimeError(f"Shape mismatch for {label}: {matrix.shape}")
    if genes.duplicated().any() or metadata["cell_id"].duplicated().any():
        raise RuntimeError(f"Duplicated genes or cells in {label}")
    return matrix, genes, metadata


def donor_capped_indices(metadata: pd.DataFrame, state_offset: int) -> np.ndarray:
    selected: list[int] = []
    for donor_rank, (_, frame) in enumerate(
        metadata.groupby("Patient", sort=True), start=1
    ):
        idx = frame.index.to_numpy(dtype=int)
        if len(idx) > DONOR_CAP:
            rng = np.random.default_rng(SEED + state_offset + donor_rank)
            idx = np.sort(rng.choice(idx, size=DONOR_CAP, replace=False))
        selected.extend(idx.tolist())
    return np.asarray(sorted(selected), dtype=int)


def gene_metrics(
    matrix: sparse.csc_matrix, metadata: pd.DataFrame
) -> pd.DataFrame:
    detected = matrix > 0
    cell_fraction = np.asarray(detected.mean(axis=1)).ravel()
    patient_count = np.zeros(matrix.shape[0], dtype=np.int32)
    for _, frame in metadata.groupby("Patient", sort=True):
        patient_count += np.asarray(detected[:, frame.index].sum(axis=1)).ravel() > 0

    dense = matrix.toarray().astype(np.float64)
    library_size = dense.sum(axis=0)
    if np.any(library_size <= 0):
        raise RuntimeError("Zero-library cell passed into HVG selection")
    lognorm = np.log1p(dense * (1e4 / library_size)[None, :])
    mean = lognorm.mean(axis=1)
    variance = lognorm.var(axis=1, ddof=1)

    valid = (
        (cell_fraction >= MIN_CELL_FRACTION)
        & (patient_count >= MIN_PATIENTS)
        & np.isfinite(mean)
        & np.isfinite(variance)
    )
    hvg_score = np.full(matrix.shape[0], np.nan)
    valid_idx = np.flatnonzero(valid)
    mean_rank = pd.Series(mean[valid_idx]).rank(method="first")
    bins = pd.qcut(mean_rank, q=min(20, len(valid_idx)), labels=False, duplicates="drop")
    log_variance = np.log1p(variance[valid_idx])
    for bin_id in np.unique(bins):
        in_bin = np.flatnonzero(np.asarray(bins) == bin_id)
        values = log_variance[in_bin]
        sd = values.std(ddof=1)
        z = np.zeros_like(values) if not np.isfinite(sd) or sd == 0 else (values - values.mean()) / sd
        hvg_score[valid_idx[in_bin]] = z

    return pd.DataFrame(
        {
            "cell_fraction": cell_fraction,
            "patient_count": patient_count,
            "lognorm_mean": mean,
            "lognorm_variance": variance,
            "hvg_score": hvg_score,
            "passes_detection_filter": valid,
        }
    )


def write_mtx_gz(matrix: sparse.spmatrix, path: Path) -> None:
    temporary = path.with_suffix("")
    mmwrite(temporary, matrix, field="integer", symmetry="general")
    with temporary.open("rb") as src, gzip.open(path, "wb", compresslevel=6) as dst:
        for chunk in iter(lambda: src.read(1024 * 1024), b""):
            dst.write(chunk)
    temporary.unlink()


programs = pd.read_csv(
    RESOURCE / "myeloid_receiver_programs_frozen.tsv", sep="\t", dtype=str
)
programs = programs.loc[programs["program"].isin(PRIMARY_PROGRAMS)].copy()
if set(programs["program"].unique()) != set(PRIMARY_PROGRAMS):
    raise RuntimeError("A frozen primary program is missing")

loaded: dict[str, dict[str, object]] = {}
reference_genes: pd.Series | None = None
for offset, (label, annotation) in enumerate(STATES.items(), start=1000):
    matrix, genes, metadata = read_state(label)
    if reference_genes is None:
        reference_genes = genes
    elif not reference_genes.equals(genes):
        raise RuntimeError("State matrices do not share the same ordered gene universe")
    selected = donor_capped_indices(metadata, offset)
    capped_matrix = matrix[:, selected].tocsc()
    capped_metadata = metadata.iloc[selected].reset_index(drop=True)
    metrics = gene_metrics(capped_matrix, capped_metadata)
    metrics.insert(0, "gene_symbol", genes.to_numpy())

    technical = metrics["gene_symbol"].str.match(
        re.compile(r"^(MT-|RPL\d|RPS\d|HB[ABDEGQZ])", flags=re.IGNORECASE)
    )
    eligible_hvg = metrics["passes_detection_filter"] & ~technical
    hvg_idx = (
        metrics.loc[eligible_hvg]
        .sort_values(["hvg_score", "gene_symbol"], ascending=[False, True])
        .head(N_HVG_PER_STATE)
        .index
    )
    hvg_genes = set(metrics.loc[hvg_idx, "gene_symbol"])
    if len(hvg_genes) != N_HVG_PER_STATE:
        raise RuntimeError(f"Insufficient HVGs for {label}")
    loaded[label] = {
        "matrix": capped_matrix,
        "genes": genes,
        "metadata": capped_metadata,
        "metrics": metrics,
        "hvg_genes": hvg_genes,
        "selected_original_indices": selected,
        "annotation": annotation,
    }

assert reference_genes is not None
all_gene_set = set(reference_genes)
missing_forced = sorted((set(AXIS_GENES) | set(MECHANISM_GUARD_GENES)) - all_gene_set)
if missing_forced:
    raise RuntimeError(f"Forced genes absent from expression matrix: {missing_forced}")

sender_metrics = loaded["sender_CTHRC1_CAF"]["metrics"].set_index("gene_symbol")
receiver_metrics = loaded["receiver_CD14_monocyte"]["metrics"].set_index("gene_symbol")
detected_either = (
    sender_metrics["passes_detection_filter"]
    | receiver_metrics["passes_detection_filter"]
)
eligible_program_genes = set(programs["gene_symbol"]) & set(
    detected_either.index[detected_either]
)

universe = (
    loaded["sender_CTHRC1_CAF"]["hvg_genes"]
    | loaded["receiver_CD14_monocyte"]["hvg_genes"]
    | eligible_program_genes
    | set(AXIS_GENES)
    | set(MECHANISM_GUARD_GENES)
)
universe = sorted(universe)

program_membership = (
    programs.groupby("gene_symbol")["program"]
    .apply(lambda x: ";".join(sorted(set(x))))
    .to_dict()
)
universe_table = pd.DataFrame({"gene_symbol": universe})
universe_table["sender_hvg"] = universe_table["gene_symbol"].isin(
    loaded["sender_CTHRC1_CAF"]["hvg_genes"]
)
universe_table["receiver_hvg"] = universe_table["gene_symbol"].isin(
    loaded["receiver_CD14_monocyte"]["hvg_genes"]
)
universe_table["primary_programs"] = universe_table["gene_symbol"].map(
    program_membership
).fillna("")
universe_table["forced_axis"] = universe_table["gene_symbol"].isin(AXIS_GENES)
universe_table["forced_mechanism_guard"] = universe_table["gene_symbol"].isin(
    MECHANISM_GUARD_GENES
)
for label in STATES:
    prefix = "sender" if label.startswith("sender") else "receiver"
    metrics = loaded[label]["metrics"].set_index("gene_symbol").loc[universe]
    for col in [
        "cell_fraction",
        "patient_count",
        "lognorm_mean",
        "lognorm_variance",
        "hvg_score",
        "passes_detection_filter",
    ]:
        universe_table[f"{prefix}_{col}"] = metrics[col].to_numpy()

universe_path = RESOURCE / "virtual_knockout_sanity_gene_universe.tsv"
universe_table.to_csv(universe_path, sep="\t", index=False)

filtered_programs = programs.loc[
    programs["gene_symbol"].isin(universe)
    & ~programs["gene_symbol"].isin(CIRCULARITY_EXCLUSION)
].copy()
filtered_program_path = RESOURCE / "virtual_knockout_sanity_programs_axis_excluded.tsv"
filtered_programs.to_csv(filtered_program_path, sep="\t", index=False)

state_audit: dict[str, object] = {}
universe_position = pd.Index(reference_genes).get_indexer(universe)
if np.any(universe_position < 0):
    raise RuntimeError("Universe contains genes absent from input matrix")
for label, annotation in STATES.items():
    item = loaded[label]
    selected_matrix = item["matrix"][universe_position, :].tocsc()
    prefix = OUT / f"GSE206785_{label}_sanity_counts"
    matrix_path = Path(str(prefix) + ".mtx.gz")
    genes_path = Path(str(prefix) + "_genes.tsv")
    metadata_path = Path(str(prefix) + "_metadata.tsv")
    write_mtx_gz(selected_matrix, matrix_path)
    pd.Series(universe, name="gene_symbol").to_csv(genes_path, sep="\t", index=False)
    item["metadata"].to_csv(metadata_path, sep="\t", index=False)
    patient_sizes = item["metadata"].groupby("Patient").size()
    state_audit[label] = {
        "annotation": annotation,
        "cells_before_cap": int(read_state(label)[0].shape[1]),
        "cells_after_cap": int(selected_matrix.shape[1]),
        "patients": int(item["metadata"]["Patient"].nunique()),
        "maximum_cells_per_patient_after_cap": int(patient_sizes.max()),
        "minimum_cells_per_patient_after_cap": int(patient_sizes.min()),
        "matrix_nonzero": int(selected_matrix.nnz),
        "matrix_sha256": sha256_file(matrix_path),
        "metadata_sha256": sha256_file(metadata_path),
    }

audit = {
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "dataset": "GSE206785",
    "purpose": "frozen compact GRN universe for scTenifoldKnk sanity calibration",
    "selection_is_outcome_blind": True,
    "seed": SEED,
    "donor_cap_unique_cells_without_replacement": DONOR_CAP,
    "hvg_per_state": N_HVG_PER_STATE,
    "hvg_method": "log1p(CP10K) variance z-scored within 20 expression-rank bins",
    "hvg_technical_exclusion_regex": "^(MT-|RPL\\d|RPS\\d|HB[ABDEGQZ])",
    "detection_filter": {
        "minimum_cell_fraction_within_donor_capped_state": MIN_CELL_FRACTION,
        "minimum_patients_detected": MIN_PATIENTS,
    },
    "primary_programs": PRIMARY_PROGRAMS,
    "axis_genes": AXIS_GENES,
    "mechanism_guard_genes": MECHANISM_GUARD_GENES,
    "program_circularity_exclusion": sorted(CIRCULARITY_EXCLUSION),
    "universe_genes": len(universe),
    "eligible_detected_program_genes": len(eligible_program_genes),
    "axis_excluded_program_rows": int(len(filtered_programs)),
    "axis_excluded_program_unique_genes": int(filtered_programs["gene_symbol"].nunique()),
    "universe_sha256": sha256_file(universe_path),
    "program_file_sha256": sha256_file(filtered_program_path),
    "states": state_audit,
    "warning": "This donor-capped layer is a low-cost sanity screen, not the final donor-equal bootstrap analysis.",
}
audit_path = AUDIT / "GSE206785_virtual_knockout_sanity_universe.json"
audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(audit, ensure_ascii=False, indent=2))
