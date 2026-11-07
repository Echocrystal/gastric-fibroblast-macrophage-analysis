"""Reproduce cached spatial endpoint aggregation for standalone R2 review.

Run from any directory with Python, numpy, pandas, matplotlib and threadpoolctl.
This does not rerun raw spot-level models or modify source analysis outputs.
On the first run it freezes the listed source inputs in inputs/ and scripts/.
Subsequent runs use those local frozen copies.
"""
from __future__ import annotations

import hashlib
import json
import platform
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import threadpoolctl
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[3]
SOURCE = PROJECT / "analysis/results/spatial/SPON2_ITGB2_candidate_context_v1"
SEED = 2026082654
DRAWS = 9999
DATASETS = ["GSE270678_main_ST", "GSE270678_ST_R_sensitivity", "GSE251950_primary"]
LABELS = ["GSE270678 main ST", "GSE270678 ST_R", "GSE251950 primary"]
ENDPOINTS = ["SPON2_ITGB2_receiver_composite_r1", "SPON2_ITGB2_r1"]
ENDPOINT_LABELS = {ENDPOINTS[0]: "ITGB2 + receiver composite", ENDPOINTS[1]: "ITGB2 alone"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def freeze_inputs() -> list[dict]:
    records = []
    files = {f"inputs/{name}": SOURCE / name for name in [
        "SPON2_ITGB2_section_effects.tsv", "SPON2_ITGB2_donor_effects.tsv",
        "SPON2_ITGB2_meta_summary.tsv", "SPON2_ITGB2_gene_endpoint_coverage.tsv",
        "SPON2_ITGB2_candidate_axis_calibration.tsv", "SPON2_ITGB2_analysis_audit.json",
    ]}
    files.update({
        "scripts/original_54_spatial_candidate_context.py": PROJECT / "scripts/54_analyze_SPON2_ITGB2_spatial_candidate_context.py",
        "scripts/original_35_spatial_chain.py": PROJECT / "scripts/35_analyze_GSE270678_PLAT_spatial_chain.py",
        "inputs/original_analysis_lock.md": PROJECT / "analysis/design/SPON2_ITGB2_SPATIAL_CANDIDATE_CONTEXT_LOCK_20260826.md",
    })
    for relative, source in files.items():
        dest = ROOT / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            shutil.copy2(source, dest)
        records.append({"source_path": str(source), "frozen_path": relative,
                        "bytes": dest.stat().st_size, "sha256": sha256(dest),
                        "matches_current_source": sha256(dest) == sha256(source) if source.exists() else None})
    pd.DataFrame(records).to_csv(ROOT / "input_provenance_sha256.tsv", sep="\t", index=False)
    return records


def fisher_pool(values) -> float:
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    return float(np.tanh(np.arctanh(np.clip(values, -0.999999, 0.999999)).mean()))


def exact_signflip(values) -> float:
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    z = np.arctanh(np.clip(values, -0.999999, 0.999999))
    observed = abs(float(z.mean()))
    total = 1 << len(z)
    exceeds = 0
    bits = np.arange(len(z), dtype=np.uint64)
    for first in range(0, total, 65536):
        masks = np.arange(first, min(first + 65536, total), dtype=np.uint64)[:, None]
        signs = (((masks >> bits) & 1).astype(float) * 2.0) - 1.0
        exceeds += int(np.sum(np.abs(signs @ z / len(z)) >= observed - 1e-15))
    return exceeds / total


def bh_adjust(values):
    values = np.asarray(values, float)
    order = np.argsort(values)
    q = values[order] * len(values) / np.arange(1, len(values) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    result = np.empty_like(q)
    result[order] = np.minimum(q, 1.0)
    return result


def aggregate(section):
    # Keep every original endpoint and original group iteration order. Bootstrap
    # RNG is sequential in the original code; omitting groups changes the CIs.
    accepted = section.loc[section.section_gate & (section.n_spots >= 100)].copy()
    donor = accepted.groupby(
        ["dataset", "donor_id", "endpoint", "category", "axis", "ring"], as_index=False
    ).agg(rho=("rho", fisher_pool), n_sections=("section_id", "nunique"), median_spots=("n_spots", "median"))
    rng = np.random.default_rng(SEED)
    rows = []
    for keys, d in donor.groupby(["dataset", "endpoint", "category", "axis", "ring"], sort=False, dropna=False):
        x = d.rho.to_numpy(float)
        z = np.arctanh(np.clip(x, -0.999999, 0.999999))
        sampled = rng.choice(z, size=(DRAWS, len(z)), replace=True)
        lo, hi = np.quantile(np.tanh(sampled.mean(axis=1)), [0.025, 0.975])
        rows.append(dict(zip(["dataset", "endpoint", "category", "axis", "ring"], keys)) | {
            "n_donors": len(x), "pooled_rho": fisher_pool(x), "ci_low": lo, "ci_high": hi,
            "positive_donors": int((x > 0).sum()), "negative_donors": int((x < 0).sum()),
            "exact_signflip_p": exact_signflip(x),
        })
    meta = pd.DataFrame(rows)
    meta["candidate_BH_FDR"] = np.nan
    meta["supporting_BH_FDR"] = np.nan
    for dataset, d in meta.loc[meta.category.eq("candidate_calibration")].groupby("dataset"):
        meta.loc[d.index, "candidate_BH_FDR"] = bh_adjust(d.exact_signflip_p)
        meta.loc[d.index, "candidate_effect_rank"] = d.pooled_rho.rank(method="min", ascending=False)
    for keys, d in meta.loc[meta.category.isin(["supporting", "receiver_function"])].groupby(["dataset", "category"]):
        meta.loc[d.index, "supporting_BH_FDR"] = bh_adjust(d.exact_signflip_p)
    return donor, meta, accepted


def verify(donor, meta):
    old_d = pd.read_csv(ROOT / "inputs/SPON2_ITGB2_donor_effects.tsv", sep="\t", keep_default_na=False)
    old_m = pd.read_csv(ROOT / "inputs/SPON2_ITGB2_meta_summary.tsv", sep="\t", keep_default_na=False)
    d_keys = ["dataset", "donor_id", "endpoint", "category", "axis", "ring"]
    m_keys = ["dataset", "endpoint", "category", "axis", "ring"]
    d = donor.merge(old_d, on=d_keys, suffixes=("_rerun", "_source"), validate="one_to_one")
    m = meta.merge(old_m, on=m_keys, suffixes=("_rerun", "_source"), validate="one_to_one")
    assert len(d) == len(donor) == len(old_d)
    assert len(m) == len(meta) == len(old_m)
    checks = {"donor_rows": len(d), "meta_rows": len(m), "absolute_tolerance": 1e-12}
    for name in ["rho", "n_sections", "median_spots"]:
        a, b = d[f"{name}_rerun"].astype(float), d[f"{name}_source"].astype(float)
        checks[f"donor_{name}_max_abs_difference"] = float(np.max(np.abs(a-b)))
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-12)
    for name in ["n_donors", "pooled_rho", "ci_low", "ci_high", "positive_donors", "negative_donors", "exact_signflip_p", "candidate_BH_FDR", "supporting_BH_FDR", "candidate_effect_rank"]:
        a = pd.to_numeric(m[f"{name}_rerun"], errors="coerce").to_numpy(float)
        b = pd.to_numeric(m[f"{name}_source"], errors="coerce").to_numpy(float)
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-12, equal_nan=True)
        delta = np.abs(a-b)
        checks[f"meta_{name}_max_abs_difference"] = float(np.nanmax(delta)) if np.isfinite(delta).any() else 0.0
    checks["all_checks_passed"] = True
    return checks


def compare(donor, meta, accepted):
    rows = []
    for dataset in DATASETS:
        for endpoint in ENDPOINTS:
            x = meta.loc[meta.dataset.eq(dataset) & meta.endpoint.eq(endpoint)].iloc[0]
            if endpoint == ENDPOINTS[0]:
                qrow = meta.loc[meta.dataset.eq(dataset) & meta.endpoint.eq("candidate_SPON2_ITGB2_r1")].iloc[0]
                assert abs(float(qrow.pooled_rho) - float(x.pooled_rho)) < 1e-14
                q, q_family, q_tests = qrow.candidate_BH_FDR, "Four receptor-receiver composite candidate axes within this dataset", 4
            else:
                q, q_family, q_tests = x.supporting_BH_FDR, "Three supporting endpoints within this dataset", 3
            passed = accepted.loc[accepted.dataset.eq(dataset) & accepted.endpoint.eq(endpoint)]
            rows.append({"dataset": dataset,
                "cohort_role": "Primary spatial cohort" if dataset == DATASETS[0] else "Same-donor repeated-section sensitivity" if dataset == DATASETS[1] else "Independent secondary spatial projection",
                "endpoint": ENDPOINT_LABELS[endpoint], "source_endpoint_id": endpoint,
                "n_sections": passed.section_id.nunique(), "n_donors_or_patients": int(x.n_donors),
                "pooled_adjusted_rho": float(x.pooled_rho), "ci95_low": float(x.ci_low), "ci95_high": float(x.ci_high),
                "positive_donors_or_patients": int(x.positive_donors), "negative_donors_or_patients": int(x.negative_donors),
                "nominal_exact_signflip_p": float(x.exact_signflip_p),
                "original_within_dataset_BH_q": float(q), "original_BH_family": q_family, "original_BH_number_of_tests": q_tests,
                "bootstrap_draws": DRAWS, "bootstrap_seed": SEED,
                "bootstrap_policy": "Original sequential RNG across all endpoint groups; primary-distance-profile CI for composite",
            })
    table = pd.DataFrame(rows)
    table["new_exploratory_six_endpoint_BH_q"] = bh_adjust(table.nominal_exact_signflip_p)
    table["new_exploratory_BH_family"] = "Six displayed tests: two endpoints x three datasets; new descriptive review family"
    return table


def forest(table):
    plt.rcParams.update({"font.family": "Arial", "font.size": 10, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, ax = plt.subplots(figsize=(10.8, 4.5))
    colours = ["#B45D35", "#286B8A"]
    for i, (dataset, label) in enumerate(zip(DATASETS, LABELS)):
        for j, endpoint in enumerate(ENDPOINTS):
            r = table.loc[table.dataset.eq(dataset) & table.source_endpoint_id.eq(endpoint)].iloc[0]
            y = 2 - i + (0.12 if j == 0 else -0.12)
            ax.errorbar(r.pooled_adjusted_rho, y,
                        xerr=[[r.pooled_adjusted_rho-r.ci95_low], [r.ci95_high-r.pooled_adjusted_rho]],
                        fmt="o" if j == 0 else "s", color=colours[j], markersize=6,
                        capsize=3, linewidth=1.6, label=ENDPOINT_LABELS[endpoint] if i == 0 else None)
            p = r.nominal_exact_signflip_p
            p_label = f"{p:.2g}" if p < 0.001 else f"{p:.3f}"
            text = f"{r.pooled_adjusted_rho:.4f} [{r.ci95_low:.4f}, {r.ci95_high:.4f}]   P = {p_label}"
            ax.text(0.078, y, text, va="center", fontsize=9)
    ax.axvline(0, color="#707070", linestyle="--", linewidth=1)
    ax.set_yticks([2, 1, 0], ["GSE270678 main ST\n19 donors / 33 sections", "GSE270678 ST_R\n15 donors / 36 sections", "GSE251950 primary\n9 patients / 9 sections"])
    ax.set_ylim(-0.45, 2.5)
    ax.set_xlim(-0.009, 0.077)
    ax.set_xlabel("Donor-equal pooled adjusted Spearman correlation (95% bootstrap CI)")
    ax.set_title("SPON2 versus neighbouring composite and single-gene endpoints", loc="left", fontsize=12, pad=22)
    ax.text(0.078, 2.48, "Pooled rho [95% CI]; nominal exact P", fontsize=9, fontweight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=12)
    ax.legend(loc="lower right", frameon=False, bbox_to_anchor=(1.7, -0.25), ncol=2, fontsize=9)
    fig.subplots_adjust(left=0.2, right=0.65, bottom=0.22, top=0.83)
    fig.text(0.02, 0.025, "ST_R overlaps main-cohort donors. Exact P values and original correction families are listed in the comparison table.", fontsize=8, color="#444444")
    fig.savefig(ROOT / "spatial_endpoint_comparison.png", dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(ROOT / "spatial_endpoint_comparison.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_readme(table, checks):
    footnote = (
        "SPON2 expression at each measured Visium spot was related to the exact ring-1 mean of either "
        "an ITGB2-receiver composite or ITGB2 alone under identical rank residualisation. "
        "The spot-level composite was [z(ITGB2) + z(receiver score)]/2; both terms were standardised "
        "within section before spatial averaging. Covariates were log1p total UMI, same-spot broad fibroblast "
        "and epithelial scores, ring-1 broad myeloid score, same-spot frozen ECM-direction score, fixed sender "
        "score and PLAT expression. Section correlations were Fisher-z averaged within donor and then equally "
        "across donors/patients. Intervals are percentile donor-bootstrap 95% CIs (9,999 draws); nominal P values "
        "are exhaustive two-sided donor sign-flip tests. Original BH q values retain their separate families: "
        "four candidate receptor-receiver composites for the composite, and three supporting endpoints for "
        "ITGB2 alone, separately within each dataset. ITGB2-alone q values do not derive from the four-axis "
        "candidate family. The additional six-test BH column is a new exploratory family for this standalone "
        "comparison and does not replace original inference. GSE270678 ST_R comprises repeated/extended "
        "sections from overlapping donors and is not independent validation. Section eligibility required "
        "SPON2 and ITGB2 detection in at least 1% of spots and at least 100 complete endpoint/covariate spots. "
        "These mixed-cell associations cannot establish direct binding, receptor activation or causality."
    )
    (ROOT / "English_table_footnote.txt").write_text(footnote + "\n", encoding="utf-8")
    lines = ["# Standalone spatial endpoint comparison for Reviewer 2", "",
             "This review package recomputes aggregation and inference from frozen existing section-level coefficients. "
             "It does not re-fit spot-level models from raw spatial matrices and does not alter the manuscript or original tables.", "",
             "## Execution", "", "`python run_spatial_endpoint_comparison.py`", "",
             "Dependencies: Python 3, numpy, pandas, matplotlib and threadpoolctl. Paths resolve relative to this package. "
             "Source inputs and original analysis scripts are snapshotted locally on first execution; subsequent runs use the frozen copies.", "",
             "## Formula and analysis", "", footnote, "",
             "The original composite and candidate-calibration rows represent the same coefficient and P value. "
             "Their bootstrap intervals differ slightly because they were generated at different positions in the original sequential RNG stream. "
             "The comparison retains the original primary ring-1 distance-profile CI and attaches the identical-coefficient candidate-family q value. "
             "All 57 original meta-analysis groups are reconstructed in the original order with seed 2026082654 so that CI values are reproduced.", "",
             "## Verification", "",
             f"Recomputed {checks['donor_rows']} donor-endpoint rows and {checks['meta_rows']} meta-analysis rows from 1,501 section-endpoint rows. "
             "All original donor estimates, pooled estimates, CIs, exact P values and original-family q values match within absolute tolerance 1e-12. "
             "See verification.json and input_provenance_sha256.tsv.", "",
             "## Review results", "",
             "| Dataset | Endpoint | n | Adjusted rho | 95% CI | Exact P | Original-family BH q |", "|---|---|---:|---:|---|---:|---:|"]
    for r in table.itertuples(index=False):
        lines.append(f"| {r.dataset} | {r.endpoint} | {r.n_donors_or_patients} | {r.pooled_adjusted_rho:.6f} | {r.ci95_low:.6f} to {r.ci95_high:.6f} | {r.nominal_exact_signflip_p:.8g} | {r.original_within_dataset_BH_q:.8g} |")
    lines += ["", "The single-gene endpoint is not significant in the main cohort or ST_R, and is small and positive in GSE251950. "
              "The composite therefore has stronger statistical support than ITGB2 alone in the primary cohort, but this comparison does not quantify "
              "the percentage contribution of the receiver score or prove that it causes the composite result.", "",
              "## Output files", "", "- spatial_endpoint_comparison.tsv / .json: six-row comparison with separately labelled original and new exploratory q values.",
              "- spatial_endpoint_donor_results.tsv: donor/patient values for the displayed two endpoints.",
              "- spatial_endpoint_section_results.tsv: section values and original eligibility flags for those endpoints.",
              "- all_recomputed_meta_results.tsv: all original endpoint groups used to reproduce RNG order and original q families.",
              "- spatial_endpoint_comparison.png / .pdf: forest plot.",
              "- verification.json: exact numerical checks; input_provenance_sha256.tsv: source snapshots and hashes.", ""]
    (ROOT / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    provenance = freeze_inputs()
    section = pd.read_csv(ROOT / "inputs/SPON2_ITGB2_section_effects.tsv", sep="\t", keep_default_na=False)
    section["rho"] = pd.to_numeric(section.rho, errors="coerce")
    if section.section_gate.dtype != bool:
        section["section_gate"] = section.section_gate.astype(str).str.lower().eq("true")
    with threadpool_limits(limits=1):
        donor, meta, accepted = aggregate(section)
    checks = verify(donor, meta)
    table = compare(donor, meta, accepted)
    table.to_csv(ROOT / "spatial_endpoint_comparison.tsv", sep="\t", index=False, float_format="%.17g")
    presentation = pd.DataFrame({
        "Dataset": table.dataset.map(dict(zip(DATASETS, LABELS))),
        "Endpoint": table.endpoint,
        "Sections (n)": table.n_sections,
        "Donors or patients (n)": table.n_donors_or_patients,
        "Adjusted rho": table.pooled_adjusted_rho.map(lambda x: f"{x:.6f}"),
        "95% bootstrap CI": table.apply(lambda x: f"{x.ci95_low:.6f} to {x.ci95_high:.6f}", axis=1),
        "Nominal exact P": table.nominal_exact_signflip_p.map(lambda x: f"{x:.8g}"),
        "Original BH q": table.original_within_dataset_BH_q.map(lambda x: f"{x:.8g}"),
        "Original correction family": table.original_BH_family,
    })
    presentation.to_csv(ROOT / "spatial_endpoint_presentation.tsv", sep="\t", index=False)
    (ROOT / "spatial_endpoint_comparison.json").write_text(table.to_json(orient="records", indent=2, double_precision=15), encoding="utf-8")
    donor.loc[donor.endpoint.isin(ENDPOINTS)].to_csv(ROOT / "spatial_endpoint_donor_results.tsv", sep="\t", index=False, float_format="%.17g")
    section.loc[section.endpoint.isin(ENDPOINTS)].to_csv(ROOT / "spatial_endpoint_section_results.tsv", sep="\t", index=False, float_format="%.17g")
    meta.to_csv(ROOT / "all_recomputed_meta_results.tsv", sep="\t", index=False, float_format="%.17g")
    checks.update({"seed": SEED, "bootstrap_draws": DRAWS, "python": platform.python_version(),
                   "numpy": np.__version__, "pandas": pd.__version__, "matplotlib": matplotlib.__version__,
                   "threadpoolctl": threadpoolctl.__version__,
                   "runner_sha256": sha256(Path(__file__)), "raw_spot_level_refit": False,
                   "input_snapshot_count": len(provenance)})
    (ROOT / "verification.json").write_text(json.dumps(checks, indent=2), encoding="utf-8")
    (ROOT / "requirements.txt").write_text("\n".join([
        f"numpy=={np.__version__}", f"pandas=={pd.__version__}",
        f"matplotlib=={matplotlib.__version__}", f"threadpoolctl=={threadpoolctl.__version__}", ""
    ]), encoding="utf-8")
    write_readme(table, checks)
    forest(table)
    print(table[["dataset", "endpoint", "n_donors_or_patients", "pooled_adjusted_rho", "ci95_low", "ci95_high", "nominal_exact_signflip_p", "original_within_dataset_BH_q", "new_exploratory_six_endpoint_BH_q"]].to_string(index=False))
    print("All source numerical verification checks passed.")


if __name__ == "__main__":
    main()
