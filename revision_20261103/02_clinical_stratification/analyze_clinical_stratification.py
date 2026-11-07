"""Build a descriptive Kim2018 biomarker-stratified response panel.

Run with the bundled Python interpreter. All defaults resolve from this script,
so the command is independent of the shell's working directory. This script
does not refit signatures, configurations, or inferential response models.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path


GROUPS = ["Dual-low", "Fibroblast-dominant", "Macrophage-dominant", "Dual-high"]
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[3]
UPSTREAM_INPUT = PROJECT / "figures/source_data/Figure6_final/Figure6g_Kim2018_patient_data.tsv"
SAVED_INPUT = HERE / "Kim2018_patient_input_snapshot.tsv"
DEFAULT_INPUT = SAVED_INPUT if SAVED_INPUT.is_file() else UPSTREAM_INPUT


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def quantile_linear(values: list[float], probability: float) -> float | None:
    """Linear sample quantile, equivalent to R type 7 / NumPy method='linear'."""
    if not values:
        return None
    x = sorted(values)
    position = (len(x) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    return x[lower] + (position - lower) * (x[upper] - x[lower])


def format_number(value: float | None) -> str:
    return "NA" if value is None else f"{value:.2f}".rstrip("0").rstrip(".")


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fields or list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    source = args.input.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    existing_manifest_path = output / "clinical_stratification_manifest.json"
    existing_manifest = json.loads(existing_manifest_path.read_text(encoding="utf-8")) if existing_manifest_path.is_file() else {}
    original_upstream_source = existing_manifest.get("original_upstream_source_path", existing_manifest.get("input_source", str(UPSTREAM_INPUT)))
    original_upstream_sha256 = existing_manifest.get("original_upstream_source_sha256", existing_manifest.get("input_source_sha256", sha256(source)))
    with source.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        original_fields = list(reader.fieldnames or [])
        patients = list(reader)
    required = {"patient_id", "ecotype4", "response_binary", "msi_type_pentaplex", "ebv_in_situ", "pdl1_cps"}
    if not required.issubset(original_fields):
        raise ValueError(f"Missing columns: {required - set(original_fields)}")
    if len({p["patient_id"] for p in patients}) != len(patients):
        raise ValueError("Duplicate patient identifiers")
    if set(p["ecotype4"] for p in patients) != set(GROUPS):
        raise ValueError("Unexpected configuration labels")

    for p in patients:
        response = p["response_binary"].strip()
        msi = p["msi_type_pentaplex"].strip()
        ebv = p["ebv_in_situ"].strip()
        cps = p["pdl1_cps"].strip()
        if response not in {"0", "1"}:
            raise ValueError(f"Invalid response: {p['patient_id']} {response}")
        if msi not in {"", "MSI", "MSS"} or ebv not in {"", "Positive", "Negative"}:
            raise ValueError(f"Unmapped biomarker category: {p['patient_id']}")
        if cps and (not math.isfinite(float(cps)) or not 0 <= float(cps) <= 100):
            raise ValueError(f"Invalid CPS: {p['patient_id']} {cps}")
        p["msi_positive_flag"] = int(msi == "MSI")
        p["msi_negative_flag"] = int(msi == "MSS")
        p["msi_missing_flag"] = int(msi == "")
        p["ebv_positive_flag"] = int(ebv == "Positive")
        p["ebv_negative_flag"] = int(ebv == "Negative")
        p["ebv_missing_flag"] = int(ebv == "")
        p["cps_available_flag"] = int(cps != "")
        p["mss_and_ebv_negative_flag"] = int(msi == "MSS" and ebv == "Negative")

    derived_fields = [x for x in patients[0] if x not in original_fields]
    snapshot = output / "Kim2018_patient_input_snapshot.tsv"
    if source != snapshot.resolve():
        shutil.copyfile(source, snapshot)
    write_tsv(output / "Kim2018_patients_with_restriction_flags.tsv", patients, original_fields + derived_fields)

    restrictions = {
        "All patients": lambda p: True,
        "MSS only": lambda p: bool(p["msi_negative_flag"]),
        "EBV ISH-negative only": lambda p: bool(p["ebv_negative_flag"]),
        "MSS and EBV ISH-negative": lambda p: bool(p["mss_and_ebv_negative_flag"]),
    }
    restriction_rows = []
    for name, condition in restrictions.items():
        for group in GROUPS:
            selected = [p for p in patients if p["ecotype4"] == group and condition(p)]
            responders = sum(int(p["response_binary"]) for p in selected)
            restriction_rows.append({
                "restriction": name,
                "configuration": group,
                "n": len(selected),
                "responders_n": responders,
                "nonresponders_n": len(selected) - responders,
                "response_fraction": responders / len(selected) if selected else None,
            })
    write_tsv(output / "Kim2018_response_by_restriction.tsv", restriction_rows)

    summary = []
    for group in GROUPS:
        selected = [p for p in patients if p["ecotype4"] == group]
        cps = [float(p["pdl1_cps"]) for p in selected if p["cps_available_flag"]]
        restricted = [p for p in selected if p["mss_and_ebv_negative_flag"]]
        row = {
            "configuration": group,
            "n": len(selected),
            "responders_n": sum(int(p["response_binary"]) for p in selected),
            "nonresponders_n": sum(1 - int(p["response_binary"]) for p in selected),
            "MSI_positive_n": sum(p["msi_positive_flag"] for p in selected),
            "MSI_negative_MSS_n": sum(p["msi_negative_flag"] for p in selected),
            "MSI_missing_n": sum(p["msi_missing_flag"] for p in selected),
            "EBV_ISH_positive_n": sum(p["ebv_positive_flag"] for p in selected),
            "EBV_ISH_negative_n": sum(p["ebv_negative_flag"] for p in selected),
            "EBV_ISH_missing_n": sum(p["ebv_missing_flag"] for p in selected),
            "PDL1_CPS_available_n": len(cps),
            "PDL1_CPS_missing_n": len(selected) - len(cps),
            "PDL1_CPS_median": quantile_linear(cps, 0.5),
            "PDL1_CPS_q1": quantile_linear(cps, 0.25),
            "PDL1_CPS_q3": quantile_linear(cps, 0.75),
            "MSS_EBVnegative_n": len(restricted),
            "MSS_EBVnegative_responders_n": sum(int(p["response_binary"]) for p in restricted),
            "MSS_EBVnegative_nonresponders_n": sum(1 - int(p["response_binary"]) for p in restricted),
        }
        if row["MSI_positive_n"] + row["MSI_negative_MSS_n"] + row["MSI_missing_n"] != row["n"]:
            raise AssertionError("MSI counts do not sum")
        if row["EBV_ISH_positive_n"] + row["EBV_ISH_negative_n"] + row["EBV_ISH_missing_n"] != row["n"]:
            raise AssertionError("EBV counts do not sum")
        summary.append(row)
    write_tsv(output / "Kim2018_configuration_biomarker_numeric_summary.tsv", summary)

    metrics = [
        ("Patients, n", lambda r: str(r["n"])),
        ("Responders, n", lambda r: str(r["responders_n"])),
        ("Nonresponders, n", lambda r: str(r["nonresponders_n"])),
        ("MSI-positive, n", lambda r: str(r["MSI_positive_n"])),
        ("MSI-negative (MSS), n", lambda r: str(r["MSI_negative_MSS_n"])),
        ("MSI status missing, n", lambda r: str(r["MSI_missing_n"])),
        ("EBV ISH-positive, n", lambda r: str(r["EBV_ISH_positive_n"])),
        ("EBV ISH-negative, n", lambda r: str(r["EBV_ISH_negative_n"])),
        ("EBV ISH status missing, n", lambda r: str(r["EBV_ISH_missing_n"])),
        ("PD-L1 CPS available, n", lambda r: str(r["PDL1_CPS_available_n"])),
        ("PD-L1 CPS missing, n", lambda r: str(r["PDL1_CPS_missing_n"])),
        ("PD-L1 CPS, median (Q1-Q3)", lambda r: f"{format_number(r['PDL1_CPS_median'])} ({format_number(r['PDL1_CPS_q1'])}-{format_number(r['PDL1_CPS_q3'])})"),
        ("MSS and EBV ISH-negative patients, n", lambda r: str(r["MSS_EBVnegative_n"])),
        ("Responders among MSS and EBV ISH-negative patients, n", lambda r: str(r["MSS_EBVnegative_responders_n"])),
        ("Response among MSS and EBV ISH-negative patients, n/N", lambda r: f"{r['MSS_EBVnegative_responders_n']}/{r['MSS_EBVnegative_n']}"),
    ]
    panel = [{"Metric": metric, **{r["configuration"]: formatter(r) for r in summary}} for metric, formatter in metrics]
    write_tsv(output / "Supplemental_Table6_biomarker_stratified_response_panel.tsv", panel)

    checks = {
        "n_45": len(patients) == 45,
        "responders_12": sum(int(p["response_binary"]) for p in patients) == 12,
        "patient_ids_unique": len({p["patient_id"] for p in patients}) == len(patients),
        "MSI_complete": sum(r["MSI_missing_n"] for r in summary) == 0,
        "EBV_ISH_complete": sum(r["EBV_ISH_missing_n"] for r in summary) == 0,
        "CPS_available_41_missing_4": sum(r["PDL1_CPS_available_n"] for r in summary) == 41 and sum(r["PDL1_CPS_missing_n"] for r in summary) == 4,
    }
    by_group = {r["configuration"]: r for r in summary}
    checks["restricted_macrophage_dominant_1_of_2"] = (by_group["Macrophage-dominant"]["MSS_EBVnegative_responders_n"], by_group["Macrophage-dominant"]["MSS_EBVnegative_n"]) == (1, 2)
    checks["restricted_dual_high_0_of_14"] = (by_group["Dual-high"]["MSS_EBVnegative_responders_n"], by_group["Dual-high"]["MSS_EBVnegative_n"]) == (0, 14)
    if not all(checks.values()):
        raise AssertionError(checks)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis": "Descriptive biomarker-stratified tabulation of frozen Kim2018 configurations",
        "input_source": str(source),
        "input_source_sha256": sha256(source),
        "original_upstream_source_path": original_upstream_source,
        "original_upstream_source_sha256": original_upstream_sha256,
        "input_default_policy": "Prefer the saved input snapshot beside the script; otherwise use the original project data path",
        "input_snapshot": snapshot.name,
        "input_snapshot_sha256": sha256(snapshot),
        "script": Path(__file__).name,
        "script_sha256": sha256(Path(__file__)),
        "runtime": {"python": sys.version, "platform": platform.platform(), "third_party_dependencies": []},
        "definition": {
            "MSI_positive": "msi_type_pentaplex == 'MSI'",
            "MSI_negative": "msi_type_pentaplex == 'MSS' (pentaplex assay, not WES subtype)",
            "EBV_positive": "ebv_in_situ == 'Positive'",
            "EBV_negative": "ebv_in_situ == 'Negative' (ISH assay, not WES estimate)",
            "restriction": "MSS AND EBV ISH-negative; missing status is never classified as negative",
            "response": "response_binary as supplied in the frozen patient table",
            "configuration": "ecotype4 as supplied in the frozen patient table; no thresholds re-estimated",
            "CPS_quartiles": "Linear interpolation, R type 7 / NumPy method='linear'; available cases only",
        },
        "inference": {"tests_performed": 0, "correction_family": "Not applicable: descriptive counts and quantiles only", "random_seed": "Not applicable: deterministic tabulation"},
        "validation": checks,
        "outputs": {},
    }
    for name in ["Kim2018_patients_with_restriction_flags.tsv", "Kim2018_response_by_restriction.tsv", "Kim2018_configuration_biomarker_numeric_summary.tsv", "Supplemental_Table6_biomarker_stratified_response_panel.tsv"]:
        manifest["outputs"][name] = {"sha256": sha256(output / name)}
    (output / "clinical_stratification_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    footnote = (
        "Configurations and response labels were taken from the frozen Kim2018 patient-level source data; "
        "no signatures or configuration thresholds were re-estimated. MSI status was defined by the pentaplex assay "
        "(MSI versus MSS), and EBV status by in situ hybridisation. The restricted analysis included only patients "
        "classified as both MSS and EBV ISH-negative. Missing biomarker values were reported separately and were "
        "not treated as negative. PD-L1 CPS is summarised among patients with available measurements as the median "
        "and first and third quartiles, calculated by linear interpolation (R type 7). Counts and response fractions "
        "are descriptive; no additional hypothesis tests were performed. The restricted macrophage-dominant group "
        "contained only two patients, including one responder, and does not establish a biomarker-independent "
        "association or treatment benefit. CPS, combined positive score; EBV ISH, Epstein-Barr virus in situ "
        "hybridisation; MSI, microsatellite instability; MSS, microsatellite stable."
    )
    (output / "Table6_panel_footnote_EN.txt").write_text(footnote + "\n", encoding="utf-8")
    readme = [
        "# Kim2018 biomarker-stratified response panel (standalone review analysis)",
        "",
        "This deterministic tabulation uses the existing frozen patient-level data. Formal manuscripts and original supplemental tables were not modified.",
        "",
        "Run: `python analyze_clinical_stratification.py` (defaults resolve from the script location). The saved input snapshot is preferred when present; otherwise the original project data path is used. The manifest separately preserves the original upstream source path and hash. No third-party packages are required.",
        "",
        "Verified source totals: 45 patients, 12 responders; MSI and EBV ISH are complete; CPS is available for 41 patients and missing for 4.",
        "",
        "| Configuration | Original responders/patients | MSI-positive | EBV ISH-positive | MSS and EBV-negative responders/patients | CPS available | CPS median (Q1-Q3) |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for r in summary:
        readme.append(f"| {r['configuration']} | {r['responders_n']}/{r['n']} | {r['MSI_positive_n']} | {r['EBV_ISH_positive_n']} | {r['MSS_EBVnegative_responders_n']}/{r['MSS_EBVnegative_n']} | {r['PDL1_CPS_available_n']}/{r['n']} | {format_number(r['PDL1_CPS_median'])} ({format_number(r['PDL1_CPS_q1'])}-{format_number(r['PDL1_CPS_q3'])}) |")
    readme += [
        "",
        "The dual-high restricted group still contains 14 patients with zero responders; zero is the response count, not the group size. The macrophage-dominant restricted patients are PB-16-005 (responder) and PB-16-060 (nonresponder).",
        "",
        "No new inferential test or response model was fitted. The small restricted group is descriptive and cannot establish biomarker-independent association. No multiplicity family or random seed applies to this deterministic tabulation.",
        "",
        "Files:",
        "- `Supplemental_Table6_biomarker_stratified_response_panel.tsv`: publication-friendly metric rows across the four configurations.",
        "- `Kim2018_configuration_biomarker_numeric_summary.tsv`: numeric summary, one row per configuration.",
        "- `Kim2018_response_by_restriction.tsv`: original, MSS-only, EBV-negative-only and joint-restricted counts and response fractions.",
        "- `Kim2018_patients_with_restriction_flags.tsv`: all original columns plus explicit biomarker/restriction flags.",
        "- `Kim2018_patient_input_snapshot.tsv`: unchanged input snapshot.",
        "- `clinical_stratification_manifest.json`: source path, hashes, definitions, runtime and validation.",
        "- `Table6_panel_footnote_EN.txt`: English table footnote.",
        "",
        "English table footnote:",
        "",
        footnote,
    ]
    (output / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(output), "validation": checks, "summary": summary}, indent=2))


if __name__ == "__main__":
    main()
