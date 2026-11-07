"""Download per-model RNA-seq files for the GDSC2 gastric carcinoma models.

The Cell Model Passports API exposes model-level file URLs. Downloading only
the 28 gastric models avoids pulling the full ~1.2 GB merged RNA-seq archive.
"""

from pathlib import Path
import json
import time
import zipfile
import io

import pandas as pd
import requests


ROOT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026\analysis\drug_sensitivity")
RAW = ROOT / "raw_20260807"
EXPR = RAW / "cmp_expression_by_model"
EXPR.mkdir(parents=True, exist_ok=True)

GDSC = RAW / "GDSC2_fitted_dose_response_27Oct23.xlsx"
API = "https://api.cellmodelpassports.sanger.ac.uk"
session = requests.Session()
session.headers.update({"Accept": "application/vnd.api+json", "User-Agent": "Codex-GC-drug-sensitivity/1.0"})


def gastric_models():
    d = pd.read_excel(GDSC, usecols=["SANGER_MODEL_ID", "CELL_LINE_NAME", "CANCER_TYPE"])
    d = d[d["CANCER_TYPE"].eq("Gastric Carcinoma")].drop_duplicates().sort_values("SANGER_MODEL_ID")
    return d


def expression_url(model_id):
    r = session.get(f"{API}/models/{model_id}/files", timeout=60)
    r.raise_for_status()
    for entry in r.json().get("data", []):
        attrs = entry.get("attributes", {})
        meta = attrs.get("meta", {}) or {}
        url = attrs.get("url", "") or ""
        if meta.get("data_type") == "RNASeq" and url.endswith("_rnaseq.csv.zip"):
            return url
    return None


def main():
    models = gastric_models()
    manifest = []
    for row in models.itertuples(index=False):
        model_id = row.SANGER_MODEL_ID
        out = EXPR / f"{model_id}_rnaseq.csv"
        if out.exists() and out.stat().st_size > 1000:
            manifest.append({"model_id": model_id, "cell_line": row.CELL_LINE_NAME,
                             "url": "cached", "file": str(out), "status": "cached"})
            continue
        url = expression_url(model_id)
        if not url:
            manifest.append({"model_id": model_id, "cell_line": row.CELL_LINE_NAME,
                             "url": None, "file": None, "status": "missing"})
            continue
        r = session.get(url, timeout=120)
        r.raise_for_status()
        with zipfile.ZipFile(io.BytesIO(r.content)) as z:
            csvs = [n for n in z.namelist() if n.lower().endswith(".csv")]
            if len(csvs) != 1:
                raise RuntimeError(f"Unexpected archive members for {model_id}: {z.namelist()}")
            with z.open(csvs[0]) as src, open(out, "wb") as dst:
                dst.write(src.read())
        manifest.append({"model_id": model_id, "cell_line": row.CELL_LINE_NAME,
                         "url": url, "file": str(out), "status": "downloaded"})
        time.sleep(.2)
    (EXPR / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(pd.DataFrame(manifest).to_string(index=False))


if __name__ == "__main__":
    main()
