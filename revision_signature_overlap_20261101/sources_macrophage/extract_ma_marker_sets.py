"""Extract source-derived marker sets from the unmodified Ma2025 Table2.xlsx.

Uses only the Python standard library.  The study's Methods explicitly state
adjusted P <= 0.05 and average log fold change >= 0.15.  The file supplies a
longer differential-expression table, so we apply these author-stated criteria.
No new cell-level analysis, signature fitting, or ranking by overlap is done.
"""
from pathlib import Path
import csv
import json
import re
import xml.etree.ElementTree as ET
import zipfile

HERE = Path(__file__).resolve().parent
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
XLREL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"


def read_sheet(path, sheet_name):
    with zipfile.ZipFile(path) as z:
        strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            strings = [
                "".join(x.itertext())
                for x in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS)
            ]
        relations = {
            x.attrib["Id"]: x.attrib["Target"].lstrip("/")
            for x in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        }
        workbook = ET.fromstring(z.read("xl/workbook.xml"))
        for sheet in workbook.findall("m:sheets/m:sheet", NS):
            if sheet.attrib["name"] != sheet_name:
                continue
            target = relations[sheet.attrib[f"{{{XLREL}}}id"]]
            if not target.startswith("xl/"):
                target = "xl/" + target
            data = ET.fromstring(z.read(target))
            rows = []
            for row in data.findall("m:sheetData/m:row", NS):
                values = {}
                for cell in row:
                    v = cell.find("m:v", NS)
                    value = v.text if v is not None else ""
                    if cell.attrib.get("t") == "s":
                        value = strings[int(value)]
                    elif cell.attrib.get("t") == "inlineStr":
                        value = "".join(cell.find("m:is", NS).itertext())
                    values[re.sub(r"\d", "", cell.attrib["r"])] = value
                rows.append(values)
            return rows
    raise ValueError(sheet_name)


def main():
    source = HERE / "Ma2025_Table2.xlsx"
    rows = read_sheet(source, "Table S2")[2:]
    sets_path = HERE / "verified_macrophage_sets.json"
    existing = json.loads(sets_path.read_text(encoding="utf8"))
    existing = [x for x in existing if not x["id"].startswith("Ma2025_")]
    sensitivity = []
    extracted = []
    unresolved = []
    for cluster in ["Macro_C1QC", "Macro_SPP1"]:
        cluster_rows = [r for r in rows if r["A"] == cluster]
        selected = [r for r in cluster_rows if float(r["G"]) <= 0.05 and float(r["D"]) >= 0.15]
        invalid = [r for r in selected if not re.search(r"[A-Za-z]", r["B"])]
        recoverable = [r for r in selected if re.search(r"[A-Za-z]", r["B"])]
        genes = sorted({r["B"] for r in recoverable})
        base = dict(
            id=f"Ma2025_{cluster}_author_threshold",
            label=f"{cluster} positive marker set (source-derived recoverable-symbol subset)",
            role="Receiver",
            genes=genes,
            doi="10.3389/fimmu.2025.1618591",
            source_file=str(source),
            source_url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12307293/supplementaryFiles",
            extraction_rule=(
                f"Table S2; Cluster exactly {cluster}; P-Value adjusted <= 0.05 and "
                "Average Log2 Fold change >= 0.15, the thresholds explicitly stated "
                "in Methods / Single-cell RNA-seq data processing. Retain all unique "
                "Gene symbols passing both thresholds. Exclude pure-numeric source "
                "identifiers as unresolved source corruption; no gene-name reconstruction."
            ),
            reference_note=(
                "Ma et al. 2025. Source-derived published cluster-marker set, not an "
                "author-frozen standalone signature. Same discovery dataset GSE183904 "
                "as the present study: literature comparison is not independent "
                "cross-dataset validation. C1QC and SPP1 clusters kept separate. "
                "The paper mentions top 30 representative markers; the supplied Table S2 "
                "contains full differential-expression rows. Two numeric/date-corrupted "
                "identifiers per selected set could not be unambiguously recovered; "
                "metrics describe only the recoverable-symbol subset, not the intact set."
            ),
            original_table_rows=len(cluster_rows),
            selected_unique_symbols=len(genes),
            selected_raw_entries=len(selected),
            unresolved_identifiers=[r["B"] for r in invalid],
            unresolved_n=len(invalid),
            shared_discovery_dataset="GSE183904",
            set_type="source-derived positive cluster markers",
        )
        existing.append(base)
        extracted.extend({"cluster": cluster, "gene": r["B"], "avg_log2FC": r["D"], "p_adj": r["G"]} for r in recoverable)
        unresolved.extend({"cluster": cluster, "raw_identifier": r["B"], "avg_log2FC": r["D"], "p_adj": r["G"], "status": "unresolved numeric source identifier; excluded; no reconstruction"} for r in invalid)
        top = sorted(recoverable, key=lambda r: (-float(r["D"]), r["B"]))[:30]
        sensitivity.append({**base, "id": f"Ma2025_{cluster}_top30_sensitivity", "label": f"{cluster} top 30 by published avgLog2FC (analyst-defined sensitivity)", "genes": [r["B"] for r in top], "extraction_rule": base["extraction_rule"] + " Then rank descending by published avgLog2FC, lexical Gene tie-break; first 30. Analyst-defined sensitivity, not verified as the author's frozen representative panel."})
        print(cluster, "original rows", len(cluster_rows), "author thresholds unique", len(genes))
    sets_path.write_text(json.dumps(existing, indent=2), encoding="utf8")
    (HERE / "ma_top30_sensitivity_sets.json").write_text(json.dumps(sensitivity, indent=2), encoding="utf8")
    (HERE / "Ma2025_unresolved_identifiers.json").write_text(json.dumps(unresolved, indent=2), encoding="utf8")
    with (HERE / "Ma2025_author_threshold_marker_rows.tsv").open("w", encoding="utf8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cluster", "gene", "avg_log2FC", "p_adj"], delimiter="\t")
        w.writeheader()
        w.writerows(extracted)


if __name__ == "__main__":
    main()
