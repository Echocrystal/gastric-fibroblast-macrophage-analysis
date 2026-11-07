"""Integrate the audited WJGO Figure 3/6 artwork and matching legends.

The script is guarded by the SHA-256 of the polished manuscript so that it
cannot silently overwrite a later author edit.  It preserves rollback copies,
updates one canonical DOCX, embeds all seven figures, and then creates the
second manuscript copy byte-for-byte from the canonical file.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

from docx import Document

import create_wjgo_package as typography
import embed_wjgo_figures


ROOT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
PACKAGE = ROOT / "WJGOsubmit"
PRIMARY = PACKAGE / "00_Ready_to_Upload" / "Manuscript_WJGO.docx"
SECONDARY = PACKAGE / "01_Manuscript" / "Manuscript_WJGO.docx"
CANDIDATE = ROOT / "tmp" / "wjgo_figures_candidate"
SOURCE_DIR = PACKAGE / "05_Editable_Sources" / "Figures_WJGO"
CHECK_DIR = PACKAGE / "06_Checks"
BACKUP_DIR = CHECK_DIR / "backup_before_figure_integration_20260924"
REPORT = CHECK_DIR / "FIGURE_3_6_INTEGRATION_AUDIT_20260924.json"

EXPECTED_MANUSCRIPT_SHA256 = (
    "f730ef952ba4bb169e55d2390235a23fffe326ca12edd31cd04cc88548309e25"
)
EXPECTED_CANDIDATE_SHA256 = {
    "Figure_3_WJGO_candidate.png": "b3f93929f6b80c6eb43242d48395f7883b3fcd703d2cad21cc3456799882b84b",
    "Figure_3_WJGO_candidate.tiff": "d71afd13bfbd0c18fa9056c4afcabe264c726094d2f34d8618a417d02bbdee16",
    "Figure_6_WJGO_candidate.png": "c5a3a19c858f49a92330718198e70117104a86c36ed3a36506e93f60e0e51d35",
    "Figure_6_WJGO_candidate.tiff": "ad421ee003f1e5b1b0b488e11150f360769bd77e7ecdacdba4b2df2a49616446",
}

FIGURE_3_OLD_D = (
    "(d) F1 PLAT-rich fibroblast communication across C1–C5 with significance markers."
)
FIGURE_3_NEW_D = (
    "(d) F1 PLAT-rich fibroblast communication across C1–C5. Black cell outlines "
    "identify observations that simultaneously met the CellPhoneDB permutation P < 0.05, "
    "LIANA magnitude rank ≤ 0.05, and LIANA specificity rank ≤ 0.05 criteria. Right-hand "
    "bars show the number of receiver states satisfying all three criteria for each "
    "ligand–receptor pair."
)
FIGURE_3_OLD_F = (
    "(f) Covariate-adjusted associations between PLAT-rich fibroblast abundance and "
    "myeloid states with FDR values."
)
FIGURE_3_NEW_F = (
    "(f) Points show adjusted Spearman correlations between PLAT-rich fibroblast abundance "
    "and myeloid states; horizontal lines show leave-one-patient-out ranges, and the "
    "right-hand column reports Benjamini–Hochberg FDR values in scientific notation."
)

FIGURE_6_OLD = (
    "Figure 6 Bulk projection links the PLAT-rich fibroblast–C1QC+APOE+ macrophage ecology "
    "to immune composition, clinicopathology, established classifications and anti-PD-1 "
    "response. (a) TCGA-STAD sender–receiver score plane with four configurations. (b) "
    "MCP-counter populations across configurations; boxes show interquartile ranges, "
    "medians, whiskers to 1.50×IQR, and pairwise Dunn tests with BH correction. (c) "
    "Checkpoint and myeloid-marker expression. (d) Random-effects clinicopathological "
    "meta-analysis; intervals are 95% CIs and q values. (e) Associations with six "
    "gastric-cancer classifications; Cramer's V and BH-adjusted q values. (f) TCGA "
    "molecular-subtype composition within each configuration; percentages rounded to "
    "integers. (g) Kim2018 pretreatment score plane for 45 patients, including 12 "
    "responders; two-sided four-configuration Fisher's exact P = 7.24 × 10−6. (h) "
    "Bias-reduced logistic response models; points show odds ratios per 1-SD increase and "
    "intervals are 95% CIs. Model definitions, nominal P values and BH-adjusted FDR values "
    "are provided in Supplementary Table S6. (i) GDSC2 gastric-cancer cell-line AUC "
    "distributions; lower AUC indicates greater in-vitro sensitivity, pairwise comparisons "
    "for each drug, BH-adjusted q values. Significance codes: *q < 0.05, **q < 0.01, "
    "***q < 1.00 × 10−3, ****q < 1.00 × 10−4 and ns q ≥ 0.05. Abbreviations: AUC, area "
    "under the dose–response curve; BH, Benjamini–Hochberg; CI, confidence interval; FDR, "
    "false discovery rate; GDSC2, Genomics of Drug Sensitivity in Cancer version 2; IQR, "
    "interquartile range; MCP-counter, Microenvironment Cell Populations-counter; "
    "TCGA-STAD, The Cancer Genome Atlas stomach adenocarcinoma cohort."
)

FIGURE_6_NEW = (
    "Figure 6 Bulk projection links the PLAT-rich fibroblast–C1QC+APOE+ macrophage ecology "
    "to immune composition, clinicopathology, established classifications and anti-PD-1 "
    "response. (a) TCGA-STAD sender–receiver score plane with four configurations. (b) "
    "MCP-counter populations across configurations. (c) Checkpoint and myeloid-marker "
    "expression. In panels b and c, boxes show medians and interquartile ranges; whiskers "
    "extend to 1.50 times the interquartile range, with outliers not displayed. Horizontal "
    "brackets show pairwise Dunn tests with BH-adjusted q values. Only significant "
    "comparisons are annotated: a, q < 0.05; b, q < 0.01; c, q < 0.001; and d, q < "
    "0.0001. Non-significant comparisons are not annotated. (d) Random-effects "
    "clinicopathological meta-analysis; intervals are 95% CIs and q values. (e) Points "
    "show Cramer's V for associations of the four configurations with the indicated "
    "gastric-cancer classifications. Tests were Pearson χ² tests or Monte Carlo Pearson "
    "χ² tests with 99999 replicates, as appropriate for table sparsity. Letters denote "
    "BH-adjusted values across the prespecified classification tests: a, q < 0.05; b, q < "
    "0.01; c, q < 0.001; and d, q < 0.0001. Non-significant associations are not "
    "annotated. (f) TCGA molecular-subtype composition within each configuration; "
    "percentages are rounded to integers. (g) Kim2018 pretreatment score plane for 45 "
    "patients, including 12 responders; two-sided four-configuration Fisher's exact P = "
    "7.24 × 10−6. (h) Points and horizontal lines show odds ratios and 95% CIs from "
    "mean-bias-reduced continuous logistic-regression models for observed Kim2018 response. "
    "Only the nominally significant PLAT-rich fibroblast-state estimate is annotated (a, "
    "P = 0.036; two-sided unadjusted model-level P value); non-significant P values are not "
    "annotated. These analyses are exploratory single-arm associations and do not estimate "
    "treatment benefit. Model definitions, nominal P values and BH-adjusted FDR values are "
    "provided in Supplementary Table S6. (i) Violin shapes show GDSC2 gastric-cancer "
    "cell-line AUC distributions; lower AUC indicates greater in-vitro sensitivity. "
    "Internal boxes show medians and interquartile ranges, and whiskers span the observed "
    "range. Pairwise Mann–Whitney comparisons were BH-adjusted across the prespecified "
    "drug-by-pair family. No displayed pairwise comparison met q < 0.05; therefore, "
    "non-significant brackets and labels are not shown. Abbreviations: AUC, area under the "
    "dose–response curve; BH, Benjamini–Hochberg; CI, confidence interval; FDR, false "
    "discovery rate; GDSC2, Genomics of Drug Sensitivity in Cancer version 2; IQR, "
    "interquartile range; MCP-counter, Microenvironment Cell Populations-counter; "
    "TCGA-STAD, The Cancer Genome Atlas stomach adenocarcinoma cohort."
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_with_audit(source: Path, target: Path) -> dict[str, object]:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    if sha256(source) != sha256(target):
        raise AssertionError(f"Copy verification failed: {source} -> {target}")
    return {
        "source": str(source),
        "target": str(target),
        "bytes": target.stat().st_size,
        "sha256": sha256(target),
    }


def replace_exact(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise AssertionError(f"{label}: expected one old fragment, found {count}")
    return text.replace(old, new, 1)


def main() -> None:
    if sha256(PRIMARY) != EXPECTED_MANUSCRIPT_SHA256:
        raise AssertionError("Primary manuscript no longer matches the reviewed input hash")
    if sha256(SECONDARY) != EXPECTED_MANUSCRIPT_SHA256:
        raise AssertionError("Secondary manuscript no longer matches the reviewed input hash")

    for name, expected_hash in EXPECTED_CANDIDATE_SHA256.items():
        path = CANDIDATE / name
        if sha256(path) != expected_hash:
            raise AssertionError(f"Candidate hash mismatch: {name}")

    BACKUP_DIR.mkdir(parents=True, exist_ok=True)
    backup_items = []
    for source in (
        PRIMARY,
        SECONDARY,
        PACKAGE / "00_Ready_to_Upload" / "Figure_3.tiff",
        PACKAGE / "00_Ready_to_Upload" / "Figure_6.tiff",
        PACKAGE / "02_Figures" / "Figure_3.tiff",
        PACKAGE / "02_Figures" / "Figure_6.tiff",
    ):
        relative = source.relative_to(PACKAGE)
        backup_items.append(copy_with_audit(source, BACKUP_DIR / relative))

    source_copies = []
    for number in (3, 6):
        for extension in ("png", "tiff", "pdf", "svg"):
            source = CANDIDATE / f"Figure_{number}_WJGO_candidate.{extension}"
            target = SOURCE_DIR / f"Figure_{number}_WJGO.{extension}"
            source_copies.append(copy_with_audit(source, target))

    upload_copies = []
    for number in (3, 6):
        source = SOURCE_DIR / f"Figure_{number}_WJGO.tiff"
        for directory in ("00_Ready_to_Upload", "02_Figures"):
            target = PACKAGE / directory / f"Figure_{number}.tiff"
            upload_copies.append(copy_with_audit(source, target))

    document = Document(PRIMARY)
    legend_changes = []
    for paragraph_index, paragraph in enumerate(document.paragraphs):
        text = paragraph.text
        if text.startswith("Figure 3 "):
            revised = replace_exact(text, FIGURE_3_OLD_D, FIGURE_3_NEW_D, "Figure 3d")
            revised = replace_exact(revised, FIGURE_3_OLD_F, FIGURE_3_NEW_F, "Figure 3f")
            paragraph.text = revised
            typography.rebuild_paragraph_runs(
                paragraph, include_statistics=True, include_citations=True
            )
            legend_changes.append({"figure": 3, "paragraph": paragraph_index, "text": revised})
        elif text.startswith("Figure 6 "):
            if text != FIGURE_6_OLD:
                raise AssertionError("Figure 6 legend does not match the reviewed source text")
            paragraph.text = FIGURE_6_NEW
            typography.rebuild_paragraph_runs(
                paragraph, include_statistics=True, include_citations=True
            )
            legend_changes.append(
                {"figure": 6, "paragraph": paragraph_index, "text": FIGURE_6_NEW}
            )

    if [item["figure"] for item in legend_changes] != [3, 6]:
        raise AssertionError(f"Unexpected legend-change set: {legend_changes}")

    typography.force_black_and_remove_borders(document)
    document.save(PRIMARY)
    typography.sanitize_style_parts(PRIMARY)
    embed_result = embed_wjgo_figures.embed_figures(PRIMARY)
    shutil.copy2(PRIMARY, SECONDARY)

    if sha256(PRIMARY) != sha256(SECONDARY):
        raise AssertionError("The two final manuscript copies are not byte-identical")
    final_document = Document(PRIMARY)
    if len(final_document.inline_shapes) != 7:
        raise AssertionError(f"Expected seven inline figures, found {len(final_document.inline_shapes)}")

    figure_media = {
        f"Figure_{number}": {
            "ready_sha256": sha256(PACKAGE / "00_Ready_to_Upload" / f"Figure_{number}.tiff"),
            "figures_sha256": sha256(PACKAGE / "02_Figures" / f"Figure_{number}.tiff"),
        }
        for number in (3, 6)
    }
    if any(value["ready_sha256"] != value["figures_sha256"] for value in figure_media.values()):
        raise AssertionError("Independent TIFF copies are not byte-identical")

    report = {
        "status": "PASS",
        "input_manuscript_sha256": EXPECTED_MANUSCRIPT_SHA256,
        "backup_items": backup_items,
        "editable_source_copies": source_copies,
        "upload_tiff_copies": upload_copies,
        "legend_changes": legend_changes,
        "embed_result": embed_result,
        "final_manuscript": {
            "bytes": PRIMARY.stat().st_size,
            "sha256": sha256(PRIMARY),
            "secondary_sha256": sha256(SECONDARY),
            "byte_identical": PRIMARY.read_bytes() == SECONDARY.read_bytes(),
            "inline_shapes": len(final_document.inline_shapes),
        },
        "figure_media": figure_media,
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
