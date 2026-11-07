from __future__ import annotations

import copy
import hashlib
import io
import json
import math
import os
import re
import shutil
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

from lxml import etree
from PIL import Image, ImageChops, ImageDraw, ImageFont
from scipy import ndimage
import numpy as np


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
SUBMIT = PROJECT / "WJGOsubmit"
DOCX = SUBMIT / "00_Ready_to_Upload" / "Supplementary Information.docx"
CHECKS = SUBMIT / "06_Checks"
PREVIEWS = CHECKS / "supplementary_uppercase_panels_previews_20260926"
BACKUP = CHECKS / "backup_before_uppercase_supplementary_panels_20260926"
REPORT = CHECKS / "SUPPLEMENTARY_UPPERCASE_PANEL_VERIFICATION_20260926.json"
FONT = Path(r"C:\Windows\Fonts\arialbd.ttf")

SOURCES = {
    1: "supplementary/Supplementary_Figure_S1_donor_support",
    2: "ejc_revision_20260827/Supplementary_Figure_S2_PLAT_candidate_and_atlas_ejc_20260827",
    3: "main/FigureS3_patient_state_LR_programme_heatmap_v4_candidate",
    4: "main/FigureS4_cross_cohort_myeloid_mapping_v4_candidate",
    5: "supplementary/Supplementary_Figure_S5_Figure4_robustness",
    6: "ejc_revision_20260827/submission_assets/Figure_S6",
    7: "ejc_revision_20260827/submission_assets/Figure_S7",
    8: "ejc_revision_20260827/submission_assets/Figure_S8",
    9: "ejc_revision_20260827/submission_assets/Figure_S9",
    10: "supplementary/FigureS11_SPON2_ITGB2_spatial_candidate_context",
    11: "supplementary/Supplementary_Figure_S8_receiver_virtual_perturbation_b_aligned_cd_relaxed_20260829",
    12: "supplementary/Supplementary_Figure_S9_ITGB2_directional_connectivity",
    13: "supplementary/Supplementary_Figure_S10_ITGB2_diffreg_GSEA",
}
EXPECTED = {1: "ab", 2: "abcdef", 3: "a", 4: "", 5: "ab", 6: "abcd",
            7: "abcd", 8: "abcdef", 9: "abcd", 10: "abc", 11: "abcd",
            12: "ab", 13: "abc"}
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
REFERENCE = re.compile(
    r"\b(?:Figure(?:s)?|Fig\.)\s+(?:S)?\d+"
    r"(?P<panels>[a-z](?![A-Za-z])(?:\s*[-\u2013,]\s*[a-z](?![A-Za-z]))*)"
)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def panel_elements(svg: Path) -> tuple[ET.Element, list[ET.Element]]:
    root = ET.parse(svg).getroot()
    items = [e for e in root.iter() if e.tag.endswith("}text")
             and (e.text or "").strip() in "abcdef" and (e.text or "").strip()]
    return root, items


def font_size(style: str) -> float:
    match = re.search(r"font-size:\s*([0-9.]+)px", style)
    if not match:
        match = re.search(r"font:\s*700\s+([0-9.]+)px", style)
    assert match, style
    return float(match.group(1))


def patch_image(number: int, raw: bytes) -> tuple[bytes, dict, list]:
    stem = PROJECT / "figures" / SOURCES[number]
    source_png = stem.with_suffix(".png")
    source_svg = stem.with_suffix(".svg")
    assert sha(source_png.read_bytes()) == sha(raw), f"S{number}: source PNG mismatch"
    root, labels = panel_elements(source_svg)
    assert sorted(e.text.strip() for e in labels) == sorted(EXPECTED[number])
    original = Image.open(io.BytesIO(raw))
    original.load()
    assert original.mode in ("RGB", "RGBA")
    result = original.copy()
    view = [float(v) for v in root.get("viewBox").split()]
    sx, sy = original.width / view[2], original.height / view[3]
    assert abs(sx / sy - 1) < 0.005, (number, sx, sy)
    all_regions = Image.new("L", original.size, 0)
    region_draw = ImageDraw.Draw(all_regions)
    details, snippets = [], []
    for element in labels:
        lower = element.text.strip()
        style = element.get("style", "")
        x, y = float(element.get("x")) * sx, float(element.get("y")) * sy
        pixels = font_size(style) * (sx + sy) / 2
        font = ImageFont.truetype(str(FONT), round(pixels))
        glyph_box = font.getbbox(lower, anchor="ls")
        margin = max(8, math.ceil(pixels * 0.45))
        search_box = (max(0, math.floor(x + glyph_box[0]) - margin),
                      max(0, math.floor(y + glyph_box[1]) - margin),
                      min(original.width, math.ceil(x + glyph_box[2]) + margin),
                      min(original.height, math.ceil(y + glyph_box[3]) + margin))
        crop = original.crop(search_box).convert("L")
        foreground = np.asarray(crop) < 245
        components, count = ndimage.label(foreground, structure=np.ones((3, 3)))
        expected_center = (x + (glyph_box[0] + glyph_box[2]) / 2,
                           y + (glyph_box[1] + glyph_box[3]) / 2)
        candidates = []
        for component in range(1, count + 1):
            yy, xx = np.where(components == component)
            if len(xx) < pixels * 0.15:
                continue
            cx = search_box[0] + (xx.min() + xx.max()) / 2
            cy = search_box[1] + (yy.min() + yy.max()) / 2
            distance = math.hypot(cx - expected_center[0], cy - expected_center[1])
            if distance <= pixels * 0.75:
                candidates.append((distance, component, xx, yy))
        assert candidates, (number, lower, "no label component", search_box)
        distance, component, xx, yy = min(candidates, key=lambda item: item[0])
        old_box = (search_box[0] + int(xx.min()), search_box[1] + int(yy.min()),
                   search_box[0] + int(xx.max()) + 1,
                   search_box[1] + int(yy.max()) + 1)
        assert (old_box[2] - old_box[0]) <= pixels * 1.4
        assert (old_box[3] - old_box[1]) <= pixels * 1.4
        component_mask = ndimage.binary_dilation(components == component, iterations=2)
        component_mask = Image.fromarray((component_mask * 255).astype("uint8"))
        background = (255, 255, 255, 255) if original.mode == "RGBA" else (255, 255, 255)
        blank = Image.new(original.mode, component_mask.size, background)
        result.paste(blank, (search_box[0], search_box[1]), component_mask)
        new = "" if number == 3 else lower.upper()
        new_box = old_box
        if new:
            # Fractional SVG positions and antialiased text are retained on the original canvas.
            factor = 4
            high_font = ImageFont.truetype(str(FONT), round(pixels * factor))
            hb = high_font.getbbox(new, anchor="ls")
            new_box = (max(0, math.floor(x + hb[0] / factor) - 2),
                       max(0, math.floor(y + hb[1] / factor) - 2),
                       min(original.width, math.ceil(x + hb[2] / factor) + 2),
                       min(original.height, math.ceil(y + hb[3] / factor) + 2))
            mask_high = Image.new("L", ((new_box[2] - new_box[0]) * factor,
                                        (new_box[3] - new_box[1]) * factor), 0)
            ImageDraw.Draw(mask_high).text(((x - new_box[0]) * factor,
                                           (y - new_box[1]) * factor), new,
                                          fill=255, font=high_font, anchor="ls")
            mask = mask_high.resize((new_box[2] - new_box[0], new_box[3] - new_box[1]),
                                    Image.Resampling.LANCZOS)
            current = result.crop(new_box).convert("RGB")
            current_ink = Image.eval(current.convert("L"), lambda v: 255 if v < 245 else 0)
            overlay_ink = Image.eval(mask, lambda v: 255 if v > 8 else 0)
            assert not ImageChops.multiply(current_ink, overlay_ink).getbbox(), (
                number, lower, "uppercase label overlaps existing content")
            color_match = re.search(r"fill:\s*(#[0-9a-fA-F]{6})", style)
            color = color_match.group(1) if color_match else "#000000"
            text_layer = Image.new(original.mode, mask.size, color)
            result.paste(text_layer, (new_box[0], new_box[1]), mask)
        allowed_box = (min(search_box[0], new_box[0]), min(search_box[1], new_box[1]),
                       max(search_box[2], new_box[2]), max(search_box[3], new_box[3]))
        region_draw.rectangle((allowed_box[0], allowed_box[1], allowed_box[2] - 1,
                               allowed_box[3] - 1), fill=255)
        preview_box = (max(0, allowed_box[0] - 10), max(0, allowed_box[1] - 10),
                       min(original.width, allowed_box[2] + 10),
                       min(original.height, allowed_box[3] + 10))
        snippets.append((f"S{number}{lower}: {new or 'removed'}",
                         original.crop(preview_box), result.crop(preview_box)))
        details.append({"from": lower, "to": new or None, "x": x, "baseline_y": y,
                        "font_px": pixels, "old_box": old_box, "new_box": new_box,
                        "allowed_box": allowed_box})
    difference = ImageChops.difference(original, result).convert("RGB")
    difference_mask = difference.convert("L").point(lambda v: 255 if v else 0)
    outside = ImageChops.multiply(difference_mask, ImageChops.invert(all_regions))
    assert not outside.getbbox(), f"S{number}: pixels changed outside label regions"
    changed_pixels = sum(v for index, v in enumerate(difference_mask.histogram()) if index)
    result.save(PREVIEWS / f"Figure_S{number}.png")
    if not labels:
        output = raw
    else:
        buffer = io.BytesIO()
        metadata = {k: original.info[k] for k in ("dpi", "icc_profile") if k in original.info}
        result.save(buffer, format="PNG", **metadata)
        output = buffer.getvalue()
    return output, {"figure": f"S{number}", "source_svg": str(source_svg),
                    "source_png_sha256": sha(raw), "size": original.size,
                    "mode": original.mode, "labels": details,
                    "changed_pixels": changed_pixels, "outside_changed_pixels": 0}, snippets


def patch_text(xml: bytes) -> tuple[bytes, dict]:
    root = etree.fromstring(xml)
    before = copy.deepcopy(root)
    edits = []
    counters = {"caption_panel_letters": 0, "figure_reference_letters": 0,
                "inline_panel_letters": 0}
    for index, paragraph in enumerate(root.xpath("//w:body/w:p", namespaces=NS)):
        nodes = paragraph.xpath(".//w:t", namespaces=NS)
        text = "".join(n.text or "" for n in nodes)
        replacements = {}
        for match in REFERENCE.finditer(text):
            for position in range(*match.span("panels")):
                if "a" <= text[position] <= "z":
                    replacements[position] = text[position].upper()
                    counters["figure_reference_letters"] += 1
        if re.match(r"Fig\. S\d+\.", text):
            for match in re.finditer(r"\(([a-f])\)", text):
                replacements[match.start(1)] = match.group(1).upper()
                counters["caption_panel_letters"] += 1
            for match in re.finditer(r"\bpanels ([a-f]) and ([a-f])\b", text):
                for group in (1, 2):
                    replacements[match.start(group)] = match.group(group).upper()
                    counters["inline_panel_letters"] += 1
        if not replacements:
            continue
        offset = 0
        for node in nodes:
            value = node.text or ""
            node.text = "".join(replacements.get(offset + j, char)
                                for j, char in enumerate(value))
            offset += len(value)
        revised = "".join(n.text or "" for n in nodes)
        assert len(text) == len(revised)
        assert all(a == b or a.upper() == b for a, b in zip(text, revised))
        edits.append({"paragraph_index": index, "before": text, "after": revised,
                      "character_changes": len(replacements)})
    original_nodes = [[node.text or "" for node in paragraph.xpath(".//w:t", namespaces=NS)]
                      for paragraph in before.xpath("//w:body/w:p", namespaces=NS)]
    for document in (before, root):
        for node in document.xpath("//w:t", namespaces=NS):
            node.text = ""
    assert etree.tostring(before, method="c14n") == etree.tostring(root, method="c14n"), (
        "Document markup or formatting changed")
    # Patch only raw text content so ZIP XML formatting remains byte-for-byte intact.
    paragraph_pattern = re.compile(rb"<w:p(?: [^>]*)?>.*?</w:p>", re.DOTALL)
    text_pattern = re.compile(rb"(<w:t(?: [^>]*)?>)(.*?)(</w:t>)", re.DOTALL)
    raw_paragraphs = list(paragraph_pattern.finditer(xml))
    assert len(raw_paragraphs) == len(root.xpath("//w:body/w:p", namespaces=NS))
    byte_changes = []
    for edit in edits:
        match = raw_paragraphs[edit["paragraph_index"]]
        raw = match.group(0)
        nodes = original_nodes[edit["paragraph_index"]]
        chunks = list(text_pattern.finditer(raw))
        assert len(chunks) == len(nodes)
        changes = {j: b for j, (a, b) in enumerate(zip(edit["before"], edit["after"])) if a != b}
        offset = 0
        for chunk, decoded in zip(chunks, nodes):
            tokens = list(re.finditer(r"&(?:#[0-9]+|#x[0-9A-Fa-f]+|[A-Za-z]+);|.",
                                     chunk.group(2).decode("utf-8"), re.DOTALL))
            assert len(tokens) == len(decoded)
            value = chunk.group(2).decode("utf-8")
            for j, token in enumerate(tokens):
                if offset + j not in changes:
                    continue
                assert token.group(0) == decoded[j] and len(token.group(0)) == 1
                position = match.start() + chunk.start(2) + len(value[:token.start()].encode("utf-8"))
                assert xml[position:position + 1] == decoded[j].encode("ascii")
                byte_changes.append((position, changes[offset + j].encode("ascii")[0]))
            offset += len(decoded)
    output_array = bytearray(xml)
    for position, char in byte_changes:
        output_array[position] = char
    output = bytes(output_array)
    assert len(output) == len(xml)
    assert sum(a != b for a, b in zip(xml, output)) == len(byte_changes)
    revised = etree.fromstring(output)
    for paragraph in revised.xpath("//w:body/w:p", namespaces=NS):
        text = "".join(paragraph.xpath(".//w:t/text()", namespaces=NS))
        assert not REFERENCE.search(text), text
        if re.match(r"Fig\. S\d+\.", text):
            assert not re.search(r"\([a-f]\)|\bpanels [a-f] and [a-f]\b", text), text
    return output, {"counters": counters, "changed_paragraphs": len(edits), "edits": edits,
                    "formatting_and_markup_unchanged": True,
                    "raw_xml_byte_changes": len(byte_changes),
                    "raw_xml_other_bytes_identical": True}


def contact_sheet(snippets: list) -> None:
    width, cell_height, columns = 780, 115, 3
    rows = math.ceil(len(snippets) / columns)
    canvas = Image.new("RGB", (width, rows * cell_height), "white")
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype(str(FONT), 16)
    for index, (title, original, revised) in enumerate(snippets):
        x, y = (index % columns) * (width // columns), (index // columns) * cell_height
        draw.text((x + 4, y + 3), title, fill="black", font=font)
        for slot, crop in enumerate((original, revised)):
            crop = crop.convert("RGB")
            crop.thumbnail((100, 78), Image.Resampling.LANCZOS)
            canvas.paste(crop, (x + 12 + slot * 116, y + 30))
    canvas.save(PREVIEWS / "all_panel_changes_contact_sheet.png")


def main() -> None:
    PREVIEWS.mkdir(parents=True, exist_ok=True)
    BACKUP.mkdir(parents=True, exist_ok=True)
    current = DOCX.read_bytes()
    backup = BACKUP / DOCX.name
    if backup.exists():
        original = backup.read_bytes()
        if current != original:
            prior = json.loads(REPORT.read_text(encoding="utf-8"))
            assert sha(current) == prior["updated_sha256"], "DOCX changed outside this script"
    else:
        shutil.copy2(DOCX, backup)
        original = current
    replacement, image_report, snippets = {}, [], []
    with zipfile.ZipFile(io.BytesIO(original)) as archive:
        assert archive.testzip() is None
        for number in range(1, 14):
            member = f"word/media/image{number + 24}.png"
            output, details, crops = patch_image(number, archive.read(member))
            if output != archive.read(member):
                replacement[member] = output
            image_report.append(details)
            snippets.extend(crops)
        xml, text_report = patch_text(archive.read("word/document.xml"))
        replacement["word/document.xml"] = xml
        target = DOCX.with_name(DOCX.stem + ".panels-tmp.docx")
        with zipfile.ZipFile(target, "w") as output:
            for info in archive.infolist():
                output.writestr(info, replacement.get(info.filename, archive.read(info.filename)))
    with zipfile.ZipFile(target) as output, zipfile.ZipFile(io.BytesIO(original)) as archive:
        assert output.testzip() is None
        assert output.namelist() == archive.namelist()
        changed = [n for n in archive.namelist() if archive.read(n) != output.read(n)]
        assert sorted(changed) == sorted(replacement)
    assert DOCX.read_bytes() == current, "DOCX changed concurrently"
    os.replace(target, DOCX)
    contact_sheet(snippets)
    report = {"docx": str(DOCX), "backup": str(backup), "original_sha256": sha(original),
              "updated_sha256": sha(DOCX.read_bytes()), "figures": image_report,
              "text": text_report, "changed_zip_members": changed,
              "all_other_zip_members_identical": True, "zip_integrity": "passed",
              "uppercase_panel_count": 40, "removed_isolated_panel_count": 1,
              "figure_count": 13}
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"docx": str(DOCX), "report": str(REPORT), "text": text_report["counters"],
                      "changed_paragraphs": text_report["changed_paragraphs"],
                      "changed_zip_members": changed, "status": "passed"}, indent=2))


if __name__ == "__main__":
    main()
