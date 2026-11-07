"""Synchronize uppercase figure-panel references in the final WJGO manuscript.

Run after the image-only pass writes embedded/image1.png through image7.png.
Only word/document.xml and the seven existing media entries are replaced.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
from copy import copy
from collections import Counter
from pathlib import Path
from zipfile import ZipFile

from lxml import etree
from PIL import Image


SUBMISSION = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026\WJGOsubmit")
DOCUMENT = SUBMISSION / "00_Ready_to_Upload" / "Manuscript.docx"
CHECKS = SUBMISSION / "06_Checks"
IMAGES = CHECKS / "uppercase_main_panels_20260926" / "embedded"
BACKUP = CHECKS / "backup_before_uppercase_main_panels_20260926" / "Manuscript.docx"
AUDIT = CHECKS / "uppercase_main_panels_20260926" / "manuscript_audit.json"
NAMESPACES = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
TEXT_TAG = "{%s}t" % NAMESPACES["w"]

# These full strings are all contained in individual w:t nodes in the current manuscript.
BODY_REFS = {
    54: ["Figure 6i"],
    72: ["Figure 1a", "Figure 1b-d", "Figure 1e"],
    74: ["Figure 2a,b", "Figure 2c,d"],
    75: ["Figure 2e-k"],
    78: ["Figure 3a,b", "Figure 3c-e"],
    79: ["Figure 3f"],
    81: ["Figure 4a-c"],
    82: ["Figure 4d", "Figure 4e,f", "Figure 4g,h"],
    84: ["Figure 5a,b", "Figure 5c,d"],
    88: ["Figure 6a-f"],
    89: ["Figure 6g", "Figure 6h", "Figure 6i"],
    91: ["Figure 7a\u2013c"],
    92: ["Figure 7d,e", "Figure 7f"],
}
LEGENDS = {
    166: "abcde",
    168: "abcdefghijk",
    170: "abcdef",
    172: "abcdefgh",
    174: "abcd",
    176: "abcdefghi",
    178: "abcdef",
}
PANEL_PHRASES = {
    "In panels b and c": "In panels B and C",
    "In panels b, c and e": "In panels B, C and E",
    "in panels b, c and e": "in panels B, C and E",
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def substitute_once(paragraph: etree._Element, old: str, new: str) -> None:
    hits = [node for node in paragraph.iter(TEXT_TAG) if node.text and old in node.text]
    if len(hits) != 1 or hits[0].text.count(old) != 1:
        raise ValueError(f"Expected exactly one intact Word text node for {old!r}")
    hits[0].text = hits[0].text.replace(old, new)


def updated_xml(original: bytes) -> tuple[bytes, dict]:
    root = etree.fromstring(original)
    formatting_before = etree.fromstring(original)
    paragraphs = root.findall(".//w:body/w:p", NAMESPACES)
    changes = []
    for index, refs in BODY_REFS.items():
        for old in refs:
            new = old[:8] + old[8:].upper()
            substitute_once(paragraphs[index], old, new)
            changes.append({"paragraph": index, "from": old, "to": new})
    if sum(map(len, BODY_REFS.values())) != 23:
        raise AssertionError("Unexpected body reference count")

    for index, letters in LEGENDS.items():
        for letter in letters:
            old, new = f"({letter})", f"({letter.upper()})"
            substitute_once(paragraphs[index], old, new)
            changes.append({"paragraph": index, "from": old, "to": new})
    if sum(map(len, LEGENDS.values())) != 49:
        raise AssertionError("Unexpected legend panel count")

    for old, new in PANEL_PHRASES.items():
        substitute_once(paragraphs[176], old, new)
        changes.append({"paragraph": 176, "from": old, "to": new})

    result = etree.tostring(root, encoding="UTF-8", xml_declaration=True, standalone=True)
    revised = etree.fromstring(result)
    formatting_after = etree.fromstring(result)
    for formatting in (formatting_before, formatting_after):
        for node in formatting.iter(TEXT_TAG):
            node.text = ""
    if etree.tostring(formatting_before, method="c14n") != etree.tostring(formatting_after, method="c14n"):
        raise AssertionError("Word formatting or document structure changed")

    text = " ".join("".join(p.xpath(".//w:t/text()", namespaces=NAMESPACES))
                    for p in revised.findall(".//w:body/w:p", NAMESPACES))
    if re.search(r"Figure [1-7][a-k](?:[,\-\u2013][a-k])*", text):
        raise AssertionError("Lowercase main-figure reference remains")
    if ("Annotation a" not in text or
            "letters denote BH-adjusted thresholds: a, q < 0.05; b, q < 0.01; c, q < 0.001; and d, q < 0.0001" not in text):
        raise AssertionError("Figure 6 statistical notation changed unexpectedly")
    return result, {"changes": changes, "count": Counter(str(x["paragraph"]) for x in changes)}


def main() -> None:
    image_data = {}
    original_data = DOCUMENT.read_bytes()
    with ZipFile(DOCUMENT) as source:
        original_xml = source.read("word/document.xml")
        xml, audit = updated_xml(original_xml)
        for number in range(1, 8):
            name = f"word/media/image{number}.png"
            image_path = IMAGES / f"image{number}.png"
            data = image_path.read_bytes()
            with Image.open(image_path) as im:
                new_size = im.size
            from io import BytesIO
            with Image.open(BytesIO(source.read(name))) as im:
                old_size = im.size
            if old_size != new_size:
                raise AssertionError(f"Embedded Figure {number} size changed: {old_size} -> {new_size}")
            image_data[name] = data

        temporary = DOCUMENT.with_suffix(".uppercase.tmp")
        try:
            with ZipFile(temporary, "w") as destination:
                for entry in source.infolist():
                    payload = (xml if entry.filename == "word/document.xml" else
                               image_data.get(entry.filename, source.read(entry.filename)))
                    destination.writestr(copy(entry), payload)
                destination.comment = source.comment
            with ZipFile(temporary) as result:
                bad = result.testzip()
                if bad is not None:
                    raise AssertionError(f"Corrupt DOCX entry: {bad}")
                altered = {name for name in source.namelist()
                           if sha256(source.read(name)) != sha256(result.read(name))}
                expected = {"word/document.xml", *image_data}
                if altered != expected:
                    raise AssertionError(f"Unexpected changed DOCX parts: {altered ^ expected}")
        except Exception:
            temporary.unlink(missing_ok=True)
            raise

    try:
        if sha256(DOCUMENT.read_bytes()) != sha256(original_data):
            raise RuntimeError("Manuscript changed during processing; refusing to overwrite")
        BACKUP.parent.mkdir(parents=True, exist_ok=True)
        if BACKUP.exists():
            raise FileExistsError(f"Existing backup, refusing to rerun: {BACKUP}")
        shutil.copy2(DOCUMENT, BACKUP)
        temporary.replace(DOCUMENT)
    finally:
        temporary.unlink(missing_ok=True)

    audit.update({"source_sha256": sha256(original_data),
                  "result_sha256": sha256(DOCUMENT.read_bytes()),
                  "changed_parts": sorted({"word/document.xml", *image_data}),
                  "backup": str(BACKUP), "document": str(DOCUMENT)})
    AUDIT.parent.mkdir(parents=True, exist_ok=True)
    AUDIT.write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"text_changes": len(audit["changes"]),
                      "changed_parts": audit["changed_parts"],
                      "document": str(DOCUMENT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
