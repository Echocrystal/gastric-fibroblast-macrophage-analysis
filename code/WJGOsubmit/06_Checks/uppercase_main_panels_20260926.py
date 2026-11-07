"""Replace only main-figure panel identifiers in the current raster images."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from lxml import etree
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
SUBMISSION = ROOT / "WJGOsubmit"
READY = SUBMISSION / "00_Ready_to_Upload"
CHECKS = SUBMISSION / "06_Checks"
OUTPUT = CHECKS / "uppercase_main_panels_20260926"
BACKUP = CHECKS / "backup_before_uppercase_main_panels_20260926"
FONT = Path(r"C:\Windows\Fonts\arialbd.ttf")
COUNTS = {1: 5, 2: 11, 3: 6, 4: 8, 5: 4, 6: 9, 7: 6}
SOURCES = {
    **{n: ROOT / "figures" / "main" / "final_20260817" / f"Figure_{n}.svg"
       for n in (1, 2, 4, 5)},
    **{n: SUBMISSION / "05_Editable_Sources" / "Figures_WJGO" / f"Figure_{n}_WJGO.svg"
       for n in (3, 6)},
    7: ROOT / "figures" / "ejc_revision_20260827" / "submission_assets" / "Figure_7.svg",
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def labels(number):
    root = etree.parse(str(SOURCES[number])).getroot()
    viewbox = [float(v) for v in root.get("viewBox").replace(",", " ").split()]
    found = []
    for node in root.xpath("//*[local-name()='text']"):
        letter = "".join(node.itertext()).strip()
        style = node.get("style", "")
        if (len(letter) == 1 and letter in "abcdefghijk" and
                node.get("x") and node.get("y") and ("bold" in style or "700" in style)):
            match = re.search(r"font-size:\s*([\d.]+)px", style)
            if not match:
                match = re.search(r"font:\s*700\s+([\d.]+)px", style)
            if not match:
                raise ValueError(f"Unknown font size: {style}")
            found.append({"letter": letter, "x": float(node.get("x")),
                          "y": float(node.get("y")), "size": float(match.group(1)),
                          "new_x_shift_svg": -1.0 if number == 2 and letter == "b" else 0.0})
    if "".join(x["letter"] for x in found) != "abcdefghijk"[:COUNTS[number]]:
        raise AssertionError(f"Unexpected panel identifiers in Figure {number}")
    return viewbox, found


def components(mask):
    pending = mask.copy()
    result = []
    for y, x in zip(*np.nonzero(mask)):
        if not pending[y, x]:
            continue
        stack = [(int(y), int(x))]
        pending[y, x] = False
        pixels = []
        while stack:
            cy, cx = stack.pop()
            pixels.append((cy, cx))
            for ny in range(max(0, cy - 1), min(mask.shape[0], cy + 2)):
                for nx in range(max(0, cx - 1), min(mask.shape[1], cx + 2)):
                    if pending[ny, nx]:
                        pending[ny, nx] = False
                        stack.append((ny, nx))
        if len(pixels) >= 3:
            result.append(pixels)
    return result


def glyph_mask(letter, size):
    multiplier = 4
    font = ImageFont.truetype(str(FONT), round(size * multiplier))
    canvas = Image.new("L", (round(size * 2 * multiplier), round(size * 2 * multiplier)))
    origin = (round(size * .5 * multiplier), round(size * multiplier))
    ImageDraw.Draw(canvas).text(origin, letter, font=font, fill=255, anchor="ls")
    small = canvas.resize((round(size * 2), round(size * 2)), Image.Resampling.LANCZOS)
    bbox = small.getbbox()
    origin_small = (origin[0] / multiplier, origin[1] / multiplier)
    return small.crop(bbox), (bbox[0] - origin_small[0], bbox[1] - origin_small[1])


def plan_label(im, spec, viewbox):
    x = (spec["x"] - viewbox[0]) * im.width / viewbox[2]
    y = (spec["y"] - viewbox[1]) * im.height / viewbox[3]
    size = spec["size"] * im.height / viewbox[3]
    font = ImageFont.truetype(str(FONT), round(size))
    letter = spec["letter"]
    expected = font.getbbox(letter, anchor="ls")
    box = (max(0, int(x + expected[0] - 15)), max(0, int(y + expected[1] - 15)),
           min(im.width, round(x + expected[2] + 15)), min(im.height, round(y + expected[3] + 15)))
    pixels = np.array(im.crop(box))
    ink = np.any(pixels[:, :, :3] < 254, axis=2)
    groups = components(ink)
    if not groups:
        raise AssertionError(f"No glyph detected for {letter} at {x}, {y}")
    groups.sort(key=len, reverse=True)
    selected = groups[:2] if letter in "ij" else groups[:1]
    points = np.array([p for group in selected for p in group])
    bbox = (box[0] + int(points[:, 1].min()), box[1] + int(points[:, 0].min()),
            box[0] + int(points[:, 1].max()) + 1, box[1] + int(points[:, 0].max()) + 1)
    if bbox[0] <= box[0] or bbox[1] <= box[1] or bbox[2] >= box[2] or bbox[3] >= box[3]:
        raise AssertionError(f"Glyph may be clipped for {letter}: {bbox} in {box}")
    clear = (max(0, bbox[0] - 1), max(0, bbox[1] - 1),
             min(im.width, bbox[2] + 1), min(im.height, bbox[3] + 1))
    glyph_pixels = np.array(im.crop(bbox))[:, :, :3]
    darkest = tuple(int(v) for v in glyph_pixels.reshape(-1, 3).min(axis=0))
    new_mask, offset = glyph_mask(letter.upper(), size)
    new_x_shift = spec["new_x_shift_svg"] * im.width / viewbox[2]
    new_x, new_y = round(x + offset[0] + new_x_shift), round(y + offset[1])
    if new_x < 0 or new_y < 0:
        raise AssertionError(f"New glyph clipped for {letter}: {(new_x, new_y)}")
    new_box = (new_x, new_y, new_x + new_mask.width, new_y + new_mask.height)
    # The entire edit rectangle must contain only this glyph and white background.
    edit = (min(clear[0], new_box[0]), min(clear[1], new_box[1]),
            max(clear[2], new_box[2]), max(clear[3], new_box[3]))
    region = np.array(im.crop(edit))[:, :, :3]
    foreign = np.any(region < 254, axis=2)
    foreign[clear[1] - edit[1]:clear[3] - edit[1], clear[0] - edit[0]:clear[2] - edit[0]] = False
    proposed = np.zeros(foreign.shape, dtype=np.uint8)
    proposed[new_box[1] - edit[1]:new_box[3] - edit[1], new_box[0] - edit[0]:new_box[2] - edit[0]] = np.array(new_mask)
    collision = foreign & (proposed > 0)
    if np.any(collision):
        coords = [(edit[0] + int(px), edit[1] + int(py), int(proposed[py, px]))
                  for py, px in zip(*np.nonzero(collision))]
        raise AssertionError(f"New glyph overlaps neighboring ink for {letter}: {edit}, pixels {coords[:20]}")
    unselected = {tuple(p) for group in groups[len(selected):] for p in group}
    if any(clear[0] <= box[0] + px < clear[2] and clear[1] <= box[1] + py < clear[3]
           for py, px in unselected):
        raise AssertionError(f"Clearing another mark near {letter}")
    return {"letter": letter, "uppercase": letter.upper(), "svg_position": [spec["x"], spec["y"]],
            "baseline": [x, y], "font_size_px": size, "old_ink_bbox": list(bbox),
            "new_x_shift_px": new_x_shift,
            "clear_bbox": list(clear), "new_bbox": list(new_box), "edit_bbox": list(edit),
            "color": list(darkest), "component_sizes": [len(g) for g in selected]}, new_mask


def replace_labels(image, specs, viewbox):
    im = image.convert("RGB")
    before = np.array(im)
    after = im.copy()
    plans = []
    for spec in specs:
        plan, mask = plan_label(im, spec, viewbox)
        ImageDraw.Draw(after).rectangle((plan["clear_bbox"][0], plan["clear_bbox"][1],
                                        plan["clear_bbox"][2] - 1, plan["clear_bbox"][3] - 1), fill="white")
        fill = Image.new("RGB", mask.size, tuple(plan["color"]))
        after.paste(fill, tuple(plan["new_bbox"][:2]), mask)
        plans.append(plan)
    difference = np.any(np.array(after) != before, axis=2)
    allowable = np.zeros(difference.shape, dtype=bool)
    for plan in plans:
        left, top, right, bottom = plan["edit_bbox"]
        allowable[top:bottom, left:right] = True
        plan["changed_pixels"] = int(difference[top:bottom, left:right].sum())
    if np.any(difference & ~allowable):
        raise AssertionError("Pixels outside panel-label edit rectangles changed")
    return after, plans, int(difference.sum())


def contact_sheet(number, before, after, plans, kind):
    columns = 4
    tile_width, tile_height = 210, 150
    rows = (len(plans) + columns - 1) // columns
    sheet = Image.new("RGB", (columns * tile_width, rows * tile_height), "#eeeeee")
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.truetype(str(FONT), 14)
    for index, plan in enumerate(plans):
        x, y = index % columns * tile_width, index // columns * tile_height
        l, t, r, b = plan["edit_bbox"]
        box = (max(0, l - 5), max(0, t - 5), min(before.width, r + 12), min(before.height, b + 12))
        draw.text((x + 6, y + 5), f"Figure {number}{plan['uppercase']} {kind}", fill="black", font=font)
        for offset, picture in ((8, before), (108, after)):
            crop = picture.crop(box)
            crop.thumbnail((90, 104))
            sheet.paste(crop, (x + offset, y + 30))
        draw.text((x + 8, y + 134), "Before", fill="black", font=font)
        draw.text((x + 108, y + 134), "After", fill="black", font=font)
    sheet.save(OUTPUT / f"Figure_{number}_{kind}_before_after.png")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "embedded").mkdir(exist_ok=True)
    (OUTPUT / "tiff").mkdir(exist_ok=True)
    audit = {"document": str(READY / "Manuscript.docx"), "mode": "apply" if args.apply else "dry_run", "figures": []}
    with ZipFile(READY / "Manuscript.docx") as doc:
        for n in range(1, 8):
            viewbox, specs = labels(n)
            item = {"number": n, "svg_reference": str(SOURCES[n])}
            for kind in ("TIFF", "PNG"):
                path = READY / f"Figure_{n}.tiff"
                media_name = f"word/media/image{n}.png"
                data = path.read_bytes() if kind == "TIFF" else doc.read(media_name)
                with Image.open(BytesIO(data)) as original:
                    dpi = original.info.get("dpi")
                    if dpi:
                        dpi = tuple(float(value) for value in dpi)
                    info = dict(original.info)
                    size = original.size
                    mode = original.mode
                    before = original.convert("RGB")
                    after, plans, changed = replace_labels(original, specs, viewbox)
                    if mode == "RGBA":
                        after.putalpha(original.getchannel("A"))
                contact_sheet(n, before, after, plans, kind)
                destination = OUTPUT / ("tiff" if kind == "TIFF" else "embedded") / (f"Figure_{n}.tiff" if kind == "TIFF" else f"image{n}.png")
                options = {"dpi": dpi} if dpi else {}
                if info.get("icc_profile"):
                    options["icc_profile"] = info["icc_profile"]
                if kind == "TIFF":
                    options["compression"] = info.get("compression", "tiff_lzw")
                after.save(destination, **options)
                with Image.open(destination) as saved:
                    if saved.size != size or saved.mode != mode or not np.array_equal(np.array(saved), np.array(after)):
                        raise AssertionError(f"Saved image differs from prepared pixels: {destination}")
                    if dpi and any(abs(float(a) - float(b)) > .01 for a, b in zip(dpi, saved.info["dpi"])):
                        raise AssertionError(f"DPI changed: {destination}")
                if args.apply:
                    backup = BACKUP / ("tiff" if kind == "TIFF" else "embedded") / destination.name
                    backup.parent.mkdir(parents=True, exist_ok=True)
                    if backup.exists():
                        raise FileExistsError(f"Backup already exists: {backup}")
                    if kind == "TIFF":
                        shutil.copy2(path, backup)
                        shutil.copy2(destination, path)
                        shutil.copy2(destination, SUBMISSION / "02_Figures" / path.name)
                    else:
                        backup.write_bytes(data)
                item[kind] = {"source_sha256": digest(data), "result_sha256": digest(destination.read_bytes()),
                              "size": size, "dpi": dpi, "mode": mode, "changed_pixels": changed,
                              "pixels_outside_edit_boxes_unchanged": True, "labels": plans,
                              "output": str(destination)}
            audit["figures"].append(item)
            print(f"Figure {n}: {COUNTS[n]} panel labels prepared in TIFF and PNG")
    audit["total_labels_per_format"] = sum(COUNTS.values())
    (OUTPUT / "image_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(f"Completed {audit['mode']}; output: {OUTPUT}")


if __name__ == "__main__":
    main()
