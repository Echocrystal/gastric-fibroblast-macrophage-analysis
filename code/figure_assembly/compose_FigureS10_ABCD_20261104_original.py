"""Cosmetic composition of frozen S10 A-C and the approved new D.

Internal working files stay in scripts/_FigureS10_composition_work.
The user delivery folder contains exactly one 600-dpi PNG.
"""
from pathlib import Path
from io import BytesIO
import hashlib
import json
import subprocess
import sys

from PIL import Image
from pypdf import PdfReader, PdfWriter, PageObject, Transformation
from pypdf.generic import RectangleObject
from reportlab.pdfgen.canvas import Canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

sys.stdout.reconfigure(encoding='utf-8')
ROOT = Path(r'I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026')
REVIEW = ROOT/'WJGOsubmit/review/R2_supplementary_analyses_20261103'
ABC = ROOT/'figures/supplementary/FigureS11_SPON2_ITGB2_spatial_candidate_context.pdf'
D = REVIEW/'13_insertion_ready_figures_20261104/spatial/Spatial_endpoints_with_overall_INSERTION.pdf'
WORK = REVIEW/'scripts/_FigureS10_composition_work'
DELIVERY = REVIEW/'15_FigureS10_ABCD_20261104'
WORK.mkdir(exist_ok=True)
DELIVERY.mkdir(exist_ok=True)
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
before = {str(p): sha(p) for p in (ABC, D)}

pdfmetrics.registerFont(TTFont('ArialBold', r'C:/Windows/Fonts/arialbd.ttf'))

def make_overlay(width, height, draw):
    stream = BytesIO()
    canvas = Canvas(stream, pagesize=(width, height))
    draw(canvas)
    canvas.save()
    stream.seek(0)
    return PdfReader(stream).pages[0]

abc = PdfReader(ABC).pages[0]
d = PdfReader(D).pages[0]
abc_w, abc_h = map(float, (abc.mediabox.width, abc.mediabox.height))
d_w, d_h = map(float, (d.mediabox.width, d.mediabox.height))

def label_abc(c):
    # Erase only the three original lowercase panel glyphs, then use uppercase.
    for x, y, new in [(119.42681, 368.87257, 'A'),
                       (378.48634, 365.03522, 'B'),
                       (378.48634, 162.21897, 'C')]:
        c.setFillColorRGB(1, 1, 1)
        c.rect(x-0.8, y-1.5, 9.0, 11.0, fill=1, stroke=0)
        c.setFillColorRGB(0, 0, 0)
        c.setFont('ArialBold', 9)
        c.drawString(x, y, new)

abc.merge_page(make_overlay(abc_w, abc_h, label_abc))
crop_bottom = abc_h - 363.0
abc.cropbox = RectangleObject([0, crop_bottom, abc_w, abc_h])
abc_trimmed = PageObject.create_blank_page(width=abc_w, height=363)
abc_trimmed.merge_transformed_page(abc, Transformation().translate(0, -crop_bottom))

def label_d(c):
    # Move the existing small title right by 8 pt to make space for panel D.
    # All data, legend, axes, and statistical text are copied without changes.
    c.setFillColorRGB(1, 1, 1)
    c.rect(11.8, 299.9, 238, 12, fill=1, stroke=0)
    c.setFillColorRGB(0, 0, 0)
    c.setFont('ArialBold', 9)
    c.drawString(7.2, 301.7709, 'D')
    c.setFont('ArialBold', 8.5)
    c.drawString(21.0, 301.7709, 'Composite and component endpoints')

d.merge_page(make_overlay(d_w, d_h, label_d))

scale = d_w / abc_w
gap_pt = 6.0
height = d_h + gap_pt + 363*scale
combined = PageObject.create_blank_page(width=d_w, height=height)
combined.merge_page(d)
combined.merge_transformed_page(abc_trimmed, Transformation().scale(scale).translate(0, d_h+gap_pt))
writer = PdfWriter()
writer.add_page(combined)
pdf_path = WORK/'Figure_S10_ABCD_vector_internal.pdf'
with pdf_path.open('wb') as handle:
    writer.write(handle)

poppler = Path(r'C:/Users/ASUS/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe')
prefix = WORK/'Figure_S10_ABCD_600dpi'
subprocess.run([str(poppler), '-png', '-r', '600', '-singlefile', str(pdf_path), str(prefix)], check=True)
final = DELIVERY/'Figure_S10_ABCD.png'
with Image.open(prefix.with_suffix('.png')) as im:
    im.convert('RGB').save(final, dpi=(600, 600), optimize=True)
    preview = im.copy()
    preview.thumbnail((1600, 2300))
    preview.convert('RGB').save(WORK/'composition_preview.png')

after = {str(p): sha(p) for p in (ABC, D)}
assert before == after
assert [p.name for p in DELIVERY.iterdir()] == ['Figure_S10_ABCD.png']
with Image.open(final) as im:
    report = dict(delivered_file=str(final), pixels=im.size, dpi=im.info.get('dpi'),
        layout='A left, B upper right, C lower right; full-width D below',
        width_mm=d_w/72*25.4, height_mm=height/72*25.4,
        original_ABC_scale=scale, minimum_font_pt=6.5*scale,
        source_hashes=before, source_files_unchanged=True,
        scientific_data_recomputed=False, original_ABC_footer_removed=True,
        only_final_image_in_delivery_folder=True, visual_QA='pending')
    (WORK/'composition_QA.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(report, ensure_ascii=False))
