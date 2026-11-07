# Academic Figure Skill Asset Confirmation (verified against project production assets)
# Figure 3 -> figures/gen_fig3_reference_communication_layout.py -> native production run with annotation-only patch
# RULE: data loading, panel geometry, colors, and statistical values remain in the authoritative source script.

from pathlib import Path

from PIL import Image


ROOT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
SOURCE_SCRIPT = ROOT / "figures" / "gen_fig3_reference_communication_layout.py"
OUTPUT_DIR = ROOT / "tmp" / "wjgo_figures_candidate"
OUTPUT_STEM = "Figure_3_WJGO_candidate"


def replace_once(source, old, new, description):
    count = source.count(old)
    if count != 1:
        raise RuntimeError(
            f"Expected one {description} block in {SOURCE_SCRIPT}, found {count}."
        )
    return source.replace(old, new, 1)


def convert_tiff_to_rgb(path):
    with Image.open(path) as image:
        if image.mode == "RGBA":
            rgb = Image.new("RGB", image.size, "white")
            rgb.paste(image, mask=image.getchannel("A"))
        else:
            rgb = image.convert("RGB")
        rgb.save(path, format="TIFF", compression="tiff_lzw", dpi=(600, 600))


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    source = SOURCE_SCRIPT.read_text(encoding="utf-8")

    source = replace_once(
        source,
        'FIG_DIR = PROJECT / "figures" / "main"',
        'FIG_DIR = PROJECT / "tmp" / "wjgo_figures_candidate"',
        "Figure 3 output-directory",
    )
    source = replace_once(
        source,
        'OUT_DIR = PROJECT / "figures" / "data" / "Figure3_PLAT_reference_layout_v4"',
        'OUT_DIR = PROJECT / "tmp" / "wjgo_figures_candidate" / "Figure3_source_data_audit"',
        "Figure 3 audit-directory",
    )
    source = replace_once(
        source,
        'stem = "Figure3_PLAT_reference_communication_layout_v4_candidate"',
        f'stem = "{OUTPUT_STEM}"',
        "Figure 3 output-stem",
    )
    source = replace_once(
        source,
        '''for y in range(len(pair_order)):
    for x in range(5):
        if bool(star_matrix.iloc[y, x]):
            ax_dheat.text(x, y, "*", ha="center", va="center",
                          fontsize=7.0, fontweight="bold", color=BLACK, zorder=4)''',
        '''for y in range(len(pair_order)):
    for x in range(5):
        if bool(star_matrix.iloc[y, x]):
            ax_dheat.add_patch(Rectangle(
                (x - 0.47, y - 0.47), 0.94, 0.94,
                fill=False, edgecolor=BLACK, linewidth=0.82, zorder=4,
            ))''',
        "Figure 3 panel-d rank-gate annotation",
    )
    source = replace_once(
        source,
        '''    star = " *" if fdr < 0.05 else ""
    ax_feco.text(0.84, y, f"FDR {fdr:.2e}{star}",''',
        '''    mantissa, exponent = f"{fdr:.2e}".split("e")
    fdr_label = rf"{float(mantissa):.2f}$\\times 10^{{{int(exponent)}}}$"
    ax_feco.text(0.84, y, fdr_label,''',
        "Figure 3 panel-f FDR annotation",
    )
    source = replace_once(
        source,
        '''eco_order = [1, 3, 4, 5]
for y, cluster in enumerate(eco_order):''',
        '''eco_order = [1, 3, 4, 5]
ax_feco.text(0.985, 1.015, "BH FDR", transform=ax_feco.transAxes,
             ha="right", va="bottom", fontsize=4.7, color="#000000")
for y, cluster in enumerate(eco_order):''',
        "Figure 3 panel-f FDR column heading",
    )
    source = replace_once(
        source,
        '"starred_liana_cells": int(star_matrix.to_numpy().sum()),',
        '"outlined_liana_combined_gate_cells": int(star_matrix.to_numpy().sum()),',
        "Figure 3 QA rank-gate label",
    )

    namespace = {"__name__": "__main__", "__file__": str(SOURCE_SCRIPT)}
    exec(compile(source, str(SOURCE_SCRIPT), "exec"), namespace)
    convert_tiff_to_rgb(OUTPUT_DIR / f"{OUTPUT_STEM}.tiff")


if __name__ == "__main__":
    main()
