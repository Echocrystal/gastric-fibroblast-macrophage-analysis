# Academic Figure Skill Asset Confirmation (verified against project production assets)
# Figure 6 -> figures/gen_figure6_PLAT_bulk_groupwise_clinical_drug_v22_readability_20260827.py -> native production run with annotation-only patch
# RULE: data loading, panel geometry, colors, and statistical values remain in the authoritative source script.

from pathlib import Path
import sys

from PIL import Image


ROOT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
SOURCE_SCRIPT = (
    ROOT
    / "figures"
    / "gen_figure6_PLAT_bulk_groupwise_clinical_drug_v22_readability_20260827.py"
)
OUTPUT_DIR = ROOT / "tmp" / "wjgo_figures_candidate"
OUTPUT_STEM = "Figure_6_WJGO_candidate"


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
    sys.path.insert(0, str(ROOT / "figures"))
    source = SOURCE_SCRIPT.read_text(encoding="utf-8")

    source = replace_once(
        source,
        'OUT = ROOT / "figures" / "main" / "revised_20260827" / "Figure_6_readability_v22"',
        f'OUT = ROOT / "tmp" / "wjgo_figures_candidate" / "{OUTPUT_STEM}"',
        "Figure 6 output-stem",
    )
    source = replace_once(
        source,
        'SOURCE_V18 = ROOT / "figures" / "source_data" / "Figure6_v19_readability_20260827"',
        'SOURCE_V18 = ROOT / "tmp" / "wjgo_figures_candidate" / "Figure6_source_data_audit"',
        "Figure 6 audit-directory",
    )
    source = replace_once(
        source,
        '''def star(q):
    if q < 1e-4: return "****"
    if q < 1e-3: return "***"
    if q < 1e-2: return "**"
    if q < 5e-2: return "*"
    return "ns"


def bracket_label(q):
    """Reference-style label: stars for significant pairs and lowercase ns otherwise."""
    return star(q) if q < 0.05 else "ns"''',
        '''def significance_letter(q):
    """WJGO-style superscript code for significant adjusted q values only."""
    q = float(q)
    if q < 1e-4:
        return r"$^{d}$"
    if q < 1e-3:
        return r"$^{c}$"
    if q < 1e-2:
        return r"$^{b}$"
    if q < 5e-2:
        return r"$^{a}$"
    return None''',
        "Figure 6 significance-code helper",
    )
    source = replace_once(
        source,
        '''    # Asterisks have a lower visible glyph edge than lowercase letters and
    # need slightly more clearance to avoid touching the bracket cap.
    text_gap = .0008 if label == "ns" else .012
    ax.text((x1 + x2) / 2, y + height + text_gap, label, ha="center", va="bottom",
            fontsize=fontsize, color=INK, fontweight="bold", clip_on=False)''',
        '''    ax.text((x1 + x2) / 2, y + height + .004, label, ha="center", va="bottom",
            fontsize=fontsize, color=INK, fontweight="bold", clip_on=False)''',
        "Figure 6 bracket-label placement",
    )
    source = replace_once(
        source,
        '''            q_pair = float(row.iloc[0].pairwise_FDR_BH)
            x1 = centers[i] + offsets[GROUPS.index(g1)]
            x2 = centers[i] + offsets[GROUPS.index(g2)]
            draw_bracket(ax, x1, x2, top + .045 + bracket_level * bracket_step,
                         bracket_height, bracket_label(q_pair))''',
        '''            q_pair = float(row.iloc[0].pairwise_FDR_BH)
            label = significance_letter(q_pair)
            if label is None:
                continue
            x1 = centers[i] + offsets[GROUPS.index(g1)]
            x2 = centers[i] + offsets[GROUPS.index(g2)]
            draw_bracket(ax, x1, x2, top + .045 + bracket_level * bracket_step,
                         bracket_height, label)''',
        "Figure 6 panels-b/c pairwise annotations",
    )
    source = replace_once(
        source,
        '''        ax.text(.985, yy, star(row.p_BH_all_specific_tests),
                transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=4.8, fontweight="bold", color=INK)''',
        '''        label = significance_letter(row.p_BH_all_specific_tests)
        if label is not None:
            ax.text(.985, yy, label,
                    transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=4.8, fontweight="bold", color=INK)''',
        "Figure 6 panel-e adjusted-q annotations",
    )
    source = replace_once(
        source,
        '''        stars = p_stars(row["p"])
        suffix = f" {stars}" if stars else ""
        ax.text(0.985, y[i], f"OR={row['OR']:.2f}{suffix}",
                transform=ax.get_yaxis_transform(), fontsize=4.85,
                ha="right", va="center", clip_on=True)''',
        '''        marker = r" $^{a}$" if float(row["p"]) < .05 else ""
        ax.text(0.985, y[i], f"OR={row['OR']:.2f}" + marker,
                transform=ax.get_yaxis_transform(), fontsize=4.85,
                ha="right", va="center", clip_on=True)''',
        "Figure 6 panel-h unadjusted-P annotation",
    )
    source = replace_once(
        source,
        '''            q_pair = float(row.iloc[0].q_BH_across_drugs_pairs)
            draw_bracket(
                ax, GROUPS.index(g1), GROUPS.index(g2),
                bracket_base + level * bracket_step, .005,
                bracket_label(q_pair), fontsize=3.5,
            )''',
        '''            q_pair = float(row.iloc[0].q_BH_across_drugs_pairs)
            label = significance_letter(q_pair)
            if label is None:
                continue
            draw_bracket(
                ax, GROUPS.index(g1), GROUPS.index(g2),
                bracket_base + level * bracket_step, .005,
                label, fontsize=3.5,
            )''',
        "Figure 6 panel-i pairwise annotations",
    )
    source = replace_once(
        source,
        '(REVIEW / "Figure6_v22_readability_20260827_QA.json").write_text(',
        '(OUT.parent / "Figure_6_WJGO_candidate_QA.json").write_text(',
        "Figure 6 QA output",
    )

    namespace = {"__name__": "__main__", "__file__": str(SOURCE_SCRIPT)}
    exec(compile(source, str(SOURCE_SCRIPT), "exec"), namespace)
    convert_tiff_to_rgb(OUTPUT_DIR / f"{OUTPUT_STEM}.tiff")


if __name__ == "__main__":
    main()
