"""Render Supplementary Figure S7 page 2: PLAT sampling-location sensitivity."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paper_plot_style import PROJECT, clean_axis, panel_label  # noqa: E402


FIG_DATA = PROJECT / "figures" / "data"
SUMMARY = FIG_DATA / "Figure6_PLAT_core_edge_summary_20260809.tsv"
OUT_DIR = PROJECT / "figures" / "ejc_revision_20260827"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PLAT = "#8F5F6B"
PLAT_LIGHT = "#C8A5AE"
BLUE = "#456F9D"
BLUE_LIGHT = "#9CB5D3"
GREY = "#707980"
LIGHT_GREY = "#D7DDE1"
INK = "#252525"
LOCATIONS = ["Core", "In-between", "Edge"]
LOCATION_LABELS = ["Core", "In-between", "Edge"]


def add_jitter(ax, x, values, color, marker, label, seed):
    """Add deterministic small-jitter points when section-level values exist."""
    rng = np.random.default_rng(seed)
    xj = x + rng.uniform(-0.055, 0.055, len(values))
    ax.scatter(xj, values, s=15, marker=marker, color=color, alpha=0.82,
               edgecolor="white", linewidth=0.45, label=label, zorder=3)


def build_figure():
    summary = pd.read_csv(SUMMARY, sep="\t")
    summary["dataset"] = summary["dataset"].replace({
        "Main_ST": "GSE270678 Main ST",
        "ST_R_sensitivity": "GSE270678 ST_R sensitivity",
    })
    summary["location"] = pd.Categorical(summary["location"], LOCATIONS, ordered=True)
    summary = summary.sort_values(["dataset", "location"])
    summary.to_csv(FIG_DATA / "FigureS7G_core_edge_summary_source.tsv", sep="\t", index=False)

    # These are the donor-paired differences from the locked analysis.  Only
    # D26 and D27 have both Core and Edge sections in each dataset.
    paired = pd.DataFrame([
        {"dataset": "GSE270678 Main ST", "donor": "D26", "difference": 0.0786},
        {"dataset": "GSE270678 Main ST", "donor": "D27", "difference": -0.0838},
        {"dataset": "GSE270678 ST_R sensitivity", "donor": "D26", "difference": -0.0241},
        {"dataset": "GSE270678 ST_R sensitivity", "donor": "D27", "difference": -0.0678},
    ])
    paired.to_csv(FIG_DATA / "FigureS7I_core_edge_paired_source.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(7.2, 4.45), facecolor="white")
    outer = fig.add_gridspec(2, 1, height_ratios=[1.10, 0.92], hspace=0.48)
    top = outer[0].subgridspec(1, 2, width_ratios=[1.0, 0.90], wspace=0.42)
    bottom = outer[1].subgridspec(1, 2, width_ratios=[1.15, 0.85], wspace=0.48)

    # a | Section-level aggregate PLAT intensity.
    ax_a = fig.add_subplot(top[0, 0])
    panel_label(ax_a, "a", -0.16, 1.04)
    for ds, color, marker in [
        ("GSE270678 Main ST", PLAT, "o"),
        ("GSE270678 ST_R sensitivity", BLUE, "s"),
    ]:
        d = summary.loc[summary.dataset.eq(ds)].set_index("location").loc[LOCATIONS]
        ax_a.plot(range(3), d["mean_spot_log1p_CP10K"], color=color, lw=1.2,
                  marker=marker, ms=4.5, mec="white", mew=0.55,
                  label=("GSE270678 Main ST\n(33 sections; 19 donors)"
                         if ds == "GSE270678 Main ST" else
                         "GSE270678 ST_R sensitivity\n(37 sections; 16 donors)"), zorder=3)
    ax_a.set_xticks(range(3), LOCATION_LABELS)
    ax_a.set_ylabel("mean PLAT intensity\nlog1p(CP10K)", fontsize=6.2, labelpad=2)
    ax_a.set_ylim(0, 0.52)
    ax_a.tick_params(labelsize=5.4)
    ax_a.legend(frameon=False, fontsize=5.0, loc="upper left", handlelength=1.3,
                borderaxespad=0.2, labelspacing=0.45)
    clean_axis(ax_a)

    # b | Detection fraction, displayed as a companion metric.
    ax_b = fig.add_subplot(top[0, 1])
    panel_label(ax_b, "b", -0.19, 1.04)
    for ds, color, marker in [
        ("GSE270678 Main ST", PLAT, "o"),
        ("GSE270678 ST_R sensitivity", BLUE, "s"),
    ]:
        d = summary.loc[summary.dataset.eq(ds)].set_index("location").loc[LOCATIONS]
        ax_b.plot(range(3), 100 * d["mean_PLAT_detection_fraction"], color=color, lw=1.2,
                  marker=marker, ms=4.5, mec="white", mew=0.55,
                  label=("GSE270678 Main ST\n(33 sections; 19 donors)"
                         if ds == "GSE270678 Main ST" else
                         "GSE270678 ST_R sensitivity\n(37 sections; 16 donors)"), zorder=3)
    ax_b.set_xticks(range(3), LOCATION_LABELS)
    ax_b.set_ylabel("PLAT detection\n(% in-tissue spots)", fontsize=6.2, labelpad=2)
    ax_b.set_ylim(0, 43)
    ax_b.tick_params(labelsize=5.4)
    ax_b.legend(frameon=False, fontsize=5.0, loc="upper left", handlelength=1.3,
                borderaxespad=0.2, labelspacing=0.45)
    clean_axis(ax_b)

    # c | Donor-paired Core-to-Edge contrasts.
    ax_c = fig.add_subplot(bottom[0, 0])
    panel_label(ax_c, "c", -0.16, 1.04)
    groups = [("GSE270678 Main ST", PLAT), ("GSE270678 ST_R sensitivity", BLUE)]
    y_positions = {
        ("GSE270678 Main ST", "D26"): 3,
        ("GSE270678 Main ST", "D27"): 2,
        ("GSE270678 ST_R sensitivity", "D26"): 1,
        ("GSE270678 ST_R sensitivity", "D27"): 0,
    }
    for ds, color in groups:
        d = paired.loc[paired.dataset.eq(ds)].copy()
        for _, row in d.iterrows():
            y = y_positions[(ds, row.donor)]
            ax_c.hlines(y, 0, row.difference, color=color, lw=1.35, alpha=0.9, zorder=2)
            ax_c.scatter(row.difference, y, s=26, color=color,
                         edgecolor="white", linewidth=0.5, zorder=3)
    ax_c.axvline(0, color=GREY, lw=0.65, ls=(0, (3, 2)), zorder=0)
    ax_c.set_yticks([3, 2, 1, 0], ["Main ST · D26", "Main ST · D27", "ST_R · D26", "ST_R · D27"])
    ax_c.set_xlabel("Edge − Core PLAT intensity", fontsize=6.2, labelpad=2)
    ax_c.set_xlim(-0.11, 0.105)
    ax_c.set_ylim(-0.65, 3.65)
    ax_c.tick_params(labelsize=5.4)
    ax_c.legend(handles=[
        Line2D([0], [0], color=PLAT, marker="o", lw=1.3, markersize=3.6, label="Main ST"),
        Line2D([0], [0], color=BLUE, marker="s", lw=1.3, markersize=3.6, label="ST_R sensitivity"),
    ], frameon=False, fontsize=4.8, loc="upper left", handlelength=1.2,
       borderaxespad=0.2, labelspacing=0.35)
    clean_axis(ax_c)

    # d | Effect-size summary with explicit zero reference and P values.
    ax_d = fig.add_subplot(bottom[0, 1])
    panel_label(ax_d, "d", -0.20, 1.04)
    effects = [
        ("Main ST sections", 0.4357103 - 0.2518455, "P = 0.3969", PLAT),
        ("Main ST paired D26/D27", (0.0786 - 0.0838) / 2, "P = 1.00", GREY),
        ("ST_R paired D26/D27", (-0.0241 - 0.0678) / 2, "P = 0.50", BLUE),
    ]
    y = np.arange(len(effects))[::-1]
    ax_d.axvline(0, color=GREY, lw=0.65, ls=(0, (3, 2)), zorder=0)
    for yi, (label, estimate, pval, color) in zip(y, effects):
        ax_d.hlines(yi, 0, estimate, color=color, lw=1.45, zorder=2)
        ax_d.scatter(estimate, yi, s=28, color=color, edgecolor="white", linewidth=0.5, zorder=3)
        ax_d.text(0.285, yi, pval, fontsize=5.2, color=INK, va="center", ha="right")
    ax_d.set_yticks(y, [x[0] for x in effects])
    ax_d.set_xlabel("Edge − Core PLAT intensity", fontsize=5.9, labelpad=2)
    ax_d.set_xlim(-0.14, 0.30)
    ax_d.set_ylim(-0.55, 2.55)
    ax_d.tick_params(labelsize=5.1)
    clean_axis(ax_d)

    fig.subplots_adjust(left=0.08, right=0.98, top=0.97, bottom=0.10)
    base = OUT_DIR / "Supplementary_Figure_S9_core_edge_PLAT_ejc_20260827"
    outputs = []
    for suffix, dpi in [(".pdf", 300), (".svg", 300), (".png", 400), (".tiff", 600)]:
        path = base.with_suffix(suffix)
        kwargs = {"dpi": dpi, "bbox_inches": "tight", "pad_inches": 0.04}
        if suffix == ".tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(path, **kwargs)
        outputs.append(path)
    plt.close(fig)
    return tuple(outputs)


if __name__ == "__main__":
    print("\n".join(str(x) for x in build_figure()))
