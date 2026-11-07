# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) seed-stability interval plot -> assets/figures/Forest -> param inherit
# (b) matched-control strip plot -> assets/figures/PairedBoxScatter -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.

# Academic Figure Skill Typography Baseline — COPY VERBATIM
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 8,
    "figure.titlesize": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "legend.frameon": False,
})

# Academic Figure Skill Nature/Cell/Science Color Palette — COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED, GREY, BLACK = "#B2182B", "#999999", "#222222"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none", "savefig.bbox": "tight", "savefig.dpi": 300})
def save_cns_figure(fig, filename):
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

PROJECT = Path(r"I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026")
SOURCE = PROJECT / "figures/source_data/FigureS9_ITGB2_directional_connectivity"
OUT = PROJECT / "figures/supplementary/Supplementary_Figure_S9_ITGB2_directional_connectivity"
PROGRAM_ORDER = [
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
    "GOBP_PHAGOCYTOSIS",
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
    "GOBP_CELL_MATRIX_ADHESION",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "ID_INFLAMMATORY_MONO",
]
LABELS = {
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION": "Myeloid migration",
    "GOBP_PHAGOCYTOSIS": "Phagocytosis",
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX": "ECM degradation",
    "GOBP_CELL_MATRIX_ADHESION": "Matrix adhesion",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNF–NF-κB",
    "ID_INFLAMMATORY_MONO": "Inflammatory monocyte",
}

def label(ax, s):
    ax.text(-0.12, 1.03, s, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", color=BLACK)

def quiet(ax):
    ax.spines["left"].set_color("#777777"); ax.spines["bottom"].set_color("#777777")
    ax.tick_params(color="#777777", labelcolor=BLACK)

def load():
    by_seed = pd.read_csv(SOURCE / "FigureS9A_ITGB2_by_seed.tsv", sep="\t")
    all_nodes = pd.read_csv(SOURCE / "FigureS9B_ITGB2_vs_controls.tsv", sep="\t")
    summ = pd.read_csv(SOURCE / "FigureS9B_summary.tsv", sep="\t")
    return by_seed, all_nodes, summ

def draw_a(ax, by_seed):
    label(ax, "a")
    q = by_seed.copy(); q["program"] = pd.Categorical(q.program, PROGRAM_ORDER, ordered=True)
    med = q.groupby("program", observed=True).mean(numeric_only=True).sort_index()
    y = np.arange(len(PROGRAM_ORDER))
    for yi, program in enumerate(PROGRAM_ORDER):
        z = q[q.program == program].sort_values("mean_percent_change_total_abs_strength")
        ax.scatter(z.mean_percent_change_total_abs_strength, [yi-0.05, yi, yi+0.05], s=18,
                   color="#6BAED6", edgecolor="white", linewidth=0.35, zorder=2)
        ax.scatter(med.loc[program, "mean_percent_change_total_abs_strength"], yi, marker="D", s=36,
                   color=ACCENT_RED, edgecolor="white", linewidth=0.55, zorder=4)
    ax.axvline(0, color="#777777", lw=0.7, ls=(0, (3, 2)))
    ax.set_yticks(y); ax.set_yticklabels([LABELS[p] for p in PROGRAM_ORDER])
    ax.set_xlabel("KO − WT change in network strength (%)")
    ax.set_xlim(-0.34, 0.06)
    ax.text(0.99, 1.01, "negative values = connectivity loss", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.4, color="#666666")
    quiet(ax)

def draw_b(ax, all_nodes, summ):
    label(ax, "b")
    y = np.arange(len(PROGRAM_ORDER))
    for yi, program in enumerate(PROGRAM_ORDER):
        c = all_nodes[(all_nodes.role == "matched_control") & (all_nodes.program == program)]
        t = all_nodes[(all_nodes.role == "target") & (all_nodes.program == program)].iloc[0]
        ax.scatter(c.percent_change_total_abs_strength, [yi + j for j in np.linspace(-0.17, 0.17, len(c))],
                   s=15, color=GREY, alpha=0.72, edgecolor="none", zorder=2)
        ax.scatter(t.percent_change_total_abs_strength, yi, marker="D", s=38, color=ACCENT_RED,
                   edgecolor="white", linewidth=0.55, zorder=4)
    ax.axvline(0, color="#777777", lw=0.7, ls=(0, (3, 2)))
    ax.set_yticks(y); ax.set_yticklabels([LABELS[p] for p in PROGRAM_ORDER])
    ax.set_xlabel("KO − WT change in network strength (%)")
    ax.set_xlim(-7.0, 0.8)
    ax.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=GREY, ms=4, label="20 matched nodes"),
                       Line2D([0], [0], marker="D", color="none", markerfacecolor=ACCENT_RED, markeredgecolor="white", ms=5, label="ITGB2")],
              loc="lower right", fontsize=5.2, handletextpad=0.25, borderaxespad=0.2)
    quiet(ax)

def build():
    by_seed, all_nodes, summ = load()
    fig = plt.figure(figsize=(183/25.4, 95/25.4), facecolor="white")
    gs = fig.add_gridspec(1, 2, left=0.16, right=0.98, top=0.92, bottom=0.20, wspace=0.72)
    draw_a(fig.add_subplot(gs[0, 0]), by_seed)
    draw_b(fig.add_subplot(gs[0, 1]), all_nodes, summ)
    with mpl.rc_context({"savefig.bbox": None}):
        fig.savefig(f"{OUT}.pdf", bbox_inches=None, dpi=300)
        fig.savefig(f"{OUT}.png", bbox_inches=None, dpi=300)
        fig.savefig(f"{OUT}.svg", bbox_inches=None)
        fig.savefig(f"{OUT}.tiff", bbox_inches=None, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    return OUT

if __name__ == "__main__":
    print(build())
