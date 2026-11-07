# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) NES heatmap -> assets/figures/heatmap -> param inherit
# (b) running-enrichment curves -> assets/figures/LineTrend -> param inherit
# (c) leading-edge recurrence dot heatmap -> assets/figures/heatmap -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.

# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})

# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,         # TrueType font embedding
    "svg.fonttype": "none",     # editable text in SVG
    "savefig.bbox": "tight",    # trim whitespace
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

from pathlib import Path
import json
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


PROJECT = Path(r"I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026")
GSEA_DIR = PROJECT / "analysis/results/virtual_knockout/GSE183904/C1QC_APOE_macrophage_receptor_screen_v2/ITGB2_GSEA_v1"
DIFF_PATH = PROJECT / "figures/source_data/Figure7_receiver_screen_v3/Figure7B_ITGB2_diffreg_all_seeds.tsv"
PROGRAM_PATH = PROJECT / "resources/gene_sets/GSE183904_C1QC_APOE_receptor_screen_v2_scTenifoldKnk_programs_axis_excluded.tsv"
SOURCE_DIR = PROJECT / "figures/source_data/FigureS10_ITGB2_diffreg_GSEA"
OUT = PROJECT / "figures/supplementary/Supplementary_Figure_S10_ITGB2_diffreg_GSEA"
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

PROGRAM_ORDER = [
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
    "GOBP_PHAGOCYTOSIS",
    "ID_INFLAMMATORY_MONO",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "GOBP_CELL_MATRIX_ADHESION",
]
FOCAL_PROGRAMS = PROGRAM_ORDER[:3]
LABELS = {
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX": "ECM degradation",
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION": "Myeloid migration",
    "GOBP_PHAGOCYTOSIS": "Phagocytosis",
    "ID_INFLAMMATORY_MONO": "Inflammatory monocyte",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNF–NF-κB",
    "GOBP_CELL_MATRIX_ADHESION": "Matrix adhesion",
}
SEED_COLORS = {
    2026072739: "#2166AC",
    2026072740: "#1B7837",
    2026072741: "#F1A340",
}
SEED_LABELS = {
    2026072739: "seed 1",
    2026072740: "seed 2",
    2026072741: "seed 3",
}
DIVERGING_CMAP = LinearSegmentedColormap.from_list("cns_diverging", DIVERGING, N=256)


def panel_label(ax, value, x=-0.12, y=1.04):
    ax.text(x, y, value, transform=ax.transAxes, fontsize=9, fontweight="bold",
            ha="left", va="bottom", color=BLACK)


def quiet_axis(ax):
    ax.spines["left"].set_color("#777777")
    ax.spines["bottom"].set_color("#777777")
    ax.tick_params(color="#777777", labelcolor=BLACK)


def enrichment_running(scores, gene_set):
    genes = np.asarray(scores.index)
    vals = scores.to_numpy(dtype=float)
    hit = np.isin(genes, list(gene_set))
    nh = int(hit.sum())
    weights = np.abs(vals) * hit
    denom = weights.sum()
    if denom <= 0:
        weights = hit.astype(float)
        denom = weights.sum()
    running = np.cumsum(weights / denom) - np.cumsum((~hit) / (len(genes) - nh))
    return running, np.flatnonzero(hit)


def load_data():
    result = pd.read_csv(GSEA_DIR / "ITGB2_diffreg_GSEA_by_seed.tsv", sep="\t")
    leading = pd.read_csv(GSEA_DIR / "ITGB2_diffreg_GSEA_leading_edge.tsv", sep="\t")
    diff = pd.read_csv(DIFF_PATH, sep="\t")
    programmes = pd.read_csv(PROGRAM_PATH, sep="\t")
    diff["gene"] = diff["gene"].str.upper()
    programmes["gene_symbol"] = programmes["gene_symbol"].str.upper()
    gene_sets = {p: set(g.gene_symbol.dropna()) for p, g in programmes.groupby("program")}
    return result, leading, diff, gene_sets


def validate_data(result, leading, diff, gene_sets):
    expected_seeds = {2026072739, 2026072740, 2026072741}
    if set(result.seed) != expected_seeds:
        raise ValueError("GSEA seed universe differs from the three frozen Figure 7 seeds")
    if len(result) != 18 or result.program.nunique() != 6:
        raise ValueError("Expected 18 seed-programme tests across six frozen programmes")
    if result.NES.isna().any() or not np.isfinite(result.NES).all():
        raise ValueError("Non-finite NES values")
    if diff[diff.gene != "ITGB2"].groupby("seed").gene.nunique().nunique() != 1:
        raise ValueError("Ranked-gene universe is inconsistent across seeds")
    for p in PROGRAM_ORDER:
        if p not in gene_sets:
            raise ValueError(f"Missing frozen programme: {p}")
    if leading.empty:
        raise ValueError("Leading-edge table is empty")


def make_source_data(result, leading, diff, gene_sets):
    result.to_csv(SOURCE_DIR / "FigureS10a_seed_programme_NES.tsv", sep="\t", index=False)
    curve_rows = []
    for seed in sorted(result.seed.unique()):
        d = diff[(diff.seed == seed) & (diff.gene != "ITGB2")].drop_duplicates("gene")
        d = d[np.isfinite(d.Z)].sort_values("Z", ascending=False)
        scores = pd.Series(d.Z.to_numpy(float), index=d.gene.to_numpy())
        for program in FOCAL_PROGRAMS:
            running, hits = enrichment_running(scores, gene_sets[program] & set(scores.index))
            hit_set = set(hits.tolist())
            for rank, (gene, z, es) in enumerate(zip(scores.index, scores.values, running), 1):
                curve_rows.append({
                    "seed": int(seed), "program": program, "rank": rank,
                    "gene": gene, "Z": float(z), "running_ES": float(es),
                    "is_gene_set_member": (rank - 1) in hit_set,
                })
    curve = pd.DataFrame(curve_rows)
    curve.to_csv(SOURCE_DIR / "FigureS10b_running_enrichment_curves.tsv.gz", sep="\t", index=False, compression="gzip")

    recurrence = (leading.drop_duplicates(["seed", "program", "gene"])
                  .groupby(["program", "gene"]).seed.nunique()
                  .rename("leading_edge_seed_count").reset_index())
    focal = recurrence[recurrence.program.isin(FOCAL_PROGRAMS)].copy()
    cross = (focal.assign(present=focal.leading_edge_seed_count > 0)
             .groupby("gene").present.sum())
    selected = sorted(cross[cross >= 2].index)
    matrix = (focal[focal.gene.isin(selected)]
              .pivot(index="gene", columns="program", values="leading_edge_seed_count")
              .reindex(columns=FOCAL_PROGRAMS).fillna(0).astype(int))
    matrix["total"] = matrix.sum(axis=1)
    matrix = matrix.sort_values(["total"] + FOCAL_PROGRAMS, ascending=False).drop(columns="total")
    matrix.reset_index().to_csv(SOURCE_DIR / "FigureS10c_cross_programme_leading_edge.tsv", sep="\t", index=False)
    return curve, matrix


def draw_a(ax, result):
    panel_label(ax, "a", -0.20)
    seeds = sorted(result.seed.unique())
    pivot = (result.pivot(index="program", columns="seed", values="NES")
             .reindex(index=PROGRAM_ORDER, columns=seeds))
    norm = TwoSlopeNorm(vmin=-2.0, vcenter=0.0, vmax=2.0)
    im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap=DIVERGING_CMAP, norm=norm)
    for i, program in enumerate(PROGRAM_ORDER):
        for j, seed in enumerate(seeds):
            r = result[(result.program == program) & (result.seed == seed)].iloc[0]
            color = "white" if abs(r.NES) >= 1.35 else BLACK
            ax.text(j, i, f"{r.NES:.2f}", ha="center", va="center", fontsize=6.1,
                    color=color, fontweight="bold" if r.FDR_global < 0.05 else "normal")
            if r.FDR_global < 0.05:
                ax.scatter(j + 0.34, i - 0.34, s=8, color=BLACK, marker="o", zorder=4)
    ax.set_xticks(range(len(seeds)))
    ax.set_xticklabels([SEED_LABELS[s] for s in seeds])
    ax.set_yticks(range(len(PROGRAM_ORDER)))
    ax.set_yticklabels([LABELS[p] for p in PROGRAM_ORDER])
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.045, pad=0.035)
    cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
    cb.ax.set_title("NES", fontsize=6.5, pad=3)
    cb.ax.tick_params(labelsize=6, width=0.5, length=2)
    ax.text(0.0, -0.14, "dot: global BH FDR < 0.05 for the individual seed–programme test",
            transform=ax.transAxes, fontsize=5.4, color=BLACK, ha="left", va="top")


def draw_running_axis(ax, curve, result, program, add_legend=False):
    q = curve[curve.program == program]
    for row, seed in enumerate(sorted(q.seed.unique())):
        z = q[q.seed == seed].sort_values("rank")
        ax.plot(z["rank"], z["running_ES"], color=SEED_COLORS[seed], lw=1.05,
                label=SEED_LABELS[seed])
        hits = z.loc[z.is_gene_set_member, "rank"].to_numpy()
        rug_y = -0.23 - row * 0.035
        ax.vlines(hits, rug_y - 0.012, rug_y + 0.012, color=SEED_COLORS[seed],
                  lw=0.28, alpha=0.45)
    ax.axhline(0, color="#999999", lw=0.55, ls=(0, (3, 2)))
    ax.set_ylim(-0.35, 0.74)
    ax.set_xlim(1, int(q["rank"].max()))
    ax.set_title(LABELS[program], fontsize=7.2, fontweight="normal", pad=3)
    ax.set_xlabel("Ranked genes")
    if program == FOCAL_PROGRAMS[0]:
        ax.set_ylabel("Running ES")
    else:
        ax.set_yticklabels([])
        ax.set_ylabel("")
    ax.set_xticks([1, 250, 500, 774])
    ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
    quiet_axis(ax)
    if add_legend:
        ax.legend(loc="upper right", fontsize=5.2, handlelength=1.7,
                  handletextpad=0.35, borderaxespad=0.2)
    med = result[result.program == program].NES.median()
    ax.text(0.98, 0.04, f"median NES {med:.2f}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=5.3, color=BLACK)


def draw_c(ax, matrix):
    panel_label(ax, "c", -0.055, 1.02)
    genes = list(matrix.index)
    programmes = FOCAL_PROGRAMS
    ax.set_xlim(-0.6, len(programmes) - 0.4)
    ax.set_ylim(len(genes) - 0.5, -0.5)
    for yi, gene in enumerate(genes):
        for xi, program in enumerate(programmes):
            n = int(matrix.loc[gene, program])
            if n > 0:
                color = CATEGORICAL[xi]
                ax.scatter(xi, yi, s={1: 20, 2: 48, 3: 86}[n], color=color,
                           edgecolor="white", linewidth=0.5)
    ax.set_xticks(range(len(programmes)))
    ax.set_xticklabels([LABELS[p] for p in programmes])
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=6.0)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.grid(axis="y", color="#E6E6E6", lw=0.35)
    ax.set_axisbelow(True)
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor="#777777",
                      markeredgecolor="white", markersize={1: 3.4, 2: 5.0, 3: 6.6}[n],
                      label=f"{n}/3 seeds") for n in [1, 2, 3]]
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 0.995),
              fontsize=5.5, borderaxespad=0, handletextpad=0.4)


def build():
    result, leading, diff, gene_sets = load_data()
    validate_data(result, leading, diff, gene_sets)
    curve, matrix = make_source_data(result, leading, diff, gene_sets)

    fig = plt.figure(figsize=(183 / 25.4, 168 / 25.4), facecolor="white")
    outer = fig.add_gridspec(2, 12, height_ratios=[0.92, 1.28],
                             left=0.16, right=0.975, top=0.96, bottom=0.10,
                             wspace=1.25, hspace=0.48)
    ax_a = fig.add_subplot(outer[0, 0:5])
    draw_a(ax_a, result)

    b_grid = outer[0, 6:12].subgridspec(1, 3, wspace=0.40)
    b_axes = [fig.add_subplot(b_grid[0, i]) for i in range(3)]
    panel_label(b_axes[0], "b", -0.35)
    for i, (ax, program) in enumerate(zip(b_axes, FOCAL_PROGRAMS)):
        draw_running_axis(ax, curve, result, program, add_legend=(i == 2))

    ax_c = fig.add_subplot(outer[1, 0:12])
    draw_c(ax_c, matrix)

    with mpl.rc_context({"savefig.bbox": None}):
        fig.savefig(f"{OUT}.pdf", bbox_inches=None, dpi=300)
        fig.savefig(f"{OUT}.png", bbox_inches=None, dpi=300)
        fig.savefig(f"{OUT}.svg", bbox_inches=None)
        fig.savefig(f"{OUT}.tiff", bbox_inches=None, dpi=600,
                    pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    shutil.copy2(GSEA_DIR / "ITGB2_diffreg_GSEA_manifest.json",
                 SOURCE_DIR / "ITGB2_diffreg_GSEA_manifest.json")
    manifest = {
        "figure": "Supplementary Figure S10",
        "generated": "2026-08-13",
        "panels": {
            "a": "all 18 seed-programme NES values",
            "b": "running enrichment curves for three prespecified focal programmes and all seeds",
            "c": "leading-edge genes recurring across at least two focal programmes",
        },
        "excluded_from_ranking": "ITGB2 self-hit",
        "n_ranked_genes_per_seed": 774,
        "permutations_per_test": 2000,
    }
    (SOURCE_DIR / "FigureS10_source_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    return OUT


if __name__ == "__main__":
    print(build())
