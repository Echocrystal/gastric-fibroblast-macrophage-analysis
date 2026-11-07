# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# panel | chart type | production asset | inheritance mode | rationale
# a | annotated expression heatmap | heatmap/plot_composition.py | cross-type inherit | 60-gene row-z matrix plus six annotation strips requires custom compositing
# b | stability-gate small multiples | LinePlot/available line-chart assets | cross-type inherit | four prespecified metrics have distinct scientific thresholds
# c | complete signed checkpoint matrix | heatmap/plot_composition.py | cross-type inherit | signed beta values, missing VISTA coverage and adjusted estimates require a diverging matrix
# d | phenotype heatmap | heatmap/plot_composition.py | cross-type inherit | signed median z scores require a zero-centred scale
# e | paired coefficient forest | BarComparison/plot_comparison_GeneRegulatory.py | cross-type inherit | HC3 intervals require paired coefficient geometry
# f | predicted-response composition | BarComposition/plot_rewriting.py | param inherit | zero-based two-part stacked bars retain exact model-predicted counts
# g | paired pathway evidence | GroupedBarChart/plot_GroupedBarChartv1.py | cross-type inherit | paired adjustment tiers are connected lollipops
# h | cohort-specific clinicopathology forest | forest plot production conventions | param inherit | four cohort estimates and one meta-estimate are grouped by harmonized feature

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
import hashlib
import json
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.backends.backend_pdf import PdfPages

from gen_fig5_immune_extension_review import effect_matrix


ROOT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
DATA = ROOT / "figures" / "source_data" / "FigureS6_final_v1"
OUTDIR = ROOT / "figures" / "ejc_revision_20260827"
BASE_S6 = OUTDIR / "Supplementary_Figure_S6_bulk_immune_sensitivity_ejc_20260827"
BASE_S7 = OUTDIR / "Supplementary_Figure_S7_tide_clinical_sensitivity_ejc_20260827"
COHORTS = ["TCGA-STAD", "GSE62254", "GSE84437", "GSE15459"]
COHORT_LABELS = ["TCGA", "GSE62254", "GSE84437", "GSE15459"]
INK = "#26323E"
LIGHT_GREY = "#E7EBEE"
MID_GREY = "#9AA3AA"
SENDER = "#8F5F6B"
RECEIVER = "#456F9D"
ECM = "#6F8B70"
HEAT_CMAP = LinearSegmentedColormap.from_list("signed", DIVERGING, N=256)


def read(name, index_col=None):
    return pd.read_csv(DATA / name, sep="\t", index_col=index_col)


def label(ax, text, x=-0.14, y=1.10):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=9, fontweight="bold",
            ha="left", va="top", color=BLACK, clip_on=False)


def title(ax, text, y=None):
    kwargs = dict(loc="left", fontweight="bold", color=BLACK, pad=4)
    if y is not None:
        kwargs.update(y=y, pad=0)
    ax.set_title(text, **kwargs)


def clean(ax, hide_left=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if hide_left:
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    ax.tick_params(length=2.2, width=0.6, colors=INK)


def diamond(ax, center, low, high, y, height, color):
    ax.add_patch(Polygon([[low, y], [center, y + height], [high, y],
                          [center, y - height]], closed=True,
                         facecolor=color, edgecolor="white", linewidth=0.45, zorder=4))


def draw_heatmap(fig, spec):
    heat = read("Figure5b_TCGA_60gene_rowZ_heatmap.tsv.gz", index_col=0)
    anno = read("Figure5b_TCGA_sample_annotations.tsv")
    gene = read("Figure5b_TCGA_gene_annotations.tsv")
    frame = fig.add_subplot(spec, frameon=False)
    frame.set_xticks([]); frame.set_yticks([])
    for spine in frame.spines.values(): spine.set_visible(False)
    label(frame, "a", x=-0.10, y=1.05)
    grid = spec.subgridspec(6, 1, height_ratios=[0.055, 0.055, 0.055, 0.055, 0.11, 1.0], hspace=0.025)
    groups = anno["ecological_configuration"].fillna("Dual-low")
    group_colors = {"Dual-low": "#D9D9D9", "Fibroblast-dominant": "#CDACB4",
                    "Macrophage-dominant": "#8FABDA", "Dual-high": "#8F5F6B"}
    subtype_colors = {"CIN": "#91A8B8", "GS": "#B77C88", "MSI": "#AFC5A9",
                      "EBV": "#6F91B7", "Unknown": "#E6E8EA"}
    stage_colors = {"I": "#E8F0E5", "II": "#BFD2B8", "III": "#849F82", "IV": "#4F6E55",
                    "Unknown": "#E6E8EA"}
    sex_colors = {"FEMALE": "#D8A8B4", "MALE": "#7896B6", "Unknown": "#E6E8EA"}
    arrays = []
    arrays.append(np.array([mpl.colors.to_rgb(group_colors.get(x, "#D9D9D9")) for x in groups])[None, :, :])
    arrays.append(np.array([mpl.colors.to_rgb(subtype_colors.get(x, "#E6E8EA")) for x in anno["tcga2014_subtype"].fillna("Unknown")])[None, :, :])
    arrays.append(np.array([mpl.colors.to_rgb(stage_colors.get(x, "#E6E8EA")) for x in anno["stage_group_frozen"].fillna("Unknown")])[None, :, :])
    arrays.append(np.array([mpl.colors.to_rgb(sex_colors.get(x, "#E6E8EA")) for x in anno["sex"].fillna("Unknown")])[None, :, :])
    score = np.vstack([anno["sender_primary_z"], anno["receiver_specific_z"]])
    axes = []
    for i in range(4):
        ax = fig.add_subplot(grid[i, 0]); ax.imshow(arrays[i], aspect="auto"); ax.axis("off"); axes.append(ax)
    ax_s = fig.add_subplot(grid[4, 0])
    ax_s.imshow(score, aspect="auto", cmap=HEAT_CMAP, norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5))
    ax_s.axis("off"); axes.append(ax_s)
    ax_h = fig.add_subplot(grid[5, 0])
    ax_h.imshow(heat.to_numpy(), aspect="auto", interpolation="nearest", cmap=HEAT_CMAP,
                norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5), rasterized=True)
    selected = ["VSTM2A", "FST", "APOD", "FABP4", "PTGDS", "SMOC1",
                "LGMN", "CTSZ", "MS4A4A", "DAB2", "MSR1", "SERPING1"]
    yt = [i for i, g in enumerate(heat.index) if g in selected]
    ax_h.set_yticks(yt); ax_h.set_yticklabels([heat.index[i] for i in yt], fontsize=5.0)
    ax_h.set_xticks([])
    for spine in ax_h.spines.values(): spine.set_visible(False)
    boundaries = np.where(groups.to_numpy()[1:] != groups.to_numpy()[:-1])[0] + 0.5
    for xx in boundaries:
        for ax in axes + [ax_h]: ax.axvline(xx, color="white", linewidth=0.6)
    boundary = int((gene["state_block"] == "Fibroblast sender").sum()) - 0.5
    ax_h.axhline(boundary, color="white", linewidth=0.9)
    names = ["Config.", "TCGA", "Stage", "Sex", "S / R"]
    for ax, name in zip(axes, names):
        ax.text(-0.006, 0.5, name, transform=ax.transAxes, ha="right", va="center", fontsize=5.0)


def draw_gate(fig, spec):
    metrics = read("Figure5b_consensus_gate_metrics.tsv")
    frame = fig.add_subplot(spec, frameon=False)
    frame.set_xticks([]); frame.set_yticks([])
    for s in frame.spines.values(): s.set_visible(False)
    grid = spec.subgridspec(2, 2, hspace=0.62, wspace=0.38)
    fields = [
        ("PAC_0.1_0.9", "PAC", 0.25, "below"),
        ("mean_silhouette", "Mean silhouette", 0.25, "above"),
        ("minimum_cluster_fraction", "Minimum cluster fraction", 0.10, "above"),
        ("primary_vs_sensitivity_ARI", "ARI", 0.50, "above"),
    ]
    for i, (field, ttl, threshold, direction) in enumerate(fields):
        ax = fig.add_subplot(grid[i // 2, i % 2])
        ax.plot(metrics["k"], metrics[field], color="#666666", linewidth=0.8,
                marker="o", markersize=3, markerfacecolor="white")
        best = metrics[metrics["k"] == 3].iloc[0]
        ax.scatter([3], [best[field]], s=24, color=ACCENT_RED, zorder=3)
        ax.axhline(threshold, color=CATEGORICAL[0], linewidth=0.7, linestyle="--")
        ax.set_title(ttl, loc="left", fontsize=5.8, fontweight="bold", pad=2)
        ax.set_xticks(metrics["k"]); ax.set_xlabel("k", fontsize=5.2)
        ax.tick_params(labelsize=5.0)
        clean(ax)
    label(frame, "b", x=-0.12, y=1.12)


def draw_checkpoint(ax):
    cohort = read("four_cohort_checkpoint_continuous_HC3_models.tsv")
    meta = read("four_cohort_checkpoint_random_effects_meta.tsv")
    adj = read("TCGA_checkpoint_adjusted_HC3_models.tsv")
    adj = adj[adj["model_tier"] == "lineage_subtype_adjusted"]
    order = ["CD86", "SIGLEC10", "HAVCR2", "ICOS", "TIGIT", "SIRPA", "PDCD1LG2",
             "TNFRSF4", "LAG3", "CD80", "TNFRSF9", "CTLA4", "CD274", "PDCD1",
             "CD276", "CD24", "CD47", "VSIR"]
    labels = ["CD86", "SIGLEC10", "TIM-3", "ICOS", "TIGIT", "SIRPalpha", "PD-L2",
              "OX40", "LAG-3", "CD80", "4-1BB", "CTLA-4", "PD-L1", "PD-1", "B7-H3",
              "CD24", "CD47", "VISTA"]
    cp = cohort.pivot(index="gene_symbol", columns="cohort", values="beta").reindex(order)
    mm = meta.set_index("gene_symbol").reindex(order)
    aa = adj.set_index("gene_symbol").reindex(order)
    matrix = np.column_stack([cp[c].to_numpy(float) for c in COHORTS] +
                             [mm["beta_random_DL"].to_numpy(float), aa["beta"].to_numpy(float)])
    effect_matrix(ax, matrix, labels, COHORT_LABELS + ["RE meta", "TCGA full adj."],
                  mm["FDR_BH"].fillna(1).to_numpy(float), mm["I2_percent"].to_numpy(float),
                  adjusted_fdr=aa["FDR_BH"].to_numpy(float))
    label(ax, "c", x=-0.09, y=1.13)


def draw_tide_config(ax):
    d = read("Figure5d_TIDE_configuration_medians.tsv").set_index("ecotype4")
    order = ["Dual-low", "Fibroblast-dominant", "Macrophage-dominant", "Dual-high"]
    endpoints = ["TIDE", "Dysfunction", "Exclusion", "CAF", "MDSC", "TAM M2", "CD8", "CTL"]
    mat = d.reindex(order)[[x + "_z" for x in endpoints]].to_numpy().T
    im = ax.imshow(mat, aspect="auto", cmap=HEAT_CMAP,
                   norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1))
    ax.set_xticks(range(4)); ax.set_xticklabels(["Dual-low", "Fibroblast-dom.", "Macrophage-dom.", "Dual-high"],
                                                rotation=25, ha="right", fontsize=5.5)
    ax.set_yticks(range(8)); ax.set_yticklabels(endpoints, fontsize=5.7)
    for i in range(8):
        for j in range(4):
            val = mat[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.0,
                    color="white" if abs(val) > 0.58 else BLACK)
    for s in ax.spines.values(): s.set_visible(False)
    label(ax, "d", x=-0.25, y=1.13)
    return im


def page1():
    fig = plt.figure(figsize=(183 / 25.4, 235 / 25.4), facecolor="white")
    outer = fig.add_gridspec(3, 12, height_ratios=[1.28, 1.18, 0.92],
                             left=0.16, right=0.975, top=0.955, bottom=0.065,
                             hspace=0.62, wspace=0.76)
    draw_heatmap(fig, outer[0, :])
    draw_gate(fig, outer[1, :5])
    ax_c = fig.add_subplot(outer[1, 6:]); draw_checkpoint(ax_c)
    ax_d = fig.add_subplot(outer[2, 2:10]); im = draw_tide_config(ax_d)
    cax = fig.add_axes([0.89, 0.095, 0.012, 0.15])
    cb = fig.colorbar(im, cax=cax); cb.set_label("Median endpoint Z score", fontsize=5.2)
    cb.ax.tick_params(labelsize=5.0)
    return fig


def draw_tide_forest(ax):
    tide = read("TCGA_TIDE_numeric_endpoint_models_existing.tsv")
    order = ["Dysfunction", "Exclusion", "TIDE", "MDSC", "CAF", "TAM M2",
             "IFNG", "CD274", "CD8", "CTL", "MSI Score"]
    y = np.arange(len(order))[::-1]
    ax.axvline(0, color="#777777", linewidth=0.6, linestyle="--")
    for tier, off, color, marker in [("purity_adjusted", 0.13, GREY, "o"),
                                      ("lineage_subtype_adjusted", -0.13, CATEGORICAL[0], "D")]:
        d = tide[tide["model_tier"] == tier].set_index("endpoint").reindex(order)
        for yy, (_, row) in zip(y + off, d.iterrows()):
            sig = row["FDR_BH"] < 0.05
            ax.plot([row["ci_low"], row["ci_high"]], [yy, yy], color=color, linewidth=0.8)
            ax.scatter(row["beta"], yy, s=18, marker=marker, facecolor=color if sig else "white",
                       edgecolor=color, linewidth=0.8, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(order, fontsize=5.7)
    ax.set_xlim(-0.55, 0.65); ax.set_xlabel("Standardized association, beta (95% CI)", fontsize=6.2)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="white", markeredgecolor=GREY,
               label="Purity-adjusted", markersize=4),
        Line2D([0], [0], marker="D", linestyle="none", markerfacecolor="white", markeredgecolor=CATEGORICAL[0],
               label="Lineage + subtype adjusted", markersize=4),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=BLACK, markeredgecolor=BLACK,
               label="FDR < 0.05", markersize=4)],
        loc="upper left", bbox_to_anchor=(0.0, -0.18), ncol=3, fontsize=5.0,
        labelspacing=0.2, handletextpad=0.3, columnspacing=0.8, borderaxespad=0)
    clean(ax)


def draw_predicted(ax):
    d = read("TCGA_TIDE_predicted_response_by_configuration.tsv")
    order = ["Dual-low", "Macrophage-dominant", "Dual-high", "Fibroblast-dominant"]
    d = d.set_index("ecotype4").reindex(order)
    y = np.arange(4)[::-1]; pct = d["predicted_response_fraction"].to_numpy(float) * 100
    ax.barh(y, pct, color=CATEGORICAL[0], edgecolor="white", linewidth=0.5, height=0.65)
    ax.barh(y, 100 - pct, left=pct, color="#D9D9D9", edgecolor="white", linewidth=0.5, height=0.65)
    for yy, pp, (_, row) in zip(y, pct, d.iterrows()):
        ax.text(pp / 2, yy, f"{int(row['predicted_responders'])}/{int(row['n'])}", ha="center", va="center",
                fontsize=5.3, color="white", fontweight="bold")
        ax.text(101, yy, f"{pp:.1f}%", fontsize=5.2, va="center")
    ax.set_yticks(y); ax.set_yticklabels(order, fontsize=5.5)
    ax.set_xlim(0, 112); ax.set_xlabel("TIDE-predicted response (%)", fontsize=6.2)
    ax.text(0, 1.02, "4 x 2 Fisher P=0.0118", transform=ax.transAxes, fontsize=5.3)
    clean(ax)


def draw_camera(ax):
    camera = read("Figure5c_CAMERA_focal_pathways.tsv")
    order = ["Complement & coagulation", "Coagulation", "Complement activation",
             "Myeloid migration", "Positive chemotaxis", "Macrophage migration",
             "Phagocytosis regulation", "ECM assembly", "ECM organization",
             "Inflammatory response", "Hypoxia", "EMT"]
    p = camera[camera["model_tier"] == "purity_adjusted"].set_index("label").reindex(order)
    f = camera[camera["model_tier"] == "lineage_subtype_adjusted"].set_index("label").reindex(order)
    x1 = -np.log10(p["FDR_BH"].to_numpy(float)); x2 = -np.log10(f["FDR_BH"].to_numpy(float))
    y = np.arange(len(order))[::-1]
    for yy, a, b in zip(y, x1, x2): ax.plot([a, b], [yy, yy], color="#BDBDBD", linewidth=0.7)
    ax.scatter(x1, y, s=18, facecolor="white", edgecolor=GREY, linewidth=0.8, label="Purity-adjusted")
    sig = f["FDR_BH"].to_numpy(float) < 0.05
    ax.scatter(x2, y, s=20, marker="D", facecolor=np.where(sig, ACCENT_RED, "white"),
               edgecolor=ACCENT_RED, linewidth=0.8, label="Lineage + subtype adjusted")
    ax.axvline(-math.log10(0.05), color="#777777", linewidth=0.6, linestyle="--")
    ax.set_yticks(y); ax.set_yticklabels(order, fontsize=5.5)
    ax.set_xlim(0, 7.35); ax.set_xlabel("-log10(FDR), CAMERA competitive test", fontsize=6.2)
    ax.legend(loc="lower right", bbox_to_anchor=(1, 1.0), fontsize=5.0, ncol=2, borderaxespad=0)
    clean(ax)


def draw_clinic(ax):
    models = read("four_cohort_clinicopathology_HC3_models.tsv")
    models = models[models["model_tier"] == "univariable"]
    meta = read("four_cohort_clinicopathology_random_effects_meta.tsv")
    meta = meta[(meta["model_tier"] == "univariable") &
                (meta["analysis_set"] == "all_available_including_TCGA")]
    features = ["lauren_diffuse", "stage", "t_stage", "n_stage", "male", "age10"]
    labels = ["Diffuse vs intestinal", "Pathological stage / level", "T stage / level",
              "N stage / level", "Male vs female", "Age / 10 years"]
    cohort_colors = dict(zip(COHORTS, [CATEGORICAL[1], CATEGORICAL[2], CATEGORICAL[0], CATEGORICAL[3]]))
    bases = np.arange(len(features))[::-1] * 5.5
    offsets = {"TCGA-STAD": 1.5, "GSE62254": 0.5, "GSE84437": -0.5, "GSE15459": -1.5}
    ax.axvline(0, color="#777777", linewidth=0.6, linestyle="--")
    for base, feature in zip(bases, features):
        sub = models[models["feature"] == feature]
        for cohort in COHORTS:
            row = sub[sub["cohort"] == cohort]
            if row.empty: continue
            row = row.iloc[0]; yy = base + offsets[cohort]
            ax.plot([row["ci_low"], row["ci_high"]], [yy, yy], color=cohort_colors[cohort], linewidth=0.75)
            ax.scatter(row["beta"], yy, s=15, facecolor="white", edgecolor=cohort_colors[cohort], linewidth=0.75)
        mm = meta[meta["feature"] == feature].iloc[0]
        diamond(ax, mm["beta_random_DL"], mm["ci_low"], mm["ci_high"], base - 2.35, 0.25, BLACK)
        ax.axhline(base - 2.8, color="#EFEFEF", linewidth=0.5)
    ax.set_yticks(bases - 0.4); ax.set_yticklabels(labels, fontsize=5.7)
    ax.set_xlim(-0.38, 1.15); ax.set_xlabel("Difference in ecosystem score (95% CI)", fontsize=6.2)
    ax.legend(handles=[Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="white",
                              markeredgecolor=cohort_colors[c], label=c.replace("TCGA-STAD", "TCGA"), markersize=4)
                       for c in COHORTS] +
                      [Line2D([0], [0], marker="D", linestyle="none", markerfacecolor=BLACK,
                              markeredgecolor=BLACK, label="RE meta", markersize=4)],
              loc="lower right", bbox_to_anchor=(1, 1.0), ncol=5, fontsize=5.0,
              handletextpad=0.3, columnspacing=0.7, borderaxespad=0)
    clean(ax, hide_left=True)


def page2():
    fig = plt.figure(figsize=(183 / 25.4, 235 / 25.4), facecolor="white")
    outer = fig.add_gridspec(3, 12, height_ratios=[1.05, 1.0, 1.18],
                             left=0.17, right=0.975, top=0.955, bottom=0.065,
                             hspace=0.62, wspace=0.82)
    ax_e = fig.add_subplot(outer[0, :6]); draw_tide_forest(ax_e)
    label(ax_e, "a", x=-0.34, y=1.12)
    ax_f = fig.add_subplot(outer[0, 7:]); draw_predicted(ax_f)
    label(ax_f, "b", x=-0.35, y=1.12)
    ax_g = fig.add_subplot(outer[1, 1:11]); draw_camera(ax_g)
    label(ax_g, "c", x=-0.16, y=1.12)
    ax_h = fig.add_subplot(outer[2, :]); draw_clinic(ax_h)
    label(ax_h, "d", x=-0.13, y=1.12)
    return fig


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_exact(fig, stem):
    save_cns_figure(fig, str(stem))
    with mpl.rc_context({"savefig.bbox": None}):
        fig.savefig(f"{stem}.pdf", bbox_inches=None, pad_inches=0, dpi=300)
        fig.savefig(f"{stem}.png", bbox_inches=None, pad_inches=0, dpi=300)
        fig.savefig(f"{stem}.svg", bbox_inches=None, pad_inches=0, dpi=300)
        fig.savefig(f"{stem}.tiff", bbox_inches=None, pad_inches=0, dpi=600,
                    pil_kwargs={"compression": "tiff_lzw"})


if __name__ == "__main__":
    OUTDIR.mkdir(parents=True, exist_ok=True)
    s6 = page1(); s7 = page2()
    save_exact(s6, BASE_S6); save_exact(s7, BASE_S7)
    plt.close(s6); plt.close(s7)
    outputs = [Path(f"{stem}.{ext}") for stem in [BASE_S6, BASE_S7]
               for ext in ["pdf", "png", "svg", "tiff"]]
    qa = {"figures": ["Figure S6", "Figure S7"], "pages_each": 1,
          "canvas_mm_each": [183, 235],
          "outputs": {p.name: {"bytes": p.stat().st_size, "sha256": sha256(p)} for p in outputs},
          "claim_boundary": "Consensus gate failed; TIDE response is predicted; clinicopathology k varies by feature."}
    with open(OUTDIR / "FiguresS6_S7_split_ejc_20260827_QA.json", "w", encoding="utf-8") as handle:
        json.dump(qa, handle, indent=2)
    print("\n".join(str(p) for p in outputs))
