# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) single-cell UMAP -> assets/figures/PCA -> param inherit (equal-aspect adaptation)
# (b) marker-gene dot plot -> assets/figures/MarkerGeneDotPlot -> visual adapt
# (c) sender bars + directed circle -> assets/figures/GroupedBarChart + user reference -> param inherit
# (d) LR heat map with marginal summaries -> assets/figures/heatmap -> visual adapt
# (e) donor-by-candidate heat map -> assets/figures/heatmap -> visual adapt
# (f) ecological forest plot -> assets/figures/BarComparison -> cross-type param inherit
# (S4) cross-cohort heat map -> assets/figures/heatmap -> visual adapt
# RULE: no production asset is semantically compatible with the project data;
#       all panels preserve the frozen project statistics and adapt visual parameters only.

from __future__ import annotations

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
    "text.color": "#000000",
    "axes.labelcolor": "#000000",
    "axes.edgecolor": "#000000",
    "xtick.color": "#000000",
    "ytick.color": "#000000",
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
BLACK       = "#000000"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": None,
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches=None, dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches=None, dpi=300)


from pathlib import Path
import json

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
import numpy as np
import pandas as pd


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
FIG_DIR = PROJECT / "figures" / "main"
SOURCE_DIR = PROJECT / "figures" / "data" / "Figure3_PLAT_all5_raw_communication"
OUT_DIR = PROJECT / "figures" / "data" / "Figure3_PLAT_reference_layout_v4"
CCC_DIR = PROJECT / "analysis" / "results" / "communication" / "Figure3_all5_raw"
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

FILES = {
    "cells": SOURCE_DIR / "Figure3a_myeloid_umap.tsv.gz",
    "markers": SOURCE_DIR / "Figure3b_patient_balanced_marker_profile.tsv",
    "circle": SOURCE_DIR / "Figure3c_CellChat_fibroblast_myeloid_circle.tsv",
    "receiver": SOURCE_DIR / "Figure3c_CellChat_F1_all5_strength.tsv",
    "liana": SOURCE_DIR / "Figure3d_LIANA_all5_LR_matrix.tsv",
    "order": SOURCE_DIR / "Figure3d_LIANA_LR_display_order.tsv",
    "programme": SOURCE_DIR / "Figure3e_all5_donor_programme.tsv",
    "mapping": SOURCE_DIR / "Figure3f_cross_cohort_mapping.tsv",
    "ecology": SOURCE_DIR / "Figure3g_patient_ecology.tsv",
    "donor_pair": CCC_DIR / "panel_e_all5_donor_pair_potentials.tsv",
}
for path in FILES.values():
    if not path.exists():
        raise FileNotFoundError(path)

cells = pd.read_csv(FILES["cells"], sep="\t")
markers = pd.read_csv(FILES["markers"], sep="\t")
circle = pd.read_csv(FILES["circle"], sep="\t")
receiver = pd.read_csv(FILES["receiver"], sep="\t")
liana = pd.read_csv(FILES["liana"], sep="\t")
pair_order_df = pd.read_csv(FILES["order"], sep="\t").sort_values("display_order")
programme = pd.read_csv(FILES["programme"], sep="\t")
mapping = pd.read_csv(FILES["mapping"], sep="\t")
ecology = pd.read_csv(FILES["ecology"], sep="\t")
donor_pair = pd.read_csv(FILES["donor_pair"], sep="\t")


STATE_ORDER = [1, 2, 3, 4, 5]
STATE_COLORS = {1: "#D55E00", 2: "#CC79A7", 3: "#0072B2", 4: "#E69F00", 5: "#56B4E9"}
STATE_LABELS = {
    1: "C1  monocyte",
    2: "C2  SPP1/APOE",
    3: "C3  C1QC/APOE",
    4: "C4  infl. myeloid",
    5: "C5  HLA-II/C1QC",
}
MARKER_GROUPS = [
    ["FCN1", "AREG", "S100A8"],
    ["SPP1", "APOC1", "MARCO"],
    ["C1QA", "C1QC", "APOE", "LGMN"],
    ["IL1B", "CXCL8", "NFKBIA"],
    ["HLA-DRA", "HLA-DPB1", "CD74", "FCER1A"],
]
MARKER_GENES = [gene for group in MARKER_GROUPS for gene in group]
MARKER_CLUSTER = {
    gene: cluster
    for cluster, genes in enumerate(MARKER_GROUPS, start=1)
    for gene in genes
}
MARKER_Y = {}
MARKER_BOUNDS = {}
_marker_cursor = 0.0
for _cluster, _genes in reversed(list(enumerate(MARKER_GROUPS, start=1))):
    _positions = []
    for _gene in reversed(_genes):
        MARKER_Y[_gene] = _marker_cursor
        _positions.append(_marker_cursor)
        _marker_cursor += 1.0
    MARKER_BOUNDS[_cluster] = (min(_positions) - 0.43, max(_positions) + 0.43)
    if _cluster != 1:
        _marker_cursor += 0.65
MARKER_YMAX = max(MARKER_Y.values())
MARKER_BOTTOM_TO_TOP = sorted(MARKER_GENES, key=MARKER_Y.get)

SENDER_ORDER = ["F1_PLATrich", "F2_MuralControl", "F3_AntigenPresentation"]
SENDER_LABELS = {
    "F1_PLATrich": "F1\nPLAT-rich",
    "F2_MuralControl": "F2\nMural-like",
    "F3_AntigenPresentation": "F3\nHLA-II",
}
SENDER_COLORS = {"F1_PLATrich": ACCENT_RED, "F2_MuralControl": "#6F8EA5", "F3_AntigenPresentation": "#9AA8A0"}

QUERY_ORDER = ["M02_Mono_CD14", "M04_Macro_C3", "M06_Macro_SLC40A1", "M11_Macro_PLA2G7"]
QUERY_LABELS = {
    "M02_Mono_CD14": "M02  CD14+ mono.",
    "M04_Macro_C3": "M04  C3+ mac.",
    "M06_Macro_SLC40A1": "M06  SLC40A1+ mac.",
    "M11_Macro_PLA2G7": "M11  PLA2G7+ mac.",
}
REFERENCE_ORDER = [
    "C1_FCN1_AREG_inflammatory_monocyte",
    "C3_C1QC_APOE_macrophage",
    "C4_IL1B_CXCL8_inflammatory_myeloid",
    "C5_HLAII_C1QC_macrophage",
]


def panel_label(ax, label, x=-0.12, y=1.06):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="left", va="top")


def clean_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=2.5, width=0.6, colors="#000000")
    ax.xaxis.label.set_color("#000000")
    ax.yaxis.label.set_color("#000000")
    ax.spines["left"].set_color("#000000")
    ax.spines["bottom"].set_color("#000000")


def export_tsv(frame, name):
    path = OUT_DIR / name
    frame.to_csv(path, sep="\t", index=False)
    return path


# Data validation and derived matrices.
if len(cells) != 9002 or cells.Patient_ID.nunique() != 26:
    raise ValueError("Expected 9,002 myeloid cells from 26 patients")
if len(circle) != 15 or set(circle.source) != set(SENDER_ORDER) or set(circle.cluster) != set(STATE_ORDER):
    raise ValueError("The complete 3 x 5 directed communication matrix is required")

sender_summary = (
    circle.groupby("source", observed=True)
    .agg(aggregate_strength=("aggregate_probability", "sum"), supported_lr_pairs=("supported_lr_pairs", "sum"))
    .reindex(SENDER_ORDER)
    .reset_index()
)
receiver_summary = (
    circle.groupby("cluster", observed=True)
    .agg(incoming_strength=("aggregate_probability", "sum"), supported_lr_pairs=("supported_lr_pairs", "sum"))
    .reindex(STATE_ORDER)
    .reset_index()
)
pair_order = pair_order_df.pair.astype(str).tolist()
mag_matrix = (
    liana.assign(mag=lambda x: -np.log10(x.magnitude_rank.clip(lower=1e-8)))
    .pivot(index="pair", columns="cluster", values="mag")
    .reindex(index=pair_order, columns=STATE_ORDER)
)
strict_matrix = (
    liana.assign(strict=lambda x: x.magnitude_rank.le(0.05) & x.specificity_rank.le(0.05))
    .pivot(index="pair", columns="cluster", values="strict")
    .reindex(index=pair_order, columns=STATE_ORDER)
    .astype("boolean").fillna(False).astype(bool)
)
p_matrix = (
    liana.pivot(index="pair", columns="cluster", values="cellphone_pvals")
    .reindex(index=pair_order, columns=STATE_ORDER)
)
star_matrix = strict_matrix & p_matrix.lt(0.05).fillna(False)
strict_counts = star_matrix.sum(axis=1).astype(int)
receiver_plot = receiver.set_index("cluster").reindex(STATE_ORDER)

focus_pairs = []
for prefix in ["APP", "NRG1", "CXCL14", "SPON2"]:
    matches = [pair for pair in pair_order if pair.startswith(prefix)]
    if len(matches) != 1:
        raise ValueError(f"Expected one display pair beginning with {prefix}: {matches}")
    focus_pairs.append(matches[0])
focus_donor = donor_pair[
    donor_pair.target.eq("M3_C1QC_APOE_target") & donor_pair.pair.isin(focus_pairs)
].copy()
if focus_donor.Patient_ID.nunique() != 22 or len(focus_donor) != 88:
    raise ValueError("Expected four candidate estimates in each of 22 F1-to-C3 patients")
focus_donor["display_value"] = focus_donor.lr_expression_potential.where(focus_donor.expression_gate_support)
focus_donor["within_pair_z"] = focus_donor.groupby("pair", observed=True).display_value.transform(
    lambda x: (x - x.mean()) / x.std(ddof=0) if x.std(ddof=0) > 0 else 0.0
)
focus_support = (
    focus_donor.groupby("pair", observed=True)
    .agg(
        evaluable_patients=("Patient_ID", "nunique"),
        supported_patients=("expression_gate_support", "sum"),
        median_expression_potential=("lr_expression_potential", "median"),
    )
    .reindex(focus_pairs).reset_index()
)
donor_order = (
    focus_donor.pivot(index="Patient_ID", columns="pair", values="within_pair_z")
    .reindex(columns=focus_pairs)
    .assign(_n=lambda x: x.notna().sum(axis=1), _mean=lambda x: x.mean(axis=1, skipna=True))
    .sort_values(["_n", "_mean"], ascending=[False, False]).index.tolist()
)
donor_codes = {donor: f"P{i + 1:02d}" for i, donor in enumerate(donor_order)}
focus_donor["patient_code"] = focus_donor.Patient_ID.map(donor_codes)

map_matrix = (
    mapping[mapping.query_state.isin(QUERY_ORDER) & mapping.reference_state.isin(REFERENCE_ORDER)]
    .pivot(index="query_state", columns="reference_state", values="spearman_rho")
    .reindex(index=QUERY_ORDER, columns=REFERENCE_ORDER)
)
eco_plot = ecology.set_index("cluster").reindex([1, 3, 4, 5])

export_tsv(sender_summary, "Figure3c_sender_communication_summary.tsv")
export_tsv(receiver_summary, "Figure3c_receiver_communication_summary.tsv")
export_tsv(circle, "Figure3c_fibroblast_myeloid_circle.tsv")
export_tsv(liana, "Figure3d_LIANA_LR_matrix.tsv")
export_tsv(pair_order_df, "Figure3d_LIANA_LR_display_order.tsv")
export_tsv(focus_donor, "Figure3e_C3_candidate_patient_recurrence.tsv")
export_tsv(focus_support, "Figure3e_C3_candidate_support_summary.tsv")
export_tsv(programme, "FigureS3_patient_state_programme.tsv")


# Main Figure 3: asymmetric mixed-modality layout, 183 mm x 225 mm.
fig = plt.figure(figsize=(183 / 25.4, 225 / 25.4), facecolor="white")
outer = fig.add_gridspec(
    3, 1, height_ratios=[1.28, 0.92, 1.58],
    left=0.110, right=0.985, top=0.982, bottom=0.052, hspace=0.18,
)
top = outer[0].subgridspec(1, 2, width_ratios=[0.94, 1.06], wspace=0.18)
agrid = top[0, 0].subgridspec(1, 2, width_ratios=[0.78, 0.22], wspace=0.01)
cgrid = outer[1].subgridspec(1, 2, width_ratios=[1.55, 1.25], wspace=0.12)
cbargrid = cgrid[0, 0].subgridspec(1, 2, width_ratios=[1.04, 0.96], wspace=0.06)
lower = outer[2].subgridspec(1, 2, width_ratios=[1.48, 0.96], wspace=0.12)
rightstack = lower[0, 1].subgridspec(2, 1, height_ratios=[0.75, 0.25], hspace=0.38)


# a: all five myeloid states; UMAP coordinates retain equal scaling.
ax_a = fig.add_subplot(agrid[0, 0])
panel_label(ax_a, "a", -0.08, 1.04)
for cluster in [2, 5, 4, 1, 3]:
    block = cells[cells.cluster.eq(cluster)]
    ax_a.scatter(block.UMAP1, block.UMAP2, s=1.7, color=STATE_COLORS[cluster],
                 alpha=0.48, linewidths=0, rasterized=True)
    centroid = block[["UMAP1", "UMAP2"]].median().to_numpy(float)
    label = ax_a.text(centroid[0], centroid[1], f"C{cluster}",
                      ha="center", va="center", fontsize=6.2,
                      fontweight="bold", color="#000000", zorder=6)
    label.set_path_effects([pe.withStroke(linewidth=1.8, foreground="white")])
ax_a.set_aspect("equal", adjustable="box")
ax_a.set_anchor("W")
ax_a.margins(x=0.018, y=0.018)
ax_a.set_xticks([]); ax_a.set_yticks([])
for spine in ax_a.spines.values():
    spine.set_visible(False)
handles = [Line2D([0], [0], marker="o", linestyle="", markersize=3.6,
                  markerfacecolor=STATE_COLORS[c], markeredgecolor="none",
                  label=STATE_LABELS[c]) for c in STATE_ORDER]
ax_aleg = fig.add_subplot(agrid[0, 1]); ax_aleg.axis("off")
ax_aleg.legend(handles=handles, ncol=1, loc="center left", bbox_to_anchor=(-0.12, 0.50),
               handletextpad=0.34, labelspacing=0.48, fontsize=5.1)


# b: patient-balanced marker profile.
ax_b = fig.add_subplot(top[0, 1])
panel_label(ax_b, "b", -0.08, 1.04)
dot_cmap = LinearSegmentedColormap.from_list("marker_z", DIVERGING)
for row in markers.itertuples(index=False):
    x = STATE_ORDER.index(int(row.cluster)); y = MARKER_Y[row.gene_symbol]
    detect = float(row.patient_median_detection)
    zvalue = float(np.clip(row.expression_z_across_states, -1.8, 1.8))
    ax_b.scatter(x, y, s=7 + 68 * detect**1.3,
                 facecolor=dot_cmap((zvalue + 1.8) / 3.6),
                 edgecolor="white", linewidth=0.30, zorder=3)
for cluster in STATE_ORDER:
    low, high = MARKER_BOUNDS[cluster]
    ax_b.axhspan(low, high, color=STATE_COLORS[cluster], alpha=0.035, zorder=0)
    ax_b.text(-0.68, (low + high) / 2, f"C{cluster}", ha="center", va="center",
              fontsize=5.2, fontweight="bold", color=STATE_COLORS[cluster], zorder=4)
ax_b.set_xlim(-0.84, 4.62); ax_b.set_ylim(-0.62, MARKER_YMAX + 0.62)
ax_b.set_xticks(range(5), [f"C{x}" for x in STATE_ORDER])
for tick, cluster in zip(ax_b.get_xticklabels(), STATE_ORDER):
    tick.set_color(STATE_COLORS[cluster]); tick.set_fontweight("bold")
ax_b.set_yticks(
    [MARKER_Y[gene] for gene in MARKER_BOTTOM_TO_TOP],
    MARKER_BOTTOM_TO_TOP,
    fontstyle="italic",
)
for tick, gene in zip(ax_b.get_yticklabels(), MARKER_BOTTOM_TO_TOP):
    tick.set_color(STATE_COLORS[MARKER_CLUSTER[gene]])
ax_b.tick_params(length=0, pad=1.2)
ax_b.grid(axis="x", color="#EDF0F2", linewidth=0.40, zorder=0)
for spine in ax_b.spines.values():
    spine.set_visible(False)
for lower_cluster in [5, 4, 3, 2]:
    boundary = (MARKER_BOUNDS[lower_cluster][1] + MARKER_BOUNDS[lower_cluster - 1][0]) / 2
    ax_b.axhline(boundary, color="#CED5DA", linewidth=0.50)
for index, detect in enumerate([0.25, 0.50, 0.75, 1.00]):
    xpos = 0.08 + index * 0.115
    ax_b.scatter(xpos, -0.10, s=7 + 68 * detect**1.3, transform=ax_b.transAxes,
                 facecolor="#AAB3BA", edgecolor="white", linewidth=0.3, clip_on=False)
    ax_b.text(xpos, -0.17, f"{int(100 * detect)}", transform=ax_b.transAxes,
              ha="center", va="top", fontsize=5.0, color="#000000", clip_on=False)
cax_b = ax_b.inset_axes([0.59, -0.135, 0.34, 0.030])
cb_b = fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(-1.8, 1.8), cmap=dot_cmap),
                    cax=cax_b, orientation="horizontal")
cb_b.set_ticks([-1.5, 0, 1.5]); cb_b.ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f")); cb_b.ax.tick_params(labelsize=5, length=1, pad=0.4)
cb_b.outline.set_linewidth(0.35)


# c: complete outgoing sender and incoming receiver summaries with the 3 x 5 network.
ax_cnum = fig.add_subplot(cbargrid[0, 0])
panel_label(ax_cnum, "c", -0.30, 1.05)
bar_positions = np.array([0, 1, 2, 4, 5, 6, 7, 8], dtype=float)
bar_labels = ["F1 →", "F2 →", "F3 →", "→ C1", "→ C2", "→ C3", "→ C4", "→ C5"]
bar_colors = [SENDER_COLORS[s] for s in SENDER_ORDER] + [STATE_COLORS[c] for c in STATE_ORDER]
count_values = np.r_[
    sender_summary.supported_lr_pairs.to_numpy(float),
    receiver_summary.supported_lr_pairs.to_numpy(float),
]
strength_values = np.r_[
    sender_summary.aggregate_strength.to_numpy(float),
    receiver_summary.incoming_strength.to_numpy(float),
]
bars_num = ax_cnum.barh(bar_positions, count_values, height=0.66,
                        color=bar_colors, alpha=0.84, edgecolor="white", linewidth=0.50)
for y, bar in zip(bar_positions, bars_num):
    ax_cnum.text(bar.get_width() + 15, y, f"{int(bar.get_width())}",
                 ha="left", va="center", fontsize=5.0, color="#000000")
ax_cnum.set_xlim(0, 860); ax_cnum.set_ylim(8.65, -0.65)
ax_cnum.set_yticks(bar_positions, bar_labels)
for tick, color in zip(ax_cnum.get_yticklabels(), bar_colors):
    tick.set_color(color); tick.set_fontweight("bold")
ax_cnum.axhline(3.0, color="#D7DCE0", linewidth=0.55)
ax_cnum.set_xlabel("Supported LR pairs", labelpad=1)
ax_cnum.tick_params(axis="y", length=0, pad=2); clean_axis(ax_cnum)

ax_cstrength = fig.add_subplot(cbargrid[0, 1])
bars_strength = ax_cstrength.barh(bar_positions, strength_values, height=0.66,
                                  color=bar_colors, alpha=0.84, edgecolor="white", linewidth=0.50)
for y, bar in zip(bar_positions, bars_strength):
    ax_cstrength.text(bar.get_width() + 0.25, y, f"{bar.get_width():.2f}",
                      ha="left", va="center", fontsize=5.0, color="#000000")
ax_cstrength.set_xlim(0, 16.0); ax_cstrength.set_ylim(8.65, -0.65)
ax_cstrength.set_yticks(bar_positions, [])
ax_cstrength.axhline(3.0, color="#D7DCE0", linewidth=0.55)
ax_cstrength.set_xlabel("Aggregate strength", labelpad=1)
ax_cstrength.tick_params(axis="y", length=0); clean_axis(ax_cstrength)

ax_cnet = fig.add_subplot(cgrid[0, 1])
sender_angles = {"F1_PLATrich": np.deg2rad(145), "F2_MuralControl": np.deg2rad(180),
                 "F3_AntigenPresentation": np.deg2rad(215)}
receiver_angles = {cluster: np.deg2rad(angle) for cluster, angle in zip(STATE_ORDER, [58, 29, 0, -29, -58])}
node_positions = {sender: np.array([np.cos(angle), np.sin(angle)]) for sender, angle in sender_angles.items()}
node_positions.update({f"C{cluster}": np.array([np.cos(angle), np.sin(angle)]) for cluster, angle in receiver_angles.items()})
for _position in node_positions.values():
    _position[1] -= 0.08
edge_max = float(circle.aggregate_probability.max())
for source in ["F3_AntigenPresentation", "F2_MuralControl", "F1_PLATrich"]:
    source_index = SENDER_ORDER.index(source)
    for row in circle[circle.source.eq(source)].sort_values("cluster").itertuples(index=False):
        start = node_positions[source]; end = node_positions[f"C{int(row.cluster)}"]
        width = 0.35 + 2.20 * np.sqrt(float(row.aggregate_probability) / edge_max)
        curvature = 0.025 * (int(row.cluster) - 3) + 0.015 * (source_index - 1)
        ax_cnet.add_patch(FancyArrowPatch(
            start, end, arrowstyle="-|>", mutation_scale=5.0,
            connectionstyle=f"arc3,rad={curvature:.3f}", linewidth=width,
            color=SENDER_COLORS[source],
            alpha={"F1_PLATrich": 0.58, "F2_MuralControl": 0.31, "F3_AntigenPresentation": 0.22}[source],
            shrinkA=7.3, shrinkB=7.3, zorder=1 + source_index,
        ))
for source in SENDER_ORDER:
    x, y = node_positions[source]
    ax_cnet.add_patch(Circle((x, y), 0.090, facecolor=SENDER_COLORS[source],
                             edgecolor="white", linewidth=0.75, zorder=6))
    sender_dy = {"F1_PLATrich": 0.10, "F2_MuralControl": 0.00, "F3_AntigenPresentation": -0.10}[source]
    ax_cnet.text(x - 0.13, y + sender_dy, SENDER_LABELS[source],
                 ha="right", va="center", fontsize=5.2, color=SENDER_COLORS[source],
                 fontweight="bold" if source == "F1_PLATrich" else "normal")
for cluster in STATE_ORDER:
    x, y = node_positions[f"C{cluster}"]
    ax_cnet.add_patch(Circle((x, y), 0.090, facecolor=STATE_COLORS[cluster],
                             edgecolor="white", linewidth=0.75, zorder=6))
    receiver_dy = {1: 0.10, 2: 0.05, 3: 0.00, 4: -0.05, 5: -0.10}[cluster]
    ax_cnet.text(x + 0.14, y + receiver_dy, f"C{cluster}",
                 ha="left", va="center", fontsize=5.7,
                 color=STATE_COLORS[cluster], fontweight="bold")
ax_cnet.set_aspect("equal", adjustable="box")
ax_cnet.set_xlim(-1.42, 1.35); ax_cnet.set_ylim(-1.22, 1.22); ax_cnet.axis("off")


# d: LIANA LR heat map with strict-state recurrence at right.
dgrid = lower[0, 0].subgridspec(1, 3, width_ratios=[0.06, 0.78, 0.16], wspace=0.03)
ax_dspacer = fig.add_subplot(dgrid[0, 0]); ax_dspacer.axis("off")
ax_dheat = fig.add_subplot(dgrid[0, 1])
panel_label(ax_dheat, "d", -0.14, 1.04)
mag_max = max(4.0, float(np.nanquantile(mag_matrix.to_numpy(float), 0.98)))
lr_cmap = LinearSegmentedColormap.from_list("liana_magnitude", SEQUENTIAL)
lr_cmap.set_bad("#EEF0F2")
im_d = ax_dheat.imshow(np.ma.masked_invalid(mag_matrix.to_numpy(float)), cmap=lr_cmap,
                       vmin=0, vmax=mag_max, interpolation="nearest", aspect="auto")
ax_dheat.set_xticks(range(5), [f"C{x}" for x in STATE_ORDER])
for tick, cluster in zip(ax_dheat.get_xticklabels(), STATE_ORDER):
    tick.set_color(STATE_COLORS[cluster]); tick.set_fontweight("bold")
ax_dheat.set_yticks(range(len(pair_order)), pair_order)
for tick, pair in zip(ax_dheat.get_yticklabels(), pair_order):
    tick.set_fontstyle("italic")
    if pair.startswith("SPON2"):
        tick.set_color(ACCENT_RED); tick.set_fontweight("bold")
    elif pair.startswith("PLAT"):
        tick.set_color("#7A8288")
ax_dheat.tick_params(length=0, pad=0.8)
for y in range(len(pair_order) + 1):
    ax_dheat.axhline(y - 0.5, color="white", linewidth=0.55)
for x in range(6):
    ax_dheat.axvline(x - 0.5, color="white", linewidth=0.55)
for y in range(len(pair_order)):
    for x in range(5):
        if bool(star_matrix.iloc[y, x]):
            ax_dheat.text(x, y, "*", ha="center", va="center",
                          fontsize=7.0, fontweight="bold", color=BLACK, zorder=4)
for spine in ax_dheat.spines.values():
    spine.set_visible(False)
cax_d = ax_dheat.inset_axes([0.12, -0.075, 0.50, 0.030])
cb_d = fig.colorbar(im_d, cax=cax_d, orientation="horizontal")
cb_d.set_ticks([0, round(mag_max / 2, 2), round(mag_max, 2)])
cb_d.ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
cb_d.ax.tick_params(labelsize=5, length=1, pad=0.4)
cb_d.set_label(r"$-\log_{10}$(magnitude rank)", fontsize=5.2, labelpad=0)
cb_d.outline.set_linewidth(0.35)

ax_dright = fig.add_subplot(dgrid[0, 2])
ypos = np.arange(len(pair_order))
ax_dright.barh(ypos, strict_counts.to_numpy(), height=0.62,
               color="#AAB3BA", edgecolor="white", linewidth=0.4)
ax_dright.set_ylim(len(pair_order) - 0.5, -0.5); ax_dright.set_xlim(0, 5.3)
ax_dright.set_xticks([0, 2, 4]); ax_dright.set_yticks([])
ax_dright.set_xlabel("States", labelpad=1); ax_dright.tick_params(axis="x", labelsize=5, length=2)
clean_axis(ax_dright)


# e: patient recurrence of four strict F1-to-C3 LR candidates.
egrid = rightstack[0, 0].subgridspec(
    2, 2, height_ratios=[0.24, 0.76], width_ratios=[0.92, 0.08],
    hspace=0.04, wspace=0.20,
)
ax_etop = fig.add_subplot(egrid[0, 0])
panel_label(ax_etop, "e", -0.18, 1.10)
support_counts = focus_support.set_index("pair").reindex(focus_pairs).supported_patients.astype(int).to_numpy()
evaluable_counts = focus_support.set_index("pair").reindex(focus_pairs).evaluable_patients.astype(int).to_numpy()
e_colors = ["#6F8EA5", "#7A9CB5", "#8AA8BD", ACCENT_RED]
ax_etop.bar(np.arange(4), support_counts, width=0.68, color=e_colors, alpha=0.84)
for x, (supported, evaluable) in enumerate(zip(support_counts, evaluable_counts)):
    ax_etop.text(x, supported + 0.55, f"{supported}/{evaluable}",
                 ha="center", va="bottom", fontsize=5.0, color="#000000")
ax_etop.set_xlim(-0.5, 3.5); ax_etop.set_ylim(0, 25)
ax_etop.set_ylabel("Patients", labelpad=1); ax_etop.set_xticks([])
ax_etop.spines["bottom"].set_visible(False); ax_etop.tick_params(axis="y", labelsize=5, length=2)

ax_eheat = fig.add_subplot(egrid[1, 0])
z_matrix = (
    focus_donor.pivot(index="Patient_ID", columns="pair", values="within_pair_z")
    .reindex(index=donor_order, columns=focus_pairs)
)
e_cmap = LinearSegmentedColormap.from_list("donor_z", DIVERGING)
e_cmap.set_bad("#E3E6E8")
im_e = ax_eheat.imshow(np.ma.masked_invalid(z_matrix.to_numpy(float)), cmap=e_cmap,
                       vmin=-2, vmax=2, interpolation="nearest", aspect="auto")
ax_eheat.set_xticks(range(4), focus_pairs, rotation=45, ha="right", rotation_mode="anchor")
for tick in ax_eheat.get_xticklabels():
    tick.set_fontstyle("italic")
ax_eheat.set_yticks(range(len(donor_order)), [donor_codes[x] for x in donor_order])
ax_eheat.tick_params(length=0, pad=0.3, labelsize=5)
for y in range(len(donor_order) + 1):
    ax_eheat.axhline(y - 0.5, color="white", linewidth=0.35)
for x in range(5):
    ax_eheat.axvline(x - 0.5, color="white", linewidth=0.50)
for spine in ax_eheat.spines.values():
    spine.set_visible(False)
sp_idx = next(i for i, pair in enumerate(focus_pairs) if pair.startswith("SPON2"))
ax_eheat.add_patch(Rectangle((sp_idx - 0.47, -0.47), 0.94, len(donor_order) - 0.06,
                             fill=False, edgecolor=ACCENT_RED, linewidth=0.80, clip_on=False))
cax_e = fig.add_subplot(egrid[1, 1])
cb_e = fig.colorbar(im_e, cax=cax_e, orientation="vertical")
cb_e.set_ticks([-2, 0, 2]); cb_e.ax.tick_params(labelsize=5, length=1, pad=1)
cb_e.set_label("Within-pair Z score", fontsize=5.2, labelpad=1)
cb_e.outline.set_linewidth(0.35)
ax_eblank = fig.add_subplot(egrid[0, 1]); ax_eblank.axis("off")


# f: patient-level ecological association for recurrent states.
ax_feco = fig.add_subplot(rightstack[1, 0])
panel_label(ax_feco, "f", -0.18, 1.16)
eco_order = [1, 3, 4, 5]
for y, cluster in enumerate(eco_order):
    row = eco_plot.loc[cluster]
    ax_feco.hlines(y, float(row.loo_adjusted_rho_min), float(row.loo_adjusted_rho_max),
                   color=STATE_COLORS[cluster], linewidth=1.8, alpha=0.82)
    ax_feco.scatter(float(row.adjusted_rho), y, s=34, color=STATE_COLORS[cluster],
                    edgecolor="white", linewidth=0.60, zorder=3)
    fdr = float(row.adjusted_permutation_FDR_BH)
    star = " *" if fdr < 0.05 else ""
    ax_feco.text(0.84, y, f"FDR {fdr:.2e}{star}",
                 ha="right", va="center", fontsize=5.0, color="#000000")
ax_feco.axvline(0, color="#7F888E", linewidth=0.65, linestyle=(0, (3, 2)))
ax_feco.set_xlim(-0.48, 0.88); ax_feco.set_ylim(3.60, -0.60)
ax_feco.set_yticks(range(4), [f"C{x}" for x in eco_order])
for tick, cluster in zip(ax_feco.get_yticklabels(), eco_order):
    tick.set_color(STATE_COLORS[cluster]); tick.set_fontweight("bold")
ax_feco.set_xlabel(r"Adjusted $\rho$", labelpad=1)
ax_feco.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
ax_feco.grid(axis="x", color="#EDF0F2", linewidth=0.45, zorder=0)
ax_feco.tick_params(axis="y", length=0, labelsize=5.2)
ax_feco.tick_params(axis="x", labelsize=5.0)
clean_axis(ax_feco)


stem = "Figure3_PLAT_reference_communication_layout_v4_candidate"
base = FIG_DIR / stem
save_cns_figure(fig, str(base))
fig.savefig(str(base) + ".svg", bbox_inches=None, dpi=300)
fig.savefig(str(base) + ".tiff", bbox_inches=None, dpi=600,
            pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)


# Supplementary Figure S3: patient x C1-C5 audit heat map.
programme_wide = programme.pivot(index="Patient_ID", columns="cluster", values="lr_programme_potential").reindex(columns=STATE_ORDER)
supp_order = (
    programme_wide.assign(
        _n=programme_wide.notna().sum(axis=1),
        _c3=programme_wide[3].fillna(-np.inf),
    )
    .sort_values(["_n", "_c3"], ascending=[False, False])
    .index.tolist()
)
supp_codes = {donor: f"P{i + 1:02d}" for i, donor in enumerate(supp_order)}
supp_matrix = programme_wide.reindex(index=supp_order, columns=STATE_ORDER)

fig_s = plt.figure(figsize=(89 / 25.4, 118 / 25.4), facecolor="white")
sg = fig_s.add_gridspec(
    2, 1, height_ratios=[0.20, 0.80],
    left=0.23, right=0.93, top=0.96, bottom=0.13, hspace=0.04,
)
ax_stop = fig_s.add_subplot(sg[0, 0])
panel_label(ax_stop, "a", -0.28, 1.14)
eval_counts = programme_wide.notna().sum(axis=0).reindex(STATE_ORDER).astype(int).to_numpy()
ax_stop.bar(np.arange(5), eval_counts, width=0.70,
            color=[STATE_COLORS[c] for c in STATE_ORDER], alpha=0.82)
for x, value in enumerate(eval_counts):
    ax_stop.text(x, value + 0.6, str(value), ha="center", va="bottom", fontsize=5.2)
ax_stop.set_xlim(-0.5, 4.5); ax_stop.set_ylim(0, 25)
ax_stop.set_ylabel("Patients", labelpad=1); ax_stop.set_xticks([])
ax_stop.spines["bottom"].set_visible(False); ax_stop.tick_params(axis="y", labelsize=5, length=2)

ax_sh = fig_s.add_subplot(sg[1, 0])
supp_cmap = LinearSegmentedColormap.from_list("supp_programme", SEQUENTIAL)
supp_cmap.set_bad("#E3E6E8")
vmin = float(np.nanquantile(supp_matrix.to_numpy(float), 0.02))
vmax = float(np.nanquantile(supp_matrix.to_numpy(float), 0.98))
im_s = ax_sh.imshow(np.ma.masked_invalid(supp_matrix.to_numpy(float)), cmap=supp_cmap,
                    vmin=vmin, vmax=vmax, interpolation="nearest", aspect="auto")
ax_sh.set_xticks(range(5), [f"C{x}" for x in STATE_ORDER])
for tick, cluster in zip(ax_sh.get_xticklabels(), STATE_ORDER):
    tick.set_color(STATE_COLORS[cluster]); tick.set_fontweight("bold")
ax_sh.set_yticks(range(len(supp_order)), [supp_codes[x] for x in supp_order])
ax_sh.tick_params(length=0, pad=1.5, labelsize=5)
for y in range(len(supp_order) + 1):
    ax_sh.axhline(y - 0.5, color="white", linewidth=0.35)
for x in range(6):
    ax_sh.axvline(x - 0.5, color="white", linewidth=0.45)
for spine in ax_sh.spines.values():
    spine.set_visible(False)
cax_s = ax_sh.inset_axes([0.18, -0.095, 0.64, 0.025])
cb_s = fig_s.colorbar(im_s, cax=cax_s, orientation="horizontal")
cb_s.ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
cb_s.ax.tick_params(labelsize=5, length=1, pad=0.4)
cb_s.set_label("LR programme potential", fontsize=5.2, labelpad=0)
cb_s.outline.set_linewidth(0.35)

supp_stem = "FigureS3_patient_state_LR_programme_heatmap_v4_candidate"
supp_base = FIG_DIR / supp_stem
save_cns_figure(fig_s, str(supp_base))
fig_s.savefig(str(supp_base) + ".svg", bbox_inches=None, dpi=300)
fig_s.savefig(str(supp_base) + ".tiff", bbox_inches=None, dpi=600,
              pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig_s)


# Supplementary Figure S4: independent cross-cohort myeloid-state mapping.
# M04 C3+ denotes expression of the complement gene C3 in the query cohort;
# it is not the C3 cluster label used for the reference cohort columns.
fig_s4 = plt.figure(figsize=(89 / 25.4, 72 / 25.4), facecolor="white")
ax_s4 = fig_s4.add_axes([0.40, 0.24, 0.52, 0.68])
s4_cmap = LinearSegmentedColormap.from_list("mapping_rho", DIVERGING)
im_s4 = ax_s4.imshow(map_matrix.to_numpy(float), cmap=s4_cmap,
                     vmin=-0.8, vmax=0.8, interpolation="nearest", aspect="auto")
ax_s4.set_xticks(range(len(REFERENCE_ORDER)), ["C1", "C3", "C4", "C5"])
for tick, cluster in zip(ax_s4.get_xticklabels(), [1, 3, 4, 5]):
    tick.set_color(STATE_COLORS[cluster]); tick.set_fontweight("bold")
ax_s4.set_yticks(range(len(QUERY_ORDER)), [QUERY_LABELS[x] for x in QUERY_ORDER])
ax_s4.tick_params(length=0, pad=2, labelsize=5.5)
for y in range(map_matrix.shape[0]):
    row_values = map_matrix.iloc[y].to_numpy(float)
    row_best = int(np.nanargmax(row_values))
    for x, value in enumerate(row_values):
        ax_s4.text(x, y, f"{value:.2f}", ha="center", va="center",
                   fontsize=5.4, color="white" if abs(value) >= 0.50 else BLACK,
                   fontweight="bold" if x == row_best else "normal")
for y in range(map_matrix.shape[0] + 1):
    ax_s4.axhline(y - 0.5, color="white", linewidth=0.55)
for x in range(map_matrix.shape[1] + 1):
    ax_s4.axvline(x - 0.5, color="white", linewidth=0.55)
for spine in ax_s4.spines.values():
    spine.set_visible(False)
cax_s4 = fig_s4.add_axes([0.46, 0.10, 0.40, 0.035])
cb_s4 = fig_s4.colorbar(im_s4, cax=cax_s4, orientation="horizontal")
cb_s4.set_ticks([-0.8, 0, 0.8]); cb_s4.ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f")); cb_s4.ax.tick_params(labelsize=5, length=1, pad=0.4)
cb_s4.set_label("Spearman ρ", fontsize=5.4, labelpad=0)
cb_s4.outline.set_linewidth(0.35)

s4_stem = "FigureS4_cross_cohort_myeloid_mapping_v4_candidate"
s4_base = FIG_DIR / s4_stem
save_cns_figure(fig_s4, str(s4_base))
fig_s4.savefig(str(s4_base) + ".svg", bbox_inches=None, dpi=300)
fig_s4.savefig(str(s4_base) + ".tiff", bbox_inches=None, dpi=600,
               pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig_s4)


manifest = pd.DataFrame([
    ["a", "Figure3a_myeloid_umap.tsv.gz", "cell (descriptive)", "myeloid-state overview"],
    ["b", "Figure3b_patient_balanced_marker_profile.tsv", "patient-balanced state", "state annotation"],
    ["c", "Figure3c_sender_communication_summary.tsv", "donor-balanced pooled cells", "sender breadth and strength"],
    ["c", "Figure3c_receiver_communication_summary.tsv", "donor-balanced pooled cells", "receiver breadth and strength"],
    ["c", "Figure3c_fibroblast_myeloid_circle.tsv", "donor-balanced pooled cells", "complete fibroblast-to-myeloid network"],
    ["d", "Figure3d_LIANA_LR_matrix.tsv", "donor-balanced pooled cells", "specific LR landscape"],
    ["d", "Figure3d_LIANA_LR_display_order.tsv", "ligand-receptor pair", "best-magnitude-rank display order"],
    ["e", "Figure3e_C3_candidate_patient_recurrence.tsv", "patient", "four strict F1-to-C3 candidates"],
    ["f", "Figure3f_patient_ecology.tsv", "patient", "sender-receiver ecology"],
    ["S3", "FigureS3_patient_state_programme.tsv", "patient", "all-five-state LR-programme audit"],
    ["S4", "FigureS4_cross_cohort_mapping.tsv", "donor-equal state centroid", "independent state mapping"],
], columns=["panel", "source_file", "statistical_unit", "purpose"])
export_tsv(ecology, "Figure3f_patient_ecology.tsv")
export_tsv(mapping, "FigureS4_cross_cohort_mapping.tsv")
export_tsv(manifest, "Figure3_reference_layout_source_data_manifest.tsv")

qa = {
    "figure": stem,
    "canvas_mm": [183, 225],
    "cells": int(len(cells)),
    "patients": int(cells.Patient_ID.nunique()),
    "umap_equal_aspect": True,
    "circle_edges": int(len(circle)),
    "sender_summary": sender_summary.to_dict(orient="records"),
    "receiver_summary": receiver_summary.to_dict(orient="records"),
    "f1_receiver_supported_pairs": {f"C{c}": int(receiver_plot.loc[c, "supported_lr_pairs"]) for c in STATE_ORDER},
    "liana_display_pairs": int(len(pair_order)),
    "liana_ordering": "best magnitude rank across C1-C5, ascending",
    "starred_liana_cells": int(star_matrix.to_numpy().sum()),
    "supplement_S4_mapping_query_labels": QUERY_LABELS,
    "focus_candidate_support": focus_support.to_dict(orient="records"),
    "focus_candidate_patients": int(focus_donor.Patient_ID.nunique()),
    "supplement_patient_state_estimates": int(len(programme)),
    "outputs": [str(base) + suffix for suffix in [".pdf", ".png", ".svg", ".tiff"]],
    "supplement_S3_outputs": [str(supp_base) + suffix for suffix in [".pdf", ".png", ".svg", ".tiff"]],
    "supplement_S4_outputs": [str(s4_base) + suffix for suffix in [".pdf", ".png", ".svg", ".tiff"]],
}
(OUT_DIR / "Figure3_reference_layout_numeric_QA.json").write_text(
    json.dumps(qa, indent=2), encoding="utf-8"
)
print(json.dumps(qa, indent=2))
