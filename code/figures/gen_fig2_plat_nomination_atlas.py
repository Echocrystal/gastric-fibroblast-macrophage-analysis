# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) four-set intersection -> project-native R/ggvenn production panel -> native run
# (b) ranked candidate scatter -> BarComparison cross-type -> param inherit
# (c) patient-state scatter -> project-native Figure 2 production panel -> param inherit
# (d) paired patient slope plot -> PairedBoxScatter -> param inherit
# (e) single-cell cluster UMAP -> PCA cross-type + project-native Figure 2 production panel -> param inherit
# (f) single-cell lineage UMAP -> PCA cross-type + project-native Figure 2 production panel -> param inherit
# (g) broad-lineage marker dot plot -> MarkerGeneDotPlot -> param inherit
# (h) fibroblast-state UMAP -> PCA cross-type + project-native Figure 2 production panel -> param inherit
# (i) PLAT feature UMAP -> PCA cross-type + project-native Figure 2 production panel -> param inherit
# (j) fibroblast-state marker dot plot -> MarkerGeneDotPlot -> param inherit
# (k) patient-by-state heatmap -> heatmap/composition_heatmap -> param inherit
# RULE: all panels use complete frozen source tables; no cell or patient rows are subsampled here.

from __future__ import annotations

from pathlib import Path
import hashlib
import subprocess
import sys

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import matplotlib.patheffects as pe
from scipy.cluster.hierarchy import linkage, leaves_list


# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
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
BLACK       = "#222222"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})


def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)


PROJECT = Path(r"I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026")
SOURCE = PROJECT / "figures" / "source_data" / "Figure2"
CANDIDATE_SOURCE = PROJECT / "figures" / "data" / "Figure1_GSE206785_unbiased"
R_EXECUTABLE = Path(r"C:/Program Files/R/R-4.4.0/bin/Rscript.exe")
GGVENN_SCRIPT = PROJECT / "figures" / "gen_fig2a_ggvenn.R"
GGVENN_PNG = (PROJECT / "figures" / "panels" / "Figure2_PLAT_nomination" /
              "Figure2A_ggvenn.png")

INK = "#000000"
TEXT_GREY = "#000000"
LIGHT_GREY = "#D9DEE2"
PLAT = "#9E3D4F"
FIB = "#9D6875"

LINEAGE_COLORS = {
    "Epithelial": "#D59A91",
    "Endothelial": "#2A9D78",
    "Fibroblast": FIB,
    "Mural": "#98974A",
    "T_cell": "#66A7D5",
    "NK": "#D59628",
    "B_cell": "#76549B",
    "Plasma": "#C678A4",
    "Myeloid": "#2878A8",
    "Mast": "#8C6D45",
    "Unresolved": "#A9ADB1",
}

LINEAGE_DISPLAY = {
    "Epithelial": "Epithelial",
    "Endothelial": "Endothelial",
    "Fibroblast": "Fibroblast",
    "Mural": "Mural",
    "T_cell": "T cell",
    "NK": "NK",
    "B_cell": "B cell",
    "Plasma": "Plasma",
    "Myeloid": "Myeloid",
    "Mast": "Mast",
    "Unresolved": "Unresolved",
}

FIB_COLORS = {
    1: PLAT,
    2: "#7E8589",
    3: "#68A9D3",
    4: "#B57A9A",
    5: "#B9B9B9",
    6: "#3B9C78",
}

FIB_LABELS = {
    1: "CXCL14+POSTN act.",
    2: "Mural-like",
    3: "Immune-associated",
    4: "Matrix",
    5: "Epithelial-like",
    6: "PI16+",
}

# Thirty explicit colours avoid matplotlib tab20/tab20b defaults. Numeric cluster
# labels provide the primary identifier; colour is a secondary cue.
CLUSTER_COLORS = [
    "#6F8FB8", "#C57A82", "#7FA37C", "#C99A5E", "#8D79A8", "#6E9D9A",
    "#B77D5B", "#8A96A3", "#A87A95", "#8FAE6D", "#7099C0", "#D08A73",
    "#6F9B83", "#BEA15E", "#927DAE", "#709C9A", "#B36E60", "#8997A5",
    "#A67893", "#89A96F", "#5D86AE", "#C36D78", "#669071", "#B78C4F",
    "#7E6A9C", "#5F8D8B", "#A9634B", "#74838E", "#966A84", "#75945E",
]

MAT_CLASS_COLORS = {
    "ECM regulator": PLAT,
    "ECM glycoprotein": "#8B75A0",
    "Proteoglycan": "#6F9277",
    "Collagen": "#6E91B3",
}

DOT_CMAP = LinearSegmentedColormap.from_list(
    "scaled_expression", ["#477AA6", "#FFFFFF", "#B65C6A"]
)
DOT_VMAX = 1.5
DOT_MIN_SIZE = 3.0
DOT_MAX_ADDITION = 38.0

PLAT_CMAP = LinearSegmentedColormap.from_list(
    "plat_expression", ["#F4EFEA", "#E6C7BE", "#D79782", "#B85E61", PLAT]
)

PATIENT_CMAP = LinearSegmentedColormap.from_list(
    "patient_fraction_rose",
    ["#FAF9F7", "#EEE3E1", "#D8B7BA", "#B97C88", "#8D485B", "#5F2E40"],
)


def dot_size(detected_fraction: float | np.ndarray) -> float | np.ndarray:
    fraction = np.clip(detected_fraction, 0, 1)
    return DOT_MIN_SIZE + DOT_MAX_ADDITION * fraction ** 1.15


def stable_jitter(label: str, amplitude: float = 0.18) -> float:
    digest = hashlib.sha256(label.encode("utf-8")).digest()
    unit = int.from_bytes(digest[:8], "big") / (2**64 - 1)
    return amplitude * (2 * unit - 1)


def panel_label(ax, label: str, x: float = -0.10, y: float = 1.05) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=9.2, fontweight="bold",
            ha="left", va="top", color=INK, clip_on=False)


def clean_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=2.5, width=0.6, colors="#000000")
    ax.xaxis.label.set_color("#000000")
    ax.yaxis.label.set_color("#000000")
    ax.spines["left"].set_color("#000000")
    ax.spines["bottom"].set_color("#000000")


def umap_axes(ax, label_axes: bool = True) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    if label_axes:
        ax.set_xlabel("UMAP 1", fontsize=5.5, color=TEXT_GREY, labelpad=2)
        ax.set_ylabel("UMAP 2", fontsize=5.5, color=TEXT_GREY, labelpad=3)


def padded_limits(values: np.ndarray, fraction: float = 0.04) -> tuple[float, float]:
    lo, hi = float(np.min(values)), float(np.max(values))
    pad = (hi - lo) * fraction
    return lo - pad, hi + pad


def whitened_text(ax, x, y, text, **kwargs):
    artist = ax.text(x, y, text, **kwargs)
    artist.set_path_effects([pe.withStroke(linewidth=1.6, foreground="white")])
    return artist


def save_exact_canvas(fig, stem: str) -> None:
    out = PROJECT / "figures" / "main"
    out.mkdir(parents=True, exist_ok=True)
    with mpl.rc_context({"savefig.bbox": None, "savefig.pad_inches": 0}):
        fig.savefig(out / f"{stem}.pdf", bbox_inches=None)
        fig.savefig(out / f"{stem}.svg", bbox_inches=None)
        fig.savefig(out / f"{stem}.png", dpi=400, bbox_inches=None)
        fig.savefig(out / f"{stem}.tiff", dpi=600, bbox_inches=None,
                    pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def load_and_validate():
    exclusive = pd.read_csv(
        CANDIDATE_SOURCE / "panel_a_four_set_activity_exclusive_regions.tsv", sep="\t"
    )
    ranking = pd.read_csv(
        CANDIDATE_SOURCE / "panel_b_activity_gated_ligand_ranking.tsv", sep="\t"
    ).sort_values("scan_rank_activity")
    atlas = pd.read_csv(SOURCE / "Figure2_panels_ab_global_atlas_cells.tsv.gz", sep="\t")
    broad_markers = pd.read_csv(
        SOURCE / "Figure2_panel_c_broad_lineage_canonical_markers.tsv", sep="\t"
    )
    fib = pd.read_csv(SOURCE / "Figure2_panels_dg_fibroblast_cells.tsv.gz", sep="\t")
    fib_markers = pd.read_csv(
        SOURCE / "Figure2_panel_h_fibroblast_state_markers.tsv", sep="\t"
    )
    patient_state = pd.read_csv(
        SOURCE / "Figure2_panel_i_fibroblast_state_by_patient.tsv", sep="\t"
    )
    localization = pd.read_csv(
        SOURCE / "Figure2_panel_j_GSE206785_PLAT_state_localization.tsv", sep="\t"
    )
    paired = pd.read_csv(
        SOURCE / "Figure2_panel_k_GSE206785_patient_paired_PLAT.tsv", sep="\t"
    )

    if len(ranking) != 32 or ranking["scan_rank_activity"].nunique() != 32:
        raise RuntimeError("Expected the frozen 32-ligand activity-gated ranking")
    plat_row = ranking.loc[ranking["ligand"].eq("PLAT")].iloc[0]
    if int(plat_row["scan_rank_activity"]) != 5:
        raise RuntimeError("PLAT must remain rank 5/32 in the frozen ranking")
    four_way = int(exclusive.loc[exclusive["membership_code"].astype(str).eq("1111"),
                                 "exclusive_gene_count"].iloc[0])
    if four_way != 32:
        raise RuntimeError("Expected an exact four-way intersection of 32 ligands")
    if len(atlas) != 109788 or len(fib) != 8904:
        raise RuntimeError("Figure 2 cell-universe mismatch")
    if patient_state["Patient_ID"].nunique() != 26:
        raise RuntimeError("Expected 26 patients in the GSE183904 recurrence panel")

    return exclusive, ranking, atlas, broad_markers, fib, fib_markers, patient_state, localization, paired


def set_sizes_from_exclusive(exclusive: pd.DataFrame) -> dict[str, int]:
    codes = exclusive["membership_code"].astype(str).str.zfill(4)
    counts = exclusive["exclusive_gene_count"].astype(int).to_numpy()
    names = ["Matrisome", "LR prior", "Fibroblast expression", "NicheNet activity"]
    return {name: int(sum(c for code, c in zip(codes, counts) if code[i] == "1"))
            for i, name in enumerate(names)}


def render_native_ggvenn() -> None:
    """Render panel a with the actual R ggvenn implementation."""
    if not R_EXECUTABLE.exists():
        raise RuntimeError(f"Rscript not found: {R_EXECUTABLE}")
    result = subprocess.run(
        [str(R_EXECUTABLE), str(GGVENN_SCRIPT)],
        cwd=str(PROJECT), capture_output=True, text=True,
        encoding="utf-8", errors="replace",
    )
    if result.returncode != 0:
        raise RuntimeError(
            "Native ggvenn render failed\n"
            f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
        )
    if not GGVENN_PNG.exists() or GGVENN_PNG.stat().st_size == 0:
        raise RuntimeError("Native ggvenn PNG was not produced")


def draw_candidate_intersection(ax, exclusive: pd.DataFrame) -> None:
    panel_label(ax, "a", -0.06, 1.04)
    image = plt.imread(GGVENN_PNG)
    # Crop only the device-level white margin introduced by ggplot; the
    # ggvenn geometry, labels and all 15 exact region counts remain unchanged.
    rgb = image[..., :3]
    content = np.any(rgb < 0.985, axis=2)
    yy, xx = np.where(content)
    pad = 28
    y0 = max(0, int(yy.min()) - pad)
    y1 = min(image.shape[0], int(yy.max()) + pad + 1)
    x0 = max(0, int(xx.min()) - pad)
    x1 = min(image.shape[1], int(xx.max()) + pad + 1)
    image = image[y0:y1, x0:x1]
    ax.imshow(image, interpolation="lanczos")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_candidate_ranking(ax, ranking: pd.DataFrame) -> None:
    panel_label(ax, "b", -0.105, 1.11)
    x = ranking["scan_rank_activity"].to_numpy(dtype=int)
    y = ranking["evidence_score_activity"].to_numpy(dtype=float)
    colors = [MAT_CLASS_COLORS.get(c, "#999999") for c in ranking["matrisome_class"]]
    ax.plot(x, y, color="#C8CDD1", linewidth=0.7, zorder=1)
    ax.scatter(x, y, s=14, color=colors, alpha=0.72,
               edgecolor="white", linewidth=0.35, zorder=3)
    plat = ranking[ranking["ligand"].eq("PLAT")].iloc[0]
    ax.scatter([plat["scan_rank_activity"]], [plat["evidence_score_activity"]],
               s=42, marker="D", color=PLAT, edgecolor="white", linewidth=0.6, zorder=5)
    label_offsets = {
        "LUM": (0, 7),
        "MMP2": (-6, -9),
        "VCAN": (-2, 7),
        "TIMP1": (4, -9),
        "PLAT": (6, 8),
    }
    for row in ranking.head(5).itertuples(index=False):
        dx, dy = label_offsets[row.ligand]
        ax.annotate(row.ligand,
                xy=(row.scan_rank_activity, row.evidence_score_activity),
                xytext=(dx, dy), textcoords="offset points",
                ha="center", va="center", fontsize=5.2,
                color=PLAT if row.ligand == "PLAT" else INK,
                fontstyle="italic", fontweight="bold" if row.ligand == "PLAT" else "normal")
    ax.text(float(plat["scan_rank_activity"]) + 1.2,
            float(plat["evidence_score_activity"]) - 0.002,
            "5/32", fontsize=5.0, color=PLAT, va="center")
    ax.set_xlim(0.2, 32.8); ax.set_ylim(0.12, 1.03)
    ax.set_xticks([1, 5, 10, 15, 20, 25, 32])
    ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
    ax.set_xlabel("Candidate rank", labelpad=1.5)
    ax.set_ylabel("Activity-aware evidence score", labelpad=2)
    ax.grid(axis="y", color="#ECEFF1", linewidth=0.45)
    clean_axis(ax)
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=3.4,
                      markerfacecolor=color, markeredgecolor="none", label=label)
               for label, color in MAT_CLASS_COLORS.items()]
    ax.legend(handles=handles, ncol=2, loc="lower left", bbox_to_anchor=(0.01, 0.01),
              fontsize=5.0, columnspacing=0.65, handletextpad=0.2, borderaxespad=0)


def main() -> None:
    (exclusive, ranking, atlas, broad_markers, fib, fib_markers,
     patient_state, localization, paired) = load_and_validate()
    render_native_ggvenn()

    mpl.rcParams.update({
        "font.size": 6.2,
        "axes.titlesize": 6.6,
        "axes.labelsize": 6.0,
        "xtick.labelsize": 5.2,
        "ytick.labelsize": 5.2,
        "legend.fontsize": 5.0,
    })

    fig = plt.figure(figsize=(7.20, 9.45), facecolor="white")
    gs = fig.add_gridspec(
        5, 13,
        height_ratios=[1.52, 1.38, 2.05, 2.08, 1.72],
        hspace=0.62, wspace=0.72,
        left=0.105, right=0.985, bottom=0.050, top=0.980,
    )

    # a-b — outcome-blind candidate nomination
    ax_a = fig.add_subplot(gs[0, 0:5])
    pos_a = ax_a.get_position()
    # The native ggvenn geometry is taller than the former four-petal icon.
    # Use the otherwise empty a-to-c gutter so its labels remain legible at
    # print size without moving or shrinking panel b.
    ax_a.set_position([
        pos_a.x0 - 0.004, pos_a.y0 - 0.044,
        pos_a.width + 0.008, pos_a.height + 0.044,
    ])
    draw_candidate_intersection(ax_a, exclusive)
    ax_b = fig.add_subplot(gs[0, 5:13])
    pos_b = ax_b.get_position()
    ax_b.set_position([pos_b.x0 + 0.012, pos_b.y0,
                       pos_b.width - 0.012, pos_b.height])
    draw_candidate_ranking(ax_b, ranking)

    # c — GSE206785 patient-level state localisation
    ax = fig.add_subplot(gs[1, 0:6])
    pos = ax.get_position()
    ax.set_position([pos.x0 + 0.020, pos.y0, pos.width - 0.020, pos.height])
    panel_label(ax, "c", -0.08, 1.08)
    state_order = ["F08_Activated_VSTM2A", "F13_Activated_CTHRC1", "F14_Activated_ADAM28",
                   "F09_Quiescent_IGF1", "F06_Quiescent_"]
    state_display = ["F08 VSTM2A act.", "F13 CTHRC1 act.", "F14 ADAM28 act.",
                     "F09 IGF1 quiet.", "F06 quiescent"]
    y_positions = np.arange(5)[::-1]
    n_by_state = localization.groupby("Annotation")["Patient"].nunique().to_dict()
    for yy, state in zip(y_positions, state_order):
        d = localization[localization["Annotation"].eq(state)].copy()
        active = "Activated" in state
        color = PLAT if active else "#AEB4B8"
        offsets = np.array([stable_jitter(f"{state}|{p}", 0.16) for p in d["Patient"].astype(str)])
        ax.scatter(d["patient_state_mean_raw"], yy + offsets, s=14,
                   color=color, alpha=0.72, edgecolor="white", linewidth=0.35, zorder=3)
        median = float(d["patient_state_mean_raw"].median())
        ax.plot([median, median], [yy - 0.23, yy + 0.23], color=INK,
                linewidth=1.5, solid_capstyle="round", zorder=5)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([f"{label}  (n={n_by_state[state]})"
                        for state, label in zip(state_order, state_display)])
    ax.set_xlim(-0.25, max(7.0, float(localization["patient_state_mean_raw"].max()) * 1.07))
    ax.set_xlabel("Patient-state PLAT mean raw count", labelpad=1.5)
    ax.grid(axis="x", color="#EBEEF0", linewidth=0.45)
    clean_axis(ax); ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, labelsize=5.0, pad=1)

    # d — GSE206785 paired patient comparison; old 72-ligand FDR intentionally omitted
    ax = fig.add_subplot(gs[1, 6:13])
    panel_label(ax, "d", -0.07, 1.08)
    paired = paired.sort_values("F13_minus_other")
    for row in paired.itertuples(index=False):
        higher = row.F13_minus_other > 0
        color = PLAT if higher else "#7D858B"
        ax.plot([0, 1], [row.F13_mean_raw, row.other_fibroblast_mean_raw],
                color=color, alpha=0.48, linewidth=0.75, zorder=1)
        ax.scatter([0, 1], [row.F13_mean_raw, row.other_fibroblast_mean_raw],
                   s=14, color=color, edgecolor="white", linewidth=0.35, zorder=3)
    ax.set_xlim(-0.30, 1.30)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["F13 CTHRC1\nactivated", "Other\nfibroblasts"])
    ax.set_ylabel("Patient mean raw count", labelpad=2)
    p_value = float(paired["wilcoxon_p"].iloc[0])
    median_delta = float(paired["median_difference"].iloc[0])
    positive_pairs = int(paired["positive_pairs"].iloc[0])
    paired_n = int(paired["paired_patients"].iloc[0])
    p_mantissa, p_exponent = f"{p_value:.2e}".split("e")
    ax.text(0.98, 0.97,
            rf"{positive_pairs}/{paired_n} higher" "\n"
            rf"median $\Delta$ = {median_delta:.2f}" "\n"
            rf"$P$ = {float(p_mantissa):.2f}$\times 10^{{{int(p_exponent)}}}$",
            transform=ax.transAxes, fontsize=5.3, color=INK, ha="right", va="top")
    clean_axis(ax)

    # e — GSE183904 genuine global unsupervised clusters
    ax = fig.add_subplot(gs[2, 0:4])
    panel_label(ax, "e", -0.12, 1.06)
    for cluster in range(1, 31):
        d = atlas[atlas["global_cluster"].eq(cluster)]
        ax.scatter(d["UMAP1"], d["UMAP2"], s=0.28, color=CLUSTER_COLORS[cluster - 1],
                   alpha=0.54, linewidths=0, rasterized=True)
        x_med, y_med = d[["UMAP1", "UMAP2"]].median().to_numpy()
        whitened_text(ax, x_med, y_med, str(cluster), fontsize=5.0, fontweight="bold",
                      color=INK, ha="center", va="center", zorder=8)
    xlim = padded_limits(atlas["UMAP1"].to_numpy())
    ylim = padded_limits(atlas["UMAP2"].to_numpy(), fraction=0.08)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); umap_axes(ax)

    # f — exact same coordinates with broad annotations
    ax = fig.add_subplot(gs[2, 4:8])
    panel_label(ax, "f", -0.12, 1.06)
    atlas = atlas.copy()
    atlas["plot_label"] = atlas["coarse_label"].replace({"Ambiguous": "Unresolved", "Unknown": "Unresolved"})
    annotation_order = list(LINEAGE_COLORS)
    for label in annotation_order:
        d = atlas[atlas["plot_label"].eq(label)]
        if not d.empty:
            ax.scatter(d["UMAP1"], d["UMAP2"], s=0.28, color=LINEAGE_COLORS[label],
                       alpha=0.58 if label != "Unresolved" else 0.34,
                       linewidths=0, rasterized=True)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); umap_axes(ax)
    handles = [Line2D([0], [0], marker="o", linestyle="", markersize=3.0,
                      markerfacecolor=LINEAGE_COLORS[x], markeredgecolor="none",
                      label=LINEAGE_DISPLAY[x]) for x in annotation_order]
    ax.legend(handles=handles, ncol=3, loc="lower center", bbox_to_anchor=(0.50, -0.385),
              columnspacing=0.65, handletextpad=0.20, borderaxespad=0, fontsize=5.0)

    # g — broad-lineage canonical-marker dot plot
    ax = fig.add_subplot(gs[2, 8:13])
    panel_label(ax, "g", -0.09, 1.06)
    lineage_order_raw = ["Epithelial", "Endothelial", "Fibroblast", "Mural", "T_cell",
                         "NK", "B_cell", "Plasma", "Myeloid", "Mast"]
    lineage_order_display = ["Epithelial", "Endothelial", "Fibroblast", "Mural", "T cell",
                             "NK", "B cell", "Plasma", "Myeloid", "Mast"]
    marker_order, marker_groups = [], []
    for lineage in lineage_order_raw:
        genes = broad_markers.loc[broad_markers["marker_group"].eq(lineage),
                                  "gene_symbol"].drop_duplicates().tolist()
        marker_order.extend(genes); marker_groups.extend([lineage] * len(genes))
    for yi, lineage in enumerate(lineage_order_raw):
        d = broad_markers[broad_markers["cell_type"].eq(lineage)].set_index("gene_symbol").reindex(marker_order)
        for xi, row in enumerate(d.itertuples()):
            z = float(np.clip(row.scaled_mean, -DOT_VMAX, DOT_VMAX))
            ax.scatter(xi, yi, s=dot_size(float(row.detected_fraction)),
                       color=DOT_CMAP((z + DOT_VMAX) / (2 * DOT_VMAX)),
                       edgecolor="white", linewidth=0.25)
    ax.set_xlim(-0.7, len(marker_order) - 0.3); ax.set_ylim(9.6, -0.7)
    ax.set_xticks(range(len(marker_order)))
    ax.set_xticklabels(marker_order, fontsize=5.0, rotation=58, ha="right",
                       rotation_mode="anchor", fontstyle="italic")
    ax.set_yticks(range(10)); ax.set_yticklabels(lineage_order_display, fontsize=5.1)
    ax.tick_params(length=0, pad=1)
    for spine in ax.spines.values(): spine.set_visible(False)
    ax.grid(color="#EDF0F2", linewidth=0.32)
    for xx, group in enumerate(marker_groups):
        ax.add_patch(Rectangle((xx - 0.45, -0.62), 0.90, 0.14,
                               facecolor=LINEAGE_COLORS[group], edgecolor="none", clip_on=False))

    # h — independently re-embedded fibroblast states
    ax = fig.add_subplot(gs[3, 0:4])
    panel_label(ax, "h", -0.12, 1.06)
    for cluster in range(1, 7):
        d = fib[fib["cluster"].eq(cluster)]
        ax.scatter(d["UMAP1"], d["UMAP2"], s=0.90, color=FIB_COLORS[cluster],
                   alpha=0.62, linewidths=0, rasterized=True)
        x_med, y_med = d[["UMAP1", "UMAP2"]].median().to_numpy()
        whitened_text(ax, x_med, y_med, str(cluster), fontsize=5.8, fontweight="bold",
                      color=INK, ha="center", va="center")
    fib_xlim = padded_limits(fib["UMAP1"].to_numpy())
    fib_ylim = padded_limits(fib["UMAP2"].to_numpy(), fraction=0.085)
    ax.set_xlim(*fib_xlim); ax.set_ylim(*fib_ylim); umap_axes(ax)
    fib_handles = [Line2D([0], [0], marker="o", linestyle="", markersize=3.0,
                          markerfacecolor=FIB_COLORS[c], markeredgecolor="none",
                          label=f"{c}  {FIB_LABELS[c]}") for c in range(1, 7)]
    ax.legend(handles=fib_handles, ncol=2, loc="lower center", bbox_to_anchor=(0.50, -0.34),
              handletextpad=0.25, columnspacing=0.55, borderaxespad=0, fontsize=5.0)

    # i — PLAT expression on identical fibroblast coordinates
    ax = fig.add_subplot(gs[3, 4:8])
    panel_label(ax, "i", -0.12, 1.06)
    ax.scatter(fib["UMAP1"], fib["UMAP2"], s=0.52, color="#DADDE0", alpha=0.24,
               linewidths=0, rasterized=True)
    positive = fib[fib["expr_PLAT"] > 0].sort_values("expr_PLAT")
    ax.scatter(positive["UMAP1"], positive["UMAP2"], s=0.78,
               c=np.clip(positive["expr_PLAT"], 0, 5), cmap=PLAT_CMAP, norm=Normalize(0, 5),
               alpha=0.74, linewidths=0, rasterized=True)
    ax.set_xlim(*fib_xlim); ax.set_ylim(*fib_ylim); umap_axes(ax)
    ax.text(0.03, 0.98, "PLAT", transform=ax.transAxes, ha="left", va="top",
            fontsize=5.8, fontstyle="italic", fontweight="bold", color=PLAT)
    cax = ax.inset_axes([0.14, -0.23, 0.72, 0.035])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(0, 5), cmap=PLAT_CMAP),
                      cax=cax, orientation="horizontal")
    cb.set_ticks([0, 1, 2, 3, 4, 5]); cb.set_ticklabels(["0", "1", "2", "3", "4", "≥5"])
    cb.ax.tick_params(labelsize=5.0, length=1.5, pad=0.8)
    cb.set_label("log1p(CP10K)", fontsize=5.0, labelpad=0.6)
    cb.outline.set_linewidth(0.35)

    # j — state marker dot plot; CXCL14 and POSTN replace two redundant FeaturePlots
    ax = fig.add_subplot(gs[3, 8:13])
    panel_label(ax, "j", -0.09, 1.06)
    gene_order = ["CXCL14", "POSTN", "PLAT", "TCF21", "RGS5", "ACTA2", "SRGN", "CCL5",
                  "FAP", "KRT19", "PI16", "CFD", "CTHRC1", "COL1A1"]
    marker_groups = [1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 6, 6, 6]
    for yi, cluster in enumerate(range(1, 7)):
        d = fib_markers[fib_markers["cluster"].eq(cluster)].set_index("gene_symbol").reindex(gene_order)
        for xi, row in enumerate(d.itertuples()):
            z = float(np.clip(row.scaled_mean, -DOT_VMAX, DOT_VMAX))
            ax.scatter(xi, yi, s=dot_size(float(row.detected_fraction)),
                       color=DOT_CMAP((z + DOT_VMAX) / (2 * DOT_VMAX)),
                       edgecolor=PLAT if row.Index == "PLAT" else "white",
                       linewidth=0.55 if row.Index == "PLAT" else 0.30)
    ax.set_xlim(-0.7, len(gene_order) - 0.3); ax.set_ylim(5.6, -0.7)
    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, fontsize=5.0, rotation=58, ha="right",
                       rotation_mode="anchor", fontstyle="italic")
    ax.set_yticks(range(6)); ax.set_yticklabels([f"{c}  {FIB_LABELS[c]}" for c in range(1, 7)], fontsize=5.0)
    ax.tick_params(length=0, pad=1)
    for spine in ax.spines.values(): spine.set_visible(False)
    ax.grid(color="#EDF0F2", linewidth=0.34)
    for xx, group in enumerate(marker_groups):
        ax.add_patch(Rectangle((xx - 0.45, -0.62), 0.90, 0.14,
                               facecolor=FIB_COLORS[group], edgecolor="none", clip_on=False))
    size_handles = [Line2D([0], [0], marker="o", linestyle="", markersize=np.sqrt(dot_size(v)),
                           markerfacecolor="#9AA2A8", markeredgecolor="white",
                           markeredgewidth=0.35, label=f"{int(v * 100)}%")
                    for v in (0.25, 0.50, 0.75, 1.00)]
    size_leg = ax.legend(handles=size_handles, title="Detected", ncol=4, loc="upper left",
                         bbox_to_anchor=(-0.01, -0.235), borderaxespad=0,
                         columnspacing=0.45, handletextpad=0.10, fontsize=5.0,
                         title_fontsize=5.0)
    ax.add_artist(size_leg)
    cax_dot = ax.inset_axes([0.60, -0.205, 0.38, 0.035])
    cb_dot = fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(-DOT_VMAX, DOT_VMAX), cmap=DOT_CMAP),
                          cax=cax_dot, orientation="horizontal")
    cb_dot.set_ticks([-1.5, 0, 1.5]); cb_dot.ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f")); cb_dot.ax.tick_params(labelsize=5.0, length=1.3, pad=0.6)
    cb_dot.set_label("Scaled mean expression", fontsize=5.0, labelpad=0.4)
    cb_dot.outline.set_linewidth(0.35)

    # k — patient recurrence of the six GSE183904 fibroblast states
    k_sub = gs[4, 1:12].subgridspec(2, 1, height_ratios=[0.21, 1.0], hspace=0.04)
    ax_cov = fig.add_subplot(k_sub[0, 0])
    ax = fig.add_subplot(k_sub[1, 0])
    panel_label(ax_cov, "k", -0.055, 1.30)
    fraction_matrix = patient_state.pivot(index="cluster", columns="Patient_ID",
                                          values="within_patient_fraction").reindex(range(1, 7))
    hellinger = np.sqrt(fraction_matrix.T.to_numpy(dtype=float))
    order = leaves_list(linkage(hellinger, method="average", metric="euclidean", optimal_ordering=True))
    fraction_matrix = fraction_matrix.iloc[:, order]
    totals = (patient_state[["Patient_ID", "patient_fibroblast_cells"]]
              .drop_duplicates("Patient_ID").set_index("Patient_ID")
              .reindex(fraction_matrix.columns)["patient_fibroblast_cells"].to_numpy(dtype=float))
    bar_colors = np.where(totals < 20, "#B9BEC2", "#B77A87")
    ax_cov.bar(np.arange(len(totals)), np.log10(totals), width=0.82,
               color=bar_colors, linewidth=0)
    ax_cov.set_xlim(-0.5, len(totals) - 0.5)
    ax_cov.set_ylim(0, max(3.55, float(np.log10(totals).max()) * 1.05))
    ax_cov.set_yticks([1, 3]); ax_cov.set_yticklabels(["10", "10³"], fontsize=5.0)
    ax_cov.set_ylabel("Cells", fontsize=5.0, labelpad=1)
    ax_cov.set_xticks([]); ax_cov.tick_params(length=0, pad=1)
    for spine in ax_cov.spines.values(): spine.set_visible(False)
    im = ax.imshow(fraction_matrix.to_numpy(dtype=float), aspect="auto", interpolation="nearest",
                   cmap=PATIENT_CMAP, vmin=0, vmax=1)
    ax.set_yticks(range(6)); ax.set_yticklabels([f"{c}  {FIB_LABELS[c]}" for c in range(1, 7)])
    ax.set_xticks([]); ax.set_xlabel("Patients (n=26)", labelpad=2)
    # Keep the state labels clear of the narrow colour strip immediately to
    # the left of the heatmap body.
    ax.tick_params(length=0, pad=5.5)
    for spine in ax.spines.values(): spine.set_visible(False)
    for yv, color in enumerate([FIB_COLORS[c] for c in range(1, 7)]):
        ax.add_patch(Rectangle((-0.75, yv - 0.42), 0.16, 0.84,
                               facecolor=color, edgecolor="none", clip_on=False))
    ax.text(1.012, 0.5 / 6, "24/26", transform=ax.transAxes,
            fontsize=5.0, color=PLAT, ha="left", va="center")
    cax_k = ax.inset_axes([0.80, -0.23, 0.19, 0.050])
    cb_k = fig.colorbar(im, cax=cax_k, orientation="horizontal")
    cb_k.set_ticks([0, 0.5, 1]); cb_k.set_ticklabels(["0%", "50%", "100%"])
    cb_k.ax.tick_params(labelsize=5.0, length=1.4, pad=0.8)
    cb_k.set_label("State fraction", fontsize=5.0, labelpad=0.4)
    cb_k.outline.set_linewidth(0.35)

    save_exact_canvas(fig, "Figure2_single_cell_atlas")


if __name__ == "__main__":
    main()
