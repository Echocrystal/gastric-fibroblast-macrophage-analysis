# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) differential-regulation volcano -> assets/figures/volcano -> param inherit
# (b) recurrent differential-regulation interval plot -> cross-type inherit from Forest -> param inherit
# (c) centred receptor network -> assets/figures/other -> param inherit
# (d) seed-level directional comparison -> assets/figures/PairedBoxScatter -> param inherit
# (e) matched-control rank calibration -> assets/figures/BarComparison -> param inherit
# (f) target-state coverage + pharmacology pipeline -> assets/figures/BarComparison -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

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

# Academic Figure Skill Nature/Cell/Science Color Palette — COPY VERBATIM
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import networkx as nx

PROJECT = Path(r"I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026")
SOURCE = PROJECT / "figures/source_data/Figure7_receiver_screen_v3"
V4_SOURCE = PROJECT / "figures/source_data/Figure7_receiver_screen_v4"
VK = PROJECT / "analysis/results/virtual_knockout/GSE183904/C1QC_APOE_macrophage_receptor_screen_v2"
DIRECTION = VK / "ITGB2_directional_network_audit"
COMPLEX = VK / "ITGAM_vs_ITGB2_audit"
META = PROJECT / "analysis/processed/virtual_knockout/GSE183904/C1QC_APOE_macrophage_receptor_screen_v2/GSE183904_C1QC_APOE_sanity_metadata.tsv"
REVISION_DIR = PROJECT / "figures" / "ejc_revision_20260827"
REVISION_DIR.mkdir(parents=True, exist_ok=True)
OUT = REVISION_DIR / "Figure7_receiver_screen_v5_ejc_20260827"

PROGRAMS = [
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
    "GOBP_PHAGOCYTOSIS",
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
    "GOBP_CELL_MATRIX_ADHESION",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "ID_INFLAMMATORY_MONO",
]
PROGRAM_LABELS = {
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION": "Myeloid migration",
    "GOBP_PHAGOCYTOSIS": "Phagocytosis",
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX": "ECM degradation",
    "GOBP_CELL_MATRIX_ADHESION": "Matrix adhesion",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNF–NF-κB",
    "ID_INFLAMMATORY_MONO": "Inflammatory mono",
}
CORE_PROGRAMS = PROGRAMS[:4]
BLUE, ROSE = CATEGORICAL[0], CATEGORICAL[1]
LIGHT_BLUE = "#B3CBE2"
LIGHT_GREY = "#D0D0D0"


def panel_label(ax, label, x=-0.11, y=1.03):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left", color=BLACK, clip_on=False)


def quiet_axis(ax):
    ax.spines["left"].set_color("#000000")
    ax.spines["bottom"].set_color("#000000")
    ax.tick_params(color="#000000", labelcolor="#000000")
    ax.xaxis.label.set_color("#000000")
    ax.yaxis.label.set_color("#000000")


def load_data():
    diff_rep = pd.read_csv(SOURCE / "Figure7B_ITGB2_diffreg_seed2026072739.tsv", sep="\t")
    recurrent = pd.read_csv(SOURCE / "Figure7C_ITGB2_recurrent_diffreg_genes.tsv", sep="\t")
    nodes = pd.read_csv(SOURCE / "Figure7D_ITGB2_network_nodes.tsv", sep="\t")
    edges = pd.read_csv(SOURCE / "Figure7D_ITGB2_network_edges.tsv", sep="\t")
    by_seed = pd.read_csv(COMPLEX / "ITGAM_ITGB2_directional_by_seed.tsv", sep="\t")
    calibration = pd.read_csv(DIRECTION / "Figure7G_ITGB2_directional_connectivity_summary.tsv", sep="\t")
    calibration.insert(0, "representative_seed", 2026072739)
    all_nodes = pd.read_csv(DIRECTION / "Figure7G_ITGB2_directional_connectivity_all_nodes.tsv", sep="\t")
    meta = pd.read_csv(META, sep="\t")
    pharmacology = pd.read_csv(V4_SOURCE / "Figure7F_ITGB2_pharmacology.tsv", sep="\t")

    if len(diff_rep) != 775 or len(nodes) != 18 or len(edges) != 23:
        raise ValueError("Figure 7 v4 core source-data integrity check failed")
    if set(by_seed.target) != {"ITGB2", "ITGAM"} or len(by_seed) != 36:
        raise ValueError("ITGB2/ITGAM directional source-data integrity check failed")
    if set(calibration.program) != set(PROGRAMS) or len(calibration) != 6:
        raise ValueError("Directional calibration source-data integrity check failed")
    if all_nodes.ko_gene.nunique() != 21 or len(all_nodes) != 126:
        raise ValueError("Matched-node source-data integrity check failed")
    if len(meta) != 706 or meta.Patient_ID.nunique() != 26:
        raise ValueError("Donor-balanced target-state metadata integrity check failed")
    if len(pharmacology) != 4:
        raise ValueError("Pharmacology source-data integrity check failed")
    coverage = pd.DataFrame({
        "target": ["ITGB2 / CD18", "ITGAM / CD11b"],
        "detection_fraction": [meta.detect_ITGB2.mean(), meta.detect_ITGAM.mean()],
        "n_cells": [len(meta), len(meta)],
        "n_patients": [meta.Patient_ID.nunique(), meta.Patient_ID.nunique()],
    })
    V4_SOURCE.mkdir(parents=True, exist_ok=True)
    coverage.to_csv(V4_SOURCE / "Figure7F_target_state_coverage.tsv", sep="\t", index=False)
    by_seed.to_csv(V4_SOURCE / "Figure7D_ITGB2_ITGAM_directional_by_seed.tsv", sep="\t", index=False)
    calibration.to_csv(V4_SOURCE / "Figure7E_ITGB2_matched_rank_summary.tsv", sep="\t", index=False)
    return dict(diff=diff_rep, recurrent=recurrent, nodes=nodes, edges=edges,
                by_seed=by_seed, calibration=calibration, all_nodes=all_nodes,
                coverage=coverage, pharmacology=pharmacology)


def draw_a(ax, diff):
    panel_label(ax, "a", -0.12)
    q = diff.copy()
    q["neglogq"] = -np.log10(q["p.adj"].clip(lower=np.nextafter(0, 1)))
    q["class"] = "background"
    q.loc[(q.Z.abs() >= 2) & (q["p.adj"] < 0.01), "class"] = "strict hit"
    q.loc[q.gene == "ITGB2", "class"] = "perturbed target"
    q["plot_neglogq"] = q["neglogq"].clip(upper=28)
    colors = {"background": "#C9C9C9", "strict hit": ROSE, "perturbed target": BLACK}
    for cls, z in [("background", 1), ("strict hit", 3), ("perturbed target", 4)]:
        s = q[q["class"] == cls]
        ax.scatter(s.Z, s.plot_neglogq, s=8 if cls == "background" else 21,
                   color=colors[cls], alpha=0.55 if cls == "background" else 0.96,
                   edgecolor="none", zorder=z, rasterized=True)
    for x in (-2, 2):
        ax.axvline(x, color="#999999", lw=0.65, ls=(0, (3, 2)))
    ax.axhline(2, color="#999999", lw=0.65, ls=(0, (3, 2)))
    offsets = {"SLC7A11": (5, 3), "SLC11A1": (5, -9), "SPP1": (28, -2),
               "VCAN": (28, -15), "S100A9": (28, 12)}
    hits = q[q["class"] == "strict hit"].sort_values("neglogq", ascending=False)
    for _, r in hits.iterrows():
        ax.annotate(r.gene, (r.Z, r.plot_neglogq), xytext=offsets.get(r.gene, (4, 4)),
                    textcoords="offset points", fontsize=5.2, color=ROSE, ha="left", va="center",
                    fontstyle="italic",
                    arrowprops=dict(arrowstyle="-", lw=0.35, color="#C8A0A8", shrinkA=1.5, shrinkB=2))
    it = q[q.gene == "ITGB2"].iloc[0]
    ax.annotate("ITGB2", (it.Z, 28), xytext=(-22, -8), textcoords="offset points",
                fontsize=5.2, color=BLACK, ha="right", va="bottom", fontstyle="italic")
    ax.set_xlabel("Differential-regulation Z score")
    ax.set_ylabel(r"−log$_{10}$(adjusted $q$)")
    ax.set_xlim(min(-3.2, q.Z.min() - 0.2), max(3.0, q.Z.max() + 0.2))
    ax.set_ylim(0, 30)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=ROSE, ms=4, label=r"|Z| ≥ 2, $q$ < 0.01"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BLACK, ms=4, label="perturbed target"),
    ], loc="upper left", fontsize=5.3, handletextpad=0.3, borderaxespad=0.2)
    quiet_axis(ax)


def draw_b(ax, recurrent):
    panel_label(ax, "b", -0.16)
    q = recurrent[recurrent.q05_seeds == 3].sort_values(["median_Z", "gene"]).copy()
    y = np.arange(len(q))
    xerr = np.vstack([q.median_Z - q.min_Z, q.max_Z - q.median_Z])
    cols = [ROSE if z >= 2 else BLUE for z in q.median_Z]
    ax.errorbar(q.median_Z, y, xerr=xerr, fmt="none", ecolor="#8E8E8E",
                elinewidth=0.8, capsize=1.5, zorder=1)
    ax.scatter(q.median_Z, y, s=22, c=cols, edgecolor="white", linewidth=0.45, zorder=3)
    ax.axvline(2, color="#999999", lw=0.65, ls=(0, (3, 2)))
    ax.set_yticks(y)
    ax.set_yticklabels(q.gene)
    for tick in ax.get_yticklabels():
        tick.set_fontstyle("italic")
    ax.set_xlabel("Median Z score across 3 seeds")
    ax.set_xlim(1.55, max(2.48, q.max_Z.max() + 0.08))
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
    quiet_axis(ax)


def draw_c(ax, nodes, edges):
    panel_label(ax, "c", -0.06)
    G = nx.Graph()
    G.add_nodes_from(nodes.gene.tolist())
    for _, r in edges.iterrows():
        G.add_edge(r["from"], r["to"], weight=float(r.weight))
    pos = nx.spring_layout(G, seed=12, k=0.95, iterations=250, weight="weight")
    xy = np.array(list(pos.values()))
    xy = xy - xy.mean(0)
    xy = xy / max(abs(xy).max(), 1e-9) * 0.91
    pos = {g: tuple(xy[i]) for i, g in enumerate(pos)}
    for _, r in edges.iterrows():
        nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[(r["from"], r["to"])],
                               width=0.6 + 2.2 * float(r.weight), edge_color="#B7B7B7", alpha=0.62)
    node_map = nodes.set_index("gene")
    class_specs = [
        ("direct_network_neighbour", LIGHT_BLUE, 105),
        ("differential_regulation_hit", ROSE, 170),
        ("perturbed_target", BLACK, 410),
    ]
    for cls, color, size in class_specs:
        genes = [g for g in G.nodes if node_map.loc[g, "node_class"] == cls]
        nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=genes, node_color=color,
                               node_size=[size] * len(genes), edgecolors="white", linewidths=0.7)
    # Keep strict-hit labels outside the node cloud so the five gene symbols
    # remain individually legible at journal export size.
    hit_offsets = {"SLC7A11": (-28, 10), "SLC11A1": (18, -18), "SPP1": (2, 20),
                   "VCAN": (-22, -16), "S100A9": (22, 10)}
    hit_genes = set(nodes.loc[nodes.node_class == "differential_regulation_hit", "gene"])
    for g, (x, y0) in pos.items():
        if g in hit_genes:
            dx, dy = hit_offsets.get(g, (0, 10))
            ax.annotate(g, (x, y0), xytext=(dx, dy), textcoords="offset points",
                        ha="center", va="center", fontsize=5.1, color=ROSE, fontweight="bold",
                        fontstyle="italic",
                        arrowprops=dict(arrowstyle="-", lw=0.35, color="#C8A0A8"))
        else:
            ax.text(x, y0, g, ha="center", va="center", fontsize=5.0 if g != "ITGB2" else 5.8,
                    color="white" if g == "ITGB2" else BLACK,
                    fontweight="bold" if g == "ITGB2" else "normal", fontstyle="italic")
    ax.set_xlim(-1.02, 0.60)
    ax.set_ylim(-0.62, 0.57)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=ROSE, markeredgecolor="white", ms=5, label="representative-seed strict hit"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=LIGHT_BLUE, markeredgecolor="white", ms=5, label="direct WT neighbour"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BLACK, markeredgecolor="white", ms=5, label="ITGB2"),
    ], loc="lower center", bbox_to_anchor=(0.5, -0.155), ncol=3, fontsize=5.0,
       handletextpad=0.25, columnspacing=0.55)


def draw_d(ax, by_seed):
    panel_label(ax, "d", -0.14)
    order = PROGRAMS[::-1]
    specs = {"ITGB2": (ROSE, "o", -0.10), "ITGAM": (BLUE, "s", 0.10)}
    for yi, program in enumerate(order):
        for target, (color, marker, offset) in specs.items():
            q = by_seed[(by_seed.target == target) & (by_seed.program == program)].sort_values("seed")
            vals = q.percent_change_total_abs_strength.to_numpy()
            jitter = np.array([-0.035, 0.0, 0.035])
            ax.scatter(vals, yi + offset + jitter, s=18, marker=marker, facecolor="white",
                       edgecolor=color, linewidth=0.8, zorder=3)
            ax.scatter(np.median(vals), yi + offset, s=38, marker="D", color=color,
                       edgecolor="white", linewidth=0.5, zorder=4)
    ax.axvline(0, color="#777777", lw=0.65)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([PROGRAM_LABELS[p] for p in order])
    ax.set_xlabel("Change in programme network strength (%)")
    ax.set_xlim(-0.31, 0.01)
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
    ax.legend(handles=[
        Line2D([0], [0], marker="D", color="none", markerfacecolor=ROSE, markeredgecolor="white", ms=5, label="ITGB2"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=BLUE, markeredgecolor="white", ms=5, label="ITGAM"),
        Line2D([0], [0], marker="o", color="#777777", markerfacecolor="white", ms=4, lw=0, label="seed"),
    ], loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=3, fontsize=5.2,
       handletextpad=0.25, columnspacing=0.6, borderaxespad=0)
    quiet_axis(ax)


def draw_e(ax, calibration):
    panel_label(ax, "e", -0.12)
    order = PROGRAMS[::-1]
    q = calibration.set_index("program").loc[order].reset_index()
    q["rank"] = q.controls_at_least_as_negative + 1
    y = np.arange(len(q))
    ax.axvspan(0.5, 1.5, color="#F4EDCB", alpha=0.85, zorder=0)
    ax.axvline(1.5, color="#A49B70", lw=0.6, ls=(0, (3, 2)))
    ax.scatter(q["rank"], y, s=42, marker="D", color=ROSE, edgecolor="white", linewidth=0.6, zorder=3)
    for _, r in q.iterrows():
        ax.text(r["rank"] + 0.45, y[r.name], f'{int(r["rank"])}/21', ha="left", va="center",
                fontsize=5.4, color=BLACK)
    ax.set_yticks(y)
    ax.set_yticklabels([PROGRAM_LABELS[p] for p in order])
    ax.set_xlabel("Representative seed 2026072739\nrank among target + 20 matched nodes")
    ax.set_xlim(0.4, 22.2)
    ax.set_xticks([1, 5, 10, 15, 21])
    quiet_axis(ax)


def draw_f(ax, coverage, pharmacology):
    panel_label(ax, "f", -0.12)
    ax.axis("off")
    top = ax.inset_axes([0.02, 0.48, 0.96, 0.50])
    bottom = ax.inset_axes([0.02, 0.02, 0.96, 0.32])

    stages = {"Direct binding": 0, "Phase 1": 1, "Phase 2": 2, "Approval": 3}
    agents = ["Direct binder", "MLN-01", "Erlizumab", "MLNM-2201"]
    display = {
        "Unnamed direct-binding molecule": "Direct binder",
        "MLN-01": "MLN-01",
        "Erlizumab / rhuMAb CD18": "Erlizumab",
        "MLNM-2201 / LDP-01": "MLNM-2201",
    }
    p = pharmacology.copy()
    p["agent_short"] = p.agent.map(display)
    p["x"] = p.development_stage.map(stages)
    p["y"] = p.agent_short.map({a: i for i, a in enumerate(agents[::-1])})
    for _, r in p.iterrows():
        top.plot([0, r.x], [r.y, r.y], color=LIGHT_GREY, lw=1.2, zorder=1)
        top.scatter(r.x, r.y, s=35, color=ROSE if r.x >= 1 else BLUE,
                    edgecolor="white", linewidth=0.6, zorder=3)
    binding = p.loc[p.agent_short == "Direct binder"].iloc[0]
    binding_label = (
        r"$K_{\mathrm{d(app)}}$ "
        + f"{float(binding['value']):.2f} {binding['unit']}"
    )
    top.text(0.05, binding["y"] + 0.16, binding_label,
             fontsize=5.2, color=BLUE, ha="left", va="bottom")
    top.set_yticks(range(len(agents)))
    top.set_yticklabels(agents[::-1], fontsize=5.5)
    top.set_xticks(range(4))
    top.set_xticklabels(["Binding", "Phase 1", "Phase 2", "Approval"], fontsize=5.5)
    top.set_xlim(-0.08, 3.08)
    top.set_ylim(-0.45, 3.45)
    top.spines["left"].set_visible(False)
    top.tick_params(axis="y", length=0)
    quiet_axis(top)

    c = coverage.sort_values("detection_fraction")
    y = np.arange(len(c))
    colors = [BLUE if "ITGAM" in t else ROSE for t in c.target]
    bottom.hlines(y, 0, c.detection_fraction * 100, color=LIGHT_GREY, lw=2.0)
    bottom.scatter(c.detection_fraction * 100, y, s=43, color=colors,
                   edgecolor="white", linewidth=0.6, zorder=3)
    for yi, (_, r) in enumerate(c.iterrows()):
        bottom.text(r.detection_fraction * 100 + 2.2, yi, f"{100*r.detection_fraction:.2f}%",
                    va="center", ha="left", fontsize=5.5, color=BLACK)
    bottom.set_yticks(y)
    target_labels = [r"$\it{ITGB2}$ / CD18" if "ITGB2" in value else r"$\it{ITGAM}$ / CD11b"
                     for value in c.target]
    bottom.set_yticklabels(target_labels, fontsize=5.5)
    bottom.set_xlim(0, 100)
    bottom.set_xlabel(r"Detection in donor-balanced $\it{C1QC}^{+}\it{APOE}^{+}$ macrophages"
                      "\n706 cells, 26 patients", fontsize=5.8)
    bottom.spines["left"].set_visible(False)
    bottom.tick_params(axis="y", length=0)
    quiet_axis(bottom)


def validate_render_data(d):
    strict = ((d["diff"].Z.abs() >= 2) & (d["diff"]["p.adj"] < 0.01)).sum()
    if strict < 5:
        raise ValueError(f"Volcano signal too sparse: {strict} strict downstream hits")
    directional_range = d["by_seed"].percent_change_total_abs_strength.max() - d["by_seed"].percent_change_total_abs_strength.min()
    if directional_range <= 0.05 * abs(d["by_seed"].percent_change_total_abs_strength.mean()):
        raise ValueError("Directional programme panel lacks visual signal")
    if not (d["coverage"].detection_fraction.between(0, 1).all()):
        raise ValueError("Coverage fractions outside [0,1]")


def build_figure():
    d = load_data()
    validate_render_data(d)
    fig = plt.figure(figsize=(183 / 25.4, 195 / 25.4), facecolor="white")
    outer = fig.add_gridspec(
        3, 12, height_ratios=[1.10, 1.16, 1.05],
        left=0.140, right=0.985, top=0.975, bottom=0.075,
        wspace=1.35, hspace=0.42,
    )
    ax_a = fig.add_subplot(outer[0, 0:7])
    ax_b = fig.add_subplot(outer[0, 8:12])
    ax_c = fig.add_subplot(outer[1, 0:6])
    ax_d = fig.add_subplot(outer[1, 7:12])
    ax_e = fig.add_subplot(outer[2, 0:6])
    ax_f = fig.add_subplot(outer[2, 7:12])
    draw_a(ax_a, d["diff"])
    draw_b(ax_b, d["recurrent"])
    draw_c(ax_c, d["nodes"], d["edges"])
    draw_d(ax_d, d["by_seed"])
    draw_e(ax_e, d["calibration"])
    draw_f(ax_f, d["coverage"], d["pharmacology"])
    with mpl.rc_context({"savefig.bbox": None}):
        fig.savefig(f"{OUT}.pdf", bbox_inches=None, dpi=300)
        fig.savefig(f"{OUT}.png", bbox_inches=None, dpi=300)
        fig.savefig(f"{OUT}.svg", bbox_inches=None)
        fig.savefig(f"{OUT}.tiff", bbox_inches=None, dpi=600,
                    pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    return OUT


if __name__ == "__main__":
    print(build_figure())
