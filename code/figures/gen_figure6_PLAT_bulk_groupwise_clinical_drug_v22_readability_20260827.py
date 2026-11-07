# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) configuration scatter -> existing validated groupwise script -> param inherit
# (b) grouped boxplot -> existing validated groupwise script -> param inherit
# (c) grouped boxplot -> existing validated groupwise script -> param inherit
# (d) clinicopathology forest -> cross-type inherit from forest/bar assets -> param inherit
# (e) subtype-association lollipop -> cross-type inherit from forest/bar assets -> param inherit
# (f) stacked subtype composition -> existing validated contingency data -> param inherit
# (g) Kim2018 observed response -> existing validated groupwise script -> param inherit
# (h) continuous response models -> existing validated groupwise script -> param inherit
# (i) grouped violin -> existing validated GDSC2 violin implementation -> param inherit
# RULE: native-run assets were not available for these composite panels; the script preserves the
# validated statistical inputs and applies only the confirmed layout/typography adaptations.

from pathlib import Path
import hashlib
import json
import importlib.util
import math
import shutil
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 7,
    "axes.titlesize": 7,
    "axes.labelsize": 6.5,
    "xtick.labelsize": 5.2,
    "ytick.labelsize": 5.2,
    "legend.fontsize": 5.5,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.dpi": 300,
    "text.color": "#000000",
    "axes.labelcolor": "#000000",
    "axes.edgecolor": "#000000",
    "xtick.color": "#000000",
    "ytick.color": "#000000",
})

ROOT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
REVIEW = ROOT / "analysis" / "reports" / "Figure5_group_comparison_review_20260806"
OUT = ROOT / "figures" / "main" / "revised_20260827" / "Figure_6_readability_v22"
CANONICAL = ROOT / "figures" / "main" / "Figure6_PLAT_bulk_groupwise_clinical_drug"
LEGACY_CANONICAL = ROOT / "figures" / "main" / "Figure6_PLAT_bulk_groupwise_clinical_chemo"
CHEMO = ROOT / "analysis" / "results" / "chemo_response_projection_v1"
DRUG_DATA = ROOT / "analysis" / "drug_sensitivity" / "results_20260807"
SOURCE_V18 = ROOT / "figures" / "source_data" / "Figure6_v19_readability_20260827"
OUT.parent.mkdir(parents=True, exist_ok=True)
SOURCE_V18.mkdir(parents=True, exist_ok=True)

GROUPS = ["Dual-low", "Fibroblast-dominant", "Macrophage-dominant", "Dual-high"]
SHORT = ["DL", "FD", "MD", "DH"]
GROUP_COLORS = {
    "Dual-low": "#AAB7C4",
    "Fibroblast-dominant": "#CDACB4",
    "Macrophage-dominant": "#8FABDA",
    "Dual-high": "#9FB9A1",
}
GROUP_DISPLAY = {
    "Dual-low": "Fibroblast-low / Macrophage-low",
    "Fibroblast-dominant": "Fibroblast-high / Macrophage-low",
    "Macrophage-dominant": "Fibroblast-low / Macrophage-high",
    "Dual-high": "Fibroblast-high / Macrophage-high",
}
GROUP_DISPLAY_COMPACT = {
    "Dual-low": "F-low/M-low",
    "Fibroblast-dominant": "F-high/M-low",
    "Macrophage-dominant": "F-low/M-high",
    "Dual-high": "F-high/M-high",
}
SENDER = "#8F5F6B"
RECEIVER = "#456F9D"
ECM = "#6F8B70"
INK = "#000000"
GREY = "#9AA3AA"
LIGHT_GREY = "#E7EBEE"
RESPONDER = "#456F9D"

spec = importlib.util.spec_from_file_location(
    "fig5v2", ROOT / "figures" / "gen_figure5_final_v2.py"
)
fig5v2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fig5v2)


def read_review(name):
    return pd.read_csv(REVIEW / name, sep="\t")


def clean(ax, hide_left=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if hide_left:
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    ax.tick_params(length=2, width=0.55, colors=INK)
    ax.xaxis.label.set_color("#000000")
    ax.yaxis.label.set_color("#000000")
    if ax.spines["left"].get_visible():
        ax.spines["left"].set_color("#000000")
    if ax.spines["bottom"].get_visible():
        ax.spines["bottom"].set_color("#000000")


def panel_label(ax, label, x=-0.08, y=1.08):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=9,
            fontweight="bold", ha="left", va="top", color=INK, clip_on=False)


def format_probability(value):
    value = float(value)
    if value >= 0.01:
        return f"{value:.2f}"
    mantissa, exponent = f"{value:.2e}".split("e")
    return rf"{float(mantissa):.2f}$\times 10^{{{int(exponent)}}}$"


def load_score():
    return pd.read_csv(
        ROOT / "analysis" / "results" / "PLAT_axis" / "TCGA_STAD" /
        "frozen_state_projection" /
        "TCGA_STAD_PLAT_ecosystem_scores_clinical.tsv.gz", sep="\t"
    )


def draw_a(ax):
    d = load_score()
    d = d[d.ecotype4.isin(GROUPS)].copy()
    xcut = d.sender_primary_z.median()
    ycut = d.receiver_specific_z.median()
    xmin, xmax, ymin, ymax = -3.7, 3.4, -3.7, 3.4
    fx = (xcut - xmin) / (xmax - xmin)
    fy = (ycut - ymin) / (ymax - ymin)
    ax.axvspan(xmin, xcut, ymin=0, ymax=fy, color="#F2F4F6", zorder=0)
    ax.axvspan(xcut, xmax, ymin=0, ymax=fy, color="#F1E8EC", alpha=.55, zorder=0)
    ax.axvspan(xmin, xcut, ymin=fy, ymax=1, color="#E8EEF7", alpha=.6, zorder=0)
    ax.axvspan(xcut, xmax, ymin=fy, ymax=1, color="#E6F0E7", alpha=.7, zorder=0)
    for g in GROUPS:
        q = d[d.ecotype4.eq(g)]
        ax.scatter(q.sender_primary_z, q.receiver_specific_z, s=7, alpha=.52,
                   color=GROUP_COLORS[g], edgecolors="white", linewidth=.15,
                   rasterized=True)
    ax.axvline(xcut, color=SENDER, ls="--", lw=.85)
    ax.axhline(ycut, color=RECEIVER, ls="--", lw=.85)
    positions = {
        "Dual-low": (.03, .05, "left", "bottom"),
        "Fibroblast-dominant": (.97, .05, "right", "bottom"),
        "Macrophage-dominant": (.03, .94, "left", "top"),
        "Dual-high": (.97, .94, "right", "top"),
    }
    for g, (xx, yy, ha, va) in positions.items():
        name = {"Dual-low": "Fibroblast-low\nMacrophage-low",
                "Fibroblast-dominant": "Fibroblast-high\nMacrophage-low",
                "Macrophage-dominant": "Fibroblast-low\nMacrophage-high",
                "Dual-high": "Fibroblast-high\nMacrophage-high"}[g]
        ax.text(xx, yy, f"{name}\nn={int(d.ecotype4.eq(g).sum())}",
                transform=ax.transAxes, ha=ha, va=va, fontsize=5.4,
                fontweight="bold", color=INK,
                bbox=dict(facecolor="white", edgecolor="none", alpha=.82, pad=1.1))
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    ax.set_xticks([-3, 0, 3]); ax.set_yticks([-3, 0, 3])
    ax.set_xlabel(r"$\it{PLAT}$-rich fibroblast-state score (z)", color=SENDER, labelpad=1)
    ax.set_ylabel(r"$\it{C1QC}^{+}\it{APOE}^{+}$ macrophage-state score (z)", color=RECEIVER, labelpad=1)
    ax.set_box_aspect(1)
    clean(ax)


def load_mcp():
    p = (ROOT / "analysis" / "results" / "PLAT_axis" / "external_validation" /
         "immune_biomarker_extension_v1" / "four_cohort_MCPcounter_leaveout_sample_scores.tsv.gz")
    d = pd.read_csv(p, sep="\t")
    d = d[d.cohort.eq("TCGA-STAD")].pivot(index="sample_id", columns="population", values="score_z").reset_index()
    return d.merge(load_score()[["sample_submitter_id", "ecotype4"]],
                   left_on="sample_id", right_on="sample_submitter_id", how="inner")


def load_checkpoint():
    p = (ROOT / "analysis" / "results" / "PLAT_axis" / "external_validation" /
         "immune_biomarker_extension_v1" / "four_cohort_checkpoint_sample_expression.tsv.gz")
    d = pd.read_csv(p, sep="\t")
    return d[d.cohort.eq("TCGA-STAD")].pivot(
        index=["sample_id", "ecotype4"], columns="display_label", values="expression_z"
    ).reset_index()


def star(q):
    if q < 1e-4: return "****"
    if q < 1e-3: return "***"
    if q < 1e-2: return "**"
    if q < 5e-2: return "*"
    return "ns"


def bracket_label(q):
    """Reference-style label: stars for significant pairs and lowercase ns otherwise."""
    return star(q) if q < 0.05 else "ns"


def compact_letters(pairwise, endpoint):
    """Compact-letter display from pairwise Dunn BH q values.

    Groups sharing a letter are not significantly different; groups with no
    shared letter have pairwise BH q < 0.05. Maximal nonsignificant cliques
    are used so the annotations remain compact for four groups.
    """
    p = pairwise[pairwise.endpoint.eq(endpoint)]
    nonsig = {(g, g) for g in GROUPS}
    for row in p.itertuples(index=False):
        if float(row.pairwise_FDR_BH) >= 0.05:
            nonsig.add((row.group1, row.group2)); nonsig.add((row.group2, row.group1))
    cliques = []
    for size in range(1, len(GROUPS) + 1):
        for subset in combinations(GROUPS, size):
            if all((a, b) in nonsig for a, b in combinations(subset, 2)):
                cliques.append(set(subset))
    maximal = [c for c in cliques if not any(c < d for d in cliques)]
    letters = {g: [] for g in GROUPS}
    for idx, clique in enumerate(maximal):
        letter = chr(ord("a") + idx)
        for g in clique:
            letters[g].append(letter)
    return letters


def draw_bracket(ax, x1, x2, y, height, label, fontsize=3.6):
    """Reference-style significance bracket over two grouped boxes."""
    ax.plot([x1, x1, x2, x2], [y, y + height, y + height, y],
            color=INK, lw=.45, clip_on=False, solid_capstyle="butt")
    # Asterisks have a lower visible glyph edge than lowercase letters and
    # need slightly more clearance to avoid touching the bracket cap.
    text_gap = .0008 if label == "ns" else .012
    ax.text((x1 + x2) / 2, y + height + text_gap, label, ha="center", va="bottom",
            fontsize=fontsize, color=INK, fontweight="bold", clip_on=False)


def draw_group_strip(ax, data, endpoints, summary_name, pairwise_name, ylabel, labels):
    """Single complete grouped boxplot; endpoint names are x-axis labels, not subtitles."""
    tests = pd.read_csv(REVIEW / summary_name.replace("group_summary", "global_tests"), sep="\t").set_index("endpoint")
    pairwise = pd.read_csv(REVIEW / pairwise_name, sep="\t")
    # Match the supplied reference figure: show each of the six pairwise
    # comparisons as a separate horizontal bracket at a different height.
    # The bracket labels are pairwise Dunn BH q values, including NS. pairs.
    all_pairs = list(combinations(GROUPS, 2))
    centers = np.arange(len(endpoints), dtype=float)
    offsets = np.array([-.28, -.093, .093, .28])
    width = .16
    allvalues = []
    displayed_ymins = []
    displayed_ymaxs = []
    for i, ep in enumerate(endpoints):
        allv = []
        for j, g in enumerate(GROUPS):
            v = pd.to_numeric(data.loc[data.ecotype4.eq(g), ep], errors="coerce").dropna().to_numpy(float)
            allv.append(v); allvalues.append(v)
            # Individual jitter points are intentionally omitted in v9. The
            # sample-level values remain in the boxplot/statistical input;
            # removing the tiny points improves legibility of the six
            # reference-style significance brackets.
        bp = ax.boxplot(allv, positions=centers[i] + offsets, widths=width,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color=INK, lw=.7),
                        whiskerprops=dict(color=GREY, lw=.55),
                        capprops=dict(color=GREY, lw=.55),
                        boxprops=dict(lw=.5))
        for patch, g in zip(bp["boxes"], GROUPS):
            patch.set_facecolor(GROUP_COLORS[g]); patch.set_alpha(.85); patch.set_edgecolor("white")
        # Use the actual whisker extrema drawn by matplotlib. The previous
        # implementation used max(raw values), which includes hidden outliers
        # and incorrectly pushed the brackets far above the visible boxes.
        whisker_y = [float(y) for whisker in bp["whiskers"] for y in whisker.get_ydata()]
        displayed_ymins.append(min(whisker_y))
        displayed_ymaxs.append(max(whisker_y))
        top = max(whisker_y)
        # Shorter brackets are placed lower and longer brackets higher, as in
        # the reference boxplot. This keeps vertical crossings minimal.
        pair_order = sorted(all_pairs,
                            key=lambda pair: GROUPS.index(pair[1]) - GROUPS.index(pair[0]))
        endpoint_range = max(max(whisker_y) - min(whisker_y), 1.0)
        # With the jitter layer removed, start the annotation closer to the
        # whisker tops. Keep a generous but shorter ladder so the six labels
        # remain separated without recreating the previous collisions.
        bracket_step = max(endpoint_range * .075, .30)
        bracket_height = max(endpoint_range * .022, .065)
        for bracket_level, (g1, g2) in enumerate(pair_order):
            row = pairwise[((pairwise.endpoint == ep) &
                            (pairwise.group1 == g1) & (pairwise.group2 == g2)) |
                           ((pairwise.endpoint == ep) &
                            (pairwise.group1 == g2) & (pairwise.group2 == g1))]
            if len(row) == 0:
                continue
            q_pair = float(row.iloc[0].pairwise_FDR_BH)
            x1 = centers[i] + offsets[GROUPS.index(g1)]
            x2 = centers[i] + offsets[GROUPS.index(g2)]
            draw_bracket(ax, x1, x2, top + .045 + bracket_level * bracket_step,
                         bracket_height, bracket_label(q_pair))
    # Match matplotlib's visible boxplot extent rather than hidden outliers.
    ymin = min(displayed_ymins); ymax = max(displayed_ymaxs)
    # Six reference-style pairwise brackets need explicit headroom. The
    # endpoint-aware ladder above determines the actual top; this conservative
    # padding prevents the highest `NS.`/star label from being clipped.
    pad = max((ymax - ymin) * .30, .1)
    ax.set_xlim(-.65, centers[-1] + .65); ax.set_ylim(ymin - pad * .2, ymax + pad)
    ax.set_xticks(centers); ax.set_xticklabels(labels, rotation=42, ha="right", fontsize=5.15)
    ax.set_ylabel(ylabel, labelpad=1)
    clean(ax)


def draw_d_compact(ax):
    """Compact two-tier CAMERA pathway comparison; retains pathway-level estimand."""
    camera = pd.read_csv(ROOT / "figures" / "source_data" / "Figure5_final_v2" /
                         "Figure5c_CAMERA_focal_pathways.tsv", sep="\t")
    order = ["Complement & coagulation", "Complement activation", "Positive chemotaxis",
             "Myeloid migration", "Macrophage migration", "Phagocytosis regulation",
             "ECM assembly", "ECM organization"]
    camera = camera[camera.label.isin(order)].copy()
    p = camera[camera.model_tier.eq("purity_adjusted")].set_index("label").reindex(order)
    f = camera[camera.model_tier.eq("lineage_subtype_adjusted")].set_index("label").reindex(order)
    y = np.arange(len(order))[::-1]
    x1 = -np.log10(p.FDR_BH.to_numpy(float)); x2 = -np.log10(f.FDR_BH.to_numpy(float))
    for yy, a, b in zip(y, x1, x2):
        ax.plot([a, b], [yy, yy], color="#C9CED2", lw=.65, zorder=1)
    ax.scatter(x1, y + .10, s=14, facecolor="white", edgecolor=SENDER, lw=.7, label="Purity")
    sig = f.FDR_BH.to_numpy(float) < .05
    ax.scatter(x2, y - .10, s=15, marker="D", facecolor=np.where(sig, RECEIVER, "white"),
               edgecolor=RECEIVER, lw=.7, label="Lineage + subtype")
    ax.axvline(-math.log10(.05), color="#777777", lw=.55, ls="--")
    display = ["Complement/coag.", "Complement act.", "Chemotaxis", "Myeloid migr.",
               "Macrophage migr.", "Phagocytosis", "ECM assembly", "ECM organization"]
    ax.set_yticks(y); ax.set_yticklabels(display, fontsize=4.4)
    ax.tick_params(axis="y", pad=2)
    # The pathway estimand is defined in the caption; suppress an x-axis title
    # here so it cannot collide with the checkpoint boxplot in panel c.
    ax.set_xlim(0, 7.35); ax.set_xlabel("")
    # Move the two-tier model legend into the whitespace above panel d and
    # increase separation so the marker meanings are readable at print size.
    ax.legend(loc="lower center", bbox_to_anchor=(0.50, 1.035), ncol=2,
              fontsize=4.3, handletextpad=.55, columnspacing=1.35,
              labelspacing=.45, borderaxespad=0, frameon=False)
    clean(ax)


def draw_e_compact(ax):
    meta = pd.read_csv(ROOT / "figures" / "source_data" / "Figure5_final_v2" /
                       "Figure5e_clinicopathology_meta_source.tsv", sep="\t")
    meta = meta[(meta.model_tier == "univariable") &
                (meta.analysis_set == "all_available_including_TCGA")].copy()
    order = ["lauren_diffuse", "stage", "t_stage", "n_stage", "male", "age10"]
    meta = meta.set_index("feature").reindex(order).reset_index()
    labels = ["Diffuse vs intestinal", "Pathological stage / level", "T stage / level",
              "N stage / level", "Male vs female", "Age / 10 years"]
    y = np.arange(len(labels))[::-1]
    ax.axvline(0, color="#777777", lw=.55, ls="--")
    for i, row in meta.iterrows():
        yy = y[i]; sig = row.p_BH_within_analysis_tier < .05
        ax.plot([row.ci_low, row.ci_high], [yy, yy], color=GREY, lw=.9)
        ax.scatter(row.beta_random_DL, yy, s=21, facecolor=ECM if sig else "white",
                   edgecolor=ECM, lw=.7, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=5.0)
    ax.set_xlim(-.26, 1.02); ax.set_xlabel("Difference in ecosystem score (95% CI)", fontsize=5.6, labelpad=1)
    clean(ax)


def draw_f_compact(ax):
    """State-specific overall-survival forest from a traceable source TSV."""
    src = (ROOT / "figures" / "source_data" / "Figure6_v15_state_specific_OS" /
           "Figure6f_state_specific_univariable_OS.tsv")
    d = pd.read_csv(src, sep="\t")
    cohort_order = ["TCGA-STAD", "GSE62254", "GSE84437", "GSE15459", "External pooled"]
    d = d[d.cohort.isin(cohort_order)].copy()
    ybase = np.arange(len(cohort_order))[::-1]
    offsets = {"sender_primary_z": 0.105, "receiver_specific_z": -0.105}
    colors = {"sender_primary_z": SENDER, "receiver_specific_z": RECEIVER}
    markers = {"sender_primary_z": "o", "receiver_specific_z": "D"}
    labels = {
        "sender_primary_z": r"$\it{PLAT}$-rich fibroblast state",
        "receiver_specific_z": r"$\it{C1QC}^{+}\it{APOE}^{+}$ macrophage state",
    }
    ax.axvline(1, color="#777777", lw=.55, ls="--", zorder=0)
    for ci, cohort in enumerate(cohort_order):
        yy = ybase[ci]
        for exposure in ["sender_primary_z", "receiver_specific_z"]:
            row = d[(d.cohort == cohort) & (d.exposure == exposure)]
            if len(row) != 1:
                continue
            row = row.iloc[0]
            y = yy + offsets[exposure]
            ax.plot([row.ci_low, row.ci_high], [y, y], color=colors[exposure], lw=.85,
                    solid_capstyle="round", zorder=1)
            ax.scatter(row.HR_per_SD, y, s=19, marker=markers[exposure],
                       facecolor=colors[exposure], edgecolor="white", lw=.45, zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(.76, 1.85)
    ax.set_xticks([0.8, 1.2, 1.6])
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
    ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.tick_params(axis="x", labelsize=4.4)
    ax.set_yticks(ybase)
    ax.set_yticklabels(cohort_order, fontsize=4.7)
    ax.set_xlabel("Overall-survival HR per 1 SD (95% CI)", fontsize=5.35, labelpad=1)
    handles = [
        mpl.lines.Line2D([0], [0], marker="o", color=colors["sender_primary_z"],
                         markerfacecolor=colors["sender_primary_z"], markeredgecolor="white",
                         lw=.8, markersize=3.4, label=labels["sender_primary_z"]),
        mpl.lines.Line2D([0], [0], marker="D", color=colors["receiver_specific_z"],
                         markerfacecolor=colors["receiver_specific_z"], markeredgecolor="white",
                         lw=.8, markersize=3.1, label=labels["receiver_specific_z"]),
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(.50, 1.02),
              ncol=2, fontsize=4.0, handletextpad=.35, columnspacing=.8,
              borderaxespad=0, frameon=False)
    clean(ax, hide_left=True)


def draw_d_km(fig, spec):
    """Four-configuration TCGA-STAD KM curve matching panel-a quadrants."""
    ax = fig.add_subplot(spec)
    src = (ROOT / "figures" / "source_data" / "Figure6_v17_TCGA_four_group_KM" /
           "Figure6d_TCGA_STAD_four_group_KM_curves.tsv")
    stats_src = (ROOT / "figures" / "source_data" / "Figure6_v17_TCGA_four_group_KM" /
                 "Figure6d_TCGA_STAD_four_group_KM_statistics.tsv")
    curves = pd.read_csv(src, sep="\t")
    stats = pd.read_csv(stats_src, sep="\t").set_index("configuration")
    for group in GROUPS:
        v = curves[curves.configuration.eq(group)].sort_values("time_months")
        n = int(stats.loc[group, "n"])
        ax.step(v.time_months, v.survival, where="post", color=GROUP_COLORS[group],
                lw=1.15, label=f"{GROUP_DISPLAY_COMPACT[group]} (n={n})")
    global_p = float(stats.global_logrank_p.iloc[0])
    ax.text(.98, .97, f"Global log-rank P={global_p:.3f}", transform=ax.transAxes,
            ha="right", va="top", fontsize=4.7, color=INK)
    ax.set_xlim(0, 72); ax.set_ylim(0, 1.03)
    ax.set_xticks([0, 24, 48, 72]); ax.set_yticks([0, .5, 1.0])
    ax.set_xlabel("Time (months)", fontsize=5.4, labelpad=1)
    ax.set_ylabel("Overall survival probability", fontsize=5.4, labelpad=1)
    ax.legend(loc="lower left", bbox_to_anchor=(.02, .02), ncol=2,
              fontsize=4.15, handlelength=1.45, handletextpad=.35,
              columnspacing=.8, borderaxespad=0, frameon=False,
              labelspacing=.28)
    panel_label(ax, "d", x=-.14, y=1.08)
    clean(ax)


def draw_g_unified(ax):
    kim = pd.read_csv(ROOT / "figures" / "source_data" / "Figure5_final_v2" /
                      "Figure5g_Kim2018_patient_data.tsv", sep="\t")
    groups = pd.read_csv(ROOT / "figures" / "source_data" / "Figure5_final_v2" /
                         "Figure5g_Kim2018_group_summary.tsv", sep="\t").set_index("ecotype4")
    xcut = kim.sender_primary_z.median(); ycut = kim.receiver_specific_z.median()
    ax.axvline(xcut, color=SENDER, lw=.7, ls="--"); ax.axhline(ycut, color=RECEIVER, lw=.7, ls="--")
    for g in GROUPS:
        q = kim[kim.ecotype4.eq(g)]
        resp = q.response_binary.astype(bool)
        ax.scatter(q.loc[~resp, "sender_primary_z"], q.loc[~resp, "receiver_specific_z"],
                   s=21, facecolor="white", edgecolor=GROUP_COLORS[g], lw=.8)
        ax.scatter(q.loc[resp, "sender_primary_z"], q.loc[resp, "receiver_specific_z"],
                   s=25, facecolor=GROUP_COLORS[g], edgecolor="white", lw=.5)
    positions = {"Dual-low": (.02, .03, "left", "bottom"), "Fibroblast-dominant": (.98, .03, "right", "bottom"),
                 "Macrophage-dominant": (.02, .99, "left", "top"), "Dual-high": (.98, .99, "right", "top")}
    for g, row in groups.iterrows():
        xx, yy, ha, va = positions[g]
        ax.text(xx, yy, f"{GROUP_DISPLAY_COMPACT[g]}\n{int(row.responders)}/{int(row.n)} responders", transform=ax.transAxes,
                fontsize=4.4, fontweight="bold", ha=ha, va=va,
                bbox=dict(facecolor="white", edgecolor="none", alpha=.88, pad=.8))
    ax.set_xlabel(r"$\it{PLAT}$-rich fibroblast-state score (z)", fontsize=5.8, color=SENDER, labelpad=1)
    ax.set_ylabel(r"$\it{C1QC}^{+}\it{APOE}^{+}$ macrophage-state score (z)", fontsize=5.8, color=RECEIVER, labelpad=1)
    # Add deliberate coordinate padding so quadrant labels sit in the empty
    # corners rather than over the observations.
    ax.set_xlim(-3.18, 1.72)
    ax.set_ylim(-2.70, 2.35)
    ax.set_xticks([-3, -1, 1])
    ax.set_yticks([-2, 0, 1, 2])
    clean(ax)


def p_stars(p):
    """Compact significance notation for the narrow response panel."""
    p = float(p)
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < .05:
        return "*"
    return ""


def draw_h_stars(ax):
    """Compact response forest: retain OR, replace numeric P with stars."""
    models = pd.read_csv(ROOT / "figures" / "source_data" / "Figure5_final_v2" /
                         "Figure5h_Kim2018_continuous_models.tsv", sep="\t")
    labels = ["Average ecosystem", r"$\it{PLAT}$-rich fibroblast" "\nstate",
              r"$\it{C1QC}^{+}\it{APOE}^{+}$ macrophage" "\nstate",
              "Fibroblast x macrophage\nstate"]
    y = np.arange(4)[::-1]
    colors = [ECM, SENDER, RECEIVER, "#8C7D96"]
    ax.axvline(1, color="#777777", linewidth=0.6, linestyle="--")
    for i, row in models.reset_index(drop=True).iterrows():
        ax.plot([row["low"], row["high"]], [y[i], y[i]], color=GREY, linewidth=0.9)
        ax.scatter(row["OR"], y[i], s=24, facecolor=colors[i], edgecolor="white",
                   linewidth=0.5, zorder=3)
        stars = p_stars(row["p"])
        suffix = f" {stars}" if stars else ""
        ax.text(0.985, y[i], f"OR={row['OR']:.2f}{suffix}",
                transform=ax.get_yaxis_transform(), fontsize=4.85,
                ha="right", va="center", clip_on=True)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=5.0)
    ax.set_xscale("log"); ax.set_xlim(0.08, 16.0)
    ax.set_xlabel("Odds ratio for observed response\n(95% CI)", fontsize=6.0, labelpad=1)
    ax.grid(axis="x", color="#EEF1F3", linewidth=0.4)
    clean(ax, hide_left=True)



def draw_d_clinicopathology(fig, spec):
    """Cross-cohort clinicopathology meta-analysis on a common SD scale."""
    ax = fig.add_subplot(spec)
    meta = pd.read_csv(
        ROOT / "figures" / "source_data" / "Figure5_final_v2" /
        "Figure5e_clinicopathology_meta_source.tsv", sep="\t"
    )
    meta = meta[(meta.model_tier == "univariable") &
                (meta.analysis_set == "all_available_including_TCGA")].copy()
    order = ["lauren_diffuse", "stage", "t_stage", "n_stage"]
    meta = meta.set_index("feature").reindex(order).reset_index()
    meta.to_csv(SOURCE_V18 / "Figure6d_clinicopathology_meta.tsv", sep="\t", index=False)

    labels = [
        "Diffuse vs intestinal",
        "Pathological stage",
        "T stage / level",
        "N stage / level",
    ]
    y = np.arange(len(labels))[::-1]
    ax.axvline(0, color="#777777", lw=.55, ls="--", zorder=0)
    for i, row in meta.iterrows():
        yy = y[i]
        ax.plot([row.ci_low, row.ci_high], [yy, yy], color="#9AA3AA",
                lw=1.05, solid_capstyle="round", zorder=1)
        ax.scatter(row.beta_random_DL, yy, s=27, facecolor=ECM,
                   edgecolor="white", lw=.55, zorder=3)
        q = float(row.p_BH_within_analysis_tier)
        qtxt = rf"$q$ < 1.00$\times 10^{{-3}}$" if q < .001 else rf"$q$ = {format_probability(q)}"
        text_y = yy + (0.13 if i == 0 else 0.0)
        ax.text(.985, text_y, qtxt, transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=4.9, color="#000000")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5.0)
    ax.set_xlim(-.08, 1.42)
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
    ax.set_xlabel("Change in standardized ecosystem score\n(95% CI)",
                  fontsize=5.6, labelpad=1)
    panel_label(ax, "d", x=-.20, y=1.08)
    clean(ax)
    # The long x-axis title otherwise enters the checkpoint panel at print
    # size. Shift the forest right and reduce its width while retaining all CIs.
    pos = ax.get_position()
    ax.set_position([pos.x0 + .016, pos.y0, pos.width * .80, pos.height])


def draw_e_subtype_associations(ax):
    """Configuration associations with established gastric-cancer classifications."""
    src = (ROOT / "analysis" / "results" / "PLAT_axis" /
           "external_validation" / "cohort_specific_clinicopathology_v1" /
           "cohort_specific_configuration_tests.tsv")
    d = pd.read_csv(src, sep="\t")
    wanted = [
        "tcga2014_subtype", "who_histology", "acrg_subtype",
        "mlh1_ihc", "ebv_ish", "lei_molecular_subtype",
    ]
    d = d.set_index("variable_id").reindex(wanted).reset_index()
    d.to_csv(SOURCE_V18 / "Figure6e_subtype_associations.tsv", sep="\t", index=False)

    display = [
        "TCGA 2014\nsubtype", "WHO histology", "ACRG subtype",
        "MLH1-IHC", "EBV-ISH", "Lei subtype",
    ]
    cohort_colors = {
        "TCGA-STAD": SENDER,
        "GSE62254": RECEIVER,
        "GSE15459": ECM,
    }
    y = np.arange(len(display))[::-1]
    for yy, row in zip(y, d.itertuples(index=False)):
        col = cohort_colors.get(row.cohort, GREY)
        ax.plot([0, row.effect_size], [yy, yy], color="#D4D9DD", lw=.8, zorder=1)
        ax.scatter(row.effect_size, yy, s=23, facecolor=col, edgecolor="white",
                   lw=.5, zorder=3)
        ax.text(.985, yy, star(row.p_BH_all_specific_tests),
                transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=4.8, fontweight="bold", color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels(display, fontsize=5.35)
    ax.set_xlim(0, .405)
    ax.set_xticks([0, .2, .4])
    ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
    ax.set_xlabel("Association with four configurations (Cramer's V)",
                  fontsize=5.35, labelpad=1)
    clean(ax)


def draw_f_tcga_subtype_composition(ax):
    """TCGA molecular-subtype composition across the four configurations."""
    p = REVIEW / "TCGA_configuration_by_official_subtype.tsv"
    counts = pd.read_csv(p, sep="\t").set_index("ecotype4").reindex(GROUPS)
    props = counts.div(counts.sum(axis=1), axis=0) * 100
    props.reset_index().to_csv(
        SOURCE_V18 / "Figure6f_TCGA_subtype_composition_percent.tsv",
        sep="\t", index=False
    )
    subtype_colors = {
        "CIN": "#AAB7C4",
        "EBV": "#8FABDA",
        "GS": "#CDACB4",
        "MSI": "#6F8B70",
        "Unclassified": "#E7EBEE",
    }
    ypos = np.arange(len(GROUPS))[::-1]
    left = np.zeros(len(GROUPS))
    for subtype in ["CIN", "EBV", "GS", "MSI", "Unclassified"]:
        vals = props[subtype].to_numpy(float)
        bars = ax.barh(ypos, vals, left=left, height=.62,
                       color=subtype_colors[subtype], edgecolor="white",
                       linewidth=.35, label=subtype)
        for yy, lft, val in zip(ypos, left, vals):
            if val >= 14:
                ax.text(lft + val / 2, yy, f"{val:.0f}%", ha="center",
                        va="center", fontsize=4.5, color=INK)
        left += vals
    ax.set_yticks(ypos)
    ax.set_yticklabels([GROUP_DISPLAY_COMPACT[g] for g in GROUPS], fontsize=5.0)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel("TCGA molecular-subtype composition (%)", fontsize=5.35, labelpad=1)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.035), ncol=3,
              fontsize=5.0, handlelength=.80, handletextpad=.28,
              columnspacing=.60, labelspacing=.25, borderaxespad=0, frameon=False)
    clean(ax)


def draw_i_gdsc2_violin(fig, spec):
    """Exploratory GDSC2 gastric-cell-line AUC distributions."""
    response = pd.read_csv(
        DRUG_DATA / "GDSC2_selected_drug_response_long.tsv", sep="\t"
    )
    tests = pd.read_csv(
        DRUG_DATA / "GDSC2_selected_drug_group_tests.tsv", sep="\t"
    )
    pairwise = pd.read_csv(
        DRUG_DATA / "GDSC2_selected_drug_pairwise_tests.tsv", sep="\t"
    )
    drugs = ["5-Fluorouracil", "Irinotecan", "Cisplatin", "Oxaliplatin"]
    selected = response[response.DRUG_NAME.isin(drugs)].dropna(subset=["AUC"]).copy()
    selected.to_csv(SOURCE_V18 / "Figure6i_GDSC2_drug_AUC.tsv", sep="\t", index=False)
    tests[tests.drug.isin(drugs)].to_csv(
        SOURCE_V18 / "Figure6i_GDSC2_drug_tests.tsv", sep="\t", index=False
    )
    tests.to_csv(
        SOURCE_V18 / "Figure6i_GDSC2_nine_drug_BH_family.tsv", sep="\t", index=False
    )
    pairwise[pairwise.drug.isin(drugs)].to_csv(
        SOURCE_V18 / "Figure6i_GDSC2_drug_pairwise_tests.tsv", sep="\t", index=False
    )

    frame = fig.add_subplot(spec, frameon=False)
    frame.set_xticks([])
    frame.set_yticks([])
    for spine in frame.spines.values():
        spine.set_visible(False)
    panel_label(frame, "i", x=-.13, y=1.045)

    grid = spec.subgridspec(2, 2, hspace=.68, wspace=.36)
    short_names = {
        "5-Fluorouracil": "5-FU",
        "Irinotecan": "Irinotecan",
        "Cisplatin": "Cisplatin",
        "Oxaliplatin": "Oxaliplatin",
    }
    compact = {
        "Dual-low": "F-lo/M-lo",
        "Fibroblast-dominant": "F-hi/M-lo",
        "Macrophage-dominant": "F-lo/M-hi",
        "Dual-high": "F-hi/M-hi",
    }
    for idx, drug in enumerate(drugs):
        ax = fig.add_subplot(grid[idx // 2, idx % 2])
        g = selected[selected.DRUG_NAME.eq(drug)]
        values = [
            g.loc[g.configuration.eq(group), "AUC"].to_numpy(float)
            for group in GROUPS
        ]
        vp = ax.violinplot(
            values, positions=np.arange(4), widths=.78, points=120,
            bw_method=.38, showmeans=False, showextrema=False, showmedians=False
        )
        for body, group in zip(vp["bodies"], GROUPS):
            body.set_facecolor(GROUP_COLORS[group])
            body.set_edgecolor("#3F4850")
            body.set_linewidth(.65)
            body.set_alpha(.78)
        bp = ax.boxplot(
            values, positions=np.arange(4), widths=.20, patch_artist=True,
            showfliers=False, whis=(0, 100),
            medianprops=dict(color=INK, linewidth=1.20),
            whiskerprops=dict(color=INK, linewidth=.80),
            capprops=dict(color=INK, linewidth=.80),
            boxprops=dict(facecolor="white", edgecolor=INK,
                          linewidth=.85, alpha=.92),
        )
        labels = [compact[group] for group in GROUPS]
        ax.set_title(short_names[drug], loc="left", fontsize=5.7,
                     fontweight="bold", pad=1.5)
        # Display the same complete six-pair comparison structure used in
        # panels b and c. The correction family is the prespecified set of
        # drug-by-pair tests in the source table.
        drug_pairs = pairwise[pairwise.drug.eq(drug)]
        pair_order = sorted(combinations(GROUPS, 2),
                            key=lambda pair: GROUPS.index(pair[1]) - GROUPS.index(pair[0]))
        bracket_base = 1.000
        bracket_step = .020
        for level, (g1, g2) in enumerate(pair_order):
            row = drug_pairs[
                ((drug_pairs.group_a == g1) & (drug_pairs.group_b == g2)) |
                ((drug_pairs.group_a == g2) & (drug_pairs.group_b == g1))
            ]
            if row.empty:
                continue
            q_pair = float(row.iloc[0].q_BH_across_drugs_pairs)
            draw_bracket(
                ax, GROUPS.index(g1), GROUPS.index(g2),
                bracket_base + level * bracket_step, .005,
                bracket_label(q_pair), fontsize=3.5,
            )
        ax.set_ylim(.58, 1.145)
        ax.set_yticks([.6, .8, 1.0])
        ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
        if idx % 2 == 0:
            ax.set_ylabel("GDSC2 AUC", fontsize=6.1, labelpad=2)
        else:
            ax.set_yticklabels([])
        ax.set_xticks(np.arange(4))
        ax.set_xticklabels(
            labels, fontsize=4.85, rotation=45, ha="right",
            rotation_mode="anchor",
        )
        ax.tick_params(axis="x", pad=1.5)
        # Keep sample-size information visible but subordinate to the group
        # labels by drawing it as a separate, more tightly spaced line.
        for xpos, vals in enumerate(values):
            ax.text(
                xpos, -0.105, f"n={len(vals)}",
                transform=ax.get_xaxis_transform(), fontsize=3.65,
                color=INK, rotation=45, rotation_mode="anchor",
                ha="right", va="top", clip_on=False,
            )
        clean(ax)


def save(fig, path):
    fig.savefig(f"{path}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{path}.png", bbox_inches="tight", dpi=300)
    fig.savefig(f"{path}.svg", bbox_inches="tight", dpi=300)
    fig.savefig(f"{path}.tiff", bbox_inches="tight", dpi=600, pil_kwargs={"compression": "tiff_lzw"})


def _pstars(p):
    p = float(p)
    if p < 1e-4: return "****"
    if p < 1e-3: return "***"
    if p < 1e-2: return "**"
    if p < .05: return "*"
    return "ns"


def draw_chemo_box(ax, score, ylabel):
    d = pd.read_csv(CHEMO / "GSE31811_tumor_sample_scores.tsv", sep="\t")
    tests = pd.read_csv(CHEMO / "GSE31811_response_group_tests.tsv", sep="\t").set_index("score")
    vals = [d.loc[d.response.eq(g), score].dropna().to_numpy(float) for g in ["effective", "non-effective"]]
    bp = ax.boxplot(vals, positions=[0, 1], widths=.46, patch_artist=True, showfliers=False,
                    medianprops=dict(color=INK, lw=.85), whiskerprops=dict(color=GREY, lw=.6),
                    capprops=dict(color=GREY, lw=.6), boxprops=dict(lw=.55))
    for patch, col in zip(bp["boxes"], [RESPONDER, "#D8B9B5"]):
        patch.set_facecolor(col); patch.set_alpha(.82); patch.set_edgecolor("white")
    row = tests.loc[score]
    allv = np.concatenate(vals); yr = max(np.nanmax(allv) - np.nanmin(allv), .5)
    y = np.nanmax(allv) + .12 * yr; h = .035 * yr
    ax.plot([0, 0, 1, 1], [y, y+h, y+h, y], color=INK, lw=.5, clip_on=False)
    ax.text(.5, y+h+.015*yr, _pstars(row.p), ha="center", va="bottom", fontsize=5.0, fontweight="bold")
    ax.text(.98, .98, f"n={int(row.n_effective)}/{int(row.n_non_effective)}; q={row.q_BH:.3f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=4.2, color="#555555")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Effective", "Non-effective"], fontsize=4.7)
    ax.set_ylabel(ylabel, fontsize=5.5, labelpad=1)
    ax.set_ylim(np.nanmin(allv) - .10*yr, y + .19*yr)
    clean(ax)


def draw_chemo_forest(ax):
    c = pd.read_csv(CHEMO / "GSE14208_cox_models.tsv", sep="\t")
    c = c[(c.cohort == "GSE14208_all") & (c.endpoint == "TTP")].set_index("score").reindex(["sender_z", "receiver_z", "ecosystem_z"])
    labels = ["PLAT-rich fibroblast state", "C1QC+APOE+ macrophage state", "Ecosystem score"]
    y = np.arange(3)[::-1]
    ax.axvline(1, color="#777777", lw=.55, ls="--")
    for yy, lab, label in zip(y, c.index, labels):
        row = c.loc[lab]; col = SENDER if lab == "sender_z" else RECEIVER if lab == "receiver_z" else ECM
        ax.plot([row.CI95_low, row.CI95_high], [yy, yy], color=GREY, lw=.85)
        ax.scatter(row.HR_per_SD, yy, s=18, facecolor=col, edgecolor="white", lw=.45, zorder=3)
        ax.text(.98, yy, f"HR={row.HR_per_SD:.2f} {_pstars(row.p)}", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=4.15)
    ax.set_xscale("log"); ax.set_xlim(.55, 1.45)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=4.8)
    ax.set_xlabel("GSE14208 TTP HR per 1 SD", fontsize=5.3, labelpad=1)
    ax.text(.98, .02, "122 patients; CF chemotherapy", transform=ax.transAxes, ha="right", va="bottom", fontsize=4.0, color="#555555")
    clean(ax)


def draw_i_patient_chemo(fig, spec):
    frame = fig.add_subplot(spec, frameon=False)
    frame.set_xticks([]); frame.set_yticks([])
    for spine in frame.spines.values(): spine.set_visible(False)
    panel_label(frame, "i", x=-.18, y=1.08)
    grid = spec.subgridspec(2, 2, hspace=.52, wspace=.38)
    draw_chemo_box(fig.add_subplot(grid[0, 0]), "sender_z", "PLAT-rich fibroblast-state score (z)")
    draw_chemo_box(fig.add_subplot(grid[0, 1]), "receiver_z", "C1QC+APOE+ macrophage-state score (z)")
    draw_chemo_box(fig.add_subplot(grid[1, 0]), "ecosystem_z", "Ecosystem score (z)")
    draw_chemo_forest(fig.add_subplot(grid[1, 1]))


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""): h.update(b)
    return h.hexdigest()



def main():
    mcp = load_mcp()
    cp = load_checkpoint()
    fig = plt.figure(figsize=(183 / 25.4, 242 / 25.4), facecolor="white")
    gs = fig.add_gridspec(
        4, 15, height_ratios=[1.18, 1.07, 1.00, 1.06],
        left=.10, right=.99, top=.967, bottom=.058,
        hspace=.42, wspace=.82,
    )

    ax_a = fig.add_subplot(gs[0, :5])
    draw_a(ax_a)
    panel_label(ax_a, "a", x=-.035, y=1.075)
    apos = ax_a.get_position()
    ax_a.set_position([apos.x0 - .018, apos.y0, apos.width * .95, apos.height])

    immune_eps = [
        "Fibroblasts", "Endothelial cells", "Monocytic lineage",
        "Myeloid dendritic cells", "Neutrophils", "B lineage", "T cells",
        "Cytotoxic lymphocytes", "CD8 T cells", "NK cells",
    ]
    immune_lab = [
        "Fibroblasts", "Endothelial", "Monocytic", "Myeloid DC",
        "Neutrophils", "B lineage", "T cells", "Cytotoxic", "CD8 T", "NK cells",
    ]
    ax_b = fig.add_subplot(gs[0, 5:])
    draw_group_strip(
        ax_b, mcp, immune_eps,
        "TCGA_configuration_MCPcounter_group_summary.tsv",
        "TCGA_configuration_MCPcounter_pairwise_Dunn.tsv",
        "MCP-counter score", immune_lab,
    )
    panel_label(ax_b, "b", x=-.06, y=1.08)
    handles = [
        Patch(facecolor=GROUP_COLORS[g], edgecolor="none", alpha=.85,
              label=GROUP_DISPLAY[g]) for g in GROUPS
    ]
    legend_kw = dict(
        fontsize=5.25, ncol=2, loc="lower left", borderaxespad=0,
        handlelength=.90, handleheight=.72, columnspacing=.78,
        handletextpad=.34, labelspacing=.30, frameon=False,
    )
    ax_b.legend(handles=handles, bbox_to_anchor=(.02, 1.105), **legend_kw)

    cp_eps = [
        "CD86", "SIGLEC10", "TIM-3", "TIGIT", "ICOS",
        "SIRPalpha", "CD47", "PD-L1", "PD-1", "VISTA",
    ]
    ax_c = fig.add_subplot(gs[1, :10])
    draw_group_strip(
        ax_c, cp, cp_eps,
        "TCGA_configuration_checkpoint_group_summary.tsv",
        "TCGA_configuration_checkpoint_pairwise_Dunn.tsv",
        "Expression z", cp_eps,
    )
    panel_label(ax_c, "c", x=-.10, y=1.08)
    ax_c.legend(handles=handles, bbox_to_anchor=(.02, 1.055), **legend_kw)
    cpos = ax_c.get_position()
    ax_c.set_position([cpos.x0 - .045, cpos.y0, cpos.width * .99, cpos.height])

    draw_d_clinicopathology(fig, gs[1, 10:])

    ax_e = fig.add_subplot(gs[2, :5])
    draw_e_subtype_associations(ax_e)
    panel_label(ax_e, "e", x=-.13, y=1.08)
    epos = ax_e.get_position()
    ax_e.set_position([epos.x0 - .018, epos.y0, epos.width * .86, epos.height])

    ax_f = fig.add_subplot(gs[2, 5:9])
    draw_f_tcga_subtype_composition(ax_f)
    panel_label(ax_f, "f", x=-.15, y=1.08)
    fpos = ax_f.get_position()
    ax_f.set_position([fpos.x0 - .010, fpos.y0, fpos.width * .92, fpos.height])

    ax_g = fig.add_subplot(gs[3, :5])
    draw_g_unified(ax_g)
    panel_label(ax_g, "g", x=-.10, y=1.08)
    gpos = ax_g.get_position()
    ax_g.set_position([gpos.x0 - .070, gpos.y0, gpos.width * .98, gpos.height])

    ax_h = fig.add_subplot(gs[3, 5:9])
    draw_h_stars(ax_h)
    panel_label(ax_h, "h", x=-.10, y=1.08)
    hpos = ax_h.get_position()
    ax_h.set_position([hpos.x0 + .012, hpos.y0, hpos.width * .84, hpos.height])

    draw_i_gdsc2_violin(fig, gs[2:4, 9:])

    save(fig, OUT)
    plt.close(fig)

    outputs = [Path(f"{OUT}.{ext}") for ext in ["pdf", "png", "svg", "tiff"]]
    qa = {
        "figure": OUT.name,
        "canvas_mm": [183, 242],
        "archetype": "asymmetric_mixed",
        "survival_panels": "removed from the main figure at user request",
        "d_panel": "four-feature cross-cohort random-effects clinicopathology meta-analysis on a standardized ecosystem-score scale",
        "e_panel": "effect-size summary for established molecular/pathology classifications; exact tests and multiplicity are source-data traceable",
        "f_panel": "TCGA official subtype composition within the four median-defined configurations",
        "i_panel": "GDSC2 gastric cancer cell-line AUC violins with internal box/range summaries and all six pairwise BH-adjusted comparisons for each displayed drug",
        "drug_inference_unit": "gastric cancer cell line; n=28 overall, configuration counts 10/4/4/10; cisplatin has fewer measured models",
        "drug_multiplicity": "Kruskal-Wallis and Spearman test families were each BH-adjusted across nine prespecified GDSC2 drugs; four cytotoxic agents are displayed and no prespecified drug has q<0.05",
        "claim_boundary": "drug panel measures tumour-cell-intrinsic monoculture sensitivity and does not represent patient treatment benefit or stromal/immune effects",
        "b_c_preservation": "validated groupwise boxplot implementation and source data unchanged",
        "outputs": {p.name: {"bytes": p.stat().st_size, "sha256": sha(p)} for p in outputs},
    }
    (REVIEW / "Figure6_v22_readability_20260827_QA.json").write_text(
        json.dumps(qa, indent=2), encoding="utf-8"
    )
    print("\n".join(str(p) for p in outputs))


if __name__ == "__main__":
    main()
