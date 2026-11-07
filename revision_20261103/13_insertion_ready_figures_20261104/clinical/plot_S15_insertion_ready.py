# Display-only copy of the07 production script; all original plot data and styles retained.
# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (A, B, D, E, G, H) unpaired box + all patient points -> cross-type inherit from
# I:/skills/academic-figure-skill/assets/figures/PairedBoxScatter/plot_PairedBoxScatter1.R -> param inherit
# (C, F, I) correlation scatter -> cross-type inherit from
# I:/skills/academic-figure-skill/assets/figures/MarginalDensity/plot_VariableCorrelation.py -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       No native-run panels: paired-connection and log-log genomic jointplot
#       semantics are incompatible with these unpaired clinical comparisons.
#       Companion PNGs were inspected. Basic-plots contains a README only.

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

import argparse
import hashlib
import json
import platform
import shutil
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy
from scipy.stats import mannwhitneyu, spearmanr

HERE = Path(__file__).resolve().parent
UPSTREAM = HERE.parents[1] / "04_clinical_biomarker_relationships/Kim2018_biomarker_relationships_input_snapshot.tsv"
SAVED_INPUT = HERE / "Kim2018_biomarker_relationships_input_snapshot.tsv"
SCORES = [
    ("ecosystem_z", "Ecological score"),
    ("sender_primary_z", "Sender score"),
    ("receiver_specific_z", "Receiver score"),
]
BINARY = [
    ("msi_type_pentaplex", "MSI status", "MSS", "MSI", "MSS", "MSI"),
    ("ebv_in_situ", "EBV ISH status", "Negative", "Positive", "EBV−", "EBV+"),
]
JITTER_SEED = 2026110304
RESPONSE_STYLE = {
    0: {"label": "Non-responder (SD/PD)", "marker": "o", "color": CATEGORICAL[0], "facecolor": "white"},
    1: {"label": "Responder (CR/PR)", "marker": "^", "color": CATEGORICAL[3], "facecolor": CATEGORICAL[3]},
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bh(pvalues):
    values = np.asarray(pvalues, dtype=float)
    order = np.argsort(values)
    ordered = values[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.empty_like(values)
    adjusted[order] = np.minimum(1, np.minimum.accumulate(ordered[::-1])[::-1])
    return adjusted


def pformat(value):
    return f"{value:.3f}" if value >= 0.001 else f"{value:.2g}"


def analyse(data):
    rows = []
    summaries = []
    for score, label in SCORES:
        for biomarker, name, negative, positive, negative_label, positive_label in BINARY:
            neg = data.loc[data[biomarker].eq(negative), score].to_numpy(float)
            pos = data.loc[data[biomarker].eq(positive), score].to_numpy(float)
            # Positive group is the first argument: U and rank-biserial sign are
            # positive when the score is higher in biomarker-positive patients.
            test = mannwhitneyu(pos, neg, alternative="two-sided", method="asymptotic", use_continuity=True)
            rbc = 2 * test.statistic / (len(pos) * len(neg)) - 1
            rows.append({
                "score_column": score, "score_label": label, "biomarker_column": biomarker,
                "biomarker_label": name, "analysis_status": "post_hoc_exploratory",
                "test": "Two-sided Mann-Whitney U; asymptotic, tie-corrected variance, continuity correction",
                "n": len(pos) + len(neg), "missing_biomarker_n": int(data[biomarker].isna().sum()),
                "negative_group": negative, "positive_group": positive,
                "negative_n": len(neg), "positive_n": len(pos),
                "negative_median": np.median(neg), "positive_median": np.median(pos),
                "median_difference_positive_minus_negative": np.median(pos) - np.median(neg),
                "Mann_Whitney_U_positive_group": test.statistic, "rank_biserial_positive_minus_negative": rbc,
                "Spearman_rho": np.nan, "nominal_P": test.pvalue,
                "effect_orientation": "Positive values indicate a higher score in biomarker-positive than biomarker-negative patients",
                "inferential_unit": "Patient", "multiplicity_family": "All nine score-by-biomarker tests in this standalone analysis",
            })
            for category, values in [(negative, neg), (positive, pos)]:
                summaries.append({"score_column": score, "biomarker_column": biomarker, "group": category, "n": len(values),
                                  "median": np.median(values), "Q1": np.quantile(values, .25), "Q3": np.quantile(values, .75)})
        observed = data.loc[data.pdl1_cps.notna(), ["pdl1_cps", score]]
        test = spearmanr(observed.pdl1_cps, observed[score], alternative="two-sided")
        rows.append({
            "score_column": score, "score_label": label, "biomarker_column": "pdl1_cps",
            "biomarker_label": "PD-L1 CPS", "analysis_status": "post_hoc_exploratory",
            "test": "Two-sided Spearman rank correlation; average ranks for ties; asymptotic t approximation",
            "n": len(observed), "missing_biomarker_n": len(data) - len(observed),
            "negative_group": "", "positive_group": "", "negative_n": np.nan, "positive_n": np.nan,
            "negative_median": np.nan, "positive_median": np.nan,
            "median_difference_positive_minus_negative": np.nan,
            "Mann_Whitney_U_positive_group": np.nan, "rank_biserial_positive_minus_negative": np.nan,
            "Spearman_rho": test.statistic, "nominal_P": test.pvalue,
            "effect_orientation": "Positive rho indicates a higher score with higher PD-L1 CPS",
            "inferential_unit": "Patient", "multiplicity_family": "All nine score-by-biomarker tests in this standalone analysis",
        })
    stats = pd.DataFrame(rows)
    stats["BH_q_9_tests"] = bh(stats.nominal_P)
    return stats, pd.DataFrame(summaries)


def draw_binary(ax, data, score, biomarker, negative, positive, negative_label, positive_label, stat, rng):
    arrays = [data.loc[data[biomarker].eq(category), score].to_numpy(float) for category in [negative, positive]]
    colors = [GREY, GREY]
    # Inherited box width, thin charcoal boundary and all-patient dots; paired
    # connecting lines are deliberately inapplicable to independent patients.
    parts = ax.boxplot(arrays, positions=[0, 1], widths=.48, patch_artist=True,
                       showfliers=False, whis=1.5, manage_ticks=False,
                       medianprops={"color": BLACK, "linewidth": 1},
                       boxprops={"edgecolor": BLACK, "linewidth": .55},
                       whiskerprops={"color": BLACK, "linewidth": .55},
                       capprops={"color": BLACK, "linewidth": .55})
    for i, (values, color) in enumerate(zip(arrays, colors)):
        parts["boxes"][i].set_facecolor(mpl.colors.to_rgba(color, .13))
        group = data.loc[data[biomarker].eq([negative, positive][i])].copy()
        group["_jitter"] = rng.uniform(-.12, .12, len(group))
        # Draw responders last; colour/shape encodes treatment response. Both
        # boxes use neutral grey; biomarker status is defined by x-axis labels.
        for response, style in RESPONSE_STYLE.items():
            subset = group.loc[group.response_binary.eq(response)]
            ax.scatter(i + subset._jitter, subset[score], s=19 if response else 15,
                       marker=style["marker"], facecolors=style["facecolor"], edgecolors=style["color"],
                       linewidths=.7, alpha=.95 if response else .75, zorder=4 if response else 3)
    ax.set_xlim(-.48, 1.48)
    counts = [int(data.loc[data[biomarker].eq(category), "response_binary"].sum()) for category in [negative, positive]]
    ax.set_xticks([0, 1], [f"{negative_label}\n(n = {len(arrays[0])}; R = {counts[0]})", f"{positive_label}\n(n = {len(arrays[1])}; R = {counts[1]})"])
    ax.text(0, 1.045, f"Δmedian = {stat.median_difference_positive_minus_negative:+.2f}\nP = {pformat(stat.nominal_P)}; q = {pformat(stat.BH_q_9_tests)}",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6.5, linespacing=1.4)


def draw_cps(ax, data, score, stat):
    observed = data.loc[data.pdl1_cps.notna()]
    for response, style in RESPONSE_STYLE.items():
        subset = observed.loc[observed.response_binary.eq(response)]
        ax.scatter(np.log1p(subset.pdl1_cps), subset[score], s=22 if response else 17,
                   marker=style["marker"], facecolors=style["facecolor"], edgecolors=style["color"],
                   linewidths=.75, alpha=.95 if response else .75, zorder=4 if response else 3)
    ticks = np.array([0, 1, 5, 20, 100])
    ax.set_xticks(np.log1p(ticks), [str(v) for v in ticks])
    ax.set_xlim(-.16, np.log1p(100) + .18)
    ax.set_xlabel("PD-L1 CPS (log1p scale)", labelpad=3)
    ax.text(0, 1.045, f"Spearman ρ = {stat.Spearman_rho:+.2f}\nP = {pformat(stat.nominal_P)}; q = {pformat(stat.BH_q_9_tests)}",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6.5, linespacing=1.4)


def compose_figure(data, stats):
    fig, axes = plt.subplots(3, 3, figsize=(183 / 25.4, 170 / 25.4))
    fig.subplots_adjust(left=.105, right=.975, bottom=.14, top=.825, wspace=.42, hspace=.78)
    handles = [Line2D([], [], linestyle="none", marker=RESPONSE_STYLE[r]["marker"],
                      markerfacecolor=RESPONSE_STYLE[r]["facecolor"], markeredgecolor=RESPONSE_STYLE[r]["color"],
                      markeredgewidth=.8, markersize=5, label=RESPONSE_STYLE[r]["label"]) for r in [1, 0]]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.095,.928), ncol=2,
               fontsize=7, handlelength=1, handletextpad=.5, columnspacing=2, borderaxespad=0)
    for j, title in enumerate(["MSI status (pentaplex)", "EBV status (ISH)", "PD-L1 CPS\n(n = 41; 4 missing)"]):
        box = axes[0, j].get_position()
        fig.text(box.x0 + box.width / 2, .878, title, ha="center", fontsize=8, fontweight="bold")
    rng = np.random.default_rng(JITTER_SEED)
    for i, (score, label) in enumerate(SCORES):
        y = data[score].to_numpy(float)
        span = np.ptp(y)
        limits = (y.min() - .12 * span, y.max() + .12 * span)
        for j in range(3):
            ax = axes[i, j]
            stat = stats[(stats.score_column.eq(score)) & stats.biomarker_column.eq(BINARY[j][0] if j < 2 else "pdl1_cps")].iloc[0]
            ax.set_ylim(*limits)
            ax.axhline(0, color="#D0D0D0", lw=.45, ls="--", zorder=0)
            if j < 2:
                biomarker, _, negative, positive, neg_label, pos_label = BINARY[j]
                draw_binary(ax, data, score, biomarker, negative, positive, neg_label, pos_label, stat, rng)
            else:
                draw_cps(ax, data, score, stat)
            ax.text(-.22, 1.18, chr(ord("A") + i * 3 + j), transform=ax.transAxes,
                    fontsize=9, fontweight="bold", ha="left", va="bottom")
            ax.tick_params(length=2.5, pad=2)
            if j == 0:
                ax.set_ylabel(label + " (Z)", labelpad=5)
            else:
                ax.set_yticklabels([])
    return fig



ORIGINAL_FOLDER='I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026/WJGOsubmit/review/R2_supplementary_analyses_20261103/07_response_annotated_figures_20261103/figure1'

# Display-only export: retain the183mm canvas width and panel physical geometry.
# Trim only empty space after removing figure-level explanatory text.
from matplotlib.transforms import Bbox

def export_insertion(fig,prefix):
    fig.canvas.draw()
    tight=fig.get_tightbbox(fig.canvas.get_renderer())
    width=fig.get_size_inches()[0]
    bbox=Bbox.from_extents(0,tight.y0-0.025,width,tight.y1+0.025)
    for ext in ['pdf','svg','png']:
        fig.savefig(str(prefix)+'.'+ext,dpi=300,bbox_inches=bbox,pad_inches=0,facecolor='white')


def main():
    original=Path(__file__).resolve().parent/'source_data/figure1'
    data=pd.read_csv(original/'Kim2018_biomarker_relationships_input_snapshot.tsv',sep='\t')
    stats=pd.read_csv(original/'Kim2018_score_biomarker_relationships_all9.tsv',sep='\t')
    assert len(data)==45 and data.patient_id.nunique()==45
    assert len(stats)==9 and data.pdl1_cps.notna().sum()==41
    fig=compose_figure(data,stats)
    export_insertion(fig,Path(__file__).parent/'Supplemental_Figure15_insertion_ready')
    plt.close(fig)

if __name__=="__main__": main()
