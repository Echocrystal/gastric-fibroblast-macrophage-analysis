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
    fig.text(.105, .980, "Kim2018: clinical biomarkers and treatment response", ha="left", va="top", fontsize=9, fontweight="bold")
    fig.text(.105, .951, "45 patients: 12 responders, 33 non-responders; CPS observed in 41 (4 missing)", ha="left", va="top", fontsize=7, color=BLACK)
    handles = [Line2D([], [], linestyle="none", marker=RESPONSE_STYLE[r]["marker"],
                      markerfacecolor=RESPONSE_STYLE[r]["facecolor"], markeredgecolor=RESPONSE_STYLE[r]["color"],
                      markeredgewidth=.8, markersize=5, label=RESPONSE_STYLE[r]["label"]) for r in [1, 0]]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.095,.928), ncol=2,
               fontsize=7, handlelength=1, handletextpad=.5, columnspacing=2, borderaxespad=0)
    for j, title in enumerate(["MSI status (pentaplex)", "EBV status (ISH)", "PD-L1 CPS"]):
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
    fig.text(.105, .060, "P/q test score–biomarker relationships in all 45 patients (CPS: 41); one nine-test BH family.", fontsize=6.5, ha="left", color=BLACK)
    fig.text(.105, .039, "Response symbols are descriptive; no response-stratified tests. R = responder count.", fontsize=6.5, ha="left", color=BLACK)
    fig.text(.105, .018, "Binary panels: Δmedian = positive − negative. Each point represents one patient.", fontsize=6.5, ha="left", color=BLACK)
    return fig


def main():
    parser = argparse.ArgumentParser(description="Standalone post hoc score–clinical biomarker relationships and figure")
    parser.add_argument("--input", type=Path, default=SAVED_INPUT if SAVED_INPUT.exists() else UPSTREAM)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    source = args.input.resolve()
    previous_path = out / "biomarker_relationships_manifest.json"
    previous = json.loads(previous_path.read_text(encoding="utf-8")) if previous_path.exists() else {}
    origin = str(UPSTREAM.resolve())
    origin_hash = sha256(UPSTREAM) if UPSTREAM.exists() else previous.get("original_upstream_source_sha256", sha256(source))
    snapshot = out / SAVED_INPUT.name
    if snapshot.resolve() != source:
        shutil.copyfile(source, snapshot)
    data = pd.read_csv(source, sep="\t")
    assert len(data) == 45 and data.patient_id.nunique() == 45
    assert set(data.msi_type_pentaplex) == {"MSI", "MSS"}
    assert set(data.ebv_in_situ) == {"Positive", "Negative"}
    assert data.msi_type_pentaplex.eq("MSI").sum() == 6
    assert data.ebv_in_situ.eq("Positive").sum() == 5
    assert data.pdl1_cps.notna().sum() == 41
    assert np.isfinite(data[[x[0] for x in SCORES]]).all().all()
    assert data.pdl1_cps.dropna().between(0, 100).all()
    expected_response = data.best_response.map({"CR": 1, "PR": 1, "SD": 0, "PD": 0})
    assert expected_response.notna().all() and expected_response.eq(data.response_binary).all()
    assert data.response_binary.sum() == 12
    stats, summaries = analyse(data)
    original_stats = pd.read_csv(HERE / "original_9test_statistics_snapshot.tsv", sep="\t")
    assert np.allclose(stats.nominal_P, original_stats.nominal_P, atol=1e-15, rtol=0)
    assert np.allclose(stats.BH_q_9_tests, original_stats.BH_q_9_tests, atol=1e-15, rtol=0)
    stats.to_csv(out / "Kim2018_score_biomarker_relationships_all9.tsv", sep="\t", index=False)
    summaries.to_csv(out / "Kim2018_biomarker_group_score_summary.tsv", sep="\t", index=False)
    plot_data = data[["patient_id", "msi_type_pentaplex", "ebv_in_situ", "pdl1_cps", "best_response", "response_binary"] + [x[0] for x in SCORES]].copy()
    plot_data["response_label"] = plot_data.response_binary.map({0: "Non-responder (SD/PD)", 1: "Responder (CR/PR)"})
    plot_data["CPS_log1p_for_plot_only"] = np.log1p(plot_data.pdl1_cps)
    plot_data.to_csv(out / "Kim2018_biomarker_figure_source_data.tsv", sep="\t", index=False)
    fig = compose_figure(data, stats)
    # Exact physical canvas preserves the requested 183 x 170 mm. The standard
    # export baseline is included above; tight cropping is disabled for this
    # final multi-panel export so it cannot change the specified print size.
    prefix = out / "Figure1_clinical_biomarkers_response_annotated_v2"
    with mpl.rc_context({"savefig.bbox": None}):
        for suffix in ["pdf", "svg", "png"]:
            fig.savefig(str(prefix) + "." + suffix, dpi=300, bbox_inches=None, facecolor="white")
    plt.close(fig)

    legend = (
        "Figure 1. Relationships of the ecological, sender and receiver scores with established clinical biomarkers in Kim2018, annotated by anti-PD-1 response. "
        "Rows show the fixed ecological, fibroblast sender and macrophage receiver scores (within-cohort Z scores), respectively. "
        "A, D and G compare MSI-positive (n = 6) with microsatellite-stable (MSS; n = 39) tumours according to the pentaplex assay. "
        "B, E and H compare EBV in situ hybridisation-positive (n = 5) with EBV-negative (n = 40) tumours. "
        "Boxes indicate the median and first and third quartiles, whiskers extend to the most extreme values within 1.5 interquartile ranges, "
        "and all individual patients are shown. Orange filled triangles indicate responders (CR or PR; 12 of 45 patients) and blue open circles indicate non-responders "
        "(SD or PD; 33 of 45); response symbols are consistent across every panel. R in the group labels denotes the number of responders. "
        "Both biomarker-status groups use the same pale neutral-grey boxes; the x-axis labels identify the clinical groups, "
        "whereas point colour and shape identify treatment response. "
        "Δmedian denotes the positive-minus-negative group median difference. "
        "C, F and I show score relationships with PD-L1 combined positive score (CPS), available for 41 of 45 patients; four patients with missing CPS "
        "were omitted only from these three panels (11 responders and 30 non-responders have observed CPS; one responder and three non-responders have missing CPS). "
        "For visual separation of the many low CPS values, CPS is displayed on a log1p scale with tick labels "
        "in original CPS units; the scores and statistical tests use the original data. No fitted line or prediction model is shown. "
        "Binary comparisons use two-sided Mann-Whitney U tests with asymptotic tie correction and continuity correction. CPS associations use two-sided "
        "Spearman rank correlations with average ranks for ties and asymptotic P values. Nominal P values and Benjamini-Hochberg q values adjusted across "
        "all nine comparisons are annotated; these are unchanged from the previous all-patient analysis, not response-stratified tests. "
        "Response colour and shape provide an additional descriptive overlay only; no new response-stratified hypothesis tests were performed. "
        "Every point represents one independent patient. This post hoc analysis describes relationships with existing "
        "biomarkers; it does not establish independent response prediction, treatment benefit or absence of a relationship."
    )
    (out / "Figure_legend_EN.txt").write_text(legend + "\n", encoding="utf-8")
    findings = [
        "# 临床指标关系图 v2：增加应答/不应答标识（独立新图，未合并正式稿）", "",
        "新增显示：应答者为橙色实心三角，未应答者为蓝色空心圆；应答者最后绘制以减少遮挡。两类临床状态箱体统一为浅中性灰，由x轴标签区分，颜色仅用于应答状态。临床状态标签显示R应答人数；全部P/q仍为原来45/41人的score–biomarker检验，并未按应答重新检验。", "",
        "应答标签独立核对：CR 3例、PR 9例，对应response_binary=1；SD 15例、PD 18例，对应0。MSI组5/6应答，MSS组7/39；EBV阳性5/5，阴性7/40。CPS有记录者11/41应答，缺失者1/4应答。", "",
        "本次只比较已有的连续 ecological、sender、receiver score 与 MSI、EBV ISH、PD-L1 CPS 的关系，不按治疗应答选阈值，也不训练预测模型。", "",
        "MSI/EBV 两列保留全部45例：MSI阳性6例、MSS39例；EBV ISH阳性5例、阴性40例。CPS有记录41例，缺失4例，仅在CPS散点图及相关分析中排除这4例。", "",
        "9项比较统一进行BH校正。各项q均≥0.1739，没有一项达到q<0.05。", "",
        "| 连续评分 | MSI阳性−MSS中位差；P；q | EBV阳性−阴性中位差；P；q | CPS Spearman ρ；P；q |", "|---|---|---|---|",
    ]
    for score, label in SCORES:
        ss = stats[stats.score_column.eq(score)]
        mm = ss[ss.biomarker_column.eq("msi_type_pentaplex")].iloc[0]
        ee = ss[ss.biomarker_column.eq("ebv_in_situ")].iloc[0]
        cc = ss[ss.biomarker_column.eq("pdl1_cps")].iloc[0]
        findings.append(f"| {label} | {mm.median_difference_positive_minus_negative:+.3f}; {mm.nominal_P:.4f}; {mm.BH_q_9_tests:.4f} | {ee.median_difference_positive_minus_negative:+.3f}; {ee.nominal_P:.4f}; {ee.BH_q_9_tests:.4f} | {cc.Spearman_rho:+.3f}; {cc.nominal_P:.4f}; {cc.BH_q_9_tests:.4f} |")
    findings += ["", "数值上，sender score与CPS呈负相关（ρ = −0.364，nominal P = 0.0193），但9项校正后q = 0.1739；因此只能描述该趋势，不能称为校正后显著关联。Ecological score与三种指标的关联较弱。二分类阳性组只有5–6例，未显著不等于不存在关系。", "", "本图可以展示评分与临床指标的重叠程度，但不能单独回答是否独立预测抗PD-1应答；后者需要另外的、样本量足够的联合模型或外部验证。", "", "复现：使用包含matplotlib/numpy/pandas/scipy的Python执行 `analyze_and_plot_biomarker_relationships.py`。脚本默认优先读取同目录输入快照，所有输出保存于同目录。输入来源、SHA256、具体检验、软件版本与图像尺寸记录于manifest；统计全表为 `Kim2018_score_biomarker_relationships_all9.tsv`。"]
    (out / "Results_CN.md").write_text("\n".join(findings) + "\n", encoding="utf-8")
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(), "analysis_status": "post_hoc_exploratory_standalone",
        "scientific_question": "How do our three fixed continuous scores relate to MSI, EBV ISH and PD-L1 CPS?",
        "original_upstream_source_path": origin, "original_upstream_source_sha256": origin_hash,
        "runtime_input_path": str(source), "runtime_input_sha256": sha256(source),
        "input_snapshot": snapshot.name, "input_snapshot_sha256": sha256(snapshot),
        "script": Path(__file__).name, "script_sha256": sha256(Path(__file__)),
        "prior_figure_script_source": str(UPSTREAM.parent / "analyze_and_plot_biomarker_relationships.py"),
        "prior_figure_script_source_sha256": sha256(UPSTREAM.parent / "analyze_and_plot_biomarker_relationships.py") if (UPSTREAM.parent / "analyze_and_plot_biomarker_relationships.py").exists() else None,
        "prior_statistics_snapshot": "original_9test_statistics_snapshot.tsv",
        "unchanged_statistics_validation": "All 9 nominal P and BH q values agree with prior04 snapshot to absolute tolerance 1e-15",
        "scores": dict(SCORES), "n_patients": len(data), "n_MSI": 6, "n_MSS": 39, "n_EBV_positive": 5, "n_EBV_negative": 40,
        "n_CPS_available": 41, "n_CPS_missing": 4, "inferential_unit": "Patient",
        "tests": {"binary": "Two-sided Mann-Whitney U; method='asymptotic', use_continuity=True; SciPy tie-corrected variance",
                  "CPS": "Two-sided Spearman; average ranks for ties; asymptotic t approximation",
                  "multiplicity": "One BH family across all 9 score-by-biomarker tests; no tests selected by result",
                  "effect_binary": "Rank-biserial = 2*U_positive/(n_positive*n_negative)-1; median difference = median_positive-median_negative",
                  "effect_continuous": "Spearman rho on original CPS and score values"},
        "plot": {"archetype": "quantitative_grid", "physical_width_mm": 183, "physical_height_mm": 170, "raster_dpi": 300,
                 "all_data": "All45 in MSI and EBV; all41 observed CPS; no downsampling", "CPS_axis": "log1p for display only, original-unit tick labels", "jitter_seed_plot_only": JITTER_SEED,
                 "response_overlay": "Responders CR/PR: orange filled triangles, drawn last; nonresponders SD/PD: blue open circles. All patients retained.",
                 "box_colour": "Both biomarker-status groups use the same neutral grey; colours are reserved for response points",
                 "response_totals": {"CR": 3, "PR": 9, "SD": 15, "PD": 18, "CPS_observed_responders": 11, "CPS_observed_nonresponders": 30},
                 "no_models_or_cutoffs": "No response-derived cutoffs, ROC, prediction models or outcome-based point selection; no response-stratified tests"},
        "versions": {"python": platform.python_version(), "matplotlib": mpl.__version__, "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
        "files_sha256": {},
    }
    for file in sorted(out.iterdir()):
        if file.is_file() and file.name not in {previous_path.name, "QA_report.md"}:
            manifest["files_sha256"][file.name] = sha256(file)
    previous_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(stats[["score_label", "biomarker_label", "n", "nominal_P", "BH_q_9_tests", "Spearman_rho", "rank_biserial_positive_minus_negative"]].to_string(index=False))
    print("Saved", prefix)


if __name__ == "__main__":
    main()
