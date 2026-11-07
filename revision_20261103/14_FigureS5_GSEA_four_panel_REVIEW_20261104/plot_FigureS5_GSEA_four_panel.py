# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (A-D) GSEA running-ES curves -> no matching GSEA production asset; inherit the
# existing plot_four_individual_GSEA_insertion.py visual system and frozen data.
# Layout-only revision: 2 x 2, uppercase A-D; no inferential analysis is rerun.

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

import hashlib
import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ORIGINAL_SOURCE = BASE / "11_Figure4D_four_pathway_GSEA_curves_candidate_20261104" / "source_data"
INPUT = HERE / "inputs"
OUT = HERE / "figures"
INPUT.mkdir(exist_ok=True)
OUT.mkdir(exist_ok=True)

STEMS = [
    "01_ECM_degradation", "02_myeloid_leukocyte_migration",
    "03_phagocytosis", "04_cell_matrix_adhesion",
]
LABELS = ["ECM degradation", "Myeloid leukocyte migration", "Phagocytosis", "Cell–matrix adhesion"]
PROGRAMS = [
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION", "GOBP_PHAGOCYTOSIS", "GOBP_CELL_MATRIX_ADHESION",
]
EXPECTED_NES = [1.550872505461684, 1.1300687754365901, -1.434637235467538, 1.0746747627839752]
EXPECTED_Q = [0.03636363636363636, 0.25080818965517243, 0.03636363636363636, 0.3395271743871587]
STAT_NAME = "Figure4D_four_pathway_GSEA_statistics.tsv"
source_names = [STAT_NAME] + [f"{stem}_running_ES_source.tsv" for stem in STEMS]
input_hashes = []
for name in source_names:
    src = ORIGINAL_SOURCE / name
    dst = INPUT / name
    if not dst.exists():
        shutil.copy2(src, dst)
    src_hash = hashlib.sha256(src.read_bytes()).hexdigest()
    dst_hash = hashlib.sha256(dst.read_bytes()).hexdigest()
    if src_hash != dst_hash:
        raise RuntimeError(f"Frozen input differs from original: {name}")
    input_hashes.append({"file": name, "sha256": dst_hash, "bytes": dst.stat().st_size})

stats = pd.read_csv(INPUT / STAT_NAME, sep="\t")
if list(stats.program) != PROGRAMS:
    raise ValueError("Unexpected programme order")
assert np.allclose(stats.NES_source, EXPECTED_NES, atol=1e-12)
assert np.allclose(stats.BH_q_14_sets_source, EXPECTED_Q, atol=1e-12)

# A symmetric supplementary grid: all four prespecified programme views have equal size.
# Draw at the actual 183-mm final width, rather than shrinking 183-mm single panels.
fig = plt.figure(figsize=(183 / 25.4, 166 / 25.4), facecolor="white")
outer = fig.add_gridspec(2, 2, left=0, right=1, top=1, bottom=0, wspace=0.045, hspace=0.065)
records = []
for idx, (stem, label) in enumerate(zip(STEMS, LABELS)):
    data = pd.read_csv(INPUT / f"{stem}_running_ES_source.tsv", sep="\t")
    result = stats.iloc[idx]
    n = len(data)
    hit = data.in_gene_set.astype(str).str.lower().isin(["true", "1"]).to_numpy()
    ranks = data.rank_position.to_numpy(int)
    metric = data.ranking_statistic_sign_log2FC_sqrt_QLF.to_numpy(float)
    running = np.r_[0., data.running_enrichment_score_after_gene.to_numpy(float)]
    assert n == 8187 and np.array_equal(ranks, np.arange(1, n + 1))
    assert not data.gene_symbol.duplicated().any()
    assert np.all(np.isfinite(metric)) and np.all(np.diff(metric) <= 1e-12)
    assert int(hit.sum()) == int(result.tested_set_size)
    assert abs(running[-1]) < 1e-10
    weights = np.abs(metric[hit])
    increments = np.full(n, -1 / (n - hit.sum()), dtype=float)
    increments[hit] = weights / weights.sum()
    reconstructed = np.r_[0., np.cumsum(increments)]
    assert np.allclose(reconstructed, running, atol=1e-10)
    source_es = float(result.ES_source)
    computed_es = float(running[np.argmax(np.abs(running))])
    assert abs(computed_es - source_es) < 1e-10
    nes = float(result.NES_source)
    q = float(result.BH_q_14_sets_source)
    direction_color = ACCENT_RED if nes > 0 else CATEGORICAL[0]

    subfig = fig.add_subfigure(outer[idx // 2, idx % 2])
    subfig.text(0.025, 0.985, "ABCD"[idx], ha="left", va="top", fontsize=12, fontweight="bold")
    subfig.text(0.15, 0.972, label, ha="left", va="top", fontsize=7.3, fontweight="bold")
    subfig.text(0.15, 0.91, f"ES = {source_es:.3f}    NES = {nes:.3f}", fontsize=6.5, va="top")
    subfig.text(0.15, 0.867,
                f"P = {float(result.nominal_p_source):.4f}    BH q = {q:.4f}", fontsize=6.5, va="top")
    grid = subfig.add_gridspec(3, 1, left=0.15, right=0.97, top=0.812, bottom=0.175,
                              height_ratios=[3.4, 0.52, 1.05], hspace=0.09)
    ax = subfig.add_subplot(grid[0])
    ax_hits = subfig.add_subplot(grid[1], sharex=ax)
    ax_rank = subfig.add_subplot(grid[2], sharex=ax)
    x = np.arange(n + 1)
    ax.fill_between(x, running, 0, where=(running >= 0 if nes > 0 else running <= 0),
                    color=direction_color, alpha=0.12, linewidth=0)
    ax.plot(x, running, color=direction_color, linewidth=1.05)
    ax.axhline(0, color=GREY, linestyle="--", linewidth=0.5)
    ax.set_xlim(0, n)
    ax.set_ylim(-0.39, 0.47)
    ax.set_yticks([-0.2, 0., 0.2, 0.4])
    ax.set_ylabel("Running ES", fontsize=7)
    ax.tick_params(axis="both", labelsize=6, length=2.5)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)

    ax_hits.vlines(ranks[hit], 0, 1, color=BLACK, linewidth=0.28)
    ax_hits.set_ylim(0, 1)
    ax_hits.set_yticks([])
    ax_hits.set_ylabel("Set hits", rotation=0, ha="right", va="center", fontsize=6, labelpad=5)
    ax_hits.tick_params(axis="x", bottom=False, labelbottom=False)
    for spine in ax_hits.spines.values():
        spine.set_visible(False)

    ax_rank.fill_between(ranks, 0, metric, where=(metric >= 0), color=ACCENT_RED, alpha=0.52, linewidth=0)
    ax_rank.fill_between(ranks, 0, metric, where=(metric < 0), color=CATEGORICAL[0], alpha=0.58, linewidth=0)
    ax_rank.axhline(0, color=GREY, linewidth=0.5)
    ax_rank.set_ylabel("Rank metric", fontsize=6)
    ax_rank.set_xlim(0, n)
    ax_rank.set_xticks([1, 2000, 4000, 6000, n])
    ax_rank.set_xticklabels(["1", "2,000", "4,000", "6,000", "8,187"])
    ax_rank.tick_params(axis="both", labelsize=6, length=2.5, pad=2)
    ax_rank.set_xlabel("Gene rank", fontsize=7, labelpad=3)
    records.append({
        "panel": "ABCD"[idx], "programme": label, "source_rows": n,
        "all_ranked_genes_retained": True, "set_hits": int(hit.sum()),
        "source_ES": source_es, "rendered_ES": computed_es,
        "ES_absolute_error": abs(computed_es - source_es), "NES": nes,
        "nominal_P": float(result.nominal_p_source), "BH_q_14_sets": q,
        "direction": "toward C3" if nes > 0 else "toward pooled C1/C4/C5",
        "BH_q_below_0_05": q < 0.05,
    })

# Fix final physical dimensions exactly; do not trim the panel-label margins.
stem = OUT / "Supplemental_Figure_S5_GSEA_four_panel_REVIEW"
with mpl.rc_context({"savefig.bbox": None}):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=300, bbox_inches=None, facecolor="white")
    fig.savefig(f"{stem}_600dpi.tiff", dpi=600, bbox_inches=None, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)

im = Image.open(f"{stem}.png").convert("RGB")
arr = np.asarray(im)
height, width = arr.shape[:2]
densities = []
for row in range(2):
    for col in range(2):
        patch = arr[row * height // 2:(row + 1) * height // 2,
                    col * width // 2:(col + 1) * width // 2]
        densities.append(float((patch.min(axis=2) < 245).mean()))
assert all(0.015 < d < 0.5 for d in densities)
pd.DataFrame(records).to_csv(HERE / "FigureS5_GSEA_statistics.tsv", sep="\t", index=False)
pd.DataFrame(input_hashes).to_csv(HERE / "INPUT_SHA256.tsv", sep="\t", index=False)
qa = {
    "figure": "Supplemental Figure S5", "layout": "2 x 2; uppercase A-D",
    "size_mm": [183, 166], "png_size_pixels": [width, height],
    "raster_dpi": {"PNG": 300, "TIFF": 600},
    "comparison": "C3 versus pooled C1/C4/C5; 26 paired patients",
    "statistical_analysis_rerun": False, "all_four_programmes_retained": True,
    "ranked_gene_n": 8187, "multiple_testing_family_n": 14,
    "all_source_ES_checks_passed": all(r["ES_absolute_error"] < 1e-10 for r in records),
    "panel_nonwhite_fraction": densities, "programmes": records,
    "source": str(ORIGINAL_SOURCE),
}
(HERE / "figure_QA.json").write_text(json.dumps(qa, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(qa, indent=2, ensure_ascii=False))
