# Academic Figure Skill Asset Confirmation (verified against the existing Figure 4F GSEA curve candidate)
# (a) Four preranked GSEA curves → 10_phagocytosis_gsea_curve_candidate_20261104/plot_significant_Figure4F_GSEA_curves.py → visual adaptation / native Python run
#     Reuses the existing curve specification, signed ranking metric, and frozen GSEA inputs; no inferential analysis is rerun.
# RULE: This is a candidate review figure. It is not merged into the manuscript or Supplementary Information.

import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Academic Figure Skill Typography Baseline — copied from the prior GSEA curve script
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

# Academic Figure Skill Nature/Cell/Science Color Palette — copied from the prior GSEA curve script
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
DIVERGING = ["#2166AC", "#F7F7F7", "#B2182B"]
ACCENT_RED = "#B2182B"
GREY = "#999999"
BLACK = "#222222"

# Academic Figure Skill Export Baseline — copied from the prior GSEA curve script
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})


def save_figure(fig, stem):
    fig.savefig(f"{stem}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{stem}.png", bbox_inches="tight", dpi=300)
    fig.savefig(f"{stem}.svg", bbox_inches="tight", dpi=300)


HERE = Path(__file__).resolve().parent
INPUT = HERE / "inputs"
FIGURES = HERE / "figures"
INDIVIDUAL = FIGURES / "individual_curves"
SOURCE = HERE / "source_data"
INDIVIDUAL.mkdir(parents=True, exist_ok=True)
SOURCE.mkdir(parents=True, exist_ok=True)

# Match Figure 4D's order and exact programme definitions.
PROGRAM_ORDER = [
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
    "GOBP_PHAGOCYTOSIS",
    "GOBP_CELL_MATRIX_ADHESION",
]
LABELS = {
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX": "ECM degradation",
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION": "Myeloid leukocyte migration",
    "GOBP_PHAGOCYTOSIS": "Phagocytosis",
    "GOBP_CELL_MATRIX_ADHESION": "Cell–matrix adhesion",
}
STEMS = {
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX": "01_ECM_degradation",
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION": "02_myeloid_leukocyte_migration",
    "GOBP_PHAGOCYTOSIS": "03_phagocytosis",
    "GOBP_CELL_MATRIX_ADHESION": "04_cell_matrix_adhesion",
}

ranked = pd.read_csv(INPUT / "Figure4f_preranked_gene_list.tsv", sep="\t")
results = pd.read_csv(INPUT / "Figure4f_preranked_GSEA.tsv", sep="\t")
frozen = pd.read_csv(INPUT / "myeloid_receiver_programs_frozen.tsv", sep="\t")

if len(ranked) != 8187 or ranked["gene_symbol"].duplicated().any():
    raise ValueError("Expected the original 8,187-gene Figure 4F ranking with unique gene symbols")
genes = ranked["gene_symbol"].astype(str).to_numpy()
scores = ranked["ranking_statistic"].astype(float).to_numpy()
if not np.all(np.diff(scores) <= 1e-12):
    raise ValueError("The ranked list is not sorted from highest to lowest signed statistic")
if set(PROGRAM_ORDER) - set(results["program"]):
    raise ValueError("One or more Figure 4D programmes are absent from the original GSEA results")
if set(PROGRAM_ORDER) - set(frozen["program"]):
    raise ValueError("One or more Figure 4D programmes are absent from the frozen gene sets")

gene_to_position = {gene: index for index, gene in enumerate(genes)}
curves = {}
qa_rows = []
for program in PROGRAM_ORDER:
    result_rows = results.loc[results["program"].eq(program)]
    if len(result_rows) != 1:
        raise ValueError(f"Expected one saved GSEA result for {program}")
    result = result_rows.iloc[0]
    members = set(frozen.loc[frozen["program"].eq(program), "gene_symbol"].astype(str))
    positions = np.array(sorted(gene_to_position[g] for g in members if g in gene_to_position), dtype=int)
    n_genes = len(genes)
    n_hits = len(positions)
    n_misses = n_genes - n_hits
    hit_weights = np.abs(scores[positions])
    if n_hits == 0 or n_hits == n_genes or np.isclose(hit_weights.sum(), 0.0):
        raise ValueError(f"Invalid tested gene set for {program}")
    increments = np.full(n_genes, -1.0 / n_misses, dtype=float)
    increments[positions] = hit_weights / hit_weights.sum()
    running = np.r_[0.0, np.cumsum(increments)]
    if not np.isclose(running[-1], 0.0, atol=1e-10):
        raise ValueError(f"Running ES did not return to zero for {program}")
    source_es = float(result["ES"])
    computed_es = float(running.min() if source_es < 0 else running.max())
    if not np.isclose(computed_es, source_es, atol=1e-10):
        raise ValueError(f"Recomputed ES does not match source for {program}: {computed_es} vs {source_es}")
    if n_hits != int(result["tested_set_size"]):
        raise ValueError(f"Tested-set size does not match source for {program}")
    direction = "C3" if float(result["NES"]) > 0 else "pooled C1/C4/C5"
    curves[program] = {
        "program": program,
        "label": LABELS[program],
        "stem": STEMS[program],
        "result": result,
        "members": members,
        "positions": positions,
        "running": running,
        "scores": scores,
        "genes": genes,
        "tested_set_size": n_hits,
        "frozen_set_size": len(members),
        "ES_recomputed": computed_es,
        "direction": direction,
    }
    qa_rows.append({
        "Figure4D_order": PROGRAM_ORDER.index(program) + 1,
        "program": program,
        "display_label": LABELS[program],
        "frozen_set_size": len(members),
        "tested_set_size": n_hits,
        "ES_recomputed": computed_es,
        "ES_source": source_es,
        "NES_source": float(result["NES"]),
        "GSEA_enrichment_direction": f"toward {direction}",
        "nominal_p_source": float(result["nominal_p"]),
        "BH_q_14_sets_source": float(result["FDR_BH_14_programmes"]),
        "ES_absolute_error": abs(computed_es - source_es),
        "leading_edge_size": int(result["leading_edge_size"]),
        "leading_edge_genes": result["leading_edge_genes"],
    })
    pd.DataFrame({
        "rank_position": np.arange(1, n_genes + 1),
        "gene_symbol": genes,
        "ranking_statistic_sign_log2FC_sqrt_QLF": scores,
        "in_gene_set": np.isin(genes, list(members)),
        "running_enrichment_score_after_gene": running[1:],
    }).to_csv(SOURCE / f"{STEMS[program]}_running_ES_source.tsv", sep="\t", index=False)

qa_table = pd.DataFrame(qa_rows)
qa_table.to_csv(SOURCE / "Figure4D_four_pathway_GSEA_statistics.tsv", sep="\t", index=False)


def draw_panel(subfig, curve):
    result = curve["result"]
    positions = curve["positions"]
    running = curve["running"]
    n_genes = len(curve["genes"])
    rank_x = np.arange(1, n_genes + 1)
    score_x = np.arange(0, n_genes + 1)
    positive_es = float(result["NES"]) > 0
    direction_color = ACCENT_RED if positive_es else CATEGORICAL[0]

    grid = subfig.add_gridspec(
        3, 1,
        height_ratios=[2.6, 0.45, 0.9],
        hspace=0.08,
        left=0.105, right=0.985, top=0.91, bottom=0.18,
    )
    ax_es = subfig.add_subplot(grid[0])
    ax_hits = subfig.add_subplot(grid[1], sharex=ax_es)
    ax_rank = subfig.add_subplot(grid[2], sharex=ax_es)

    ax_es.fill_between(
        score_x, running, 0,
        where=(running >= 0 if positive_es else running <= 0),
        color=direction_color, alpha=0.12, linewidth=0,
    )
    ax_es.plot(score_x, running, color=direction_color, linewidth=1.3)
    ax_es.axhline(0, color=GREY, linestyle="--", linewidth=0.55)
    ax_es.set_ylabel("Running ES")
    ax_es.set_xlim(0, n_genes)
    ax_es.set_ylim(
        min(float(running.min()) - 0.035, -0.05),
        max(float(running.max()) + 0.035, 0.05),
    )
    ax_es.tick_params(axis="x", labelbottom=False)
    ax_es.text(
        0.01, 0.94,
        f"ES = {float(result['ES']):.3f}",
        transform=ax_es.transAxes,
        va="top", ha="left", fontsize=7,
    )
    ax_es.text(
        0.99, 0.94,
        f"NES = {float(result['NES']):.3f}  |  nominal p = {float(result['nominal_p']):.3g}  |  BH q (14 sets) = {float(result['FDR_BH_14_programmes']):.3f}",
        transform=ax_es.transAxes,
        va="top", ha="right", fontsize=7,
    )
    ax_hits.vlines(positions + 1, 0, 1, color=BLACK, linewidth=0.34)
    ax_hits.set_ylim(0, 1)
    ax_hits.set_yticks([])
    ax_hits.set_ylabel("Set hits", rotation=0, ha="right", va="center", labelpad=8, fontsize=7)
    ax_hits.tick_params(axis="x", bottom=False, labelbottom=False)
    for spine in ("left", "bottom", "top", "right"):
        ax_hits.spines[spine].set_visible(False)

    positive = curve["scores"] >= 0
    negative = ~positive
    ax_rank.fill_between(rank_x, 0, curve["scores"], where=positive, color=ACCENT_RED, alpha=0.52, linewidth=0)
    ax_rank.fill_between(rank_x, 0, curve["scores"], where=negative, color=CATEGORICAL[0], alpha=0.58, linewidth=0)
    ax_rank.axhline(0, color=GREY, linewidth=0.55)
    ax_rank.set_ylabel("Rank metric", fontsize=7)
    ax_rank.set_xlabel("Genes ranked by sign(log2FC) × √QLF: higher in C3  →  higher in pooled C1/C4/C5")
    ax_rank.set_xlim(0, n_genes)
    ax_rank.set_xticks([0, 2000, 4000, 6000, n_genes])
    ax_rank.set_xticklabels(["1", "2,000", "4,000", "6,000", f"{n_genes:,}"])
    direction_text = "C3-high end" if positive_es else "pooled-comparator-high end"
    subfig.suptitle(
        f"{curve['label']}  |  GSEA enrichment toward {curve['direction']}",
        x=0.09, y=0.99, ha="left", fontweight="bold", fontsize=8,
    )


# Full-size individual panels use exactly the prior curve figure dimensions and layout.
for program in PROGRAM_ORDER:
    curve = curves[program]
    fig = plt.figure(figsize=(183 / 25.4, 112 / 25.4))
    subfig = fig.subfigures(1, 1)
    draw_panel(subfig, curve)
    save_figure(fig, INDIVIDUAL / curve["stem"])
    plt.close(fig)

# Composite: four full-width 183 mm panel cells in Figure 4D order, separated by 8 mm gutters.
PANEL_MM = 183.0
GUTTER_MM = 8.0
TOTAL_WIDTH_MM = 4 * PANEL_MM + 3 * GUTTER_MM
fig = plt.figure(figsize=(TOTAL_WIDTH_MM / 25.4, 112 / 25.4))
outer = fig.add_gridspec(
    1, 7,
    width_ratios=[PANEL_MM, GUTTER_MM, PANEL_MM, GUTTER_MM, PANEL_MM, GUTTER_MM, PANEL_MM],
    wspace=0,
)
for i, program in enumerate(PROGRAM_ORDER):
    subfig = fig.add_subfigure(outer[0, i * 2])
    draw_panel(subfig, curves[program])

combined_stem = FIGURES / "Figure4D_four_pathway_GSEA_curves_candidate"
save_figure(fig, combined_stem)
plt.close(fig)

# Save exact input checksums and a machine-readable QA record.
manifest_rows = []
for path in sorted(INPUT.glob("*.tsv")):
    manifest_rows.append({"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size})
pd.DataFrame(manifest_rows).to_csv(HERE / "inputs_SHA256.tsv", sep="\t", index=False)
qa = {
    "figure": "Figure4D four-panel candidate GSEA curves",
    "panel_order": [LABELS[p] for p in PROGRAM_ORDER],
    "layout": "one horizontal row; 183 mm panel cell × 112 mm high per pathway; 8 mm gutters",
    "contrast": "C3 macrophage cluster vs pooled patient-shared comparator myeloid clusters C1/C4/C5",
    "rank_metric": "sign(log2FC) × sqrt(edgeR QLF)",
    "ranked_background_n": int(len(genes)),
    "statistical_analysis_rerun": False,
    "description": "Reconstructed running enrichment curves from the frozen ranked list and gene sets; no differential-expression model, gene-set definition, or permutation test was rerun.",
    "original_gsea_permutations": "10,000 gene-set permutations; seed 20260801 (inherited from Figure 4F source results)",
    "ES_recomputed_matches_source_within_1e-10": bool((qa_table["ES_absolute_error"] < 1e-10).all()),
    "direction_note": "Direction means gene-set enrichment toward the C3-high or pooled-comparator-high end of the signed ranking; it does not mean biological pathway activation or inhibition.",
    "source_results": qa_rows,
}
(FIGURES / "figure_QA.json").write_text(json.dumps(qa, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps(qa, indent=2, ensure_ascii=False))
