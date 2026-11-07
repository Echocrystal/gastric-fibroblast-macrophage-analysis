"""Generate revised Figure 4: intuitive patient ecology and phenotype.

Main panels show the prespecified primary association and direct patient-level
displays. Robust P-value variants and multiplicity details remain in the legend,
supplement and source data instead of occupying a redundant main panel.
"""

from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage, optimal_leaf_ordering
from scipy.spatial.distance import squareform
from scipy.stats import fisher_exact, rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from paper_plot_style import PROJECT, COLORS, clean_axis, save_figure

# Keep this figure's typography consistent after the narrower, taller redesign.
mpl.rcParams.update({
    "font.size": 7.0,
    "axes.labelsize": 7.0,
    "axes.titlesize": 7.5,
    "xtick.labelsize": 5.8,
    "ytick.labelsize": 5.8,
    "legend.fontsize": 5.8,
    "text.color": "#000000",
    "axes.labelcolor": "#000000",
    "axes.edgecolor": "#000000",
    "xtick.color": "#000000",
    "ytick.color": "#000000",
})


def panel_label(ax, label, x=-0.14, y=1.08):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=10.0,
            fontweight="bold", va="top", ha="left")


def format_p_value(value):
    value = float(value)
    if value >= 0.01:
        return f"{value:.2f}"
    mantissa, exponent = f"{value:.2e}".split("e")
    return rf"{float(mantissa):.2f}$\times 10^{{{int(exponent)}}}$"


BASE = PROJECT / "analysis/results/PLAT_axis/GSE183904"
SOURCE = PROJECT / "figures/source_data/Figure4"
SOURCE.mkdir(parents=True, exist_ok=True)

patient = pd.read_csv(
    BASE / "corrected_state_abundance_inference/GSE183904_corrected_state_abundance_patient_data.tsv",
    sep="\t",
)
features = pd.read_csv(
    BASE / "fine_state_axis/GSE183904_fine_state_axis_patient_features.tsv", sep="\t"
)
inference = pd.read_csv(
    BASE / "corrected_state_abundance_inference/GSE183904_corrected_state_abundance_inference.tsv",
    sep="\t",
)
deltas = pd.read_csv(
    BASE / "receiver_state_characterization/GSE183904_C1QC_APOE_vs_other_recurrent_paired_differences.tsv",
    sep="\t",
)
contrasts = pd.read_csv(
    BASE / "receiver_state_characterization/GSE183904_C1QC_APOE_program_contrasts.tsv",
    sep="\t",
)
state_scores = pd.read_csv(
    BASE / "receiver_state_characterization/GSE183904_myeloid_program_scores_by_patient_state.tsv",
    sep="\t",
)
coarse_scores = pd.read_csv(BASE / "GSE183904_PLAT_coarse_replication_scores.tsv.gz", sep="\t")
coarse_assoc = pd.read_csv(BASE / "GSE183904_PLAT_coarse_replication_associations.tsv", sep="\t")
receiver_de = pd.read_csv(
    BASE / "state_boundary_and_signatures/GSE183904_C1QC_APOE_vs_other_recurrent_paired_pseudobulk_DE.tsv",
    sep="\t",
)
program_genes = pd.read_csv(
    PROJECT / "resources/gene_sets/myeloid_receiver_programs_frozen.tsv", sep="\t"
)


ROSE = COLORS["sender_dark"]
ROSE_LIGHT = COLORS["dusty_rose"]
BLUE = COLORS["receiver_dark"]
BLUE_LIGHT = COLORS["powder_blue"]
STEEL = COLORS["steel_blue"]
GREEN = COLORS["ecm_dark"]
GREEN_LIGHT = COLORS["sage"]
CREAM = COLORS["cream"]
INK = "#000000"
GREY = COLORS["grey"]
LIGHT_GREY = "#E7EBEE"

PROGRAM_LABELS = {
    "GOBP_PHAGOCYTOSIS": "Phagocytosis",
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION": "Myeloid migration",
    "GOBP_CELL_MATRIX_ADHESION": "Cell–matrix adhesion",
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX": "ECM degradation",
    "ID_INFLAMMATORY_MONO": "Inflammatory monocyte",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNF–NF-κB signalling",
}
GSEA_LABELS = {
    "GOBP_CELL_MATRIX_ADHESION": "Cell–matrix adhesion",
    "GOBP_FOCAL_ADHESION_ASSEMBLY": "Focal-adhesion assembly",
    "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY": "Integrin-mediated signalling",
    "GOBP_MACROPHAGE_MIGRATION": "Macrophage migration",
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION": "Myeloid-leukocyte migration",
    "GOBP_PHAGOCYTOSIS": "Phagocytosis",
    "HALLMARK_HYPOXIA": "Hypoxia",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING": "IL-6/JAK/STAT3 signalling",
    "HALLMARK_INFLAMMATORY_RESPONSE": "Inflammatory response",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNF/NF-κB signalling",
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX": "ECM degradation",
    "REACTOME_FCGAMMA_RECEPTOR_FCGR_DEPENDENT_PHAGOCYTOSIS": "FcγR-dependent phagocytosis",
    "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS": "Integrin cell-surface interactions",
    "REACTOME_INTEGRIN_SIGNALING": "Integrin signalling",
}
CORE_PROGRAMS = [
    "GOBP_PHAGOCYTOSIS",
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
    "GOBP_CELL_MATRIX_ADHESION",
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
]
ALL_PROGRAMS = CORE_PROGRAMS + [
    "ID_INFLAMMATORY_MONO",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
]
STATE_LABELS = {
    "C1QC_APOE_macrophage": "C1QC+APOE+ macrophage",
    "HLAII_C1QC_macrophage": "HLA-II+C1QC+ macrophage",
    "SPP1_APOE_macrophage_donor_skewed": "SPP1+APOE+ macrophage",
    "FCN1_AREG_inflammatory_monocyte": "FCN1+AREG+ monocyte",
    "IL1B_CXCL8_inflammatory_myeloid": "IL1B+CXCL8+ myeloid",
}
STATE_ORDER = list(STATE_LABELS)


def rank_residual(values, covariates):
    y = rankdata(np.asarray(values, float), method="average")
    c = np.asarray(covariates, float)
    design = np.column_stack([np.ones(len(y)), c])
    return y - design @ np.linalg.lstsq(design, y, rcond=None)[0]


def abundance_residuals(df, x_col, y_col):
    cov = np.column_stack([
        df.ECM_direction.to_numpy(float),
        np.log1p(df.fibroblast_cells.to_numpy(float)),
        np.log1p(df.myeloid_cells.to_numpy(float)),
    ])
    return rank_residual(df[x_col], cov), rank_residual(df[y_col], cov)


def bh_adjust(p_values):
    """Benjamini-Hochberg adjustment with monotonicity enforcement."""
    p = np.asarray(p_values, float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def preranked_enrichment_score(scores, hit_positions):
    """Weighted preranked-GSEA enrichment score and leading-edge positions."""
    scores = np.asarray(scores, float)
    positions = np.sort(np.asarray(hit_positions, int))
    n_genes = len(scores)
    n_hits = len(positions)
    weights = np.abs(scores[positions])
    if n_hits == 0 or n_hits == n_genes or weights.sum() == 0:
        raise ValueError("GSEA requires a non-empty, non-universal weighted gene set")
    cumulative_hit = np.cumsum(weights) / weights.sum()
    misses_before = positions - np.arange(n_hits)
    running_after = cumulative_hit - misses_before / (n_genes - n_hits)
    previous_hit = np.r_[0.0, cumulative_hit[:-1]]
    running_before = previous_hit - misses_before / (n_genes - n_hits)
    positive_index = int(np.argmax(running_after))
    negative_index = int(np.argmin(running_before))
    positive_es = float(max(0.0, running_after[positive_index]))
    negative_es = float(min(0.0, running_before[negative_index]))
    if positive_es >= abs(negative_es):
        return positive_es, positions[:positive_index + 1]
    return negative_es, positions[negative_index:]


def preranked_gsea(scores, genes, gene_sets, permutations=10000, seed=20260801):
    """Deterministic gene-set-permutation GSEA for a small frozen collection."""
    scores = np.asarray(scores, float)
    genes = np.asarray(genes, str)
    gene_to_position = {gene: idx for idx, gene in enumerate(genes)}
    rng = np.random.default_rng(seed)
    result_rows = []
    leading_rows = []
    for set_index, (program, members) in enumerate(gene_sets.items()):
        tested_members = sorted(set(members) & set(gene_to_position))
        positions = np.array([gene_to_position[gene] for gene in tested_members], dtype=int)
        observed_es, leading_positions = preranked_enrichment_score(scores, positions)
        null_es = np.empty(permutations, dtype=float)
        for permutation in range(permutations):
            sampled = rng.choice(len(scores), size=len(positions), replace=False)
            null_es[permutation] = preranked_enrichment_score(scores, sampled)[0]
        if observed_es >= 0:
            same_sign = null_es[null_es >= 0]
            nominal_p = (np.count_nonzero(same_sign >= observed_es) + 1) / (len(same_sign) + 1)
            normalization = same_sign.mean()
        else:
            same_sign = null_es[null_es < 0]
            nominal_p = (np.count_nonzero(same_sign <= observed_es) + 1) / (len(same_sign) + 1)
            normalization = abs(same_sign.mean())
        nes = observed_es / normalization
        leading_genes = genes[leading_positions]
        result_rows.append({
            "program": program,
            "tested_set_size": len(positions),
            "ES": observed_es,
            "NES": nes,
            "nominal_p": nominal_p,
            "permutations": permutations,
            "ranking_metric": "sign(logFC) * sqrt(edgeR QLF)",
            "permutation_seed": seed,
            "leading_edge_size": len(leading_genes),
            "leading_edge_genes": "|".join(leading_genes),
        })
        for rank_index, position in enumerate(leading_positions, start=1):
            leading_rows.append({
                "program": program,
                "leading_edge_order": rank_index,
                "gene_symbol": genes[position],
                "rank_position": int(position + 1),
                "ranking_statistic": float(scores[position]),
            })
    results = pd.DataFrame(result_rows)
    results["FDR_BH_14_programmes"] = bh_adjust(results["nominal_p"])
    return results, pd.DataFrame(leading_rows)


def leave_one_out_full_ranks(df, x_col, y_col):
    """Match the frozen R implementation: full ranks, then delete one patient."""
    xr_full = rankdata(df[x_col].to_numpy(float), method="average")
    yr_full = rankdata(df[y_col].to_numpy(float), method="average")
    cov = np.column_stack([
        df.ECM_direction.to_numpy(float),
        np.log1p(df.fibroblast_cells.to_numpy(float)),
        np.log1p(df.myeloid_cells.to_numpy(float)),
    ])
    out = []
    for i in range(len(df)):
        keep = np.arange(len(df)) != i
        design = np.column_stack([np.ones(keep.sum()), cov[keep]])
        xr = xr_full[keep] - design @ np.linalg.lstsq(design, xr_full[keep], rcond=None)[0]
        yr = yr_full[keep] - design @ np.linalg.lstsq(design, yr_full[keep], rcond=None)[0]
        out.append({
            "dropped_sample_id": df.iloc[i].sample_id,
            "dropped_Patient_ID": df.iloc[i].Patient_ID,
            "adjusted_rho": float(np.corrcoef(xr, yr)[0, 1]),
        })
    return pd.DataFrame(out)


def descriptive_fit(ax, x, y, color, fill, seed):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    xx = np.linspace(x.min(), x.max(), 160)
    slope, intercept = np.polyfit(x, y, 1)
    rng = np.random.default_rng(seed)
    curves = np.empty((2500, len(xx)))
    for i in range(len(curves)):
        take = rng.integers(0, len(x), len(x))
        b, a = np.polyfit(x[take], y[take], 1)
        curves[i] = b * xx + a
    low, high = np.quantile(curves, [0.025, 0.975], axis=0)
    ax.fill_between(xx, low, high, color=fill, alpha=0.32, linewidth=0, zorder=1)
    ax.plot(xx, slope * xx + intercept, color=color, linewidth=1.35, zorder=3)


# ---------------------------------------------------------------------------
# Recompute and assert all panel-defining values.
# ---------------------------------------------------------------------------
raw_rho, raw_p = spearmanr(features.sender_state_fraction, features.receiver_state_fraction)
if not np.isclose(raw_rho, 0.684561473503083, atol=1e-12):
    raise RuntimeError("Raw patient-fraction association changed")

primary_x, primary_y = abundance_residuals(patient, "sender_logratio", "receiver_logratio")
primary_rho = float(np.corrcoef(primary_x, primary_y)[0, 1])
if not np.isclose(primary_rho, inference.iloc[0].adjusted_rank_residual_rho, atol=1e-10):
    raise RuntimeError("Primary adjusted association differs from frozen R output")

loo_parts = []
for name, xcol, ycol, row_idx in [
    ("Log-ratio", "sender_logratio", "receiver_logratio", 0),
    ("Arcsine fraction", "sender_asin_fraction", "receiver_asin_fraction", 1),
]:
    tab = leave_one_out_full_ranks(patient, xcol, ycol)
    tab["representation"] = name
    row = inference.iloc[row_idx]
    if not np.isclose(tab.adjusted_rho.min(), row.loo_rho_min, atol=1e-10):
        raise RuntimeError(f"LOO minimum mismatch: {name}")
    if not np.isclose(tab.adjusted_rho.max(), row.loo_rho_max, atol=1e-10):
        raise RuntimeError(f"LOO maximum mismatch: {name}")
    if not (tab.adjusted_rho > 0).all():
        raise RuntimeError(f"Non-positive LOO estimate: {name}")
    loo_parts.append(tab)
loo = pd.concat(loo_parts, ignore_index=True)

core_delta = deltas[deltas.program.isin(CORE_PROGRAMS)].copy()
if len(core_delta) != 104 or not (core_delta.delta_target_minus_other_recurrent > 0).all():
    raise RuntimeError("Expected four positive programmes in all 26 patients")
PROGRAM_DISPLAY_ORDER = (
    contrasts.set_index("program").loc[CORE_PROGRAMS]
    .sort_values("median_delta", ascending=False)
    .index.tolist()
)
if PROGRAM_DISPLAY_ORDER != [
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
    "GOBP_PHAGOCYTOSIS",
    "GOBP_CELL_MATRIX_ADHESION",
]:
    raise RuntimeError("Core programme effect-size order changed")

state_summary = (
    state_scores[state_scores.program.isin(ALL_PROGRAMS)]
    .groupby(["refined_label", "program"], as_index=False)
    .agg(
        median_patient_score=("mean_program_score", "median"),
        mean_patient_score=("mean_program_score", "mean"),
        n_patients=("Patient_ID", "nunique"),
    )
)
state_summary["state_relative_z"] = state_summary.groupby("program")[
    "median_patient_score"
].transform(lambda x: (x - x.mean()) / x.std(ddof=0))
state_matrix = state_summary.pivot(
    index="refined_label", columns="program", values="state_relative_z"
).loc[STATE_ORDER, ALL_PROGRAMS]
if state_summary.n_patients.min() != 26:
    raise RuntimeError("Every displayed myeloid state/program must cover all 26 patients")

# Molecular layer for panel e. The rule is fixed before layout: retain every
# eligible gene assigned to at least two core programmes, then add the
# highest-logFC programme-specific gene for each core programme. Ribbons encode
# only frozen gene-set membership, never regulatory direction or flow volume.
eligible_de = receiver_de[
    (receiver_de["FDR"] <= 0.05)
    & (receiver_de["logFC_target_minus_comparator"] >= 0.5)
    & (receiver_de["paired_direction_consistency"] >= 0.65)
    & (receiver_de["target_patient_detection_fraction"] >= 0.5)
    & (~receiver_de["technical_excluded"].astype(bool))
].copy()
network_all = program_genes[program_genes["program"].isin(CORE_PROGRAMS)].merge(
    eligible_de[[
        "gene_symbol", "logFC_target_minus_comparator", "FDR",
        "paired_direction_consistency", "target_patient_detection_fraction",
    ]],
    on="gene_symbol",
    how="inner",
)
programme_gene_counts = network_all.groupby("program")["gene_symbol"].nunique().reindex(CORE_PROGRAMS)
gene_programme_counts = network_all.groupby("gene_symbol")["program"].nunique()
shared_genes = set(gene_programme_counts[gene_programme_counts >= 2].index)
selected_genes = set(shared_genes)
programme_specific_genes = {}
for programme in CORE_PROGRAMS:
    candidates = network_all[
        network_all["program"].eq(programme)
        & ~network_all["gene_symbol"].isin(shared_genes)
    ].sort_values(
        ["logFC_target_minus_comparator", "FDR", "gene_symbol"],
        ascending=[False, True, True],
    )
    gene = str(candidates.iloc[0]["gene_symbol"])
    selected_genes.add(gene)
    programme_specific_genes[programme] = gene

expected_selected = {
    "ADAM9", "C2", "CCL13", "CD300A", "CD36", "CTSK", "EMILIN1",
    "EMP2", "ITGA9", "MMP12", "MMP14", "OTOA", "TREM2",
}
if selected_genes != expected_selected:
    raise RuntimeError("Frozen panel-e gene selection changed")
if programme_gene_counts.tolist() != [27, 33, 22, 14]:
    raise RuntimeError("Core programme upregulated-gene counts changed")

network_all["displayed"] = network_all["gene_symbol"].isin(selected_genes)
network_all["selection_reason"] = np.select(
    [
        network_all["gene_symbol"].isin(shared_genes),
        network_all.apply(
            lambda row: programme_specific_genes.get(row["program"]) == row["gene_symbol"],
            axis=1,
        ),
    ],
    ["member_of_at_least_two_core_programmes", "top_logFC_programme_specific_gene"],
    default="eligible_not_displayed",
)
network_display = network_all[network_all["displayed"]].copy()
if network_display["gene_symbol"].nunique() != 13 or len(network_display) != 23:
    raise RuntimeError("Unexpected panel-e network dimensions")

gene_programmes = network_all.groupby("gene_symbol")["program"].agg(
    lambda values: "|".join(p for p in CORE_PROGRAMS if p in set(values))
)
volcano_source = receiver_de.copy()
volcano_source["minus_log10_FDR"] = -np.log10(
    volcano_source["FDR"].clip(lower=np.finfo(float).tiny)
)
volcano_source["core_programme_eligible"] = volcano_source["gene_symbol"].isin(
    network_all["gene_symbol"]
)
volcano_source["displayed_in_membership_matrix"] = volcano_source["gene_symbol"].isin(
    selected_genes
)
volcano_source["n_core_programmes"] = (
    volcano_source["gene_symbol"].map(gene_programme_counts).fillna(0).astype(int)
)
volcano_source["core_programmes"] = (
    volcano_source["gene_symbol"].map(gene_programmes).fillna("")
)
volcano_plot = volcano_source.copy()
volcano_plot["minus_log10_FDR_capped"] = volcano_plot["minus_log10_FDR"].clip(upper=30)

# Functional over-representation is computed from the full DE universe, not
# from the 13 genes selected for display. This separates pathway-level evidence
# from illustrative gene membership.
de_universe = set(
    receiver_de.loc[~receiver_de["technical_excluded"].astype(bool), "gene_symbol"]
)
upregulated_universe = set(eligible_de["gene_symbol"]) & de_universe
ora_rows = []
for programme in CORE_PROGRAMS:
    tested_set = set(
        program_genes.loc[program_genes["program"].eq(programme), "gene_symbol"]
    ) & de_universe
    overlap = upregulated_universe & tested_set
    table = np.array([
        [len(overlap), len(upregulated_universe - tested_set)],
        [len(tested_set - upregulated_universe), len(de_universe - (upregulated_universe | tested_set))],
    ])
    odds_ratio, p_value = fisher_exact(table, alternative="greater")
    ora_rows.append({
        "program": programme,
        "universe_genes": len(de_universe),
        "upregulated_genes": len(upregulated_universe),
        "tested_programme_genes": len(tested_set),
        "overlap_genes": len(overlap),
        "odds_ratio": float(odds_ratio),
        "fisher_one_sided_p": float(p_value),
        "overlap_gene_symbols": "|".join(sorted(overlap)),
    })
programme_ora = pd.DataFrame(ora_rows)
programme_ora["FDR_BH_4_programmes"] = bh_adjust(programme_ora["fisher_one_sided_p"])
programme_ora["gene_ratio"] = (
    programme_ora["overlap_genes"] / programme_ora["tested_programme_genes"]
)
expected_ora = programme_ora.set_index("program").loc[CORE_PROGRAMS]
if not np.allclose(
    expected_ora["odds_ratio"],
    [1.093550, 1.678758, 1.266988, 2.030453],
    atol=1e-6,
):
    raise RuntimeError("Panel-e ORA odds ratios changed")
if not np.allclose(
    expected_ora["FDR_BH_4_programmes"],
    [0.370348, 0.035471, 0.255757, 0.046560],
    atol=1e-6,
):
    raise RuntimeError("Panel-e ORA FDR values changed")

# Full-ranked functional interpretation for panel e. All 14 prespecified
# external functional sets are shown; no top-N or significance filter is used.
ranked_de = receiver_de.loc[~receiver_de["technical_excluded"].astype(bool)].copy()
ranked_de["ranking_statistic"] = (
    np.sign(ranked_de["logFC_target_minus_comparator"])
    * np.sqrt(ranked_de["QLF"].clip(lower=0))
)
ranked_de["absolute_ranking_statistic"] = ranked_de["ranking_statistic"].abs()
ranked_de = (
    ranked_de.sort_values(
        ["absolute_ranking_statistic", "FDR", "gene_symbol"],
        ascending=[False, True, True],
    )
    .drop_duplicates("gene_symbol", keep="first")
    .sort_values(["ranking_statistic", "gene_symbol"], ascending=[False, True])
    .reset_index(drop=True)
)
if len(ranked_de) != 8187:
    raise RuntimeError("Unexpected number of non-technical genes in GSEA ranking")

functional_programs = sorted(
    program_genes.loc[
        program_genes["evidence_tier"].eq("primary_external"), "program"
    ].unique()
)
if set(functional_programs) != set(GSEA_LABELS):
    raise RuntimeError("Frozen external functional-program collection changed")
functional_gene_sets = {
    program: set(program_genes.loc[program_genes["program"].eq(program), "gene_symbol"])
    for program in functional_programs
}
gsea_results, gsea_leading_edge = preranked_gsea(
    ranked_de["ranking_statistic"].to_numpy(float),
    ranked_de["gene_symbol"].to_numpy(str),
    functional_gene_sets,
    permutations=10000,
    seed=20260801,
)
gsea_check = gsea_results.set_index("program")
expected_nes = {
    "GOBP_MACROPHAGE_MIGRATION": 1.5737658794386498,
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX": 1.550872505461684,
    "HALLMARK_HYPOXIA": -1.5509364860942672,
    "HALLMARK_INFLAMMATORY_RESPONSE": -2.1484228252371125,
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": -2.6346708443530718,
    "GOBP_PHAGOCYTOSIS": -1.434637235467538,
    "REACTOME_FCGAMMA_RECEPTOR_FCGR_DEPENDENT_PHAGOCYTOSIS": -1.4455571439444694,
}
for program, expected in expected_nes.items():
    if not np.isclose(gsea_check.loc[program, "NES"], expected, atol=1e-10):
        raise RuntimeError(f"Panel-e preranked GSEA changed: {program}")
if set(gsea_check.index[gsea_check["FDR_BH_14_programmes"] < 0.05]) != set(expected_nes):
    raise RuntimeError("Panel-e significant GSEA programme set changed")

ranked_gene_universe = set(ranked_de["gene_symbol"])
tested_sets = {
    program: members & ranked_gene_universe
    for program, members in functional_gene_sets.items()
}
n_programs = len(functional_programs)
jaccard_matrix = np.eye(n_programs, dtype=float)
jaccard_rows = []
for i, program_a in enumerate(functional_programs):
    for j, program_b in enumerate(functional_programs):
        intersection = tested_sets[program_a] & tested_sets[program_b]
        union = tested_sets[program_a] | tested_sets[program_b]
        similarity = len(intersection) / len(union)
        jaccard_matrix[i, j] = similarity
        jaccard_rows.append({
            "program_a": program_a,
            "program_b": program_b,
            "intersection_genes": len(intersection),
            "union_genes": len(union),
            "jaccard_similarity": similarity,
        })
jaccard_long = pd.DataFrame(jaccard_rows)
distance_matrix = 1.0 - jaccard_matrix
np.fill_diagonal(distance_matrix, 0.0)
condensed_distance = squareform(distance_matrix, checks=True)
cluster_tree = linkage(condensed_distance, method="average")
cluster_tree = optimal_leaf_ordering(cluster_tree, condensed_distance)
leaf_order = leaves_list(cluster_tree)
cluster_ids = fcluster(cluster_tree, t=6, criterion="maxclust")
cluster_table = pd.DataFrame({
    "program": functional_programs,
    "functional_cluster": cluster_ids,
})
cluster_table["cluster_leaf_order"] = cluster_table["program"].map(
    {functional_programs[index]: order for order, index in enumerate(leaf_order, start=1)}
)
cluster_theme = {}
for cluster_id, frame in cluster_table.groupby("functional_cluster"):
    members = set(frame["program"])
    if "HALLMARK_INFLAMMATORY_RESPONSE" in members:
        theme = "Inflammatory signalling"
    elif "GOBP_MACROPHAGE_MIGRATION" in members:
        theme = "Migration"
    elif "HALLMARK_HYPOXIA" in members:
        theme = "Hypoxia"
    elif "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX" in members:
        theme = "ECM / adhesion"
    elif "GOBP_PHAGOCYTOSIS" in members:
        theme = "Phagocytosis"
    elif "REACTOME_INTEGRIN_SIGNALING" in members:
        theme = "Integrin signalling"
    else:
        raise RuntimeError(f"Unresolved functional cluster: {sorted(members)}")
    cluster_theme[int(cluster_id)] = theme
cluster_table["functional_theme"] = cluster_table["functional_cluster"].map(cluster_theme)
gsea_results = gsea_results.merge(cluster_table, on="program", how="left", validate="one_to_one")
gsea_results["display_label"] = gsea_results["program"].map(GSEA_LABELS)
gsea_results = gsea_results.sort_values("cluster_leaf_order").reset_index(drop=True)

focus_programs = ["GOBP_PHAGOCYTOSIS", "GOBP_MYELOID_LEUKOCYTE_MIGRATION"]
coarse_focus = coarse_scores[
    (coarse_scores.variant == "primary_20_20") & coarse_scores.program.isin(focus_programs)
].copy()
assoc_focus = coarse_assoc[
    (coarse_assoc.variant == "primary_20_20")
    & (coarse_assoc.exposure == "fibroblast_PLAT")
    & coarse_assoc.program.isin(focus_programs)
].set_index("program").loc[focus_programs].reset_index()
coarse_residual_rows = []
for program in focus_programs:
    dat = coarse_focus[coarse_focus.program == program].copy().reset_index(drop=True)
    cov = np.column_stack([
        dat.ECM_direction.to_numpy(float),
        np.log1p(dat.fibroblast_cells.to_numpy(float)),
        np.log1p(dat.myeloid_cells.to_numpy(float)),
    ])
    xr = rank_residual(dat.fibroblast_PLAT_z, cov)
    yr = rank_residual(dat.program_score, cov)
    rho = float(np.corrcoef(xr, yr)[0, 1])
    frozen = assoc_focus[assoc_focus.program == program].iloc[0]
    if not np.isclose(rho, frozen.adjusted_rho, atol=1e-10):
        raise RuntimeError(f"Coarse adjusted rho mismatch: {program}")
    tmp = dat[["sample_id", "Patient_ID", "fibroblast_PLAT_z", "program_score"]].copy()
    tmp["program"] = program
    tmp["PLAT_rank_residual"] = xr
    tmp["program_rank_residual"] = yr
    tmp["adjusted_rho"] = rho
    tmp["permutation_p"] = float(frozen.adjusted_permutation_p)
    tmp["FDR_BH_12_program_family"] = float(frozen.adjusted_permutation_FDR_BH)
    coarse_residual_rows.append(tmp)
coarse_residuals = pd.concat(coarse_residual_rows, ignore_index=True)


# ---------------------------------------------------------------------------
# Current source-data package.
# ---------------------------------------------------------------------------
source_a = features[[
    "sample_id", "Patient_ID", "sender_state_fraction", "receiver_state_fraction",
    "fibroblast_cells", "myeloid_cells",
]].copy()
source_a["raw_spearman_rho"] = raw_rho
source_a["raw_spearman_p"] = raw_p

source_b = patient[["sample_id", "Patient_ID"]].copy()
source_b["sender_rank_residual"] = primary_x
source_b["receiver_rank_residual"] = primary_y
source_b["adjusted_rho"] = primary_rho
source_b["partial_t_p"] = float(inference.iloc[0].correct_partial_asymptotic_p)
source_b["HC3_p"] = float(inference.iloc[0].HC3_p)
source_b["wild_bootstrap_p"] = float(inference.iloc[0].wild_bootstrap_t_p)

tables = [
    ("4a", "Figure4a_raw_patient_ecology.tsv", source_a,
     "Raw target-state lineage fractions; one row per patient."),
    ("4b", "Figure4b_adjusted_state_covariation.tsv", source_b,
     "Primary log-ratio rank residuals and all frozen uncertainty estimates."),
    ("4c", "Figure4c_leave_one_patient_out.tsv", loo,
     "All leave-one-patient-out adjusted correlations for both representations."),
    ("4d", "Figure4d_patient_program_differences.tsv", core_delta,
     "Four target-minus-comparator programme differences in 26 patients."),
    ("4e", "Figure4e_pseudobulk_DE.tsv", volcano_source,
     "Full patient-blocked C1QC+APOE macrophage differential-expression result with programme flags."),
    ("4e", "Figure4e_core_programme_ORA.tsv", programme_ora,
     "One-sided Fisher over-representation analysis for four frozen programmes in the full DE universe."),
    ("4e", "Figure4e_gene_program_membership.tsv", network_all,
     "Eligible patient-blocked upregulated genes in four frozen programmes; displayed matrix genes are flagged."),
    ("4f", "Figure4f_preranked_gene_list.tsv", ranked_de,
     "Complete signed edgeR-QLF ranking used for the prespecified functional GSEA."),
    ("4f", "Figure4f_preranked_GSEA.tsv", gsea_results,
     "Ten-thousand-permutation preranked GSEA for all 14 frozen external functional programmes."),
    ("4f", "Figure4f_GSEA_leading_edge.tsv", gsea_leading_edge,
     "Leading-edge genes for every frozen functional programme."),
    ("4f", "Figure4f_programme_Jaccard.tsv", jaccard_long,
     "Pairwise tested-gene Jaccard similarities used for functional clustering."),
    ("4f", "Figure4f_function_clusters.tsv", cluster_table,
     "Average-linkage functional clusters and optimal leaf order."),
    ("4g-4h", "Figure4g_h_coarse_PLAT_program_residuals.tsv", coarse_residuals,
     "Patient-level adjusted residuals for coarse fibroblast PLAT and the two displayed myeloid programmes."),
]
manifest = []
for panel, filename, frame, description in tables:
    frame.to_csv(SOURCE / filename, sep="\t", index=False)
    manifest.append({"panel": panel, "file": filename, "description": description})
pd.DataFrame(manifest).to_csv(SOURCE / "Figure4_source_data_manifest.tsv", sep="\t", index=False)


# ---------------------------------------------------------------------------
# Compose the revised A-F figure.
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(179.5 / 25.4, 232 / 25.4))
outer = fig.add_gridspec(
    4, 12, height_ratios=[1.08, 0.86, 0.76, 0.86],
    left=0.105, right=0.988, top=0.975, bottom=0.060,
    wspace=1.02, hspace=0.72,
)


# a | Raw fractions: the most intuitive view of patient co-enrichment.
ax = fig.add_subplot(outer[0, 0:5])
panel_label(ax, "a", -0.17, 1.10)
x = features.sender_state_fraction.to_numpy(float) * 100
y = features.receiver_state_fraction.to_numpy(float) * 100
depth = np.sqrt(features.fibroblast_cells * features.myeloid_cells).to_numpy(float)
sizes = 24 + 92 * (np.sqrt(depth) - np.sqrt(depth).min()) / np.ptp(np.sqrt(depth))
ax.scatter(x, y, s=sizes, facecolor=BLUE_LIGHT, edgecolor=ROSE, linewidth=0.9,
           alpha=0.90, zorder=3)
descriptive_fit(ax, x, y, color="#657686", fill=BLUE_LIGHT, seed=2026080101)
ax.set_xlabel(r"$\it{PLAT}$-rich fibroblast state" "\n(% of fibroblasts)")
ax.set_ylabel(r"$\it{C1QC}^{+}\it{APOE}^{+}$ macrophage state" "\n(% of myeloid cells)")
ax.text(0.035, 0.965,
        rf"Spearman $\rho$ = {raw_rho:.2f}" "\n"
        rf"$P$ = {format_p_value(raw_p)}" "\n"
        "n = 26 patients",
        transform=ax.transAxes, ha="left", va="top", fontsize=6.4,
        bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="none", alpha=0.88))
clean_axis(ax)


# b | Prespecified adjusted log-ratio analysis; robust sensitivities are in legend/source data.
ax = fig.add_subplot(outer[0, 5:9])
panel_label(ax, "b", -0.20, 1.10)
ax.scatter(primary_x, primary_y, s=39, facecolor=BLUE, edgecolor=ROSE,
           linewidth=0.85, alpha=0.92, zorder=3)
descriptive_fit(ax, primary_x, primary_y, color="#586D80", fill=BLUE_LIGHT, seed=2026080102)
ax.axhline(0, color=LIGHT_GREY, linewidth=0.7, zorder=0)
ax.axvline(0, color=LIGHT_GREY, linewidth=0.7, zorder=0)
ax.set_xlabel("Adjusted fibroblast state")
ax.set_ylabel("Adjusted macrophage state", labelpad=1.5)
ax.text(0.035, 0.965,
        rf"adjusted $\rho$ = {primary_rho:.2f}" "\n"
        rf"partial-rank $P$ = {format_p_value(inference.iloc[0].correct_partial_asymptotic_p)}" "\n"
        "n = 26",
        transform=ax.transAxes, ha="left", va="top", fontsize=6.4,
        bbox=dict(boxstyle="round,pad=0.27", facecolor="white", edgecolor="none", alpha=0.88))
clean_axis(ax)


# c | Complete LOO direction stability.
ax = fig.add_subplot(outer[0:2, 9:12])
panel_label(ax, "c", -0.27, 1.055)
wide = loo.pivot(index="dropped_Patient_ID", columns="representation", values="adjusted_rho")
wide["mean"] = wide[["Log-ratio", "Arcsine fraction"]].mean(axis=1)
wide = wide.sort_values("mean")
yy = np.arange(len(wide))
ax.hlines(yy, wide["Log-ratio"], wide["Arcsine fraction"], color="#CDD4D9", linewidth=0.8)
ax.scatter(wide["Log-ratio"], yy, s=22, color=ROSE, edgecolor="white", linewidth=0.45,
           label="Log-ratio", zorder=3)
ax.scatter(wide["Arcsine fraction"], yy, s=22, marker="D", color=BLUE, edgecolor="white",
           linewidth=0.45, label="Arcsine fraction", zorder=3)
ax.axvline(0, color=GREY, linestyle=(0, (3, 2)), linewidth=0.75)
ax.set_yticks(yy)
ax.set_yticklabels([str(v).replace("NGCII", "") for v in wide.index], fontsize=5.8)
ax.set_xlim(-0.02, 0.67)
ax.set_xlabel(r"Adjusted $\rho$ after patient removal")
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.03, 0.99),
          fontsize=5.8, handletextpad=0.35, borderaxespad=0)
ax.tick_params(axis="y", length=0, pad=2)
clean_axis(ax)


# d | Patient-by-program heatmap with adjacent zero-origin lollipop summary.
dgrid = outer[1, 0:9].subgridspec(1, 2, width_ratios=[6.15, 1.35], wspace=0.035)
ax = fig.add_subplot(dgrid[0, 0])
panel_label(ax, "d", -0.115, 1.16)
pivot = core_delta.pivot(index="program", columns="Patient_ID",
                         values="delta_target_minus_other_recurrent").loc[PROGRAM_DISPLAY_ORDER]
pivot = pivot.loc[:, pivot.mean(axis=0).sort_values().index]
vmax = max(0.26, float(pivot.to_numpy().max()))
seq_cmap = LinearSegmentedColormap.from_list("positive_program", ["#FCFAEF", GREEN_LIGHT, "#789B80"])
im = ax.imshow(pivot.to_numpy(float), aspect="auto", cmap=seq_cmap,
               norm=Normalize(vmin=0, vmax=vmax), interpolation="nearest")
ax.set_yticks(np.arange(4))
ax.set_yticklabels([PROGRAM_LABELS[p] for p in PROGRAM_DISPLAY_ORDER], fontsize=6.0)
ax.set_xticks(np.arange(26))
ax.set_xticklabels([str(v).replace("NGCII", "") for v in pivot.columns],
                   rotation=90, fontsize=5.0)
ax.tick_params(axis="both", length=0, pad=1.8)
for edge in np.arange(0.5, 4, 1):
    ax.axhline(edge, color="white", linewidth=0.65)
for edge in np.arange(0.5, 26, 1):
    ax.axvline(edge, color="white", linewidth=0.18, alpha=0.7)
for spine in ax.spines.values():
    spine.set_visible(False)
cax = ax.inset_axes([0.70, -0.29, 0.27, 0.055])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_ticks([0, 0.15, 0.30])
cb.set_label("Target − comparator programme score", fontsize=5.5, labelpad=0.5)
cb.ax.tick_params(labelsize=5.8, length=2, pad=1)
cb.outline.set_linewidth(0.45)

axr = fig.add_subplot(dgrid[0, 1])
summary = contrasts.set_index("program").loc[PROGRAM_DISPLAY_ORDER]
yy = np.arange(4)
# The lollipop encodes one continuous quantity (median Δ). Use a single
# sequential green scale so darker colour means a larger programme difference.
lollipop_cmap = LinearSegmentedColormap.from_list(
    "programme_delta_green", ["#DCE8D8", "#9BB99B", "#4E7F58"]
)
lollipop_norm = Normalize(
    vmin=float(summary.median_delta.min()),
    vmax=float(summary.median_delta.max()),
)
colors = [lollipop_cmap(lollipop_norm(v)) for v in summary.median_delta]
axr.hlines(yy, 0, summary.median_delta, color=colors, linewidth=2.7)
axr.scatter(summary.median_delta, yy, s=38, color=colors, edgecolor="white", linewidth=0.55, zorder=3)
axr.spines["left"].set_position(("data", 0))
axr.spines["left"].set_color("#000000")
axr.spines["left"].set_linewidth(0.65)
axr.set_xlim(0, 0.225)
axr.set_ylim(-0.45, 3.45)
axr.invert_yaxis()
axr.set_yticks(yy)
axr.set_yticklabels([])
axr.set_xlabel("Median Δ", fontsize=5.8)
axr.spines["top"].set_visible(False)
axr.spines["right"].set_visible(False)
axr.tick_params(axis="x", labelsize=5.8, length=2)
axr.tick_params(axis="y", length=0)


# e–h | A taller, narrower lower block. The GSEA panel spans both rows on the
# right; the volcano and two patient-level bridges occupy the left column.
lower = outer[2:4, :].subgridspec(
    2, 2, width_ratios=[1.15, 1.0], height_ratios=[0.82, 1.0],
    # Keep the volcano panel's rightmost genes clear of the enlarged GSEA panel.
    wspace=0.68, hspace=0.72,
)

# e | Full DE result and clustered, full-ranked functional interpretation.
axv = fig.add_subplot(lower[0, 0])
panel_label(axv, "e", -0.24, 1.15)

nontechnical = ~volcano_plot["technical_excluded"].astype(bool)
target_enriched = (
    nontechnical
    & (volcano_plot["FDR"] <= 0.05)
    & (volcano_plot["logFC_target_minus_comparator"] >= 0.5)
)
comparator_enriched = (
    nontechnical
    & (volcano_plot["FDR"] <= 0.05)
    & (volcano_plot["logFC_target_minus_comparator"] <= -0.5)
)
background = ~(target_enriched | comparator_enriched)
axv.scatter(
    volcano_plot.loc[background, "logFC_target_minus_comparator"],
    volcano_plot.loc[background, "minus_log10_FDR_capped"],
    s=3.2, color="#C9CED2", alpha=0.40, edgecolors="none", label="Other", zorder=1,
)
axv.scatter(
    volcano_plot.loc[comparator_enriched, "logFC_target_minus_comparator"],
    volcano_plot.loc[comparator_enriched, "minus_log10_FDR_capped"],
    s=7.5, color=BLUE_LIGHT, alpha=0.86, edgecolors="none",
    label="Comparator-enriched", zorder=2,
)
axv.scatter(
    volcano_plot.loc[target_enriched, "logFC_target_minus_comparator"],
    volcano_plot.loc[target_enriched, "minus_log10_FDR_capped"],
    s=7.5, color=GREEN, alpha=0.83, edgecolors="none",
    label="C1QC+APOE+ enriched", zorder=2,
)
axv.axvline(-0.5, color="#B7BDC1", lw=0.65, ls=(0, (3, 2)))
axv.axvline(0.5, color="#B7BDC1", lw=0.65, ls=(0, (3, 2)))
axv.axhline(-np.log10(0.05), color="#B7BDC1", lw=0.65, ls=(0, (3, 2)))
for gene, offset in {
    "S100A12": (4, 3), "FCN1": (4, -8), "IL1B": (-22, 4),
    "C2": (4, 3), "TREM2": (-25, 6), "MMP12": (4, -7), "MMP14": (4, 5),
}.items():
    row = volcano_plot[volcano_plot["gene_symbol"].eq(gene)].iloc[0]
    axv.annotate(
        gene,
        (row["logFC_target_minus_comparator"], row["minus_log10_FDR_capped"]),
        xytext=offset, textcoords="offset points", fontsize=5.0, color=INK,
        arrowprops=dict(arrowstyle="-", color="#929AA0", lw=0.45), zorder=4,
        fontstyle="italic",
    )
axv.set_xlim(-7.5, 5.5)
axv.set_ylim(0, 30.8)
axv.set_xlabel(r"log$_2$ fold change" "\n" r"($\it{C1QC}^{+}\it{APOE}^{+}$ − other myeloid states)")
axv.set_ylabel(r"−log$_{10}$(FDR), capped at 30")
axv.legend(frameon=False, loc="upper right", bbox_to_anchor=(0.985, 0.985),
           fontsize=5.8, markerscale=1.15, handletextpad=0.25,
           borderaxespad=0.15, ncol=1)
clean_axis(axv)

# f | Clustered full-ranked functional interpretation; enlarged across two rows.
axg = fig.add_subplot(lower[:, 1])
panel_label(axg, "f", -0.22, 1.055)
gsea_plot = gsea_results.copy().reset_index(drop=True)
gsea_y = np.arange(len(gsea_plot))[::-1]
theme_fill = {
    "Inflammatory signalling": "#F7E9E8",
    "Migration": "#E9F0F5",
    "Hypoxia": "#F0EDF4",
    "ECM / adhesion": "#ECF3E9",
    "Phagocytosis": "#E8EFF6",
    "Integrin signalling": "#F0F1F2",
}
for theme, frame in gsea_plot.groupby("functional_theme", sort=False):
    indices = frame.index.to_numpy(int)
    theme_y = gsea_y[indices]
    axg.axhspan(theme_y.min() - 0.46, theme_y.max() + 0.46,
                color=theme_fill[theme], alpha=0.70, lw=0, zorder=0)
    axg.text(2.72, theme_y.mean(), theme, fontsize=5.0, color="#000000",
             ha="right", va="center", fontweight="bold")
sizes = 18 + 1.25 * gsea_plot["leading_edge_size"].to_numpy(float)
point_colors = np.where(gsea_plot["NES"] >= 0, GREEN, BLUE)
significant = gsea_plot["FDR_BH_14_programmes"].to_numpy(float) < 0.05
axg.scatter(
    gsea_plot["NES"], gsea_y, s=sizes, color=point_colors,
    edgecolor=np.where(significant, INK, "white"),
    linewidth=np.where(significant, 0.95, 0.45), zorder=3,
)
axg.axvline(0, color="#8F989E", lw=0.75, ls=(0, (3, 2)), zorder=1)
axg.set_xlim(-3.0, 2.85)
axg.set_ylim(-0.62, len(gsea_plot) - 0.38)
axg.set_xticks([-3, -2, -1, 0, 1, 2])
axg.set_yticks(gsea_y)
axg.set_yticklabels(gsea_plot["display_label"], fontsize=5.6)
axg.set_xlabel("Preranked GSEA normalized enrichment score")
axg.grid(axis="x", color="#E9ECEE", lw=0.45, zorder=0)
clean_axis(axg)


# g/h | Direct patient-level scatterplots for the coarse PLAT bridge.
gh_frame = fig.add_subplot(lower[1, 0], frameon=False)
gh_frame.set_xticks([])
gh_frame.set_yticks([])
for spine in gh_frame.spines.values():
    spine.set_visible(False)
gh_frame.set_xlabel(r"Fibroblast $\it{PLAT}$" "\n" "covariate-adjusted relative level", labelpad=23)
gh_grid = lower[1, 0].subgridspec(1, 2, wspace=0.48)
for j, (program, color, fill) in enumerate([
    ("GOBP_PHAGOCYTOSIS", BLUE, BLUE_LIGHT),
    ("GOBP_MYELOID_LEUKOCYTE_MIGRATION", GREEN, GREEN_LIGHT),
]):
    ax = fig.add_subplot(gh_grid[0, j])
    # Keep the panel letters clearly above the statistics block and plotting area.
    panel_label(ax, "g" if j == 0 else "h", -0.30, 1.28)
    dat = coarse_residuals[coarse_residuals.program == program]
    x = dat.PLAT_rank_residual.to_numpy(float)
    y = dat.program_rank_residual.to_numpy(float)
    ax.scatter(x, y, s=34, facecolor=fill, edgecolor=ROSE, linewidth=0.75, alpha=0.92, zorder=3)
    descriptive_fit(ax, x, y, color=color, fill=fill, seed=2026080110 + j)
    ax.axhline(0, color=LIGHT_GREY, linewidth=0.65)
    ax.axvline(0, color=LIGHT_GREY, linewidth=0.65)
    frozen = assoc_focus[assoc_focus.program == program].iloc[0]
    # Place the statistic summary in the gutter above the axes so it does not
    # cover patient points in the upper-left quadrant.
    ax.text(0.04, 1.035,
            rf"adjusted $\rho$ = {frozen.adjusted_rho:.2f}" "\n"
            rf"permutation $P$ = {format_p_value(frozen.adjusted_permutation_p)}" "\n"
            f"n = {int(frozen.n)}",
            transform=ax.transAxes, fontsize=5.8, color=INK, ha="left", va="bottom",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="none", alpha=0.88))
    ax.set_xlabel("")
    if j == 0:
        ax.set_ylabel("Phagocytosis\nadjusted relative level")
    else:
        ax.set_ylabel("Migration\nadjusted relative level")
    clean_axis(ax)


save_figure(fig, "Figure4_patient_state_association")
print("Revised Figure 4 generated")
print("raw rho", raw_rho, "adjusted rho", primary_rho)
print("source data", SOURCE)
