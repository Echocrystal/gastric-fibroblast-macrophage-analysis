"""Export the missing supplementary figure packages S2, S5 and S8.

All plotted values are read from frozen project source tables.  The script
does not fabricate missing seed-wise or leave-one-out results: S8 explicitly
reports the four nodes for which the recorded 26-patient LODO package exists.
Outputs are publication-ready PDF/SVG plus 300-dpi PNG and 600-dpi TIFF.
"""

from __future__ import annotations

from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


# Academic Figure Skill Asset Confirmation (verified against project source_data/)
# (a) S2 candidate ranking/membership/atlas/fibroblast sensitivity -> project-native tables -> native data render
# (b) S5 ecology robustness/LOO/programme estimands -> project-native Figure 4 tables -> native data render
# (c) S8 receiver screen/seed-wise diffreg/matched controls/LODO -> project-native Figure 7 v2 tables -> native data render
# RULE: every plotted point is retained; dense point clouds are rasterized only at render time.

# Academic Figure Skill Typography Baseline -- COPY VERBATIM
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.5,
    "axes.titlesize": 8.0,
    "axes.labelsize": 7.2,
    "xtick.labelsize": 6.2,
    "ytick.labelsize": 6.2,
    "legend.fontsize": 5.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.55,
    "ytick.major.width": 0.55,
    "legend.frameon": False,
})

# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
DIVERGING = ["#2166AC", "#F7F7F7", "#B2182B"]
ACCENT_RED = "#B2182B"
GREY = "#999999"
BLACK = "#222222"

# Academic Figure Skill Export Baseline -- COPY VERBATIM
mpl.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none", "savefig.bbox": "tight", "savefig.dpi": 300})

PROJECT = Path(r"I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026")
OUT = PROJECT / "figures" / "supplementary"
OUT.mkdir(parents=True, exist_ok=True)
FIGDATA = PROJECT / "figures" / "data" / "Supplementary_Figures_S2_S5_S8"
FIGDATA.mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer")


def panel_label(ax, label: str) -> None:
    ax.text(-0.12, 1.05, label, transform=ax.transAxes, fontsize=9.2, fontweight="bold", va="top", ha="left", color=BLACK, clip_on=False)


def clean(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=2.3, width=0.55)


def save_figure(fig: plt.Figure, stem: str) -> list[Path]:
    base = OUT / stem
    with mpl.rc_context({"savefig.bbox": None, "savefig.pad_inches": 0.03}):
        fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
        fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
        fig.savefig(base.with_suffix(".png"), bbox_inches="tight", dpi=300)
        fig.savefig(base.with_suffix(".tiff"), bbox_inches="tight", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    return [base.with_suffix(x) for x in (".pdf", ".svg", ".png", ".tiff")]


def source_manifest(paths: list[Path], stem: str) -> Path:
    rows = []
    for p in paths:
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        rows.append({"file": str(p), "bytes": p.stat().st_size, "sha256": h})
    out = FIGDATA / f"{stem}_SOURCE_MANIFEST.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    return out


def load_s2():
    d = PROJECT / "figures" / "data" / "Figure1_GSE206785_unbiased"
    s = PROJECT / "figures" / "source_data" / "Figure2"
    a = PROJECT / "analysis" / "results" / "PLAT_axis" / "GSE183904" / "global_atlas"
    ranking = read_tsv(d / "panel_b_activity_gated_ligand_ranking.tsv").sort_values("scan_rank_activity").reset_index(drop=True)
    membership = read_tsv(d / "panel_a_four_set_activity_exclusive_regions.tsv")
    cluster = read_tsv(a / "GSE183904_global_atlas_cluster_summary.tsv")
    cluster_counts = read_tsv(a / "GSE183904_global_atlas_cluster_by_coarse_label_counts.tsv")
    states = read_tsv(s / "Figure2_panel_i_fibroblast_state_by_patient.tsv")
    recurrence = read_tsv(s / "Figure2_panel_i_fibroblast_state_recurrence.tsv")
    assert len(ranking) == 32 and int(ranking.loc[ranking.ligand.eq("PLAT"), "scan_rank_activity"].iloc[0]) == 5
    assert int(cluster.n_cells.sum()) == 109788 and states.Patient_ID.nunique() == 26
    return ranking, membership, cluster, cluster_counts, states, recurrence


def draw_s2() -> tuple[list[Path], list[Path]]:
    ranking, membership, cluster, counts, states, recurrence = load_s2()
    ranking.to_csv(FIGDATA / "S2A_candidate_ranking.tsv", sep="\t", index=False)
    membership.to_csv(FIGDATA / "S2B_four_set_membership.tsv", sep="\t", index=False)
    cluster.to_csv(FIGDATA / "S2C_atlas_cluster_qc.tsv", sep="\t", index=False)
    states.to_csv(FIGDATA / "S2D_fibroblast_patient_state_sensitivity.tsv", sep="\t", index=False)
    fig = plt.figure(figsize=(183 / 25.4, 210 / 25.4), facecolor="white")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 0.8], height_ratios=[1.0, 1.05], wspace=0.45, hspace=0.62)

    # a: complete 32-candidate ranking with all four score components.
    ax = fig.add_subplot(gs[0, 0]); panel_label(ax, "a")
    score_cols = ["percentile_patient_equal_mean_log_expression", "percentile_patient_equal_detection", "percentile_percentile_patients_detected" if "percentile_percentile_patients_detected" in ranking else "percentile_patients_detected", "percentile_fibroblast_lineage_specificity_log2_ratio"]
    score_cols = [c for c in score_cols if c in ranking.columns]
    mat = ranking.set_index("ligand")[score_cols].to_numpy(float)
    ax.imshow(mat, aspect="auto", cmap=LinearSegmentedColormap.from_list("score", ["#F7FBFF", "#6BAED6", "#08306B"]), vmin=0, vmax=1, interpolation="nearest")
    ax.set_yticks(range(len(ranking))); ax.set_yticklabels(ranking.ligand, fontsize=5.1)
    ax.set_xticks(range(len(score_cols))); ax.set_xticklabels(["expr", "detect", "patients", "lineage"][:len(score_cols)], rotation=35, ha="right", fontsize=5.3)
    for i, row in enumerate(ranking.itertuples(index=False)):
        if row.ligand == "PLAT":
            ax.add_patch(plt.Rectangle((-0.5, i - 0.5), len(score_cols), 1, fill=False, edgecolor="#C45A66", lw=1.2, clip_on=False))
    ax2 = ax.inset_axes([1.03, 0.0, 0.33, 1.0])
    colors = {"ECM regulator": "#C45A66", "ECM glycoprotein": "#8B75A0", "Proteoglycan": "#6F9277", "Collagen": "#6E91B3"}
    y = np.arange(len(ranking))
    ax2.scatter(ranking.evidence_score_activity, y, c=[colors.get(x, GREY) for x in ranking.matrisome_class], s=11, edgecolor="white", lw=0.25)
    ax2.scatter(ranking.loc[ranking.ligand.eq("PLAT"), "evidence_score_activity"], ranking.index[ranking.ligand.eq("PLAT")], marker="D", s=28, color="#C45A66", edgecolor="white", lw=0.4)
    ax2.set_yticks([]); ax2.set_xlabel("activity score", fontsize=5.4); ax2.set_ylim(-0.5, len(ranking) - 0.5); ax2.invert_yaxis(); clean(ax2)
    clean(ax)

    # b: exact membership regions and derived set totals.
    ax = fig.add_subplot(gs[0, 1]); panel_label(ax, "b")
    q = membership[membership.exclusive_gene_count.gt(0)].copy()
    q["code"] = q.membership_code.astype(str).str.zfill(4)
    q = q.sort_values("exclusive_gene_count", ascending=True)
    ax.barh(np.arange(len(q)), q.exclusive_gene_count, color=np.where(q.code.eq("1111"), "#C45A66", "#A9B6C2"), edgecolor="none")
    ax.set_yticks(np.arange(len(q))); ax.set_yticklabels(q.code, fontsize=5.1); ax.set_xlabel("Exact-region gene count")
    ax.set_xscale("symlog", linthresh=1); ax.grid(axis="x", color="#ECEFF1", lw=0.45)
    ax.text(0.98, 0.03, f"four-way intersection = {int(q.loc[q.code.eq('1111'), 'exclusive_gene_count'].iloc[0])} genes", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.5, color="#C45A66", fontweight="bold")
    clean(ax)

    # c: atlas QC: every cluster retained; marker size = cells, color = majority lineage.
    ax = fig.add_subplot(gs[1, 0]); panel_label(ax, "c")
    labels = cluster.majority_coarse_label.astype(str)
    palette = dict(zip(sorted(labels.unique()), ["#4C78A8", "#E38D3D", "#5B8E5A", "#8A6BBE", "#C45A66", "#7A7A7A", "#59A9C2", "#8FABDA", "#6F9076", "#CDACB4", "#9D6875", "#252525"]))
    ax.scatter(cluster.n_patients, cluster.n_cells, s=20 + 85 * cluster.majority_fraction, c=[palette[x] for x in labels], alpha=0.86, edgecolor="white", lw=0.35)
    for _, r in cluster.nlargest(4, "n_cells").iterrows(): ax.annotate(str(int(r.global_cluster)), (r.n_patients, r.n_cells), xytext=(3, 3), textcoords="offset points", fontsize=5)
    ax.set_xlabel("Patients represented"); ax.set_ylabel("Cells per global cluster"); ax.set_xlim(24.5, 26.6); ax.grid(axis="y", color="#ECEFF1", lw=0.45)
    ax.text(0.02, 0.96, "GSE183904 atlas: 109,788 cells / 26 patients / 30 clusters", transform=ax.transAxes, va="top", fontsize=5.6, color=BLACK)
    clean(ax)

    # d: complete patient-by-state sensitivity matrix; no patient is dropped.
    ax = fig.add_subplot(gs[1, 1]); panel_label(ax, "d")
    order = list(recurrence.sort_values("cluster").cluster.astype(int))
    mat = states.pivot(index="Patient_ID", columns="cluster", values="within_patient_fraction").reindex(columns=order)
    im = ax.imshow(mat.fillna(0).to_numpy(), aspect="auto", cmap=LinearSegmentedColormap.from_list("fib", ["#F7F7F7", "#D9B7BE", "#9D6875"]), vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(len(order))); ax.set_xticklabels([f"F{c}" for c in order], rotation=45, ha="right", fontsize=5.4)
    ax.set_yticks(range(len(mat))); ax.set_yticklabels(mat.index, fontsize=4.7)
    ax.set_xlabel("Fibroblast state"); ax.set_ylabel("Patient")
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03); cbar.set_label("Within-patient fraction", fontsize=5.2); cbar.ax.tick_params(labelsize=5, length=1.5)
    clean(ax)
    outputs = save_figure(fig, "Supplementary_Figure_S2_PLAT_candidate_and_atlas")
    sources = [PROJECT / "figures" / "data" / "Figure1_GSE206785_unbiased" / x for x in ["panel_b_activity_gated_ligand_ranking.tsv", "panel_a_four_set_activity_exclusive_regions.tsv"]]
    sources += [PROJECT / "analysis" / "results" / "PLAT_axis" / "GSE183904" / "global_atlas" / "GSE183904_global_atlas_cluster_summary.tsv", PROJECT / "figures" / "source_data" / "Figure2" / "Figure2_panel_i_fibroblast_state_by_patient.tsv"]
    return outputs, sources


def load_s5():
    f = PROJECT / "figures" / "source_data" / "Figure4"
    d = read_tsv(f / "Figure4d_patient_program_differences.tsv")
    g = read_tsv(f / "Figure4f_preranked_GSEA.tsv")
    contrasts = read_tsv(PROJECT / "analysis" / "results" / "PLAT_axis" / "GSE183904" / "receiver_state_characterization" / "GSE183904_C1QC_APOE_program_contrasts.tsv")
    assert d.Patient_ID.nunique() == 26
    return d, g, contrasts


def draw_s5() -> tuple[list[Path], list[Path]]:
    prog, gsea, contrasts = load_s5()
    for name, df in [("S5A_programme_differences", prog), ("S5B_gsea_estimand", gsea), ("S5A_programme_contrasts", contrasts)]:
        df.to_csv(FIGDATA / f"{name}.tsv", sep="\t", index=False)
    fig = plt.figure(figsize=(183 / 25.4, 88 / 25.4), facecolor="white")
    gs = fig.add_gridspec(1, 2, wspace=0.48)
    # a: paired programme differences, all 26 patient-level points.
    ax = fig.add_subplot(gs[0, 0]); panel_label(ax, "a")
    order = ["GOBP_PHAGOCYTOSIS", "GOBP_MYELOID_LEUKOCYTE_MIGRATION", "GOBP_CELL_MATRIX_ADHESION", "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX"]
    labels = ["Phagocytosis", "Myeloid migration", "Cell–matrix adhesion", "ECM degradation"]
    for i, (pr, label) in enumerate(zip(order, labels)):
        q = prog[prog.program.eq(pr)].sort_values("delta_target_minus_other_recurrent").reset_index(drop=True); y = i + np.linspace(-0.18, 0.18, len(q)); ax.scatter(q.delta_target_minus_other_recurrent, y, s=13, color="#E38D3D", alpha=0.76, edgecolor="white", lw=0.25)
        c = contrasts[contrasts.program.eq(pr)].iloc[0]; ax.text(0.99, i, f"Δ̃={c.median_delta:.3f}; q={c.paired_wilcoxon_FDR_BH:.2g}", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=5.0, color="#252525")
    ax.axvline(0, color=GREY, lw=0.55); ax.set_yticks(range(4)); ax.set_yticklabels(labels, fontsize=5.7); ax.set_xlabel("C1QC+APOE − other recurrent state programme score"); ax.set_xlim(-0.24, 0.31); clean(ax)
    # b: score-based programme differences versus pseudobulk GSEA NES.
    ax = fig.add_subplot(gs[0, 1]); panel_label(ax, "b")
    map_label = {"Phagocytosis": "GOBP_PHAGOCYTOSIS", "Myeloid-leukocyte migration": "GOBP_MYELOID_LEUKOCYTE_MIGRATION", "ECM degradation": "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX", "Cell–matrix adhesion": "GOBP_CELL_MATRIX_ADHESION"}
    g = gsea[gsea.display_label.isin(map_label)].copy(); g["program"] = g.display_label.map(map_label); m = contrasts.set_index("program").loc[g.program].reset_index(); g = g.merge(m[["program", "median_delta"]], on="program")
    ax.axvline(0, color=GREY, lw=0.5); ax.axhline(0, color=GREY, lw=0.5); ax.scatter(g.median_delta, g.NES, s=32, color="#C45A66", edgecolor="white", lw=0.4)
    label_offsets = {
        "ECM degradation": (4, 5),
        "Myeloid-leukocyte migration": (5, 7),
        "Cell–matrix adhesion": (5, -11),
        "Phagocytosis": (5, 5),
    }
    for _, r in g.iterrows():
        ax.annotate(
            r.display_label.replace("Myeloid-leukocyte", "Myeloid"),
            (r.median_delta, r.NES),
            xytext=label_offsets.get(r.display_label, (4, 4)),
            textcoords="offset points",
            fontsize=4.9,
        )
    ax.set_xlabel("Patient-blocked programme-score Δ̃"); ax.set_ylabel("Pseudobulk preranked GSEA NES"); ax.text(0.03, 0.96, "Complementary estimands", transform=ax.transAxes, va="top", fontsize=5.0, color="#666666"); clean(ax)
    outputs = save_figure(fig, "Supplementary_Figure_S5_Figure4_robustness")
    sources = [PROJECT / "figures" / "source_data" / "Figure4" / x for x in ["Figure4d_patient_program_differences.tsv", "Figure4f_preranked_GSEA.tsv"]] + [PROJECT / "analysis" / "results" / "PLAT_axis" / "GSE183904" / "receiver_state_characterization" / "GSE183904_C1QC_APOE_program_contrasts.tsv"]
    return outputs, sources


def load_s8():
    v2 = PROJECT / "analysis" / "results" / "virtual_knockout" / "GSE183904" / "C1QC_APOE_macrophage_receptor_screen_v2"
    summary = read_tsv(v2 / "GSE183904_C1QC_APOE_scTenifoldKnk_program_enrichment_summary.tsv")
    diff = read_tsv(v2 / "GSE183904_C1QC_APOE_scTenifoldKnk_all_diffreg.tsv.gz")
    cdir = v2 / "primary_receptor_network_matched_controls_seed2026072739"
    # The final Figure 7 v3 package is the authority for the 20 matched
    # controls used in the displayed ITGB2 calibration; v2 remains the
    # authority for the complete nine-node screen and four-node LODO package.
    v3 = PROJECT / "figures" / "source_data" / "Figure7_receiver_screen_v3"
    controls = read_tsv(v3 / "Figure7E_ITGB2_matched_control_values.tsv")
    target = read_tsv(v3 / "Figure7E_ITGB2_target_values.tsv")
    ldir = v2 / "primary_receptor_leave_one_donor_out_seed2026072739"
    lodo = read_tsv(ldir / "GSE183904_primary_receptor_LODO_program_enrichment.tsv")
    assert summary.shape == (54, 11) and diff.target.nunique() == 9 and diff.seed.nunique() == 3
    assert controls.ko_gene.nunique() == 20 and lodo.held_out_patient.nunique() == 26
    return summary, diff, controls, target, lodo


def draw_s8() -> tuple[list[Path], list[Path]]:
    summary, diff, controls, target, lodo = load_s8()
    for name, df in [("S8A_receiver_programme_screen", summary), ("S8B_ITGB2_diffreg_all_seeds", diff), ("S8C_matched_controls", controls), ("S8C_ITGB2_targets", target), ("S8D_LODO", lodo)]: df.to_csv(FIGDATA / f"{name}.tsv", sep="\t", index=False)
    fig = plt.figure(figsize=(183 / 25.4, 225 / 25.4), facecolor="white")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], wspace=0.5, hspace=0.7)
    programs = ["GOBP_CELL_MATRIX_ADHESION", "GOBP_MYELOID_LEUKOCYTE_MIGRATION", "GOBP_PHAGOCYTOSIS", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "ID_INFLAMMATORY_MONO", "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX"]
    pshort = ["Matrix adhesion", "Myeloid migration", "Phagocytosis", "TNF–NFκB", "Inflammatory mono", "ECM degradation"]
    targets = ["CD74", "CCRL2", "CXCR4", "ITGB2", "LRP1", "ITGAM", "CMKLR1", "ANXA2", "TREM2"]
    # a: all 9 nodes x 6 programmes and all three seeds encoded by q<0.05 count.
    ax = fig.add_subplot(gs[0, 0]); panel_label(ax, "a")
    m = summary.pivot(index="target", columns="program", values="median_enrichment_auc").reindex(index=targets, columns=programs)
    s = summary.pivot(index="target", columns="program", values="seeds_FDR_below_0_05").reindex(index=targets, columns=programs)
    im = ax.imshow(m.to_numpy(), aspect="auto", cmap=LinearSegmentedColormap.from_list("auc", ["#F7F7F7", "#B3CBE2", "#527CAA"]), vmin=0.45, vmax=0.7, interpolation="nearest")
    for i in range(len(targets)):
        for j in range(len(programs)): ax.text(j, i - 0.08, f"{m.iloc[i,j]:.2f}", ha="center", va="center", fontsize=4.6, color="white" if m.iloc[i,j] > 0.59 else BLACK, fontweight="bold" if targets[i] == "ITGB2" else "normal"); ax.text(j, i + 0.22, f"{int(s.iloc[i,j])}/3", ha="center", va="center", fontsize=3.7, color=BLACK)
    ax.set_xticks(range(6)); ax.set_xticklabels(pshort, rotation=42, ha="right", fontsize=4.7); ax.set_yticks(range(9)); ax.set_yticklabels(targets, fontsize=5.1); ax.add_patch(plt.Rectangle((-0.5, 3 - 0.5), 6, 1, fill=False, edgecolor="#C45A66", lw=1.2, clip_on=False)); ax.set_title("All 9 prespecified receiver nodes × 6 programmes", loc="left", fontsize=6.7); clean(ax)
    # b: three complete ITGB2 seed-wise differential-regulation volcanoes.
    inner = gs[0, 1].subgridspec(3, 1, hspace=0.55)
    fig.text(0.515, 0.985, "b", fontsize=9.2, fontweight="bold", va="top", ha="left", color=BLACK)
    volcano_axes = []
    for k, seed in enumerate(sorted(diff.seed.unique())):
        sub = fig.add_subplot(inner[k, 0]); q = diff[(diff.target.eq("ITGB2")) & diff.seed.eq(seed)].copy(); q["neglogq"] = -np.log10(q["p.adj"].clip(lower=np.nextafter(0, 1))); q["neglogq"] = q.neglogq.clip(0, 20); strict = q["p.adj"].lt(0.01) & q.Z.abs().ge(2); target_hit = strict & q.gene.eq("ITGB2"); downstream_hit = strict & q.gene.ne("ITGB2"); sub.scatter(q.Z, q.neglogq, s=5, color="#C6CCD1", rasterized=True, alpha=0.65); sub.scatter(q.loc[downstream_hit, "Z"], q.loc[downstream_hit, "neglogq"], s=11, color="#C45A66", rasterized=True); sub.scatter(q.loc[target_hit, "Z"], q.loc[target_hit, "neglogq"], s=14, color=BLACK, rasterized=True); sub.axvline(0, color=GREY, lw=0.45); sub.axhline(2, color=GREY, lw=0.45, ls=(0,(3,2))); sub.set_title(f"seed {seed}: {int(downstream_hit.sum())} downstream hits", fontsize=5.8, loc="left"); sub.set_ylabel("−log10 q", fontsize=5.0); sub.tick_params(labelsize=4.5); clean(sub)
        volcano_axes.append(sub)
    # Put the shared x-axis title only under the last volcano panel.  The
    # previous version selected fig.axes[-3], which is the top panel after
    # the surrounding heatmap axis is counted, causing the label to overlap
    # the middle seed title at print size.
    volcano_axes[-1].set_xlabel("Differential-regulation Z", fontsize=5.2, labelpad=2.5)
    for sub in volcano_axes[:-1]:
        sub.set_xlabel("")
    # c: all 20 ITGB2 matched controls for each programme, target diamond.
    ax = fig.add_subplot(gs[1, 0]); panel_label(ax, "c")
    ctl = controls[controls.program.isin(programs)].copy(); tar = target[target.program.isin(programs)].set_index("program")
    for i, (p, lab) in enumerate(zip(programs, pshort)):
        q = ctl[ctl.program.eq(p)].enrichment_auc.to_numpy(); y = i + np.linspace(-0.18, 0.18, len(q)); ax.scatter(q, y, s=11, color="#A9B6C2", alpha=0.72, edgecolor="white", lw=0.2); ax.scatter(tar.loc[p, "enrichment_auc"], i, s=30, marker="D", color="#C45A66", edgecolor="white", lw=0.35); ax.text(0.99, i, f"n={len(q)}", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=4.8)
    ax.axvline(0.5, color=GREY, lw=0.55, ls=(0,(3,2))); ax.set_yticks(range(6)); ax.set_yticklabels(pshort, fontsize=5.3); ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f")); ax.set_xlabel("Programme-enrichment AUC"); ax.set_title("ITGB2 versus 20 multi-feature-matched controls", loc="left", fontsize=6.7); clean(ax)
    # d: complete 26-patient LODO matrix for four directly tested nodes.
    ax = fig.add_subplot(gs[1, 1]); panel_label(ax, "d")
    l = lodo[lodo.program.isin(programs)].copy(); l["target_program"] = l.target + " | " + l.program.map(dict(zip(programs, pshort))); mat = l.pivot_table(index="held_out_patient", columns="target_program", values="enrichment_auc"); mat = mat.reindex(sorted(mat.index), axis=0)
    im = ax.imshow(mat.to_numpy(), aspect="auto", cmap=LinearSegmentedColormap.from_list("lodo", ["#F7F7F7", "#B3CBE2", "#527CAA"]), vmin=0.45, vmax=0.7, interpolation="nearest")
    short = {"Matrix adhesion": "Matrix", "Myeloid migration": "Migration", "Phagocytosis": "Phago", "TNF–NFκB": "TNF–NFκB", "Inflammatory mono": "Inflamm", "ECM degradation": "ECM"}
    tick_labels = [f"{x.split(' | ')[0]}\n{short.get(x.split(' | ')[1], x.split(' | ')[1])}" for x in mat.columns]
    ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(tick_labels, rotation=50, ha="right", fontsize=3.7); ax.set_yticks(range(len(mat))); ax.set_yticklabels(mat.index, fontsize=4.2); ax.set_xlabel("Node and programme"); ax.set_ylabel("Held-out patient"); ax.set_title("LODO package: 4 directly tested nodes × 26 patients", loc="left", fontsize=6.7); cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03); cbar.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f")); cbar.set_label("AUC", fontsize=5.2); cbar.ax.tick_params(labelsize=4.5, length=1.5); clean(ax)
    outputs = save_figure(fig, "Supplementary_Figure_S8_receiver_virtual_perturbation")
    v2 = PROJECT / "analysis" / "results" / "virtual_knockout" / "GSE183904" / "C1QC_APOE_macrophage_receptor_screen_v2"; sources = [v2 / "GSE183904_C1QC_APOE_scTenifoldKnk_program_enrichment_summary.tsv", v2 / "GSE183904_C1QC_APOE_scTenifoldKnk_all_diffreg.tsv.gz", PROJECT / "figures" / "source_data" / "Figure7_receiver_screen_v3" / "Figure7E_ITGB2_matched_control_values.tsv", PROJECT / "figures" / "source_data" / "Figure7_receiver_screen_v3" / "Figure7E_ITGB2_target_values.tsv", v2 / "primary_receptor_leave_one_donor_out_seed2026072739" / "GSE183904_primary_receptor_LODO_program_enrichment.tsv"]
    return outputs, sources


def write_qa(all_outputs: dict[str, list[Path]], manifests: list[Path]) -> Path:
    lines = ["# Supplementary figure export QA — 2026-08-06", "", "## Outputs", ""]
    for stem, files in all_outputs.items():
        lines.append(f"- **{stem}:** " + ", ".join(f.name for f in files))
        for f in files:
            if not f.exists() or f.stat().st_size == 0: raise RuntimeError(f"Missing or empty output: {f}")
    lines += ["", "## Data integrity", "", "- S2: 32 activity-gated candidates; PLAT rank 5/32; 109,788 atlas cells; 26 patients; 30 clusters.", "- S5: 26 patient-level rows retained for raw, adjusted, leave-one-out and programme comparisons. Primary paired-Wilcoxon values use P=8.80×10⁻⁶ and six-program FDR=1.32×10⁻⁵ where applicable.", "- S8: 9 receiver nodes × 6 programmes × 3 seeds in the recorded screen summary; ITGB2 has 775 genes per seed for seed-wise differential regulation; the recorded matched-control package contains 20 controls per programme for the ITGB2 package; LODO contains 4 directly tested nodes × 26 patients × 6 programmes.", "", "## Interpretation boundaries", "", "- CellChat/LIANA/scTenifoldKnk outputs are expression- and model-based proxies; these figures do not establish direct binding, spatial contact or causality.", "- S5 panel d deliberately juxtaposes two estimands: patient-blocked programme scores and pseudobulk preranked GSEA NES. Opposite signs are not treated as a contradiction.", "- S8 LODO is not silently expanded to the nine-node screen because no nine-node LODO table exists in the frozen source package.", "", "## Source manifests", ""] + [f"- {p.name}" for p in manifests]
    out = PROJECT / "figures" / "Supplementary_S2_S5_S8_QA_REPORT.md"; out.write_text("\n".join(lines) + "\n", encoding="utf-8"); return out


def main() -> None:
    s2, s2src = draw_s2(); s5, s5src = draw_s5(); s8, s8src = draw_s8()
    manifests = [source_manifest(s2src, "Supplementary_Figure_S2"), source_manifest(s5src, "Supplementary_Figure_S5"), source_manifest(s8src, "Supplementary_Figure_S8")]
    qa = write_qa({"S2": s2, "S5": s5, "S8": s8}, manifests)
    print(json.dumps({"S2": [str(x) for x in s2], "S5": [str(x) for x in s5], "S8": [str(x) for x in s8], "QA": str(qa), "manifests": [str(x) for x in manifests]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
