"""Render the readability-revised spatial-chain main Figure 5.

The figure combines the independent main-ST analysis with same-donor ST_R
repeated-section sensitivity.  ST_R is never presented as an independent
cohort.  Visium spots remain mixed-cell capture units.
"""

from __future__ import annotations

import gzip
import importlib.util
import io
import json
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import EllipseCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize, PowerNorm
from matplotlib.lines import Line2D
from PIL import Image

from paper_plot_style import clean_axis, panel_label, save_figure


PROJECT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
PUBLIC = Path(r"I:\ECM-new\ECM-public")
RESULT = PROJECT / "analysis" / "results" / "spatial" / "GSE270678_PLAT_ecosystem_v1"
STR_RESULT = RESULT / "ST_R_sensitivity"
RAW_STR = PUBLIC / "data" / "raw" / "spatial" / "GSE270678" / "ST_R_public_files"
SPOT_STR = (
    PUBLIC / "data" / "processed" / "spatial" / "GSE270678" / "v4"
    / "GSE270678_ST_R_project_compartment_spot_scores_v4.tsv.gz"
)
MANIFEST = PUBLIC / "manifests" / "v4" / "GSE270678_spatial_sample_manifest.tsv"
BASE_ANALYSIS = PROJECT / "scripts" / "35_analyze_GSE270678_PLAT_spatial_chain.py"
FIG_DATA = PROJECT / "figures" / "data"
FIG_DATA.mkdir(parents=True, exist_ok=True)
# The enriched cross-cohort artwork is the promoted Figure 6. Keep the
# reproducibility switch, but make the promoted version the default so a
# normal rerun cannot silently recreate the superseded compact export.
ENRICHED_CANDIDATE = os.environ.get("FIG6_ENRICHED_CANDIDATE", "1") == "1"
GSE251950_META = (
    PROJECT / "analysis" / "results" / "spatial" / "GSE251950_PLAT_frozen_chain_v1"
    / "GSE251950_PLAT_frozen_chain_patient_meta.tsv"
)

SPEC = importlib.util.spec_from_file_location("plat_spatial_analysis", BASE_ANALYSIS)
M = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(M)

PLAT = "#9D5F6B"
SENDER = "#B77B89"
RECEIVER = "#527CAA"
ECM = "#6F9076"
ADH = "#6A91B5"
INK = "#000000"
GREY = "#7D878E"
LIGHT = "#D7DEE2"

# Panel-a sequential maps deliberately keep zero/low values neutral while making
# the mid-range visible at final print size.  The values and colour-bar limits
# are unchanged; PowerNorm below only redistributes colour contrast.
PLAT_CMAP = LinearSegmentedColormap.from_list(
    "plat_spatial", ["#D7DADD", "#E5B1B8", "#C57F8D", "#8F4F5E"]
)
SENDER_CMAP = LinearSegmentedColormap.from_list(
    "sender_spatial", ["#DDD9D8", "#E1ADB8", "#C17C8D", "#8A5261"]
)
RECEIVER_CMAP = LinearSegmentedColormap.from_list(
    "receiver_spatial", ["#D9DEE3", "#ACC7DF", "#78A1C8", "#3F709F"]
)
ECM_CMAP = LinearSegmentedColormap.from_list(
    "ecm_spatial", ["#DCE0DC", "#B8D0BC", "#86AA8D", "#52775A"]
)

PRIMARY = ["L1_PLAT_sender_same", "L2_sender_receiver_r1", "L3_receiver_ECMdeg_r0"]
PRIMARY_LABEL = {
    "L1_PLAT_sender_same": "PLAT ↔\nfibroblast state (r0)",
    "L2_sender_receiver_r1": "fibroblast ↔\nmacrophage (r1)",
    "L3_receiver_ECMdeg_r0": "macrophage ↔\nECM degradation (r0)",
}
PRIMARY_COLOR = {
    "L1_PLAT_sender_same": PLAT,
    "L2_sender_receiver_r1": RECEIVER,
    "L3_receiver_ECMdeg_r0": ECM,
}
INCREMENTAL = [
    "PLAT_incremental_receiver_r1",
    "PLAT_incremental_ECMdeg_r1",
    "PLAT_incremental_adhesion_r1",
]
INCREMENTAL_LABEL = {
    "PLAT_incremental_receiver_r1": "macrophage state\n(r1)",
    "PLAT_incremental_ECMdeg_r1": "ECM degradation\n(r1)",
    "PLAT_incremental_adhesion_r1": "cell–matrix adhesion\n(r1)",
}


def load_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def percentile(values: np.ndarray) -> np.ndarray:
    return pd.Series(np.asarray(values, float)).rank(method="average", pct=True).to_numpy() * 100.0


def choose_representative() -> tuple[str, pd.DataFrame]:
    """Choose a non-extreme ST_R section by a fully deterministic rule.

    Eligible sections must have all three primary effects in the cohort-level
    positive direction and PLAT detection at or above the ST_R median.  The
    selected section minimizes the sum of standardized distances to the three
    ST_R pooled effects.  This avoids choosing the strongest-looking map.
    """
    slide = load_tsv(STR_RESULT / "GSE270678_STR_PLAT_spatial_chain_slide_effects.tsv")
    coverage = load_tsv(STR_RESULT / "GSE270678_STR_PLAT_spatial_chain_gene_coverage.tsv")
    meta = load_tsv(STR_RESULT / "GSE270678_STR_PLAT_spatial_chain_patient_meta.tsv").set_index("endpoint")
    wide = (
        slide.loc[slide.endpoint.isin(PRIMARY)]
        .pivot(index=["sample_id", "donor_id"], columns="endpoint", values="rho")
        .reset_index()
        .merge(coverage[["sample_id", "PLAT_detected_fraction"]], on="sample_id", validate="one_to_one")
    )
    med_detect = float(coverage.PLAT_detected_fraction.median())
    distance_cols = []
    for endpoint in PRIMARY:
        sd = float(wide[endpoint].std(ddof=1))
        col = f"distance_{endpoint}"
        wide[col] = np.abs(wide[endpoint] - float(meta.loc[endpoint, "pooled_rho"])) / sd
        distance_cols.append(col)
    wide["distance"] = wide[distance_cols].sum(axis=1)
    wide["eligible_positive_chain"] = (wide[PRIMARY] > 0).all(axis=1)
    wide["eligible_PLAT_detection"] = wide.PLAT_detected_fraction >= med_detect
    eligible = wide.loc[wide.eligible_positive_chain & wide.eligible_PLAT_detection].copy()
    if eligible.empty:
        raise RuntimeError("No ST_R section passed the frozen representative-section rule")
    eligible = eligible.sort_values(["distance", "sample_id"])
    chosen = str(eligible.iloc[0].sample_id)
    wide["selected"] = wide.sample_id.eq(chosen)
    wide["ST_R_median_PLAT_detection"] = med_detect
    wide.to_csv(FIG_DATA / "Figure6A_representative_selection_source.tsv", sep="\t", index=False)
    if chosen != "D16T1_3":
        raise RuntimeError(f"Frozen selection expected D16T1_3, observed {chosen}")
    return chosen, wide


def load_image_gz(path: Path) -> Image.Image:
    with gzip.open(path, "rb") as fh:
        return Image.open(io.BytesIO(fh.read())).convert("RGB")


def prepare_representative(sample_id: str):
    manifest = load_tsv(MANIFEST)
    row = manifest.loc[manifest.sample_id.eq(sample_id)].iloc[0]
    spots = load_tsv(SPOT_STR)
    base = spots.loc[spots.sample_id.eq(sample_id) & spots.in_tissue.astype(bool)].sort_values("barcode").reset_index(drop=True)

    sender_genes, receiver_genes, deg_genes, _ = M.load_signatures()
    needed = sorted(set(sender_genes) | set(receiver_genes) | set(deg_genes) | {"PLAT"})
    h5_path = RAW_STR / row.filtered_h5_file
    # Total UMI is recomputed from the ST_R H5, matching the ST_R inference.
    import h5py
    from scipy.sparse import csc_matrix
    with h5py.File(h5_path, "r") as fh:
        mat = fh["matrix"]
        h5_barcodes = np.asarray(mat["barcodes"]).astype("U")
        sparse = csc_matrix(
            (mat["data"][...], mat["indices"][...], mat["indptr"][...]),
            shape=tuple(mat["shape"][...]),
        )
    lookup = {b: i for i, b in enumerate(h5_barcodes)}
    idx = np.array([lookup[b] for b in base.barcode], dtype=int)
    totals = np.asarray(sparse[:, idx].sum(axis=0)).ravel().astype(float)
    norm = M.read_selected_expression(h5_path, base.barcode.tolist(), needed, totals)
    sender, _, _ = M.module_score(norm, sender_genes, 10)
    receiver, _, _ = M.module_score(norm, receiver_genes, 10)
    ecm_deg, _, _ = M.module_score(norm, deg_genes, 20)

    base["PLAT_log1p_CP10K"] = norm.loc["PLAT"].to_numpy(float)
    base["sender_percentile"] = percentile(sender)
    base["receiver_percentile"] = percentile(receiver)
    base["ECMdeg_percentile"] = percentile(ecm_deg)

    with gzip.open(RAW_STR / row.scalefactors, "rt", encoding="utf-8") as fh:
        sf = json.load(fh)
    scale = float(sf["tissue_hires_scalef"])
    spot_diameter = float(sf["spot_diameter_fullres"]) * scale
    base["image_x"] = base.pxl_col.astype(float) * scale
    base["image_y"] = base.pxl_row.astype(float) * scale
    image = load_image_gz(RAW_STR / row.hires_image)
    base.to_csv(FIG_DATA / "Figure6A_D16T1_3_spatial_map_source.tsv.gz", sep="\t", index=False, compression="gzip")
    return base, image, spot_diameter


def crop_limits(d: pd.DataFrame, image: Image.Image):
    # The representative section contains one isolated in-tissue spot far
    # above the two connected capture regions.  Use exact Visium hex-neighbour
    # support for the display crop so that this single coordinate cannot pull
    # unanalysed H&E tissue into the map.  All spots remain in the statistical
    # analysis and source table.
    if {"array_row", "array_col"}.issubset(d.columns):
        keys = set(zip(d.array_row.astype(int), d.array_col.astype(int)))
        offsets = ((-1, -1), (-1, 1), (0, -2), (0, 2), (1, -1), (1, 1))
        supported = [
            any((int(r) + dr, int(c) + dc) in keys for dr, dc in offsets)
            for r, c in zip(d.array_row, d.array_col)
        ]
        crop_d = d.loc[np.asarray(supported, dtype=bool)]
    else:
        crop_d = d
    xlo, xhi = float(crop_d.image_x.min()), float(crop_d.image_x.max())
    ylo, yhi = float(crop_d.image_y.min()), float(crop_d.image_y.max())
    margin_x = max(20.0, 0.025 * (xhi - xlo))
    margin_y = max(20.0, 0.025 * (yhi - ylo))
    return (
        max(0.0, xlo - margin_x),
        min(float(image.width), xhi + margin_x),
        max(0.0, ylo - margin_y),
        min(float(image.height), yhi + margin_y),
    )


def add_spots(ax, d: pd.DataFrame, diameter: float, values: np.ndarray, cmap, norm):
    coll = EllipseCollection(
        widths=np.full(len(d), diameter * 0.93),
        heights=np.full(len(d), diameter * 0.93),
        angles=np.zeros(len(d)),
        units="xy",
        offsets=d[["image_x", "image_y"]].to_numpy(float),
        transOffset=ax.transData,
        cmap=cmap,
        norm=norm,
        edgecolors="#D9DDDC",
        linewidths=0.08,
        alpha=0.985,
        zorder=3,
    )
    coll.set_array(np.asarray(values, float))
    ax.add_collection(coll)
    return coll


def spatial_map(ax, d, image, diameter, column, cmap, norm, label, ticks, ticklabels, scale_bar=False):
    ax.imshow(image, alpha=0.64)
    coll = add_spots(ax, d, diameter, d[column].to_numpy(float), cmap, norm)
    x0, x1, y0, y1 = crop_limits(d, image)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_aspect("equal")
    ax.axis("off")
    # Keep map titles outside the tissue image so they cannot cover spots or
    # obscure the local expression pattern. The tight bounding box retains
    # this title margin in all four exports.
    ax.text(
        0.025, 1.025, label, transform=ax.transAxes, ha="left", va="bottom",
        fontsize=6.1, fontweight="bold", color=INK, clip_on=False,
    )
    cax = ax.inset_axes([0.18, -0.045, 0.64, 0.027])
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation="horizontal")
    cb.set_ticks(ticks)
    cb.set_ticklabels(ticklabels)
    cb.ax.tick_params(labelsize=4.5, pad=0.7, length=1.5, width=0.4)
    cb.outline.set_linewidth(0.35)
    if scale_bar:
        bar_um = 500
        bar_len = diameter * (bar_um / 55.0)
        bx = x0 + 0.06 * (x1 - x0)
        by = y1 - 0.07 * (y1 - y0)
        ax.plot([bx, bx + bar_len], [by, by], color="black", lw=1.2, solid_capstyle="butt", zorder=8)
        ax.text(bx + bar_len / 2, by - 0.018 * (y1 - y0), "500 µm", ha="center", va="bottom", fontsize=4.5)
    return coll


def load_meta():
    main = load_tsv(RESULT / "GSE270678_PLAT_spatial_chain_patient_meta.tsv")
    repeated = load_tsv(STR_RESULT / "GSE270678_STR_PLAT_spatial_chain_patient_meta.tsv")
    main["dataset"] = "Main ST"
    repeated["dataset"] = "ST_R sensitivity"
    independent = None
    frames = [main, repeated]
    if ENRICHED_CANDIDATE:
        independent = load_tsv(GSE251950_META).rename(
            columns={
                "n_patients": "n_donors",
                "positive_patients": "positive_donors",
                "negative_patients": "negative_donors",
            }
        )
        independent["dataset"] = "GSE251950 independent"
        frames.append(independent)
    combined = pd.concat(frames, ignore_index=True, sort=False)
    return main, repeated, independent, combined


def draw_effect_pair(ax, y: float, main_row, str_row, color: str):
    offsets = {"main": 0.115, "str": -0.115}
    ax.hlines(y + offsets["main"], main_row.ci_low, main_row.ci_high, color=color, lw=1.35, zorder=2)
    ax.scatter(main_row.pooled_rho, y + offsets["main"], s=30, color=color, edgecolor="white", linewidth=0.55, zorder=3)
    ax.hlines(y + offsets["str"], str_row.ci_low, str_row.ci_high, color=color, lw=1.1, alpha=0.64, zorder=2)
    ax.scatter(str_row.pooled_rho, y + offsets["str"], s=27, facecolor="white", edgecolor=color, linewidth=1.05, zorder=3)


def draw_effect_triplet(ax, y: float, main_row, independent_row, str_row, color: str):
    # Keep the two GSE270678 renders adjacent (Main ST above ST_R) and place
    # the independent GSE251950 cohort as the third row.
    offsets = {"main": 0.17, "str": 0.0, "independent": -0.17}
    ax.hlines(y + offsets["main"], main_row.ci_low, main_row.ci_high, color=color, lw=1.25, zorder=2)
    ax.scatter(main_row.pooled_rho, y + offsets["main"], s=28, color=color,
               edgecolor="white", linewidth=0.5, zorder=3)
    ax.hlines(y + offsets["independent"], independent_row.ci_low, independent_row.ci_high,
              color=color, lw=1.15, alpha=0.82, zorder=2)
    ax.scatter(independent_row.pooled_rho, y + offsets["independent"], s=30, marker="D",
               facecolor=color, edgecolor="white", linewidth=0.5, alpha=0.86, zorder=3)
    ax.hlines(y + offsets["str"], str_row.ci_low, str_row.ci_high,
              color=color, lw=1.0, alpha=0.55, zorder=2)
    ax.scatter(str_row.pooled_rho, y + offsets["str"], s=25, facecolor="white",
               edgecolor=color, linewidth=1.0, zorder=3)


def build_figure():
    selected, _ = choose_representative()
    rep, image, diameter = prepare_representative(selected)
    main, repeated, independent, combined = load_meta()

    primary_source = combined.loc[combined.endpoint.isin(PRIMARY)].copy()
    source_suffix = "_enriched" if ENRICHED_CANDIDATE else ""
    primary_source.to_csv(FIG_DATA / f"Figure6B_primary_chain{source_suffix}_source.tsv", sep="\t", index=False)
    ring_source = combined.loc[
        combined.endpoint.str.startswith(("L3_receiver_ECMdeg_r", "L3_receiver_adhesion_r"))
    ].copy()
    ring_source.to_csv(FIG_DATA / f"Figure6C_distance_decay{source_suffix}_source.tsv", sep="\t", index=False)
    incremental_source = combined.loc[combined.endpoint.isin(INCREMENTAL)].copy()
    incremental_source.to_csv(FIG_DATA / f"Figure6D_PLAT_incremental{source_suffix}_source.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(7.0, 5.5), facecolor="white")
    # Keep the map row visually dominant while giving b/c/d shorter plotting
    # bodies and more breathing room for the legends above the axes.
    outer = fig.add_gridspec(2, 1, height_ratios=[1.00, 0.90], hspace=0.22)

    maps = outer[0].subgridspec(1, 4, wspace=0.055)
    map_axes = [fig.add_subplot(maps[0, i]) for i in range(4)]
    panel_label(map_axes[0], "a", -0.12, 1.05)
    plat_vmax = max(0.5, float(np.quantile(rep.loc[rep.PLAT_log1p_CP10K.gt(0), "PLAT_log1p_CP10K"], 0.98)))
    spatial_map(map_axes[0], rep, image, diameter, "PLAT_log1p_CP10K", PLAT_CMAP,
                PowerNorm(gamma=0.62, vmin=0, vmax=plat_vmax),
                r"$\mathit{PLAT}$ expression", [0, plat_vmax], ["0", f"{plat_vmax:.2f}"], scale_bar=True)
    spatial_map(map_axes[1], rep, image, diameter, "sender_percentile", SENDER_CMAP,
                PowerNorm(gamma=0.72, vmin=0, vmax=100),
                "sender state\n" + r"$\mathit{CXCL14}^{+}\mathit{POSTN}^{+}$ fibroblast", [0, 50, 100], ["0", "50", "100"])
    spatial_map(map_axes[2], rep, image, diameter, "receiver_percentile", RECEIVER_CMAP,
                PowerNorm(gamma=0.72, vmin=0, vmax=100),
                "receiver state\n" + r"$\mathit{C1QC}^{+}\mathit{APOE}^{+}$ macrophage", [0, 50, 100], ["0", "50", "100"])
    spatial_map(map_axes[3], rep, image, diameter, "ECMdeg_percentile", ECM_CMAP,
                PowerNorm(gamma=0.72, vmin=0, vmax=100),
                "ECM degradation\nprogramme score", [0, 50, 100], ["0", "50", "100"])

    bottom = outer[1].subgridspec(1, 3, width_ratios=[1.16, 1.0, 1.04], wspace=0.48)

    # b | The three prespecified links in main ST and repeated-section sensitivity.
    ax_b = fig.add_subplot(bottom[0, 0])
    panel_label(ax_b, "b", -0.22, 1.06)
    y_positions = {endpoint: 2 - i for i, endpoint in enumerate(PRIMARY)}
    for endpoint in PRIMARY:
        mr = main.loc[main.endpoint.eq(endpoint)].iloc[0]
        sr = repeated.loc[repeated.endpoint.eq(endpoint)].iloc[0]
        if ENRICHED_CANDIDATE:
            ir = independent.loc[independent.endpoint.eq(endpoint)].iloc[0]
            draw_effect_triplet(ax_b, y_positions[endpoint], mr, ir, sr, PRIMARY_COLOR[endpoint])
        else:
            draw_effect_pair(ax_b, y_positions[endpoint], mr, sr, PRIMARY_COLOR[endpoint])
    ax_b.axvline(0, color=GREY, lw=0.55, ls=(0, (3, 2)))
    ax_b.set_yticks([y_positions[e] for e in PRIMARY])
    ax_b.set_yticklabels([PRIMARY_LABEL[e] for e in PRIMARY], fontsize=5.2)
    if ENRICHED_CANDIDATE:
        ax_b.set_xlim(-0.022, 0.245)
        ax_b.set_xticks([0, 0.05, 0.10, 0.15, 0.20])
    else:
        ax_b.set_xlim(-0.016, 0.132)
        ax_b.set_xticks([0, 0.04, 0.08, 0.12])
    ax_b.set_xlabel("Pooled adjusted partial\nSpearman ρ", fontsize=6.1, labelpad=2)
    ax_b.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
    ax_b.tick_params(axis="x", labelsize=5.0)
    ax_b.tick_params(axis="y", length=0, pad=2)
    ax_b.scatter([], [], s=25, color=GREY, edgecolor="white", label="GSE270678 Main ST")
    ax_b.scatter([], [], s=25, facecolor="white", edgecolor=GREY, label="GSE270678 ST_R sensitivity")
    if ENRICHED_CANDIDATE:
        ax_b.scatter([], [], s=27, marker="D", color=GREY, edgecolor="white", label="GSE251950 independent")
    ax_b.legend(loc="lower left", bbox_to_anchor=(0.0, 1.07), frameon=False,
                ncol=1,
                handletextpad=0.38, columnspacing=0.45, fontsize=5.15,
                borderaxespad=0, labelspacing=0.72, borderpad=0.15)
    clean_axis(ax_b)

    # c | Exact-shell distance profiles; same-spot to third ring.
    ax_c = fig.add_subplot(bottom[0, 1])
    panel_label(ax_c, "c", -0.20, 1.06)
    if ENRICHED_CANDIDATE:
        styles = [
            ("Main ST", "L3_receiver_ECMdeg_r", ECM, "o", "-"),
            ("Main ST", "L3_receiver_adhesion_r", ADH, "o", "-"),
            ("GSE251950 independent", "L3_receiver_ECMdeg_r", ECM, "D", "-"),
            ("GSE251950 independent", "L3_receiver_adhesion_r", ADH, "D", "-"),
        ]
    else:
        styles = [
            ("Main ST", "L3_receiver_ECMdeg_r", ECM, "o", "-"),
            ("Main ST", "L3_receiver_adhesion_r", ADH, "o", "-"),
            ("ST_R sensitivity", "L3_receiver_ECMdeg_r", ECM, "s", (0, (2, 2))),
            ("ST_R sensitivity", "L3_receiver_adhesion_r", ADH, "s", (0, (2, 2))),
        ]
    for dataset, prefix, color, marker, ls in styles:
        d = combined.loc[(combined.dataset.eq(dataset)) & combined.endpoint.str.startswith(prefix)].sort_values("ring")
        face = color if dataset in ("Main ST", "GSE251950 independent") else "white"
        ax_c.plot(d.ring, d.pooled_rho, color=color, lw=1.0, ls=ls, alpha=0.95)
        ax_c.errorbar(d.ring, d.pooled_rho,
                      yerr=[d.pooled_rho - d.ci_low, d.ci_high - d.pooled_rho],
                      fmt=marker, ms=3.8, mfc=face, mec=color, mew=0.85,
                      ecolor=color, elinewidth=0.75, capsize=1.6, zorder=3)
    ax_c.axhline(0, color=GREY, lw=0.5, ls=(0, (3, 2)))
    ax_c.set_xticks([0, 1, 2, 3])
    ax_c.set_yticks([0, 0.05, 0.10, 0.15, 0.20])
    ax_c.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
    ax_c.set_xlabel("exact hex-grid ring", fontsize=5.8, labelpad=2)
    ax_c.set_ylabel("pooled adjusted ρ", fontsize=5.8, labelpad=2)
    ax_c.tick_params(labelsize=5.0)
    endpoint_handles = [
        Line2D([0], [0], color=ECM, marker="o", lw=1.0, markersize=3.2, label="ECM degradation"),
        Line2D([0], [0], color=ADH, marker="o", lw=1.0, markersize=3.2, label="cell–matrix adhesion"),
    ]
    endpoint_legend = ax_c.legend(
        handles=endpoint_handles, loc="lower center", bbox_to_anchor=(0.5, 1.105),
        frameon=False, ncol=2, fontsize=5.05, handlelength=1.30,
        handletextpad=0.3, columnspacing=0.55, borderaxespad=0,
    )
    ax_c.add_artist(endpoint_legend)
    cohort_handles = [
        Line2D([0], [0], color=GREY, marker="o", lw=1.0, markersize=3.0,
               label="GSE270678 Main ST"),
        Line2D([0], [0], color=GREY, marker="D", lw=1.0, ls="-", markersize=3.0,
               label="GSE251950 independent"),
    ]
    ax_c.legend(
        handles=cohort_handles, loc="lower center", bbox_to_anchor=(0.5, 1.018),
        frameon=False, ncol=2, fontsize=5.0, handlelength=1.30,
        handletextpad=0.28, columnspacing=0.5, borderaxespad=0,
    )
    clean_axis(ax_c)
    # Increase the c-d gutter so panel d's long x-axis title cannot extend
    # underneath panel c at final double-column size.
    cpos = ax_c.get_position()
    ax_c.set_position([cpos.x0 - 0.032, cpos.y0, cpos.width * 0.98, cpos.height])

    # d | Incremental PLAT information after adjustment for frozen cell states.
    ax_d = fig.add_subplot(bottom[0, 2])
    panel_label(ax_d, "d", -0.20, 1.06)
    inc_y = {endpoint: 2 - i for i, endpoint in enumerate(INCREMENTAL)}
    for endpoint in INCREMENTAL:
        mr = main.loc[main.endpoint.eq(endpoint)].iloc[0]
        sr = repeated.loc[repeated.endpoint.eq(endpoint)].iloc[0]
        color = RECEIVER if "receiver" in endpoint else (ECM if "ECMdeg" in endpoint else ADH)
        if ENRICHED_CANDIDATE:
            ir = independent.loc[independent.endpoint.eq(endpoint)].iloc[0]
            draw_effect_triplet(ax_d, inc_y[endpoint], mr, ir, sr, color)
        else:
            draw_effect_pair(ax_d, inc_y[endpoint], mr, sr, color)
    ax_d.axvline(0, color=GREY, lw=0.55, ls=(0, (3, 2)))
    ax_d.set_yticks([inc_y[e] for e in INCREMENTAL])
    ax_d.set_yticklabels([INCREMENTAL_LABEL[e] for e in INCREMENTAL], fontsize=5.2)
    if ENRICHED_CANDIDATE:
        ax_d.set_xlim(-0.008, 0.074)
        ax_d.set_xticks([0, 0.02, 0.04, 0.06])
    else:
        ax_d.set_xlim(-0.008, 0.052)
        ax_d.set_xticks([0, 0.02, 0.04])
    ax_d.set_xlabel("PLAT → ring-1 endpoint\nAdjusted partial Spearman ρ", fontsize=5.9, labelpad=2)
    ax_d.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
    ax_d.tick_params(axis="x", labelsize=5.0)
    ax_d.tick_params(axis="y", length=0, pad=2)
    cohort_handles_d = [
        Line2D([0], [0], color=GREY, marker="o", lw=1.0, markersize=3.0,
               label="GSE270678 Main ST"),
        Line2D([0], [0], color=GREY, marker="o", lw=1.0, markerfacecolor="white",
               markersize=3.0, label="GSE270678 ST_R sensitivity"),
        Line2D([0], [0], color=GREY, marker="D", lw=1.0, markersize=3.0,
               label="GSE251950 independent"),
    ]
    ax_d.legend(
        handles=cohort_handles_d, loc="lower left", bbox_to_anchor=(0.0, 1.07),
        frameon=False, ncol=1, fontsize=5.05, handlelength=1.20,
        handletextpad=0.32, borderaxespad=0, labelspacing=0.55,
    )
    clean_axis(ax_d)

    dpos = ax_d.get_position()
    ax_d.set_position([dpos.x0 + 0.010, dpos.y0, dpos.width * 0.97, dpos.height])

    fig.subplots_adjust(left=0.07, right=0.985, top=0.985, bottom=0.09)
    # Standardized deliverable stem for the promoted enriched render.
    revised_dir = PROJECT / "figures" / "main" / "revised_20260827"
    revised_dir.mkdir(parents=True, exist_ok=True)
    stem = "revised_20260827/Figure_5_readability_v6"
    return save_figure(fig, stem)


if __name__ == "__main__":
    outputs = build_figure()
    print("\n".join(str(x) for x in outputs))
