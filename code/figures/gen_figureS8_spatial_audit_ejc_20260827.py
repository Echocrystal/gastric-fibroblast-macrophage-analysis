from __future__ import annotations

import gzip
import json
import os
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Polygon
from PIL import Image
from scipy.io import mmread


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paper_plot_style import PROJECT, clean_axis, panel_label  # noqa: E402


SOURCE_PROJECT = Path(
    r"I:/ECM-new/ECM_DEG_WGCNA_CoxLASSO_public_2026/ecm_immune_ecotype_public_2026"
)
SAMPLE_ROOT = SOURCE_PROJECT / "data/raw/spatial/GSE251950/samples"
MANIFEST_FILE = (
    SOURCE_PROJECT
    / "data/processed/spatial/GSE251950/GSE251950_spatial_sample_manifest.tsv"
)
SPOT_FILE = (
    PROJECT
    / "analysis/processed/spatial/GSE251950_candidate_axes/"
    "GSE251950_candidate_axis_spot_scores.tsv.gz"
)
PROGRAM_FILE = PROJECT / "resources/gene_sets/myeloid_receiver_programs_frozen.tsv"
SELECTION_DIR = PROJECT / "analysis/results/PLAT_axis/GSE251950/selection_corrected"
FORMAL_NULL_DIR = PROJECT / "analysis/results/PLAT_axis/GSE251950/formal_null"
FIG_DATA = PROJECT / "figures/data"
FIG_DATA.mkdir(parents=True, exist_ok=True)
REVISION_DIR = PROJECT / "figures" / "ejc_revision_20260827"
REVISION_DIR.mkdir(parents=True, exist_ok=True)

# The spatial figure was originally Figure 7.  The locked 2026-07-30
# architecture assigns it to Figure 6 and moves virtual perturbation to
# Figure 7.  Keep this renderer reusable so the historical Figure 7 bundle
# remains reproducible while new exports follow the current numbering.
FIGURE_ID = os.environ.get("PLAT_SPATIAL_FIGURE_ID", "FigureS8")
FIGURE_STEM = os.environ.get(
    "PLAT_SPATIAL_FIGURE_STEM", "Supplementary_Figure_S8_spatial_audit_ejc_20260827"
)
FIGURE_WIDTH_IN = float(os.environ.get("PLAT_SPATIAL_FIGURE_WIDTH_IN", str(183 / 25.4)))
FIGURE_HEIGHT_IN = float(os.environ.get("PLAT_SPATIAL_FIGURE_HEIGHT_IN", str(235 / 25.4)))

SLIDE_FILE = SELECTION_DIR / "GSE251950_PLAT_all11_selection_corrected_slide_statistics.tsv"
META_FILE = SELECTION_DIR / "GSE251950_PLAT_all11_selection_corrected_patient_meta.tsv"
MATCHED_FILE = FORMAL_NULL_DIR / "GSE251950_frozen_ECM_matched_ligand_patient_meta.tsv"

DEG = "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX"
ADH = "GOBP_CELL_MATRIX_ADHESION"
PROGRAM_FAMILY = [
    "ID_INFLAMMATORY_MONO",
    "ID_MONO_FCN1_VCAN",
    "GOBP_CELL_MATRIX_ADHESION",
    "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY",
    "GOBP_MYELOID_LEUKOCYTE_MIGRATION",
    "GOBP_PHAGOCYTOSIS",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
    "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
    "REACTOME_INTEGRIN_SIGNALING",
]
RECEPTORS = ["ANXA2", "ITGAM", "ITGB2", "LRP1", "ILK"]
COMPETITORS = ["CCL2", "MIF", "THBS1", "SPON2", "CTHRC1"]

PLAT = "#8F5F6B"
ECM = "#6F8B70"
MYELOID = "#456F9D"
INK = "#252525"
GREY = "#707980"
LIGHT_GREY = "#D7DDE1"
PALE_GREY = "#EEF1F3"
CREAM = "#F4EDCB"
SAGE = "#D6E1D3"
MIST_BLUE = "#B3CBE2"

PLAT_CMAP = LinearSegmentedColormap.from_list(
    "plat_spatial", ["#F6F2ED", "#E9D2C8", "#CDACB4", PLAT]
)
ECM_CMAP = LinearSegmentedColormap.from_list(
    "ecm_spatial", ["#F7F4E6", "#D6E1D3", "#A9C0A9", ECM]
)
ADH_CMAP = LinearSegmentedColormap.from_list(
    "adh_spatial", ["#F7F4E6", "#B3CBE2", "#8FABDA", MYELOID]
)


def save_revision_figure(fig, stem):
    outputs = []
    base = REVISION_DIR / stem
    for suffix, dpi in [(".pdf", 300), (".svg", 300), (".png", 400), (".tiff", 600)]:
        path = base.with_suffix(suffix)
        kwargs = {"dpi": dpi, "bbox_inches": "tight", "pad_inches": 0.04}
        if suffix == ".tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(path, **kwargs)
        outputs.append(path)
    plt.close(fig)
    return tuple(outputs)


def read_inputs():
    manifest = pd.read_csv(MANIFEST_FILE, sep="\t")
    manifest = manifest.loc[manifest["primary_analysis"].astype(bool)].copy()
    if len(manifest) != 9 or manifest["patient_id"].nunique() != 9:
        raise RuntimeError("Expected nine patient-independent primary sections")

    base_cols = [
        "GSM",
        "sample_title",
        "patient_id",
        "barcode",
        "array_row",
        "array_col",
        "pxl_row",
        "pxl_col",
        "total_UMI",
        "frozen_ECM_direction",
        "MARKER_Fibroblast",
        "MARKER_Myeloid",
        "lag_MARKER_Myeloid",
    ]
    spots = pd.read_csv(SPOT_FILE, sep="\t", usecols=base_cols)
    spots = spots.drop_duplicates().reset_index(drop=True)
    slide = pd.read_csv(SLIDE_FILE, sep="\t")
    meta = pd.read_csv(META_FILE, sep="\t")
    matched = pd.read_csv(MATCHED_FILE, sep="\t")
    return manifest, spots, slide, meta, matched


def choose_representative(slide: pd.DataFrame) -> pd.DataFrame:
    deg = slide.loc[slide.program.eq(DEG), ["patient_id", "n_spots", "rho", "PLAT_detected_fraction"]]
    deg = deg.rename(columns={"rho": "rho_ECM_degradation"})
    adh = slide.loc[slide.program.eq(ADH), ["patient_id", "rho"]]
    adh = adh.rename(columns={"rho": "rho_cell_matrix_adhesion"})
    out = deg.merge(adh, on="patient_id", validate="one_to_one")
    out["rank_ECM_degradation"] = out["rho_ECM_degradation"].rank(method="first")
    out["rank_cell_matrix_adhesion"] = out["rho_cell_matrix_adhesion"].rank(method="first")
    out["rank_PLAT_detection"] = out["PLAT_detected_fraction"].rank(method="first")
    out["representativeness_score"] = (
        (out["rank_ECM_degradation"] - 5).abs()
        + (out["rank_cell_matrix_adhesion"] - 5).abs()
        + (out["rank_PLAT_detection"] - 5).abs()
    )
    out = out.sort_values(
        ["representativeness_score", "n_spots", "patient_id"],
        ascending=[True, False, True],
    ).reset_index(drop=True)
    out["selected"] = False
    out.loc[0, "selected"] = True
    return out


def sample_paths(folder_id: str) -> dict[str, Path]:
    root = SAMPLE_ROOT / folder_id
    prefix = root / folder_id
    return {
        "root": root,
        "features": Path(str(prefix) + "_features.tsv.gz"),
        "barcodes": Path(str(prefix) + "_barcodes.tsv.gz"),
        "matrix": Path(str(prefix) + "_matrix.mtx.gz"),
        "lowres": Path(str(prefix) + "_tissue_lowres_image.png"),
        "scales": Path(str(prefix) + "_scalefactors_json.json"),
    }


def read_programs(matched: pd.DataFrame) -> dict[str, list[str]]:
    matched_ligands = matched["ligand"].drop_duplicates().tolist()
    excluded = set(["PLAT", *matched_ligands, *COMPETITORS, *RECEPTORS])
    p = pd.read_csv(PROGRAM_FILE, sep="\t")
    p = p.loc[p.program.isin(PROGRAM_FAMILY) & ~p.gene_symbol.isin(excluded)].copy()
    programs = {k: v.gene_symbol.drop_duplicates().tolist() for k, v in p.groupby("program")}
    if set(programs) != set(PROGRAM_FAMILY):
        raise RuntimeError("The complete 11-program family is unavailable")
    return programs


def first_ring_lag(values: np.ndarray, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    offsets = [(-1, -1), (-1, 1), (0, -2), (0, 2), (1, -1), (1, 1)]
    lookup = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(rows, cols))}
    out = np.full(len(values), np.nan, dtype=float)
    for i, (r, c) in enumerate(zip(rows, cols)):
        idx = [lookup.get((int(r + dr), int(c + dc))) for dr, dc in offsets]
        idx = [j for j in idx if j is not None]
        if idx:
            out[i] = np.nanmean(values[idx])
    return out


def reconstruct_representative_spots(
    patient_id: str,
    manifest: pd.DataFrame,
    spots: pd.DataFrame,
    matched: pd.DataFrame,
) -> tuple[pd.DataFrame, Path, float, float]:
    info = manifest.loc[manifest.patient_id.eq(patient_id)].iloc[0]
    paths = sample_paths(info.folder_id)
    base = spots.loc[spots.GSM.eq(info.GSM)].copy().sort_values("barcode").reset_index(drop=True)

    programs = read_programs(matched)
    needed = sorted(set(["PLAT", *programs[DEG], *programs[ADH]]))
    features = pd.read_csv(paths["features"], sep="\t", header=None, compression="gzip")
    barcodes = pd.read_csv(paths["barcodes"], sep="\t", header=None, compression="gzip")[0].astype(str)
    symbol_to_rows: dict[str, list[int]] = {}
    for i, symbol in enumerate(features[1].astype(str)):
        if symbol in needed:
            symbol_to_rows.setdefault(symbol, []).append(i)

    with gzip.open(paths["matrix"], "rb") as fh:
        # SciPy <1.15 does not expose the ``spmatrix`` keyword.  Matrix Market
        # input is sparse here in every supported runtime, so the compatible
        # call followed by ``tocsr`` preserves the same matrix representation.
        raw = mmread(fh).tocsr()
    barcode_index = pd.Series(np.arange(len(barcodes)), index=barcodes).to_dict()
    col_idx = np.array([barcode_index[x] for x in base.barcode], dtype=int)
    denom = np.maximum(base.total_UMI.to_numpy(float), 1.0)

    expression: dict[str, np.ndarray] = {}
    for gene in needed:
        rows = symbol_to_rows.get(gene, [])
        if not rows:
            continue
        counts = np.asarray(raw[rows, :][:, col_idx].sum(axis=0)).ravel()
        expression[gene] = np.log1p(counts * 1e4 / denom)

    if "PLAT" not in expression:
        raise RuntimeError(f"PLAT unavailable in representative section {patient_id}")

    def score(program_id: str) -> np.ndarray:
        genes = [g for g in programs[program_id] if g in expression]
        z_rows = []
        for gene in genes:
            x = expression[gene]
            sd = np.std(x, ddof=1)
            if np.isfinite(sd) and sd > 0:
                z_rows.append((x - np.mean(x)) / sd)
        if len(z_rows) < 5:
            raise RuntimeError(f"Insufficient informative genes for {program_id}")
        return np.mean(np.vstack(z_rows), axis=0)

    deg_score = score(DEG)
    adh_score = score(ADH)
    base["PLAT_log1p_CP10K"] = expression["PLAT"]
    base["neighbor_ECM_degradation_score"] = first_ring_lag(
        deg_score, base.array_row.to_numpy(), base.array_col.to_numpy()
    )
    base["neighbor_cell_matrix_adhesion_score"] = first_ring_lag(
        adh_score, base.array_row.to_numpy(), base.array_col.to_numpy()
    )
    base["neighbor_ECM_degradation_percentile"] = (
        base["neighbor_ECM_degradation_score"].rank(method="average", pct=True) * 100
    )
    base["neighbor_cell_matrix_adhesion_percentile"] = (
        base["neighbor_cell_matrix_adhesion_score"].rank(method="average", pct=True) * 100
    )

    with open(paths["scales"], "r", encoding="utf-8") as fh:
        scales = json.load(fh)
    lowres_scale = float(scales["tissue_lowres_scalef"])
    spot_diameter = float(scales["spot_diameter_fullres"]) * lowres_scale
    base["lowres_x"] = base.pxl_col * lowres_scale
    base["lowres_y"] = base.pxl_row * lowres_scale
    return base, paths["lowres"], lowres_scale, spot_diameter


def crop_limits(d: pd.DataFrame, image: Image.Image, pad_fraction: float = 0.045):
    x0, x1 = d.lowres_x.min(), d.lowres_x.max()
    y0, y1 = d.lowres_y.min(), d.lowres_y.max()
    pad = pad_fraction * max(x1 - x0, y1 - y0)
    return (
        max(0, x0 - pad),
        min(image.size[0], x1 + pad),
        max(0, y0 - pad),
        min(image.size[1], y1 + pad),
    )


def draw_hne(ax, image: Image.Image, d: pd.DataFrame, spot_size=1.2, alpha=0.55):
    ax.imshow(image)
    ax.scatter(
        d.lowres_x,
        d.lowres_y,
        s=spot_size,
        facecolors="none",
        edgecolors="#64717A",
        linewidths=0.20,
        alpha=alpha,
    )
    x0, x1, y0, y1 = crop_limits(d, image)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_aspect("equal")
    ax.axis("off")


def draw_spatial_map(
    ax,
    image: Image.Image,
    d: pd.DataFrame,
    value: str,
    cmap,
    norm,
    title: str,
    cbar_ticks,
    cbar_ticklabels,
    spot_size=8,
    zero_is_background=False,
):
    ax.imshow(image, alpha=0.34)
    values = d[value].to_numpy(float)
    ok = np.isfinite(values)
    if zero_is_background:
        background = ok & (values <= 0)
        ax.scatter(
            d.loc[background, "lowres_x"],
            d.loc[background, "lowres_y"],
            s=max(1.5, spot_size * 0.62),
            color="#E8EAEB",
            linewidths=0,
            alpha=0.52,
            zorder=1,
        )
        ok = ok & (values > 0)
    ax.scatter(
        d.loc[ok, "lowres_x"],
        d.loc[ok, "lowres_y"],
        c=values[ok],
        cmap=cmap,
        norm=norm,
        s=spot_size,
        linewidths=0,
        alpha=0.95,
        zorder=2,
    )
    x0, x1, y0, y1 = crop_limits(d, image)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=6.1, pad=1.5, color=INK, fontweight="bold")
    cax = ax.inset_axes([0.14, -0.05, 0.72, 0.027])
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation="horizontal")
    cb.set_ticks(cbar_ticks)
    cb.set_ticklabels(cbar_ticklabels)
    cb.ax.tick_params(labelsize=4.5, length=1.6, pad=0.8, width=0.4)
    cb.outline.set_linewidth(0.35)


def short_program(program: str) -> str:
    lookup = {
        DEG: "ECM degradation",
        ADH: "Cell–matrix adhesion",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNF–NF-κB",
        "REACTOME_INTEGRIN_SIGNALING": "Integrin signaling",
        "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY": "Integrin-mediated signaling",
        "HALLMARK_INFLAMMATORY_RESPONSE": "Inflammatory response",
        "ID_INFLAMMATORY_MONO": "Inflammatory monocyte",
        "GOBP_MYELOID_LEUKOCYTE_MIGRATION": "Myeloid migration",
        "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS": "Integrin surface interactions",
        "GOBP_PHAGOCYTOSIS": "Phagocytosis",
        "ID_MONO_FCN1_VCAN": "FCN1/VCAN monocyte",
    }
    return lookup.get(program, program)


def p_text(value: float) -> str:
    if value < 0.01:
        return f"{value:.2e}"
    return f"{value:.2f}"


def draw_patient_effect(ax, slide, meta, program, title, color, letter):
    d = slide.loc[slide.program.eq(program)].sort_values("rho").reset_index(drop=True)
    m = meta.loc[meta.program.eq(program)].iloc[0]
    y = np.arange(len(d), 0, -1)
    ax.hlines(y, 0, d.rho, color=LIGHT_GREY, lw=1.1, zorder=1)
    ax.scatter(d.rho, y, s=19, color=color, edgecolor="white", linewidth=0.45, zorder=3)
    ax.axvline(0, color=GREY, lw=0.55, ls=(0, (3, 2)), zorder=0)
    pooled_y = -0.35
    ax.hlines(pooled_y, m.ci_low, m.ci_high, color=color, lw=1.5, zorder=3)
    half_h = 0.22
    diamond = np.array(
        [[m.pooled_rho, pooled_y + half_h], [m.ci_high, pooled_y],
         [m.pooled_rho, pooled_y - half_h], [m.ci_low, pooled_y]]
    )
    ax.add_patch(Polygon(diamond, closed=True, facecolor=color, edgecolor="white", lw=0.45))
    ax.set_yticks([*y, pooled_y])
    ax.set_yticklabels([*d.patient_id, "Pooled"], fontsize=6.0)
    ax.set_ylim(-0.8, len(d) + 0.8)
    ax.set_xlim(-0.018, 0.125)
    ax.set_xticks([0, 0.04, 0.08, 0.12])
    ax.set_xlabel("Adjusted partial Spearman ρ", fontsize=6.5, labelpad=2)
    ax.tick_params(axis="x", labelsize=5.7)
    ax.tick_params(axis="y", length=0, pad=2)
    clean_axis(ax)
    panel_label(ax, letter, -0.16, 1.06)
    ax.set_title(title, loc="left", pad=4, fontsize=7.4, fontweight="bold")
    ax.text(
        0.98,
        0.96,
        f"{int(m.positive_patients)}/9 positive\n"
        f"pooled ρ = {m.pooled_rho:.2f}\n"
        f"maxT-FWER P = {p_text(m.patient_signflip_maxT_FWER_p)}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.7,
        color=INK,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.82, pad=1.4),
    )


def source_data_exports(manifest, spots, selection, rep, slide, meta, matched):
    counts = spots.groupby(["GSM", "patient_id"], as_index=False).size().rename(columns={"size": "analysed_spots"})
    a = manifest.merge(counts, on=["GSM", "patient_id"], how="left", validate="one_to_one")
    a.to_csv(FIG_DATA / f"{FIGURE_ID}A_spatial_manifest_source.tsv", sep="\t", index=False)
    selection.to_csv(FIG_DATA / f"{FIGURE_ID}B_representative_selection_source.tsv", sep="\t", index=False)
    map_cols = [
        "GSM", "sample_title", "patient_id", "barcode", "array_row", "array_col",
        "pxl_row", "pxl_col", "lowres_x", "lowres_y", "total_UMI",
        "PLAT_log1p_CP10K", "neighbor_ECM_degradation_score",
        "neighbor_ECM_degradation_percentile", "neighbor_cell_matrix_adhesion_score",
        "neighbor_cell_matrix_adhesion_percentile",
    ]
    rep[map_cols].to_csv(
        FIG_DATA / f"{FIGURE_ID}B_{rep.patient_id.iloc[0]}_spatial_map_source.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    slide.loc[slide.program.isin([DEG, ADH])].to_csv(
        FIG_DATA / f"{FIGURE_ID}C_D_patient_effects_source.tsv", sep="\t", index=False
    )
    meta.to_csv(FIG_DATA / f"{FIGURE_ID}E_all11_programs_source.tsv", sep="\t", index=False)
    matched.loc[matched.program.isin([DEG, ADH])].to_csv(
        FIG_DATA / f"{FIGURE_ID}F_matched_ligands_source.tsv", sep="\t", index=False
    )


def build_figure():
    manifest, spots, slide, meta, matched = read_inputs()
    selection = choose_representative(slide)
    representative = selection.loc[selection.selected, "patient_id"].iloc[0]
    if representative != "GC2":
        raise RuntimeError(f"Frozen representative rule expected GC2, observed {representative}")
    rep, rep_image_path, _, _ = reconstruct_representative_spots(
        representative, manifest, spots, matched
    )
    source_data_exports(manifest, spots, selection, rep, slide, meta, matched)

    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN), facecolor="white")
    outer = fig.add_gridspec(3, 1, height_ratios=[1.52, 0.88, 1.08], hspace=0.37)
    top = outer[0].subgridspec(1, 2, width_ratios=[0.44, 0.56], wspace=0.12)
    middle = outer[1].subgridspec(1, 2, wspace=0.34)
    bottom = outer[2].subgridspec(1, 2, width_ratios=[1.08, 0.92], wspace=0.35)

    # a | Nine patient-independent primary sections.
    ax_a = fig.add_subplot(top[0, 0])
    ax_a.axis("off")
    panel_label(ax_a, "a", -0.055, 1.04)
    ax_a.set_title("Nine primary gastric tumours", loc="left", pad=4, fontsize=7.7, fontweight="bold")
    for i, row in enumerate(manifest.sort_values("patient_id").itertuples(index=False)):
        rr, cc = divmod(i, 3)
        sub = ax_a.inset_axes([0.005 + cc * 0.333, 0.675 - rr * 0.325, 0.305, 0.275])
        d = spots.loc[spots.GSM.eq(row.GSM)].copy()
        paths = sample_paths(row.folder_id)
        with open(paths["scales"], "r", encoding="utf-8") as fh:
            scale = float(json.load(fh)["tissue_lowres_scalef"])
        d["lowres_x"] = d.pxl_col * scale
        d["lowres_y"] = d.pxl_row * scale
        image = Image.open(paths["lowres"])
        draw_hne(sub, image, d, spot_size=0.55, alpha=0.28)
        sub.text(
            0.02, 0.98, f"{row.patient_id}  n={len(d):,}", transform=sub.transAxes,
            ha="left", va="top", fontsize=4.7, fontweight="bold", color=INK,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, pad=0.6),
        )

    # b | Objective representative section.
    ax_b = fig.add_subplot(top[0, 1])
    ax_b.axis("off")
    panel_label(ax_b, "b", -0.045, 1.04)
    ax_b.set_title(
        "Representative section selected by prespecified median-rank rule",
        loc="left", pad=4, fontsize=7.7, fontweight="bold",
    )
    image = Image.open(rep_image_path)
    b_positions = [(0.00, 0.56), (0.51, 0.56), (0.00, 0.11), (0.51, 0.11)]
    b_axes = [ax_b.inset_axes([x, y, 0.48, 0.36]) for x, y in b_positions]
    draw_hne(b_axes[0], image, rep, spot_size=1.0, alpha=0.34)
    b_axes[0].set_title("GC2 H&E · 4,128 profiled spots", fontsize=6.1, pad=1.5, color=INK, fontweight="bold")
    plat_vmax = float(np.nanquantile(rep.loc[rep.PLAT_log1p_CP10K.gt(0), "PLAT_log1p_CP10K"], 0.98))
    draw_spatial_map(
        b_axes[1], image, rep, "PLAT_log1p_CP10K", PLAT_CMAP,
        Normalize(0, plat_vmax), "PLAT expression", [0, plat_vmax], ["0", f"{plat_vmax:.2f}"], 7,
        zero_is_background=True,
    )
    draw_spatial_map(
        b_axes[2], image, rep, "neighbor_ECM_degradation_percentile", ECM_CMAP,
        Normalize(0, 100), "Neighbouring ECM degradation", [0, 50, 100], ["0", "50", "100"], 7,
    )
    draw_spatial_map(
        b_axes[3], image, rep, "neighbor_cell_matrix_adhesion_percentile", ADH_CMAP,
        Normalize(0, 100), "Neighbouring cell–matrix adhesion", [0, 50, 100], ["0", "50", "100"], 7,
    )

    # c-d | Patient-level effects.
    ax_c = fig.add_subplot(middle[0, 0])
    draw_patient_effect(ax_c, slide, meta, DEG, "PLAT and neighbouring ECM degradation", ECM, "c")
    ax_d = fig.add_subplot(middle[0, 1])
    draw_patient_effect(ax_d, slide, meta, ADH, "PLAT and neighbouring cell–matrix adhesion", MYELOID, "d")

    # e | Complete searched programme family with maxT correction.
    ax_e = fig.add_subplot(bottom[0, 0])
    panel_label(ax_e, "e", -0.16, 1.04)
    ax_e.set_title("All 11 searched spatial programmes", loc="left", pad=4, fontsize=7.7, fontweight="bold")
    d = meta.sort_values("pooled_rho", ascending=True).reset_index(drop=True)
    y = np.arange(len(d))
    sig = d.maxT_FWER_gate_0_05.astype(bool).to_numpy()
    colors = np.where(sig, ECM, "#A9B6C2")
    ax_e.hlines(y, d.ci_low, d.ci_high, color=colors, lw=1.15, zorder=2)
    ax_e.scatter(
        d.pooled_rho, y, s=np.where(sig, 24, 17), c=colors,
        edgecolor="white", linewidth=0.45, zorder=3,
    )
    ax_e.axvline(0, color=GREY, lw=0.55, ls=(0, (3, 2)), zorder=0)
    ax_e.set_yticks(y)
    ax_e.set_yticklabels([short_program(x) for x in d.program], fontsize=5.5)
    ax_e.set_xlim(-0.04, 0.135)
    ax_e.set_xticks([-0.02, 0.02, 0.06, 0.10])
    ax_e.set_xlabel("Patient-equal pooled ρ (pointwise 95% CI)", fontsize=6.3, labelpad=2)
    ax_e.tick_params(axis="x", labelsize=5.6)
    ax_e.tick_params(axis="y", length=0, pad=2)
    for yi, row in enumerate(d.itertuples(index=False)):
        ax_e.text(
            0.132, yi, p_text(row.patient_signflip_maxT_FWER_p),
            ha="right", va="center", fontsize=5.1,
            color=ECM if row.maxT_FWER_gate_0_05 else INK,
            fontweight="bold" if row.maxT_FWER_gate_0_05 else "normal",
        )
    ax_e.text(0.132, len(d) - 0.1, "maxT-FWER P", ha="right", va="bottom", fontsize=5.0, color=INK)
    clean_axis(ax_e)

    # f | Common-base matched-ligand calibration.
    ax_f = fig.add_subplot(bottom[0, 1])
    panel_label(ax_f, "f", -0.16, 1.04)
    ax_f.set_title("PLAT within a shared ECM-ligand programme", loc="left", pad=4, fontsize=7.7, fontweight="bold")
    x = matched.loc[matched.program.eq(DEG), ["ligand", "pooled_rho", "rank_descending"]].rename(
        columns={"pooled_rho": "rho_deg", "rank_descending": "rank_deg"}
    )
    ydata = matched.loc[matched.program.eq(ADH), ["ligand", "pooled_rho", "rank_descending"]].rename(
        columns={"pooled_rho": "rho_adh", "rank_descending": "rank_adh"}
    )
    xy = x.merge(ydata, on="ligand", validate="one_to_one")
    controls = xy.loc[~xy.ligand.eq("PLAT")]
    plat = xy.loc[xy.ligand.eq("PLAT")].iloc[0]
    ax_f.scatter(
        controls.rho_deg, controls.rho_adh, s=22, color="#B8C5CE",
        edgecolor="white", linewidth=0.45, alpha=0.94, zorder=2,
    )
    ax_f.axvline(plat.rho_deg, color=PLAT, lw=0.7, ls=(0, (3, 2)), alpha=0.8)
    ax_f.axhline(plat.rho_adh, color=PLAT, lw=0.7, ls=(0, (3, 2)), alpha=0.8)
    ax_f.scatter(
        [plat.rho_deg], [plat.rho_adh], s=54, marker="D", color=PLAT,
        edgecolor="white", linewidth=0.65, zorder=4,
    )
    ax_f.annotate(
        "PLAT\nrank 12/21 on both axes",
        (plat.rho_deg, plat.rho_adh), xytext=(8, -9), textcoords="offset points",
        fontsize=5.7, color=PLAT, fontweight="bold", ha="left", va="top",
        arrowprops=dict(arrowstyle="-", color=PLAT, lw=0.55),
    )
    for ligand in ["TIMP2", "COL4A1", "COL6A1"]:
        row = xy.loc[xy.ligand.eq(ligand)].iloc[0]
        ax_f.text(row.rho_deg + 0.002, row.rho_adh + 0.0015, ligand, fontsize=4.8, color=INK)
    ax_f.set_xlim(-0.07, 0.145)
    ax_f.set_ylim(-0.02, 0.115)
    ax_f.set_xlabel("ECM degradation pooled ρ", fontsize=6.3, labelpad=2)
    ax_f.set_ylabel("Cell–matrix adhesion pooled ρ", fontsize=6.3, labelpad=2)
    ax_f.tick_params(labelsize=5.6)
    clean_axis(ax_f)

    fig.subplots_adjust(left=0.075, right=0.985, top=0.985, bottom=0.055)
    outputs = save_revision_figure(fig, FIGURE_STEM)
    return outputs


if __name__ == "__main__":
    outputs = build_figure()
    print("\n".join(str(x) for x in outputs))
