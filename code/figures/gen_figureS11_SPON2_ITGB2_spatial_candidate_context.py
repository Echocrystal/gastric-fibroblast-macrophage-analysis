# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Candidate-axis forest plot -> cross-type inherit -> param inherit
# (b) Exact-ring distance profile -> cross-type inherit -> param inherit
# (c) Receiver-function forest plot -> cross-type inherit -> param inherit
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

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026")
DATA = ROOT / "analysis" / "results" / "spatial" / "SPON2_ITGB2_candidate_context_v1"
OUT = ROOT / "figures" / "supplementary" / "FigureS11_SPON2_ITGB2_spatial_candidate_context"
OUT.parent.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(DATA / "SPON2_ITGB2_meta_summary.tsv", sep="\t")

def fmt_dataset(x):
    return {"GSE270678_main_ST": "GSE270678 main ST", "GSE270678_ST_R_sensitivity": "GSE270678 ST_R", "GSE251950_primary": "GSE251950"}[x]

def forest(ax, df, title, color_map, ylabels, xlabel="Pooled Spearman ρ"):
    df = df.copy().sort_values(["dataset", "axis"], ascending=[True, True])
    rows = []
    for dataset in ["GSE270678_main_ST", "GSE270678_ST_R_sensitivity", "GSE251950_primary"]:
        sub = df[df.dataset.eq(dataset)]
        for _, r in sub.iterrows():
            rows.append(r)
    y = np.arange(len(rows))[::-1]
    for yi, r in zip(y, rows):
        col = color_map.get(r.axis, CATEGORICAL[0])
        ax.plot([r.ci_low, r.ci_high], [yi, yi], color=col, lw=1.0, solid_capstyle="round")
        ax.scatter(r.pooled_rho, yi, s=22, color=col, edgecolor="white", linewidth=0.4, zorder=3)
    labels = [f"{fmt_dataset(r.dataset)} | {r.axis.replace('_', '–')}" for r in rows]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.axvline(0, color=GREY, lw=0.6, ls="--")
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left", fontweight="bold", pad=5)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.07, 0.09)
    ax.grid(False)

fig = plt.figure(figsize=(183/25.4, 142/25.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], height_ratios=[1.1, 1], wspace=0.55, hspace=0.75)
ax_a = fig.add_subplot(gs[:, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, 1])

candidate = meta[meta.category.eq("candidate_calibration")].copy()
forest(ax_a, candidate, "Candidate-axis calibration", {
    "SPON2_ITGB2": ACCENT_RED, "APP_CD74": CATEGORICAL[0], "NRG1_HLA_DPB1": CATEGORICAL[2], "CXCL14_CXCR4": CATEGORICAL[3]
}, None)
ax_a.text(-0.07, 1.02, "a", transform=ax_a.transAxes, fontsize=9, fontweight="bold", va="bottom")

distance = meta[(meta.category.eq("distance_profile")) & (meta.axis.eq("SPON2_ITGB2"))].copy()
for dataset, col in [("GSE270678_main_ST", CATEGORICAL[0]), ("GSE270678_ST_R_sensitivity", CATEGORICAL[2]), ("GSE251950_primary", ACCENT_RED)]:
    d = distance[distance.dataset.eq(dataset)].sort_values("ring")
    ax_b.plot(d.ring, d.pooled_rho, marker="o", ms=3.8, lw=1.2, color=col, label=fmt_dataset(dataset))
    ax_b.fill_between(d.ring.to_numpy(), d.ci_low.to_numpy(), d.ci_high.to_numpy(), color=col, alpha=0.12, linewidth=0)
ax_b.axhline(0, color=GREY, lw=0.6, ls="--")
ax_b.set_xticks([0,1,2,3]); ax_b.set_xlabel("Exact Visium ring"); ax_b.set_ylabel("Pooled Spearman ρ")
ax_b.set_title("SPON2–ITGB2 receiver context", loc="left", fontweight="bold", pad=5)
ax_b.legend(loc="upper right", fontsize=6.5, handlelength=1.6)
ax_b.text(-0.16, 1.02, "b", transform=ax_b.transAxes, fontsize=9, fontweight="bold", va="bottom")

func = meta[(meta.category.eq("receiver_function")) & (meta.ring.eq(1)) & (meta.dataset.eq("GSE270678_main_ST"))].copy()
func["label"] = func.endpoint.str.replace("ITGB2_", "", regex=False).str.replace("_r1", "", regex=False).str.replace("_", " ", regex=False)
func = func.sort_values("pooled_rho")
y = np.arange(len(func))
ax_c.hlines(y, func.ci_low, func.ci_high, color=CATEGORICAL[0], lw=1.0)
ax_c.scatter(func.pooled_rho, y, s=24, color=CATEGORICAL[0], edgecolor="white", linewidth=0.4, zorder=3)
ax_c.axvline(0, color=GREY, lw=0.6, ls="--")
ax_c.set_yticks(y); ax_c.set_yticklabels(func.label.str.title(), fontsize=7)
ax_c.set_xlabel("Pooled Spearman ρ"); ax_c.set_title("ITGB2 receiver-function links", loc="left", fontweight="bold", pad=5)
ax_c.text(-0.16, 1.02, "c", transform=ax_c.transAxes, fontsize=9, fontweight="bold", va="bottom")

fig.text(0.01, 0.005, "Points and lines show pooled estimates and 95% bootstrap confidence intervals; donors or patients were weighted equally.", fontsize=6.5)
save_cns_figure(fig, str(OUT))
fig.savefig(str(OUT.with_suffix(".svg")), bbox_inches="tight")
plt.close(fig)
print(OUT)
