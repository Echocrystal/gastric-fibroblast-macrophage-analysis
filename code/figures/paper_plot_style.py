from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

PROJECT = Path(r"I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026")
FIG_DIR = PROJECT / "figures" / "main"
FIG_DIR.mkdir(parents=True, exist_ok=True)

COLORS = {
    "navy": "#2F5597",
    "blue": "#4C78A8",
    "cyan": "#59A9C2",
    "orange": "#E38D3D",
    "red": "#C94C4C",
    "green": "#5B8E5A",
    "purple": "#8A6BBE",
    "grey": "#7A7A7A",
    "light_grey": "#D9D9D9",
    "dark": "#000000",
    # Low-saturation manuscript palette derived from the user-selected scheme.
    "dusty_rose": "#CDACB4",
    "peach": "#E9D2C8",
    "cream": "#F4EDCB",
    "sage": "#D6E1D3",
    "powder_blue": "#B3CBE2",
    "steel_blue": "#8FABDA",
    "plat": "#C45A66",
    "sender_dark": "#9D6875",
    "receiver_dark": "#527CAA",
    "ecm_dark": "#6F9076",
}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8.2,
    "axes.labelsize": 8.5,
    "axes.titlesize": 9,
    "xtick.labelsize": 7.4,
    "ytick.labelsize": 7.4,
    "legend.fontsize": 7.2,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.04,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "text.color": "#000000",
    "axes.labelcolor": "#000000",
    "axes.edgecolor": "#000000",
    "xtick.color": "#000000",
    "ytick.color": "#000000",
})


def panel_label(ax, label, x=-0.14, y=1.08):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="top", ha="left")


def clean_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=3, width=0.6, colors="#000000")
    ax.xaxis.label.set_color("#000000")
    ax.yaxis.label.set_color("#000000")
    ax.spines["left"].set_color("#000000")
    ax.spines["bottom"].set_color("#000000")


def save_figure(fig, stem):
    pdf = FIG_DIR / f"{stem}.pdf"
    png = FIG_DIR / f"{stem}.png"
    svg = FIG_DIR / f"{stem}.svg"
    tiff = FIG_DIR / f"{stem}.tiff"
    fig.savefig(pdf)
    fig.savefig(svg)
    fig.savefig(png, dpi=400)
    fig.savefig(tiff, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    return pdf, png, svg, tiff
