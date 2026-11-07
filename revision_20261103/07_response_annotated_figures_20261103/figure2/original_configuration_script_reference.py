# ASSET CONFIRMATION TABLE
# A | response + Wilson CI | GroupedBarChart/plot_GroupedBarChartv1.py inspected;
#   | semantic mismatch (model AUROC grouped bars): inherit uncluttered axes/CI.
# B,C | binary status composition | BarComposition/plot_rewriting.py inspected;
#   | semantic mismatch (unstacked before/after proportions): custom stacked bars.
# D | unpaired individual CPS + median/IQR | PairedBoxScatter1.R inspected;
#   | semantic mismatch (paired before/after): retain individual points, no lines.
# All companion PNGs inspected before construction. Author-specified 183 mm
# four-panel overview supersedes asset layout and the skill confirmation step.

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
import math
import platform
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parents[1] / "02_clinical_stratification/Kim2018_patient_input_snapshot.tsv"
GROUPS = ["Dual-low", "Fibroblast-dominant", "Macrophage-dominant", "Dual-high"]
DISPLAY = ["Dual-low", "Fibroblast-\ndominant", "Macrophage-\ndominant", "Dual-high"]
# Existing Figure 6 group colours; preserve roles in panels A and D.
GROUP_COLORS = dict(zip(GROUPS, ["#AAB7C4", "#CDACB4", "#8FABDA", "#9FB9A1"]))
POSITIVE, NEGATIVE, MISSING = "#B45D35", "#E3E8EC", "#999999"
WIDTH_MM, HEIGHT_MM = 183, 165


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def wilson(k, n):
    z = 1.959963984540054
    p = k/n
    denom = 1 + z*z/n
    center = (p + z*z/(2*n))/denom
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/denom
    return 0.0 if k == 0 else max(0, center-half), 1.0 if k == n else min(1, center+half)


def read_input():
    dest = ROOT / "Kim2018_patient_input_snapshot.tsv"
    if not dest.exists():
        shutil.copy2(SOURCE, dest)
    d = pd.read_csv(dest, sep="\t")
    assert len(d) == 45 and d.patient_id.nunique() == 45
    assert int(d.response_binary.sum()) == 12
    assert set(d.ecotype4) == set(GROUPS)
    assert int(d.msi_type_pentaplex.eq("MSI").sum()) == 6
    assert int(d.ebv_in_situ.eq("Positive").sum()) == 5
    assert d.pdl1_cps.notna().sum() == 41
    assert d.pdl1_cps.isna().sum() == 4
    assert d.pdl1_cps.dropna().between(0, 100).all()
    return d, dest


def make_summary(d):
    rows = []
    for group in GROUPS:
        g = d.loc[d.ecotype4.eq(group)]
        n, k = len(g), int(g.response_binary.sum())
        lo, hi = wilson(k,n)
        cps = g.pdl1_cps.dropna().to_numpy(float)
        rows.append({"configuration": group, "patients_n": n,
                     "response_n": k, "nonresponse_n": n-k, "response_fraction": k/n,
                     "response_Wilson95_low": lo, "response_Wilson95_high": hi,
                     "MSI_positive_n": int(g.msi_type_pentaplex.eq("MSI").sum()),
                     "MSI_negative_MSS_n": int(g.msi_type_pentaplex.eq("MSS").sum()),
                     "MSI_missing_n": int((~g.msi_type_pentaplex.isin(["MSI","MSS"])).sum()),
                     "EBV_ISH_positive_n": int(g.ebv_in_situ.eq("Positive").sum()),
                     "EBV_ISH_negative_n": int(g.ebv_in_situ.eq("Negative").sum()),
                     "EBV_ISH_missing_n": int((~g.ebv_in_situ.isin(["Positive","Negative"])).sum()),
                     "CPS_available_n": len(cps), "CPS_missing_n": n-len(cps),
                     "CPS_median": float(np.median(cps)),
                     "CPS_q1": float(np.quantile(cps,.25,method="linear")),
                     "CPS_q3": float(np.quantile(cps,.75,method="linear"))})
    s = pd.DataFrame(rows)
    s.to_csv(ROOT / "configuration_panel_numeric_summary.tsv", sep="\t", index=False, float_format="%.17g")
    return s


def base_axes(ax, letter, title, xlabel):
    ax.set_title(title, loc="left", fontweight="bold", pad=19)
    ax.text(-.28, 1.135, letter, transform=ax.transAxes, fontweight="bold", fontsize=10)
    ax.set_yticks([0,1,2,3],DISPLAY)
    ax.set_ylim(3.5,-.65)
    ax.tick_params(axis="y",length=0,pad=7)
    ax.set_xlim(-2,102)
    ax.set_xticks([0,25,50,75,100])
    ax.set_xlabel(xlabel, labelpad=6)
    ax.spines["left"].set_visible(False)


def response(ax,s):
    base_axes(ax,"A","Observed anti-PD-1 response","Response (%)")
    for i,r in enumerate(s.itertuples(index=False)):
        value = r.response_fraction*100
        ax.errorbar(value,i,xerr=[[value-r.response_Wilson95_low*100],[r.response_Wilson95_high*100-value]],
                    fmt="o",ms=5.5,mfc=GROUP_COLORS[r.configuration],mec=BLACK,mew=.5,
                    ecolor=BLACK,elinewidth=.8,capsize=2.5,zorder=3)
        ax.text(50,i-.32,f"{r.response_n}/{r.patients_n} responders",ha="center",va="bottom",fontsize=7)
    ax.text(0,1.035,"Points and Wilson 95% confidence intervals",transform=ax.transAxes,fontsize=6.7,color="#555555")


def composition(ax,s,letter,kind):
    is_msi = kind == "MSI"
    title = "MSI status (pentaplex)" if is_msi else "EBV status (in situ hybridisation)"
    base_axes(ax,letter,title,"Patients within configuration (%)")
    cols = ["MSI_positive_n","MSI_negative_MSS_n","MSI_missing_n"] if is_msi else ["EBV_ISH_positive_n","EBV_ISH_negative_n","EBV_ISH_missing_n"]
    labels = ["MSI", "MSS", "Missing"] if is_msi else ["Positive", "Negative", "Missing"]
    for i,r in enumerate(s.to_dict(orient="records")):
        start = 0
        for col,colr in zip(cols,[POSITIVE,NEGATIVE,MISSING]):
            value = r[col]/r["patients_n"]*100
            ax.barh(i,value,left=start,height=.36,color=colr,edgecolor="white",linewidth=.5)
            start += value
        text = f"{labels[0]} {r[cols[0]]}/{r['patients_n']}   |   {labels[1]} {r[cols[1]]}/{r['patients_n']}"
        ax.text(50,i-.27,text,ha="center",va="bottom",fontsize=6.6)
    handles = [Patch(facecolor=POSITIVE,edgecolor="none",label=labels[0]), Patch(facecolor=NEGATIVE,edgecolor="none",label=labels[1])]
    ax.legend(handles=handles,loc="lower left",bbox_to_anchor=(-.02,1.005),ncol=2,fontsize=6.7,
              handlelength=1,handletextpad=.4,columnspacing=1.2,borderaxespad=0)
    missing = int(s[cols[2]].sum())
    ax.text(1,-.23,f"Missing status: {missing}/45",transform=ax.transAxes,ha="right",va="top",fontsize=6.5,color="#555555")


def offsets_for_identical(values):
    offsets = np.zeros(len(values),float)
    for value in np.unique(values):
        ix = np.where(values == value)[0]
        offsets[ix] = np.linspace(-.19,.19,len(ix)) if len(ix)>1 else 0
    return offsets


def cps(ax,s,d):
    base_axes(ax,"D","PD-L1 combined positive score","PD-L1 CPS")
    plotted = []
    for i,r in enumerate(s.itertuples(index=False)):
        g = d.loc[d.ecotype4.eq(r.configuration) & d.pdl1_cps.notna()].sort_values(["pdl1_cps","patient_id"])
        values = g.pdl1_cps.to_numpy(float)
        yy = i+offsets_for_identical(values)
        ax.scatter(values,yy,s=16,c=GROUP_COLORS[r.configuration],edgecolors=BLACK,linewidths=.35,zorder=3)
        # Put summary just below points so zero-heavy samples remain countable.
        ysum = i+.30
        ax.plot([r.CPS_q1,r.CPS_q3],[ysum,ysum],color=BLACK,lw=1.8,zorder=4)
        ax.plot([r.CPS_median,r.CPS_median],[ysum-.045,ysum+.045],color=BLACK,lw=1.2,zorder=4)
        ax.text(50,i-.32,f"Available {r.CPS_available_n}/{r.patients_n}; missing {r.CPS_missing_n}",ha="center",va="bottom",fontsize=6.5)
        for p,y in zip(g.itertuples(index=False),yy):
            plotted.append({"patient_id":p.patient_id,"configuration":r.configuration,"CPS":p.pdl1_cps,"plot_y":y})
    ax.text(0,1.035,"All available patients; black bars: median and IQR",transform=ax.transAxes,fontsize=6.7,color="#555555")
    ax.text(1,-.23,"41 observed CPS values; 4 missing",transform=ax.transAxes,ha="right",va="top",fontsize=6.5,color="#555555")
    pd.DataFrame(plotted).to_csv(ROOT / "CPS_individual_plot_source.tsv",sep="\t",index=False,float_format="%.17g")
    assert len(plotted) == 41 and len({x["patient_id"] for x in plotted}) == 41


def draw(d,s):
    fig,axes = plt.subplots(2,2,figsize=(WIDTH_MM/25.4,HEIGHT_MM/25.4))
    # Fixed canvas keeps the physical output width at exactly 183 mm.
    fig.subplots_adjust(left=.16,right=.98,bottom=.13,top=.84,wspace=.61,hspace=.72)
    response(axes[0,0],s)
    composition(axes[0,1],s,"B","MSI")
    composition(axes[1,0],s,"C","EBV")
    cps(axes[1,1],s,d)
    fig.suptitle("Kim2018: response and established biomarkers across ecological configurations",x=.05,y=.976,ha="left",fontsize=9,fontweight="bold")
    fig.text(.05,.939,"Frozen configurations; 45 pretreatment tumours, 12 responders",fontsize=7.5,color="#555555")
    fig.text(.05,.022,"Descriptive distributions in a single-arm cohort; these panels do not establish biomarker-independent prediction.",fontsize=6.5,color="#555555")
    with mpl.rc_context({"savefig.bbox":None}):
        for suffix in ["pdf","svg","png"]:
            # Disable tight cropping to retain the exact physical journal width.
            fig.savefig(ROOT / f"Kim2018_configuration_biomarker_overview.{suffix}",dpi=300,bbox_inches=None,facecolor="white")
    plt.close(fig)


def document(d,s,input_file):
    legend = (
        "Figure. Established biomarker distributions across ecological configurations in the Kim2018 anti-PD-1 cohort. "
        "A, Observed response proportions in the four frozen, outcome-independent median-defined configurations; "
        "points show responders/total patients and whiskers show descriptive Wilson 95% binomial confidence intervals. "
        "B, Pentaplex MSI status, classified as MSI or microsatellite stable (MSS). C, EBV status by in situ hybridisation, "
        "classified as positive or negative. Bar labels report counts and configuration-specific denominators; MSI "
        "and EBV-ISH status were complete for all 45 patients. D, Individual available PD-L1 combined positive scores "
        "(CPS; n = 41) with black median and interquartile-range summaries. Identical CPS values are offset vertically "
        "for visibility without changing their values. Four missing CPS measurements are explicitly counted and are "
        "not classified as zero or negative. No between-group hypothesis tests were performed. These descriptive "
        "distributions contextualise confounding by established biomarkers in a small single-arm cohort and do not "
        "establish independent predictive performance or treatment benefit. Configuration colours follow the original Figure 6."
    )
    (ROOT / "Figure_legend_EN.txt").write_text(legend+"\n",encoding="utf-8")
    manifest = {
        "source":str(SOURCE),"snapshot":input_file.name,"snapshot_sha256":sha256(input_file),
        "source_sha256":sha256(SOURCE),"source_equals_snapshot":sha256(SOURCE)==sha256(input_file),
        "script_sha256":sha256(Path(__file__)),
        "n_patients":len(d),"responders":int(d.response_binary.sum()),"MSI_positive":int(d.msi_type_pentaplex.eq("MSI").sum()),
        "EBV_ISH_positive":int(d.ebv_in_situ.eq("Positive").sum()),"CPS_observed":int(d.pdl1_cps.notna().sum()),"CPS_missing":int(d.pdl1_cps.isna().sum()),
        "configuration_order":GROUPS,"group_colors_from_original_Figure6":GROUP_COLORS,
        "interval":"Wilson two-sided 95% score interval without continuity correction",
        "hypothesis_tests":0,"CPS_quantiles":"linear interpolation; R type 7 / numpy method=linear",
        "input_rows_excluded":0,"missing_policy":"Retained in group denominators; no CPS point imputed for four missing measurements",
        "width_mm":WIDTH_MM,"height_mm":HEIGHT_MM,"dpi":300,"python":platform.python_version(),
        "numpy":np.__version__,"pandas":pd.__version__,"matplotlib":mpl.__version__,
    }
    (ROOT / "figure_provenance_QA.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    qa = ["# Figure QA: configuration overview", "",
          "Scientific message: response distributions coexist with uneven MSI, EBV-ISH and CPS distributions across the four configurations.",
          "Archetype: quantitative grid; A=response, B=MSI, C=EBV, D=CPS. Each panel has a distinct endpoint.",
          "Author requested candidates directly; no additional figure-plan approval was needed.", "",
          "- All 45 unique patient records included; fixed original group order and identity retained.",
          "- Responses=12, MSI-positive=6, EBV-ISH-positive=5, CPS observed=41, CPS missing=4.",
          "- MSI and EBV use the measured pentaplex/ISH assays, not WES estimates.",
          "- No tests or performance metrics calculated. Wilson 95% intervals are descriptive binomial intervals.",
          "- All 41 available CPS values plotted. Four missing values explicitly labelled; no imputation or negative classification.",
          "- All bars originate at zero. Response and binary-composition ranges include 0–100%; CPS remains raw 0–100.",
          "- Typography baseline copied; font floor 6.5 pt; Arial; top/right spines removed; legends outside data.",
          "- Semantic status palette and existing Figure 6 phenotype colours, plus direct labels for colour redundancy.",
          "- PDF/SVG preserve vector content and fonts; PNG preview at 300 dpi; exact canvas 183 × 165 mm.",
          "- Asset scripts and companion PNGs inspected. They used incompatible paired/model structures, so their visual primitives were inherited without simulated data or inappropriate statistical overlays.",
          "- Visual render inspected after layout repair: headings and notes separated; no data clipping; direct count labels legible. Repeated CPS zero values retain all individual points with deterministic vertical offsets, although dense zero stacks partially overlap.",
          "- PDF page dimensions and PNG resolution checked: 183 × 165 mm, 300 dpi. No raw data values are changed by plotting offsets.", ""]
    (ROOT / "QA_report.md").write_text("\n".join(qa),encoding="utf-8")
    (ROOT / "README.md").write_text(
        "# Kim2018 configuration overview\n\nRun `python draw_configuration_overview.py` in an environment with numpy, pandas and matplotlib. "
        "The script uses the saved input snapshot; if absent, it copies the sibling 02_clinical_stratification snapshot. "
        "No formal manuscript, existing table or previous analysis output is edited.\n\n"+legend+"\n",encoding="utf-8")


def main():
    d,input_file = read_input()
    s = make_summary(d)
    draw(d,s)
    document(d,s,input_file)
    print(s.to_string(index=False))
    print("Wrote four-panel overview: PDF, SVG, PNG; no hypothesis tests performed.")


if __name__ == "__main__":
    main()
