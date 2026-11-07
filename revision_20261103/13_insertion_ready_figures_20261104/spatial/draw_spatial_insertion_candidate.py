# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Forest → original 01/run_spatial_endpoint_comparison.py → param inherit
# Native grouped-bar assets are semantically incompatible with adjusted-rho CIs.
# Preserve original forest colours, marker logic and grouping; all numeric inputs frozen.
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

from pathlib import Path
import hashlib,json,sys
import matplotlib
matplotlib.use('Agg')
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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.stdout.reconfigure(encoding='utf-8')
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
PREVIOUS=ROOT/'12_SPON2_ITGB2_independent_cohort_pooling_REVIEW_20261104'
ENDPOINTS=['SPON2_ITGB2_receiver_composite_r1','SPON2_receiver_r1','SPON2_ITGB2_r1']
LABELS=['ITGB2 + receiver composite','Receiver score alone','ITGB2 alone']
COLORS=['#B45D35','#1B7837','#286B8A']
MARKERS=['o','^','s']
SHA=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    original_path=PREVIOUS/'original_cohort_results_3endpoints.tsv'
    pooled_path=PREVIOUS/'pooled_results_all_schemes.tsv'
    original=pd.read_csv(original_path,sep='\t')
    pooled=pd.read_csv(pooled_path,sep='\t')
    assert len(original)==9
    individual=pooled.loc[pooled.weighting.eq('individual_equal')].copy()
    assert len(individual)==3
    cohorts=['GSE270678_main_ST','GSE270678_ST_R_sensitivity','GSE251950_primary','Overall_28_individuals']
    display=['GSE270678 main ST\n19 donors','GSE270678 ST_R\n15 donors','GSE251950\n9 patients','Overall\n28 individuals']
    rows=[]
    for i,dataset in enumerate(cohorts):
        for j,endpoint in enumerate(ENDPOINTS):
            if i<3:
                r=original.loc[original.dataset.eq(dataset)&original.source_endpoint_id.eq(endpoint)].iloc[0]
                value=float(r.rho);p=float(r.nominal_two_sided_p);q=float(r.original_BH_q)
                family=r.original_family;ci_policy='Original 9,999-draw percentile bootstrap'
                n=int(r.n_individuals)
            else:
                r=individual.loc[individual.source_endpoint_id.eq(endpoint)].iloc[0]
                value=float(r.pooled_adjusted_rho);p=float(r.nominal_exhaustive_two_sided_p);q=float(r.new_BH_q)
                family=r.new_BH_family;ci_policy='New 99,999-draw within-cohort percentile bootstrap'
                n=int(r.n_individuals)
            rows.append(dict(display_dataset=dataset,source_endpoint_id=endpoint,endpoint=LABELS[j],
                n_individuals=n,rho=value,ci95_low=float(r.ci95_low),ci95_high=float(r.ci95_high),
                P=p,q=q,q_family=family,CI_policy=ci_policy,
                y=10.5-i*3.15+[.65,0,-.65][j],colour=COLORS[j],marker=MARKERS[j],
                pooled_weighting='individual_equal' if i==3 else 'original_within_cohort_individual_equal'))
    source=pd.DataFrame(rows)
    source.to_csv(HERE/'Spatial_endpoint_insertion_source_data.tsv',sep='\t',index=False,float_format='%.17g')
    # Fixed canvas preserves 183mm width and legible print fonts; caption external.
    fig=plt.figure(figsize=(183/25.4,112/25.4),facecolor='white')
    ax=fig.add_axes([.225,.14,.305,.70])
    stat=fig.add_axes([.565,.14,.425,.70],sharey=ax)
    ax.set_ylim(-.15,11.8);ax.set_xlim(-.006,.085);ax.axvline(0,color='#777777',lw=.8,ls='--',zorder=0)
    for row in source.itertuples(index=False):
        ax.errorbar(row.rho,row.y,xerr=[[row.rho-row.ci95_low],[row.ci95_high-row.rho]],
            fmt=row.marker,color=row.colour,markersize=4.0,capsize=2.2,lw=1.15)
        stat.text(0,row.y,f'{row.rho:.4f} [{row.ci95_low:.4f}, {row.ci95_high:.4f}]',va='center',fontsize=6.7)
        stat.text(.68,row.y,f'{row.P:.3g}',va='center',fontsize=6.7)
        stat.text(.88,row.y,f'{row.q:.3g}',va='center',fontsize=6.7)
    for y in [8.925,5.775,2.625]:
        ax.axhline(y,color='#DEE2E5',lw=.5,zorder=0)
        stat.axhline(y,color='#DEE2E5',lw=.5,zorder=0)
    ax.set_yticks([10.5,7.35,4.2,1.05],display,fontsize=7)
    ax.tick_params(axis='y',length=0,pad=9)
    ax.set_xticks([0,.02,.04,.06,.08]);ax.tick_params(axis='x',labelsize=6.8)
    ax.set_xlabel('Adjusted correlation',fontsize=7,labelpad=5)
    ax.spines[['top','right','left']].set_visible(False)
    stat.axis('off');stat.set_xlim(0,1.08)
    stat.text(0,11.7,'ρ [95% CI]',fontsize=7,fontweight='bold',va='bottom')
    stat.text(.68,11.7,'P',fontsize=7,fontweight='bold',va='bottom')
    stat.text(.88,11.7,'q',fontsize=7,fontweight='bold',va='bottom')
    fig.text(.025,.97,'Composite and component endpoints',fontsize=8.5,fontweight='bold',va='top')
    handles=[Line2D([],[],marker=m,color=c,lw=1.1,ms=4,label=l) for m,c,l in zip(MARKERS,COLORS,LABELS)]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.02,.935),ncol=3,fontsize=6.9,
        handlelength=1.8,handletextpad=.5,columnspacing=1.5,frameon=False)
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    outside=[]
    for text in fig.findobj(mpl.text.Text):
        if not text.get_visible() or not text.get_text():continue
        box=text.get_window_extent(renderer)
        if box.x0<-.5 or box.y0<-.5 or box.x1>fig.bbox.width+.5 or box.y1>fig.bbox.height+.5:
            outside.append(text.get_text())
    assert not outside,outside
    for ext in ['pdf','svg','png']:
        # Export exactly at designed dimensions, no automatic canvas expansion.
        with mpl.rc_context({'savefig.bbox':None}):
            fig.savefig(HERE/f'Spatial_endpoints_with_overall_INSERTION.{ext}',dpi=300,bbox_inches=None,facecolor='white')
    plt.close(fig)
    legend = 'Supplemental Figure 10D. SPON2 associations with the neighbouring ITGB2-plus-receiver composite, receiver score alone and ITGB2 alone in GSE270678 main ST, ST_R and GSE251950. Estimates are adjusted pooled correlations with 95% bootstrap confidence intervals. Overall combines 19 main-ST donors and nine GSE251950 patients with equal weight per individual; overlapping ST_R donors are excluded. P values are nominal; q values follow the correction families specified in Supplemental Method 21.'
    (HERE/'Figure_legend_EN.txt').write_text(legend+'\n',encoding='utf-8')
    explanation=(
        'Pooled equal individuals：每个纳入的人贡献相同权重。Main ST19人占19/28≈67.9%，'
        'GSE2519509人占9/28≈32.1%。图中改为Overall (28 individuals)，图注解释每人等权。\n\n'
        'Pooled equal cohorts：两个队列各占50%权重，因此9人队列中每个人的权重大于19人队列中的个人。'
        '它是本次追加汇总中的另一种权重分析；按用户要求不单独作图。原12目录Excel/TSV/脚本继续完整保留。\n\n'
        '本图与原始01森林图同为按队列分组的横向并列CI图，新加入receiver-alone与Overall。'
        '整图长标题和图底叙述全部移除，保留小标题、数据图例、队列及样本数、坐标和数值列。'
        '必要解释置于独立英文图注。未将图片插入正式SI或替换原图。\n'
    )
    (HERE/'Weighting_explanation_CN.txt').write_text(explanation,encoding='utf-8')
    qa=dict(status='PASS_pending_visual_review',source_hashes={str(original_path):SHA(original_path),str(pooled_path):SHA(pooled_path)},
        plotted_rows=12,original_cohort_rows=9,new_individual_equal_rows=3,equal_cohort_rows_displayed=0,
        numerical_values_recomputed=False,new_statistical_tests=False,ST_R_in_overall=False,
        dimensions_mm=[183,112],minimum_font_pt=6.7,all_text_within_canvas=True,
        full_title_removed=True,narrative_footer_removed=True,panel_title='Composite and component endpoints',
        colour_shape_mapping=dict(zip(LABELS,list(zip(COLORS,MARKERS)))),
        styles_from='01/run_spatial_endpoint_comparison.py; original forest PNG inspected',
        removed_content='Figure-wide long title and narrative footer; weighting explanations are external.',
        retained='Panel title, data legend, cohort labels/n, axes and rho/CI/P/q values')
    (HERE/'figure_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(qa,ensure_ascii=False))

if __name__=='__main__':main()
