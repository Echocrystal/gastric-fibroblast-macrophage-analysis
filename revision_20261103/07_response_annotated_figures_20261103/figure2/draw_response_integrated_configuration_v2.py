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

# Asset audit: original 04 configuration script is copied alongside this script.
# heatmap/plot_comparison.py and companion inspected: semantic mismatch with a
# binary patient annotation matrix; inherit direct cell alignment and labels.
# GroupedBarChart/plot_GroupedBarChartv1.py and companion inspected: model AUROC
# bars are semantically incompatible; inherit clean CI display, use binomial
# proportions with Wilson intervals from the frozen patient records.
# CPS inherits all-patient point and median/IQR display from the 04 figure.
# User requested response integration; no response model or new test is fitted.
from pathlib import Path
import hashlib
import json
import math
import platform
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SOURCE=ROOT/'02_clinical_stratification/Kim2018_patients_with_restriction_flags.tsv'
SNAPSHOT=HERE/'Kim2018_input_snapshot.tsv'
UPSTREAM_SCRIPT=ROOT/'04_clinical_biomarker_relationships/configuration_overview/draw_configuration_overview.py'
FIGURE='Figure2_configuration_response_integrated_v2'
GROUPS=['Dual-low','Fibroblast-dominant','Macrophage-dominant','Dual-high']
DISPLAY=['Dual-low','Fibroblast-\ndominant','Macrophage-\ndominant','Dual-high']
GROUP_COLORS=['#AAB7C4','#CDACB4','#8FABDA','#9FB9A1']
RESP=CATEGORICAL[3]
NONRESP=CATEGORICAL[0]
YES='#5C6268'
NO='#EDF0F2'
W_MM,H_MM=183,178

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def wilson(k,n):
    z=1.959963984540054
    pp=k/n; den=1+z*z/n
    mid=(pp+z*z/(2*n))/den
    half=z*math.sqrt(pp*(1-pp)/n+z*z/(4*n*n))/den
    return max(0.,mid-half) if k else 0., min(1.,mid+half) if k<n else 1.

def load():
    if not SNAPSHOT.exists(): shutil.copy2(SOURCE,SNAPSHOT)
    if not (HERE/'original_configuration_script_reference.py').exists():
        shutil.copy2(UPSTREAM_SCRIPT,HERE/'original_configuration_script_reference.py')
    d=pd.read_csv(SNAPSHOT,sep='\t')
    assert len(d)==45 and d.patient_id.nunique()==45
    assert d.best_response.isin(['CR','PR','SD','PD']).all()
    assert (d.best_response.isin(['CR','PR']).astype(int)==d.response_binary).all()
    assert d.response_binary.sum()==12
    d['MSI_positive']=d.msi_type_pentaplex.eq('MSI')
    d['EBV_positive']=d.ebv_in_situ.eq('Positive')
    d['retained']=d.msi_type_pentaplex.eq('MSS') & d.ebv_in_situ.eq('Negative')
    assert d.msi_type_pentaplex.isin(['MSI','MSS']).all()
    assert d.ebv_in_situ.isin(['Positive','Negative']).all()
    assert (d.retained.astype(int)==d.mss_and_ebv_negative_flag).all()
    assert d.retained.sum()==34 and d.loc[d.retained,'response_binary'].sum()==2
    assert d.pdl1_cps.notna().sum()==41
    return d

def summarise(d):
    counts=[]
    for group in GROUPS:
        for pop,subset in [('All patients',d.ecotype4.eq(group)),('MSS and EBV-negative',d.ecotype4.eq(group)&d.retained)]:
            g=d.loc[subset]; n=len(g);k=int(g.response_binary.sum());lo,hi=wilson(k,n)
            counts.append(dict(configuration=group,population=pop,n=n,responders=k,nonresponders=n-k,
                proportion=k/n,Wilson95_low=lo,Wilson95_high=hi))
    s=pd.DataFrame(counts)
    s.to_csv(HERE/'response_full_vs_restricted.tsv',sep='\t',index=False,float_format='%.17g')
    return s

def matrix(ax,d):
    position=0.; allrecords=[]
    for gi,group in enumerate(GROUPS):
        g=d.loc[d.ecotype4.eq(group)].sort_values(
            ['response_binary','MSI_positive','EBV_positive','pdl1_cps','patient_id'],
            ascending=[False,False,False,True,True],na_position='last')
        n=len(g);start=position
        for rr in g.itertuples(index=False):
            x=position+.5
            ax.add_patch(Rectangle((position,-.36),1,.72,facecolor='#F7F8F9',edgecolor='white',linewidth=.45))
            ax.scatter(x,0,marker='^' if rr.response_binary else 'o',s=15 if rr.response_binary else 11,
                facecolors=RESP if rr.response_binary else 'white',edgecolors=RESP if rr.response_binary else NONRESP,
                linewidths=.65,zorder=4)
            for yy,value in [(1,rr.MSI_positive),(2,rr.EBV_positive),(3,rr.retained)]:
                ax.add_patch(Rectangle((position,yy-.36),1,.72,facecolor=YES if value else NO,
                    edgecolor='white',linewidth=.45))
            cps='–' if pd.isna(rr.pdl1_cps) else f'{rr.pdl1_cps:g}'
            ax.text(x,4,cps,ha='center',va='center',fontsize=5.5,color=BLACK)
            item=rr._asdict();item.update(matrix_column=len(allrecords)+1,plot_x=x)
            allrecords.append(item); position+=1
        center=start+n/2
        ax.add_patch(Rectangle((start,-1.04),n,.40,facecolor=GROUP_COLORS[gi],edgecolor='none'))
        short=DISPLAY[gi].replace('\n','')
        ax.text(center,-1.68,short,ha='center',va='center',fontsize=6.4,fontweight='bold')
        ax.text(center,-1.30,f'n = {n}',ha='center',va='center',fontsize=6)
        m=g.loc[g.MSI_positive];e=g.loc[g.EBV_positive]
        mt=f"MSI+ {len(m)}/{n}; R {int(m.response_binary.sum())}/{len(m)}" if len(m) else f'MSI+ 0/{n}'
        et=f"EBV+ {len(e)}/{n}; R {int(e.response_binary.sum())}/{len(e)}" if len(e) else f'EBV+ 0/{n}'
        ax.text(center,5.0,mt,ha='center',va='center',fontsize=5.9)
        ax.text(center,5.55,et,ha='center',va='center',fontsize=5.9)
        position+=1.1
    ax.set_xlim(-.3,position-1.1+.3);ax.set_ylim(6.0,-2.10)
    ax.set_yticks([0,1,2,3,4],['Response','MSI+','EBV+','MSS & EBV−','CPS'])
    ax.tick_params(left=False,bottom=False,labelbottom=False,labelsize=7,pad=5)
    for spine in ax.spines.values():spine.set_visible(False)
    pd.DataFrame(allrecords).to_csv(HERE/'patient_matrix_source.tsv',sep='\t',index=False,float_format='%.17g')
    assert len(allrecords)==45 and len(set(r['patient_id'] for r in allrecords))==45
    return len(allrecords)

def response_panel(ax,s):
    for i,g in enumerate(GROUPS):
        for pop,yy,marker,fill in [('All patients',i-.11,'o',BLACK),('MSS and EBV-negative',i+.17,'D','white')]:
            r=s.loc[s.configuration.eq(g)&s.population.eq(pop)].iloc[0]
            pp=r.proportion*100
            ax.errorbar(pp,yy,xerr=[[pp-r.Wilson95_low*100],[r.Wilson95_high*100-pp]],fmt=marker,
                markersize=4.3,mfc=fill,mec=BLACK,mew=.7,color=BLACK,elinewidth=.7,capsize=2,zorder=3)
            # Counts appear at fixed positions, avoiding inference from symbol size.
            ax.text(107,yy,f'{r.responders}/{r.n}',ha='left',va='center',fontsize=7,clip_on=False)
    ax.set_yticks(range(4),DISPLAY);ax.tick_params(axis='y',length=0,pad=6)
    ax.set_ylim(3.6,-.65);ax.set_xlim(-3,102);ax.set_xticks([0,25,50,75,100])
    ax.set_xlabel('Response (%)',labelpad=4)
    ax.spines['left'].set_visible(False)
    ax.text(107,-.47,'R / n',ha='left',va='center',fontsize=6.7,fontweight='bold',clip_on=False)
    ax.legend(handles=[Line2D([],[],marker='o',ms=4,color=BLACK,lw=.6,label='All patients'),
        Line2D([],[],marker='D',ms=4,color=BLACK,mfc='white',lw=.6,label='MSS & EBV−')],
        loc='upper left',bbox_to_anchor=(-.02,1.14),fontsize=6.4,ncol=2,columnspacing=1.1,
        handlelength=1.2,handletextpad=.4,borderaxespad=0)

def cps_panel(ax,d):
    records=[]
    for i,g in enumerate(GROUPS):
        group=d.loc[d.ecotype4.eq(g)]
        observed=group.loc[group.pdl1_cps.notna()].sort_values(['pdl1_cps','response_binary','patient_id'])
        values=observed.pdl1_cps.to_numpy(float)
        offsets=np.zeros(len(observed))
        for v in np.unique(values):
            ix=np.where(values==v)[0]
            offsets[ix]=np.linspace(-.20,.20,len(ix)) if len(ix)>1 else 0
        ys=i+offsets
        for response in [0,1]:
            ix=observed.response_binary.to_numpy()==response
            ax.scatter(values[ix],ys[ix],s=18 if response else 13,marker='^' if response else 'o',
                facecolors=RESP if response else 'white',edgecolors=RESP if response else NONRESP,linewidths=.6,zorder=3+response)
        q1,med,q3=np.quantile(values,[.25,.5,.75])
        ax.plot([q1,q3],[i+.32,i+.32],lw=1.6,color=BLACK)
        ax.plot([med,med],[i+.27,i+.37],lw=1.0,color=BLACK)
        missing=group.loc[group.pdl1_cps.isna()]
        # Account for missing responses explicitly without imputing CPS.
        ax.text(47,i-.33,f"Missing: {int(missing.response_binary.sum())} R / {len(missing)-int(missing.response_binary.sum())} NR",
            ha='center',va='center',fontsize=5.8,color='#555555')
        for row,yy in zip(observed.itertuples(index=False),ys):
            records.append(dict(patient_id=row.patient_id,configuration=g,CPS=row.pdl1_cps,
                response_binary=row.response_binary,best_response=row.best_response,plot_y=yy))
    ax.set_ylim(3.6,-.65);ax.set_xlim(-3,103);ax.set_xticks([0,25,50,75,100])
    ax.set_yticks(range(4),['']*4);ax.tick_params(axis='y',length=0)
    ax.spines['left'].set_visible(False);ax.set_xlabel('PD-L1 CPS',labelpad=4)
    ax.text(.01,1.105,'Points: R / NR; black: median and IQR',transform=ax.transAxes,fontsize=6.2)
    pd.DataFrame(records).to_csv(HERE/'CPS_response_source.tsv',sep='\t',index=False,float_format='%.17g')
    assert len(records)==41 and sum(r['response_binary'] for r in records)==11

def draw(d,s):
    fig=plt.figure(figsize=(W_MM/25.4,H_MM/25.4))
    fig.text(.04,.975,'Kim2018: response and clinical biomarkers within configurations',fontsize=9,fontweight='bold',va='top')
    fig.text(.04,.943,'45 patients: 12 responders (R), 33 non-responders (NR); frozen configurations',fontsize=7,va='top')
    fig.legend(handles=[Line2D([],[],marker='^',color=RESP,mfc=RESP,lw=0,ms=5,label='Responder (CR/PR)'),
        Line2D([],[],marker='o',color=NONRESP,mfc='white',lw=0,ms=4.5,label='Non-responder (SD/PD)'),
        Patch(facecolor=YES,label='Yes'),Patch(facecolor=NO,label='No')],
        loc='upper left',bbox_to_anchor=(.035,.924),ncol=4,fontsize=6.5,columnspacing=1.5,
        handlelength=1,handletextpad=.5,borderaxespad=0)
    fig.text(.04,.875,'A',fontsize=10,fontweight='bold')
    fig.text(.075,.875,'Patient-level map',fontsize=8,fontweight='bold')
    fig.text(.40,.875,'Each column is one patient; dark/light cells = yes/no',fontsize=6.5)
    top=fig.add_axes([.16,.565,.815,.295]);matrix(top,d)
    fig.text(.16,.542,'R = responders within the indicated biomarker-positive subgroup; CPS “–” = missing.',fontsize=6.1)
    fig.text(.04,.492,'B',fontsize=10,fontweight='bold')
    fig.text(.075,.492,'Full vs restricted response',fontsize=8,fontweight='bold')
    fig.text(.62,.492,'C',fontsize=10,fontweight='bold')
    fig.text(.655,.492,'CPS by response',fontsize=8,fontweight='bold')
    ax1=fig.add_axes([.205,.125,.325,.305]);response_panel(ax1,s)
    ax2=fig.add_axes([.665,.125,.31,.305]);cps_panel(ax2,d)
    fig.text(.04,.060,'Restriction: MSS and EBV ISH-negative, 34 patients / 2 responders. Whiskers in B: Wilson 95% intervals.',fontsize=6.2)
    fig.text(.04,.032,'Descriptive display; the restricted group is nested within all patients. No new hypothesis tests or predictive model.',fontsize=6.2)
    with mpl.rc_context({'savefig.bbox':None}):
        for ext in ['png','pdf','svg']:
            fig.savefig(HERE/f'{FIGURE}.{ext}',dpi=300,bbox_inches=None,facecolor='white')
    plt.close(fig)

def main():
    d=load();s=summarise(d);draw(d,s)
    legend=(
        'Figure. Treatment response and established biomarkers within the four frozen ecological configurations in Kim2018. '
        'A, Patient-level annotation map: each column represents the same patient across rows; columns are grouped by configuration '
        'and ordered within groups by response, MSI, EBV, CPS and patient identifier. Orange filled triangles denote responders '
        '(CR/PR); blue open circles denote non-responders (SD/PD). Dark and light cells indicate yes and no for the row-labelled '
        'MSI-positive, EBV ISH-positive and joint MSS/EBV-negative criteria. The CPS row reports observed numeric values, with an '
        'en dash for missing values. Under each configuration, biomarker-positive counts are shown relative to all group members, '
        'followed by responders/biomarker-positive patients where that subgroup exists. B, Observed response proportions in all '
        'patients (filled circles) and the nested restriction to MSS and EBV ISH-negative patients (open diamonds), with descriptive '
        'Wilson 95% binomial confidence intervals and responders/total counts. The restricted population contains 34 patients and '
        'two responders. C, All 41 available CPS values, using the same response symbols as A; black bars show medians and '
        'interquartile ranges. Deterministic vertical offsets separate identical values without changing CPS. Missing CPS counts '
        'are given separately for responders and non-responders; no missing value is imputed. All panels retain the original '
        'configuration order. These are descriptive distributions, with no new hypothesis tests; they do not establish '
        'biomarker-independent response prediction or treatment benefit.'
    )
    (HERE/'Figure_legend_EN.txt').write_text(legend+'\n',encoding='utf-8')
    inputs=[SOURCE,SNAPSHOT,UPSTREAM_SCRIPT]
    outputs=[p for p in HERE.iterdir() if p.is_file() and p.name not in ['figure2_manifest_QA.json','QA_report.md']]
    manifest={'input_sources':[{ 'path':str(p),'sha256':sha(p)} for p in inputs],
        'outputs':[{'file':p.name,'sha256':sha(p)} for p in sorted(outputs)],
        'patients_all':45,'responders_all':12,'patients_restricted':34,'responders_restricted':2,
        'matrix_columns':45,'CPS_observed':41,'CPS_missing':4,'new_hypothesis_tests':0,
        'response_definition':'CR/PR=1; SD/PD=0; frozen upstream labels retained',
        'restriction':'MSS by pentaplex AND EBV ISH-negative',
        'configuration_order':GROUPS,'dimensions_mm':[W_MM,H_MM],
        'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
        'matplotlib':mpl.__version__,'visual_QA':'pending'}
    (HERE/'figure2_manifest_QA.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'figure':str(HERE/(FIGURE+'.png')),'patients':len(d),'restricted_n':int(d.retained.sum())}))

if __name__=='__main__':main()
