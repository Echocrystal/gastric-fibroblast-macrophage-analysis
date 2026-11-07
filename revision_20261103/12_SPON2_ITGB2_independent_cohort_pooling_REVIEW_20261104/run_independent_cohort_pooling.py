"""Additional post hoc pooling of frozen, adjusted spatial donor coefficients.

Run with Python3,numpy,pandas,matplotlib. Never modifies previous outputs.
The statistical specification is frozen in ANALYSIS_PROTOCOL_BEFORE_POOLING.md.
"""
from __future__ import annotations
import hashlib, itertools, json, platform, shutil, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding='utf-8')
HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'01_spatial_endpoint'
COHORTS=['GSE270678_main_ST','GSE251950_primary']
ST_R='GSE270678_ST_R_sensitivity'
ENDPOINTS=['SPON2_ITGB2_receiver_composite_r1','SPON2_receiver_r1','SPON2_ITGB2_r1']
LABELS=dict(zip(ENDPOINTS,['ITGB2 + receiver composite','Receiver score alone','ITGB2 alone']))
SEED=2026110412
DRAWS=99999
TOLERANCE=1e-15
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def snapshot():
    mapping={f'inputs/{n}':OLD/'inputs'/n for n in [
        'SPON2_ITGB2_donor_effects.tsv','SPON2_ITGB2_meta_summary.tsv',
        'SPON2_ITGB2_section_effects.tsv','SPON2_ITGB2_gene_endpoint_coverage.tsv']}
    mapping['scripts/original_54_spatial_candidate_context.py']=OLD/'scripts/original_54_spatial_candidate_context.py'
    mapping['scripts/original_35_spatial_chain.py']=OLD/'scripts/original_35_spatial_chain.py'
    provenance=[]
    for rel,source in mapping.items():
        target=HERE/rel;target.parent.mkdir(parents=True,exist_ok=True)
        if not target.exists():shutil.copy2(source,target)
        assert sha(target)==sha(source),'Frozen input changed: '+rel
        provenance.append(dict(source=str(source),snapshot=rel,sha256=sha(target),bytes=target.stat().st_size))
    pd.DataFrame(provenance).to_csv(HERE/'input_provenance_sha256.tsv',sep='\t',index=False)
    return provenance

def bh(p):
    p=np.asarray(p,float);order=np.argsort(p)
    ranked=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out=np.empty_like(p);out[order]=np.minimum(ranked,1);return out

def signed_sums(a):
    sums=np.zeros(1)
    for x in a:sums=np.concatenate([sums-x,sums+x])
    return sums

def signflip_exact(z,weights):
    z=np.asarray(z,float);weights=np.asarray(weights,float)
    assert len(z)==len(weights) and np.isfinite(z).all()
    assert np.all(weights>0) and abs(weights.sum()-1)<1e-12
    v=z*weights;observed=abs(float(v.sum()));n=len(v);total=1<<n
    if observed<=TOLERANCE:return dict(p=1.,tail_count=total,total_assignments=total,observed=observed)
    cut=n//2;left=signed_sums(v[:cut]);right=np.sort(signed_sums(v[cut:]))
    upper=right.size-np.searchsorted(right,observed-TOLERANCE-left,side='left')
    lower=np.searchsorted(right,-observed+TOLERANCE-left,side='right')
    count=int(upper.sum(dtype=np.int64)+lower.sum(dtype=np.int64))
    assert 0<=count<=total
    return dict(p=count/total,tail_count=count,total_assignments=total,observed=observed)

def brute(z,w):
    null=np.array([np.dot(s,np.asarray(z)*w) for s in itertools.product([-1.,1.],repeat=len(z))])
    observed=abs(float(np.dot(z,w)))
    return int((np.abs(null)>=observed-TOLERANCE).sum())

def validate_algorithm():
    rng=np.random.default_rng(60412);checks=[]
    cases=[np.array([0.,0.,0.]),np.array([.1,-.1,.2,-.2]),np.array([.03,.04,.05])]
    cases.extend(rng.normal(size=n)*.1 for n in [1,2,5,8,12,16])
    for z in cases:
        for kind in ['equal','unequal']:
            w=np.ones(len(z)) if kind=='equal' else np.arange(1,len(z)+1,dtype=float)
            w=w/w.sum();got=signflip_exact(z,w);expected=brute(z,w)
            assert got['tail_count']==expected,(z,w,got,expected)
            checks.append(dict(n=len(z),weight_kind=kind,count=expected,passed=True))
    return checks

def inputs_and_verify():
    donor=pd.read_csv(HERE/'inputs/SPON2_ITGB2_donor_effects.tsv',sep='\t',keep_default_na=False)
    meta=pd.read_csv(HERE/'inputs/SPON2_ITGB2_meta_summary.tsv',sep='\t',keep_default_na=False)
    section=pd.read_csv(HERE/'inputs/SPON2_ITGB2_section_effects.tsv',sep='\t',keep_default_na=False)
    subset=donor.loc[donor.endpoint.isin(ENDPOINTS)].copy()
    assert not subset.duplicated(['dataset','donor_id','endpoint']).any()
    selected=subset.loc[subset.dataset.isin(COHORTS)]
    expected={COHORTS[0]:19,COHORTS[1]:9,ST_R:15}
    checks=[];cohort_rows=[]
    for dataset in COHORTS+[ST_R]:
        for endpoint in ENDPOINTS:
            ds=subset.loc[subset.dataset.eq(dataset)&subset.endpoint.eq(endpoint)].sort_values('donor_id')
            assert len(ds)==expected[dataset] and np.isfinite(ds.rho).all()
            gate=section.section_gate.astype(str).str.lower().eq('true')
            secs=section.loc[section.dataset.eq(dataset)&section.endpoint.eq(endpoint)&gate&(section.n_spots>=100)].copy()
            secs['rho']=pd.to_numeric(secs.rho)
            derived=secs.groupby('donor_id').rho.apply(lambda v:np.tanh(np.arctanh(np.clip(v,-.999999,.999999)).mean()))
            delta=float(np.max(np.abs(ds.set_index('donor_id').rho-derived)))
            assert delta<1e-12
            source=meta.loc[meta.dataset.eq(dataset)&meta.endpoint.eq(endpoint)].iloc[0]
            z=np.arctanh(np.clip(ds.rho.to_numpy(),-.999999,.999999));r=float(np.tanh(z.mean()))
            exact=signflip_exact(z,np.ones(len(z))/len(z))
            assert abs(r-float(source.pooled_rho))<1e-12
            assert abs(exact['p']-float(source.exact_signflip_p))<1e-12
            if endpoint==ENDPOINTS[0]:
                qrow=meta.loc[meta.dataset.eq(dataset)&meta.endpoint.eq('candidate_SPON2_ITGB2_r1')].iloc[0]
                q=float(qrow.candidate_BH_FDR);family='Original four candidate-axis tests within cohort';tests=4
            else:q=float(source.supporting_BH_FDR);family='Original three supporting-endpoint tests within cohort';tests=3
            cohort_rows.append(dict(dataset=dataset,endpoint=LABELS[endpoint],source_endpoint_id=endpoint,
                n_individuals=len(ds),n_sections=secs.section_id.nunique(),rho=r,
                ci95_low=float(source.ci_low),ci95_high=float(source.ci_high),
                nominal_two_sided_p=exact['p'],original_BH_q=q,original_family=family,original_family_tests=tests,
                included_in_new_pool=dataset in COHORTS,positive_individuals=int((ds.rho>0).sum()),negative_individuals=int((ds.rho<0).sum())))
            checks.append(dict(dataset=dataset,endpoint=endpoint,donor_max_abs_difference=delta,
                cohort_rho_difference=abs(r-float(source.pooled_rho)),p_difference=abs(exact['p']-float(source.exact_signflip_p)),passed=True))
    cohort_results=pd.DataFrame(cohort_rows)
    wide=selected.pivot(index=['dataset','donor_id'],columns='endpoint',values='rho').reindex(columns=ENDPOINTS).sort_index()
    ordered_keys=[(c,person) for c in COHORTS for person in sorted(wide.loc[c].index)]
    wide=wide.reindex(pd.MultiIndex.from_tuples(ordered_keys,names=wide.index.names))
    assert wide.shape==(28,3) and np.isfinite(wide.values).all()
    assert set(wide.loc[COHORTS[0]].index).isdisjoint(set(wide.loc[COHORTS[1]].index))
    return wide,cohort_results,checks

def compute_pool(wide):
    rng=np.random.default_rng(SEED);cohort_boot=[];cohort_z=[];sizes=[]
    # All endpoint vectors are resampled together; same draw under both schemes.
    for dataset in COHORTS:
        z=np.arctanh(np.clip(wide.loc[dataset].to_numpy(float),-.999999,.999999))
        sizes.append(len(z));cohort_z.append(z)
        boot=np.empty((DRAWS,len(ENDPOINTS)))
        for start in range(0,DRAWS,10000):
            end=min(start+10000,DRAWS);idx=rng.integers(0,len(z),size=(end-start,len(z)))
            boot[start:end]=z[idx].mean(axis=1)
        cohort_boot.append(boot)
    z=np.vstack(cohort_z);assert sizes==[19,9]
    records=[];lopo=[]
    for scheme in ['individual_equal','cohort_equal']:
        shares=np.array(sizes,dtype=float)/sum(sizes) if scheme=='individual_equal' else np.array([.5,.5])
        weights=np.concatenate([np.repeat(shares[i]/sizes[i],sizes[i]) for i in range(2)])
        boot=sum(shares[i]*cohort_boot[i] for i in range(2))
        for j,endpoint in enumerate(ENDPOINTS):
            exact=signflip_exact(z[:,j],weights);estimate=float(np.tanh(weights@z[:,j]))
            lo,hi=np.quantile(np.tanh(boot[:,j]),[.025,.975],method='linear')
            omissions=[]
            for pos,key in enumerate(wide.index):
                keep=np.arange(len(z))!=pos
                if scheme=='individual_equal':w=np.ones(keep.sum())/keep.sum()
                else:
                    remaining=np.array([len(wide.loc[c])-(key[0]==c) for c in COHORTS])
                    w=np.array([.5/remaining[COHORTS.index(k[0])] for k in wide.index[keep]])
                value=float(np.tanh(w@z[keep,j]));omissions.append(value)
                lopo.append(dict(weighting=scheme,source_endpoint_id=endpoint,omitted_dataset=key[0],omitted_individual=key[1],pooled_rho=value))
            records.append(dict(weighting=scheme,endpoint=LABELS[endpoint],source_endpoint_id=endpoint,
                n_individuals=28,n_independent_cohorts=2,main_weight=float(shares[0]),external_weight=float(shares[1]),
                pooled_adjusted_rho=estimate,ci95_low=float(lo),ci95_high=float(hi),
                nominal_exhaustive_two_sided_p=exact['p'],tail_assignments=exact['tail_count'],total_assignments=exact['total_assignments'],
                positive_individuals=int((z[:,j]>0).sum()),negative_individuals=int((z[:,j]<0).sum()),
                lopo_rho_min=min(omissions),lopo_rho_max=max(omissions),
                bootstrap_draws=DRAWS,bootstrap_seed=SEED,
                CI_method='Pointwise percentile bootstrap; individuals resampled within each of two fixed cohorts',
                P_method='Exhaustive individual sign-flip on weighted Fisher-z; independent zero-symmetric-null assumption'))
    results=pd.DataFrame(records)
    results['new_BH_q']=np.nan
    for scheme,group in results.groupby('weighting'):
        results.loc[group.index,'new_BH_q']=bh(group.nominal_exhaustive_two_sided_p)
    results['new_BH_family']='Three pooled endpoints; separate family for each weighting scheme; post hoc exploratory'
    return results,pd.DataFrame(lopo)

def forest(results,cohort):
    plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'axes.titlesize':8,
        'xtick.labelsize':6.5,'ytick.labelsize':7,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,3,figsize=(183/25.4,105/25.4),sharey=True)
    for j,(ax,endpoint) in enumerate(zip(axes,ENDPOINTS)):
        for k,dataset in enumerate(COHORTS):
            r=cohort.loc[cohort.dataset.eq(dataset)&cohort.source_endpoint_id.eq(endpoint)].iloc[0]
            ax.errorbar(r.rho,3-k,xerr=[[r.rho-r.ci95_low],[r.ci95_high-r.rho]],fmt='o',color='#607687',ms=4,capsize=2,lw=1)
        for k,scheme in enumerate(['individual_equal','cohort_equal']):
            r=results.loc[results.weighting.eq(scheme)&results.source_endpoint_id.eq(endpoint)].iloc[0]
            ax.errorbar(r.pooled_adjusted_rho,1-k,
                xerr=[[r.pooled_adjusted_rho-r.ci95_low],[r.ci95_high-r.pooled_adjusted_rho]],
                fmt='D' if k==0 else 's',mfc='#B45D35' if k==0 else 'white',mec='#B45D35',ecolor='#B45D35',ms=4,capsize=2,lw=1.2)
            p_label=f'{r.nominal_exhaustive_two_sided_p:.3g}'
            q_label=f'{r.new_BH_q:.3g}'
            ax.text(.02, .1 if k==0 else .01, f'{"Individuals" if k==0 else "Cohorts"}: P={p_label}; q={q_label}',
                transform=ax.transAxes,fontsize=6.1,va='bottom')
        ax.axvline(0,color='#777777',ls='--',lw=.8,zorder=0)
        ax.axhline(1.5,color='#D8DDE1',lw=.6,zorder=0)
        ax.set_title(LABELS[endpoint],pad=8)
        ax.set_xlim(-.012,.087);ax.set_xticks([0,.04,.08]);ax.set_ylim(-1.2,3.6)
        ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0)
        ax.set_xlabel('Adjusted association (95% CI)')
    axes[0].set_yticks([3,2,1,0],['Main ST (19 people)','GSE251950 (9 people)','Pooled: equal individuals','Pooled: equal cohorts'])
    fig.suptitle('SPON2 versus ring-1 endpoints: two independent cohorts',x=.025,ha='left',y=.98,fontsize=10)
    fig.text(.025,.025,'Post hoc review analysis. ST_R excluded from pooling because donors overlap. Original figures remain unchanged.',fontsize=6.5)
    fig.subplots_adjust(left=.24,right=.985,bottom=.17,top=.85,wspace=.22)
    out=HERE/'figures';out.mkdir(exist_ok=True)
    for ext in ['pdf','svg','png']:
        fig.savefig(out/f'Independent_cohort_pooling_REVIEW.{ext}',dpi=300,facecolor='white')
    plt.close(fig)

def main():
    assert (HERE/'ANALYSIS_PROTOCOL_BEFORE_POOLING.md').exists()
    protected={str(p):sha(p) for folder in [OLD,HERE.parent/'09_figures_tables_integration_REVIEW_20261103']
               for p in folder.rglob('*') if p.is_file()}
    provenance=snapshot();algorithm=validate_algorithm();wide,cohort,source_checks=inputs_and_verify()
    results,lopo=compute_pool(wide)
    for filename,df in [('pooled_results_all_schemes.tsv',results),('original_cohort_results_3endpoints.tsv',cohort),('leave_one_individual_out_effects.tsv',lopo)]:
        df.to_csv(HERE/filename,sep='\t',index=False,float_format='%.17g')
    wide.reset_index().to_csv(HERE/'included_28_individuals_3endpoints.tsv',sep='\t',index=False,float_format='%.17g')
    (HERE/'pooled_results_all_schemes.json').write_text(results.to_json(orient='records',indent=2,double_precision=15),encoding='utf-8')
    forest(results,cohort)
    assert all(sha(p)==value for p,value in protected.items())
    audit=dict(status='PASS',analysis_status='New post hoc exploratory review analysis; no formal integration',
        n_individuals=28,cohort_sizes=dict(zip(COHORTS,[19,9])),ST_R_included=False,
        raw_spot_refit=False,source_checks=source_checks,algorithm_bruteforce_checks=algorithm,
        seed=SEED,bootstrap_draws=DRAWS,exact_assignments=2**28,float_tolerance=TOLERANCE,
        input_provenance=provenance,source_script_sha256=sha(__file__),protocol_sha256=sha(HERE/'ANALYSIS_PROTOCOL_BEFORE_POOLING.md'),
        preserved_previous_files=protected,preserved_previous_files_unchanged=True,
        python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__)
    (HERE/'pooling_verification.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    (HERE/'requirements_analysis.txt').write_text(f'numpy=={np.__version__}\npandas=={pd.__version__}\nmatplotlib=={matplotlib.__version__}\n',encoding='utf-8')
    print(results[['weighting','endpoint','pooled_adjusted_rho','ci95_low','ci95_high','nominal_exhaustive_two_sided_p','new_BH_q','positive_individuals']].to_string(index=False))
    print('PASS: source coefficients/P reproduced, exact algorithm verified, prior files unchanged.')

if __name__=='__main__':main()
