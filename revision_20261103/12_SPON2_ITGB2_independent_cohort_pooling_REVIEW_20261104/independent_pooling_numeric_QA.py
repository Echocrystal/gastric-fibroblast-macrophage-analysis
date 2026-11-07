"""Independent numerical audit. Does not import the production runner.

Only this script and independent_pooling_numeric_QA.json are written.
Exact P values use a different formulation: positive subset sums and symmetry,
with long-double accumulation, rather than signed-sum enumeration.
"""
from pathlib import Path
from collections import Counter
import csv, hashlib, itertools, json, math, platform, sys
import numpy as np

sys.stdout.reconfigure(encoding='utf-8')
BASE=Path(__file__).resolve().parent
COHORTS=('GSE270678_main_ST','GSE251950_primary')
ENDPOINTS=('SPON2_ITGB2_receiver_composite_r1','SPON2_receiver_r1','SPON2_ITGB2_r1')
STR='GSE270678_ST_R_sensitivity'
TOL=1e-15
SEED=2026110412
DRAWS=99999

def read(name):
    with (BASE/name).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def zval(r):return math.atanh(max(-0.999999,min(0.999999,float(r))))
def partial_subsets(v):
    # Explicit binary masks independently enumerate subset sums in [0,sum(v)].
    bits=(np.arange(1<<len(v),dtype=np.uint64)[:,None]>>np.arange(len(v),dtype=np.uint64))&1
    return np.sum(bits.astype(np.longdouble)*np.asarray(v,dtype=np.longdouble),axis=1,dtype=np.longdouble)
def exact(z,weights):
    # Sign flips of v have the same distribution as sign flips of |v|.
    # A positive tail corresponds to a subset sum <= (sum|v|-|sum v|)/2.
    v=np.asarray(z,dtype=np.longdouble)*np.asarray(weights,dtype=np.longdouble)
    obs=abs(np.sum(v,dtype=np.longdouble));total=1<<len(v)
    if obs<=TOL:return {'tail_assignments':total,'total_assignments':total,'p':1.0}
    a=np.abs(v);cut=len(v)//2
    left=partial_subsets(a[:cut]);right=np.sort(partial_subsets(a[cut:]))
    threshold=(np.sum(a,dtype=np.longdouble)-obs+np.longdouble(TOL))/2
    positive_count=sum(int(x) for x in np.searchsorted(right,threshold-left,side='right'))
    tails=2*positive_count
    assert tails<=total
    return {'tail_assignments':tails,'total_assignments':total,'p':tails/total}
def bh_manual(p):
    # Direct definition q_(i)=min_{j>=i} m*p_(j)/j, capped at one.
    order=sorted(range(len(p)),key=p.__getitem__);out=[None]*len(p)
    for i,index in enumerate(order):
        out[index]=min(1.0,min(len(p)*p[order[j]]/(j+1) for j in range(i,len(p))))
    return out

inputs=read('inputs/SPON2_ITGB2_donor_effects.tsv')
source_meta=read('inputs/SPON2_ITGB2_meta_summary.tsv')
included=read('included_28_individuals_3endpoints.tsv')
reported=read('pooled_results_all_schemes.tsv')
original=read('original_cohort_results_3endpoints.tsv')
lopo_reported=read('leave_one_individual_out_effects.tsv')
protected_files=[BASE/'run_independent_cohort_pooling.py',BASE/'pooled_results_all_schemes.tsv',BASE/'pooled_results_all_schemes.json',BASE/'original_cohort_results_3endpoints.tsv',BASE/'leave_one_individual_out_effects.tsv']
protected_before={str(p.name):digest(p) for p in protected_files}

checks={};issues=[]
def record(name,passed,detail):
    checks[name]={'passed':bool(passed),'detail':detail}
    if not passed:issues.append(name)

selected=[r for r in inputs if r['dataset'] in COHORTS and r['endpoint'] in ENDPOINTS]
rowkeys=[(r['dataset'],r['donor_id'],r['endpoint']) for r in selected]
identities=sorted({(r['dataset'],r['donor_id']) for r in selected},key=lambda t:(COHORTS.index(t[0]),t[1]))
sizes=Counter(k[0] for k in identities)
record('pooled_source_identity',len(selected)==84 and len(set(rowkeys))==84 and len(identities)==28 and dict(sizes)==dict(zip(COHORTS,(19,9))),{'endpoint_rows':len(selected),'unique_persons':len(identities),'cohort_sizes':dict(sizes),'ST_R_included':any(k[0]==STR for k in identities)})
record('included_table_identity',len(included)==28 and {(r['dataset'],r['donor_id']) for r in included}==set(identities),{'rows':len(included),'cohorts':dict(Counter(r['dataset'] for r in included))})
lookup={(r['dataset'],r['donor_id'],r['endpoint']):float(r['rho']) for r in selected}
z=np.array([[zval(lookup[c,d,e]) for e in ENDPOINTS] for c,d in identities],dtype=float)
input_table_maxdiff=max(abs(float(r[e])-lookup[r['dataset'],r['donor_id'],e]) for r in included for e in ENDPOINTS)
record('included_numeric_values',input_table_maxdiff<1e-14,{'max_abs_difference':input_table_maxdiff})
sourcepaths={}
for filename in ('SPON2_ITGB2_donor_effects.tsv','SPON2_ITGB2_meta_summary.tsv','SPON2_ITGB2_section_effects.tsv','SPON2_ITGB2_gene_endpoint_coverage.tsv'):
    local=BASE/'inputs'/filename;old=BASE.parent/'01_spatial_endpoint/inputs'/filename
    sourcepaths[filename]={'snapshot_sha256':digest(local),'previous_snapshot_sha256':digest(old),'identical':digest(local)==digest(old)}
record('input_snapshots_unchanged',all(v['identical'] for v in sourcepaths.values()),sourcepaths)

# Small-vector external validation of the independent subset-sum algorithm.
algorithm=[]
for values in ([0,0,0],[.1,-.1,.2,-.2],[.03,.04,.05],[.01,-.02,.03,-.04,.05],list(np.random.default_rng(68411).normal(size=14))):
    for weight_kind in ('equal','unequal'):
        values=np.array(values,dtype=float)
        weights=np.ones(len(values)) if weight_kind=='equal' else np.arange(1,len(values)+1,dtype=float)
        weights/=sum(weights)
        obs=abs(math.fsum(float(x*y) for x,y in zip(values,weights)))
        brute=sum(abs(math.fsum(float(a*b*s) for a,b,s in zip(values,weights,signs)))>=obs-TOL for signs in itertools.product((-1,1),repeat=len(values)))
        result=exact(values,weights)
        algorithm.append({'n':len(values),'weights':weight_kind,'tail':result['tail_assignments'],'brute':brute,'passed':result['tail_assignments']==brute})
record('independent_MITM_small_vector_validation',all(r['passed'] for r in algorithm),algorithm)

# Independently repeat stratified bootstrap using full-size rather than chunked draws.
rng=np.random.default_rng(SEED);boot_by_cohort=[]
for c in COHORTS:
    x=z[[i for i,k in enumerate(identities) if k[0]==c]]
    indices=rng.integers(len(x),size=(DRAWS,len(x)))
    boot_by_cohort.append(np.mean(x[indices,:],axis=1))

result_rows=[];lopo_checks=[]
for scheme in ('individual_equal','cohort_equal'):
    shares=[19/28,9/28] if scheme=='individual_equal' else [.5,.5]
    weights=np.array([shares[COHORTS.index(c)]/sizes[c] for c,d in identities])
    estimates=[math.tanh(math.fsum(float(weights[i]*z[i,j]) for i in range(28))) for j in range(3)]
    exact_results=[exact(z[:,j],weights) for j in range(3)]
    q=bh_manual([r['p'] for r in exact_results])
    combined_boot=shares[0]*boot_by_cohort[0]+shares[1]*boot_by_cohort[1]
    for j,e in enumerate(ENDPOINTS):
        saved=next(r for r in reported if r['weighting']==scheme and r['source_endpoint_id']==e)
        interval=np.quantile(np.tanh(combined_boot[:,j]),[.025,.975],method='linear')
        item={'weighting':scheme,'endpoint':e,'independent_rho':estimates[j],'reported_rho':float(saved['pooled_adjusted_rho']),'rho_difference':abs(estimates[j]-float(saved['pooled_adjusted_rho'])),'independent_tail_assignments':exact_results[j]['tail_assignments'],'reported_tail_assignments':int(saved['tail_assignments']),'total_assignments':exact_results[j]['total_assignments'],'independent_P':exact_results[j]['p'],'P_difference':abs(exact_results[j]['p']-float(saved['nominal_exhaustive_two_sided_p'])),'independent_BH3_q':q[j],'q_difference':abs(q[j]-float(saved['new_BH_q'])),'independent_CI95':[float(v) for v in interval],'CI_max_abs_difference':max(abs(interval[0]-float(saved['ci95_low'])),abs(interval[1]-float(saved['ci95_high'])))}
        item['passed']=item['rho_difference']<1e-14 and item['independent_tail_assignments']==item['reported_tail_assignments'] and item['P_difference']<1e-16 and item['q_difference']<1e-16 and item['CI_max_abs_difference']<1e-14
        result_rows.append(item)

        correct_lopo=[]
        for i,key in enumerate(identities):
            remaining=[k for k in range(28) if k!=i]
            weights_lopo=np.ones(27)/27 if scheme=='individual_equal' else np.array([.5/(sizes[identities[k][0]]-int(identities[k][0]==key[0])) for k in remaining])
            value=math.tanh(math.fsum(float(weights_lopo[pos]*z[k,j]) for pos,k in enumerate(remaining)))
            correct_lopo.append(value)
            r=next(r for r in lopo_reported if r['weighting']==scheme and r['source_endpoint_id']==e and r['omitted_dataset']==key[0] and r['omitted_individual']==key[1])
            if abs(value-float(r['pooled_rho']))>1e-12:
                lopo_checks.append({'weighting':scheme,'endpoint':e,'omitted_dataset':key[0],'omitted_individual':key[1],'correct_rho':value,'reported_rho':float(r['pooled_rho']),'difference':abs(value-float(r['pooled_rho']))})
        item['correct_lopo_min']=min(correct_lopo);item['correct_lopo_max']=max(correct_lopo)
        item['reported_lopo_min']=float(saved['lopo_rho_min']);item['reported_lopo_max']=float(saved['lopo_rho_max'])
record('pooled_point_estimates_exact_tails_P_BH3_and_CIs',all(r['passed'] for r in result_rows),result_rows)

# Original cohort rows must appear exactly once, and preserve their own correction family.
original_checks=[]
for c in (*COHORTS,STR):
    for e in ENDPOINTS:
        rr=[r for r in original if r['dataset']==c and r['source_endpoint_id']==e]
        frozen=next(r for r in source_meta if r['dataset']==c and r['endpoint']==e)
        people=[r for r in inputs if r['dataset']==c and r['endpoint']==e]
        vals=[zval(r['rho']) for r in people]
        independent_rho=math.tanh(math.fsum(vals)/len(vals))
        test=exact(vals,np.repeat(1/len(vals),len(vals)))
        qsource=next(r for r in source_meta if r['dataset']==c and r['endpoint']=='candidate_SPON2_ITGB2_r1')['candidate_BH_FDR'] if e==ENDPOINTS[0] else frozen['supporting_BH_FDR']
        r=rr[0] if rr else {}
        ok=len(rr)==1 and len(people)==(19 if c==COHORTS[0] else 9 if c==COHORTS[1] else 15)
        ok=ok and abs(independent_rho-float(r['rho']))<1e-14 and abs(test['p']-float(r['nominal_two_sided_p']))<1e-14 and abs(float(r['original_BH_q'])-float(qsource))<1e-14
        ok=ok and abs(float(r['ci95_low'])-float(frozen['ci_low']))<1e-14 and abs(float(r['ci95_high'])-float(frozen['ci_high']))<1e-14
        ok=ok and r['included_in_new_pool'].lower()==str(c in COHORTS).lower()
        original_checks.append({'dataset':c,'endpoint':e,'rows':len(rr),'n_people':len(people),'independent_P':test['p'],'passed':ok})
record('original_cohort_rows_preserved_and_not_reused_as_people',len(original)==9 and all(r['passed'] for r in original_checks),original_checks)

range_checks=[{'weighting':r['weighting'],'endpoint':r['endpoint'],'correct_min':r['correct_lopo_min'],'reported_min':r['reported_lopo_min'],'correct_max':r['correct_lopo_max'],'reported_max':r['reported_lopo_max'],'passed':abs(r['correct_lopo_min']-r['reported_lopo_min'])<1e-12 and abs(r['correct_lopo_max']-r['reported_lopo_max'])<1e-12} for r in result_rows]
record('LOPO_labels_and_weights',not lopo_checks and all(r['passed'] for r in range_checks),{'mismatched_rows':len(lopo_checks),'first_mismatches':lopo_checks[:10],'summary_range_checks':range_checks,'interpretation':'Production z is stacked main19 then external9 while wide.index sorts external9 then main19. LOPO key/position becomes misaligned; individual-equal omitted-person labels are wrong and cohort-equal reweighting is wrong. Main pooled point estimates/CI/P/q unaffected.' if lopo_checks else 'All omitted-person values and six summary ranges match independent correctly ordered rows.'})
record('production_outputs_unchanged_by_this_audit',all(digest(p)==protected_before[p.name] for p in protected_files),protected_before)
out={'status':'PASS' if not issues else ('PRIMARY_NUMERICS_PASS_WITH_LOPO_DEFECT' if issues==['LOPO_labels_and_weights'] else 'FAIL'),'independent_implementation':'No import of runner; csv + math.fsum for estimates; longdouble positive subset-sum MITM + sign symmetry; direct-definition BH3; full-array stratified bootstrap.','script_sha256':digest(__file__),'python':platform.python_version(),'numpy':np.__version__,'checks':checks,'issues':issues}
(BASE/'independent_pooling_numeric_QA.json').write_text(json.dumps(out,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)),encoding='utf-8')
print(json.dumps({'status':out['status'],'issues':issues,'six_pooled_rows':[{k:r[k] for k in ('weighting','endpoint','independent_rho','independent_tail_assignments','independent_P','independent_BH3_q','CI_max_abs_difference','passed','correct_lopo_min','correct_lopo_max')} for r in result_rows],'LOPO_mismatches':len(lopo_checks)},ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
