"""Post hoc R2-5 sensitivity; writes exclusively below this script's directory.

Run with project runtime py_fig1, which supplies h5py/scipy/numpy/pandas.
Frozen small source snapshots plus extracted raw union-set counts support reruns.
No original project output is modified. No clinical outcomes are used.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import shutil
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.stats import wilcoxon

HERE = Path(__file__).resolve().parent
INP = HERE / "inputs"
RESULT = HERE / "results"
for p in [INP, RESULT]: p.mkdir(exist_ok=True)
AXIS = {"PLAT", "ANXA2", "ITGAM", "ITGB2", "LRP1", "ILK"}
GO = "GOBP_PHAGOCYTOSIS"
REACTOME = "REACTOME_FCGAMMA_RECEPTOR_FCGR_DEPENDENT_PHAGOCYTOSIS"
SEED = 2026110305

def sha256(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()

def bh(p):
    p=np.asarray(p,float); order=np.argsort(p); q=np.empty(len(p))
    q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return q

def snapshot(project, snapshots_only=False):
    paths = {
      "frozen_programs.tsv":"resources/gene_sets/myeloid_receiver_programs_frozen.tsv",
      "myeloid_metadata.tsv.gz":"analysis/processed/GSE183904_fine_states/GSE183904_Myeloid_fine_state_cell_metadata.tsv.gz",
      "atlas_metadata.tsv.gz":"analysis/processed/GSE183904_global_atlas/GSE183904_global_atlas_cells.tsv.gz",
      "atlas_genes.tsv":"analysis/processed/GSE183904_global_atlas/GSE183904_global_atlas_genes.tsv",
      "receiver_pseudobulk_DE.tsv":"analysis/results/PLAT_axis/GSE183904/state_boundary_and_signatures/GSE183904_C1QC_APOE_vs_other_recurrent_paired_pseudobulk_DE.tsv",
      "original_gsea.tsv":"figures/source_data/Figure4/Figure4f_preranked_GSEA.tsv",
      "original_program_contrasts.tsv":"analysis/results/PLAT_axis/GSE183904/receiver_state_characterization/GSE183904_C1QC_APOE_program_contrasts.tsv",
      "original_patient_differences.tsv":"analysis/results/PLAT_axis/GSE183904/receiver_state_characterization/GSE183904_C1QC_APOE_vs_other_recurrent_paired_differences.tsv",
      "original_scoring.R":"scripts/37_characterize_GSE183904_C1QC_APOE_receiver_state.R",
      "original_freeze_programs.R":"scripts/08_freeze_myeloid_receiver_programs.R",
      "original_Figure4.py":"figures/gen_fig4_patient_association.py",
    }
    if snapshots_only:
        manifest_file=INP/"source_provenance_sha256.tsv"
        if not manifest_file.is_file():
            raise FileNotFoundError("--snapshots-only requires inputs/source_provenance_sha256.tsv")
        saved=pd.read_csv(manifest_file,sep="\t")
        for name in paths:
            # Basenames stay valid after copying this folder to another machine.
            row=saved[saved.snapshot_path.map(lambda p:Path(str(p)).name).eq(name)]
            if len(row)!=1:
                raise ValueError(f"Expected one original provenance entry for snapshot {name}")
            local=INP/name
            if not local.is_file() or sha256(local)!=row.sha256.iloc[0]:
                raise ValueError(f"Saved snapshot missing or SHA256 mismatch: {name}")
        cache_manifest=INP/"count_cache_sha256.tsv"
        if not cache_manifest.is_file():
            raise FileNotFoundError("--snapshots-only requires inputs/count_cache_sha256.tsv")
        cache=pd.read_csv(cache_manifest,sep="\t")
        expected={"union_phagocytosis_raw_counts.npz","scoring_cells.tsv","scoring_genes.tsv"}
        if set(cache.filename)!=expected:
            raise ValueError("Count-cache manifest must name exactly the NPZ, cell metadata and gene metadata")
        for row in cache.itertuples():
            local=INP/row.filename
            if not local.is_file() or sha256(local)!=row.sha256:
                raise ValueError(f"Count cache missing or SHA256 mismatch: {row.filename}")
        print(f"Snapshot-only inputs verified: {len(paths)} source snapshots and {len(cache)} count-cache files; original project files are not accessed",flush=True)
        return saved.to_dict("records")
    rows=[]
    for target, rel in paths.items():
        src=project/rel; dest=INP/target
        if not dest.exists(): shutil.copy2(src,dest)
        assert sha256(dest)==sha256(src), target
        rows.append({"source_path":str(src),"snapshot_path":str(dest),"bytes":src.stat().st_size,"sha256":sha256(src)})
    return rows

def load_counts(project, members, provenance):
    cache=INP/"union_phagocytosis_raw_counts.npz"
    cellfile=INP/"scoring_cells.tsv"
    genefile=INP/"scoring_genes.tsv"
    if cache.exists():
        return sparse.load_npz(cache),pd.read_csv(cellfile,sep="\t"),pd.read_csv(genefile,sep="\t").gene_symbol.tolist()
    genes=pd.read_csv(INP/"atlas_genes.tsv",sep="\t").gene_symbol.astype(str).tolist()
    selected=sorted(set(genes)&members)
    gi=np.array([genes.index(g) for g in selected])
    my=pd.read_csv(INP/"myeloid_metadata.tsv.gz",sep="\t")
    atlas=pd.read_csv(INP/"atlas_metadata.tsv.gz",sep="\t")
    parts=[];cells=[]
    for sample, cm in my.groupby("sample_id",sort=True):
        raw=project/"analysis/processed/GSE183904_global_atlas/raw_h5"/f"{sample}_counts_csc.h5"
        with h5py.File(raw) as f:
            ids=[x.decode() if isinstance(x,bytes) else str(x) for x in f["cell_id"][:]]
            mat=sparse.csc_matrix((f["data"][:],f["indices"][:],f["indptr"][:]),shape=tuple(f["shape_genes_by_cells"][:]))
        index={s:i for i,s in enumerate(ids)}
        cols=np.array([index[s] for s in cm.cell_id])
        lib=np.asarray(mat[:,cols].sum(axis=0)).ravel()
        assert np.array_equal(lib,cm.sum_counts.to_numpy()),sample
        cc=cm[["cell_id","sample_id","Patient_ID","cluster"]].copy()
        cc["library_umi"]=lib
        parts.append(mat[gi,:][:,cols].tocsr()); cells.append(cc)
        provenance.append({"source_path":str(raw),"snapshot_path":"selected count entries in union_phagocytosis_raw_counts.npz; full library_umi in scoring_cells.tsv","bytes":raw.stat().st_size,"sha256":sha256(raw)})
        print("Extracted",sample,len(cc),flush=True)
    matrix=sparse.hstack(parts,format="csr")
    cell=pd.concat(cells,ignore_index=True)
    assert len(cell)==9002 and cell.Patient_ID.nunique()==26
    sparse.save_npz(cache,matrix)
    cell.to_csv(cellfile,sep="\t",index=False)
    pd.DataFrame({"gene_symbol":selected}).to_csv(genefile,sep="\t",index=False)
    return matrix,cell,selected

def enrichment_score(scores, positions):
    positions=np.sort(np.asarray(positions,int)); n=len(scores); k=len(positions)
    weights=np.abs(scores[positions]); cumul=np.cumsum(weights)/weights.sum()
    miss=(positions-np.arange(k))/(n-k)
    after=cumul-miss; before=np.r_[0,cumul[:-1]]-miss
    ip=int(np.argmax(after)); im=int(np.argmin(before))
    pe=max(0.,after[ip]); ne=min(0.,before[im])
    return (float(pe),positions[:ip+1]) if pe>=abs(ne) else (float(ne),positions[im:])

def gsea(ranked, sets, seed):
    scores=ranked.ranking_statistic.to_numpy(float); genes=ranked.gene_symbol.to_numpy(str)
    pos={s:i for i,s in enumerate(genes)}; rng=np.random.default_rng(seed); rows=[]
    for name,members in sorted(sets.items()):
        idx=np.array([pos[g] for g in sorted(set(members)&set(pos))])
        es,leading=enrichment_score(scores,idx)
        null=np.array([enrichment_score(scores,rng.choice(len(scores),size=len(idx),replace=False))[0] for _ in range(10000)])
        ss=null[null>=0] if es>=0 else null[null<0]
        p=(1+np.sum(ss>=es if es>=0 else ss<=es))/(1+len(ss))
        rows.append({"program":name,"ranked_background_n":len(scores),"gene_n":len(idx),"ES":es,"NES":es/abs(ss.mean()),"nominal_p":p,"leading_edge_genes":"|".join(genes[leading]),"permutation_seed":seed,"permutations":10000})
    d=pd.DataFrame(rows);d["BH_q_14_sets"]=bh(d.nominal_p)
    return d

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--project-root",type=Path,default=Path(r"I:\ECM-new\ECM_CFD_CXCL12_BCELL_GC_2026"))
    ap.add_argument("--snapshots-only",action="store_true",help="Verify saved source snapshots and cached counts; never read the original project")
    args=ap.parse_args()
    provenance=snapshot(args.project_root,args.snapshots_only)
    frozen=pd.read_csv(INP/"frozen_programs.tsv",sep="\t")
    sets={p:set(d.gene_symbol) for p,d in frozen[frozen.evidence_tier.eq("primary_external")].groupby("program")}
    de=pd.read_csv(INP/"receiver_pseudobulk_DE.tsv",sep="\t")
    if de.technical_excluded.dtype!=bool: de["technical_excluded"]=de.technical_excluded.astype(str).str.lower().eq("true")
    ranked=de[~de.technical_excluded].copy()
    ranked["ranking_statistic"]=np.sign(ranked.logFC_target_minus_comparator)*np.sqrt(ranked.QLF.clip(lower=0))
    ranked=ranked.sort_values(["ranking_statistic","gene_symbol"],ascending=[False,True]).reset_index(drop=True)
    assert len(ranked)==8187
    universe=set(ranked.gene_symbol)-AXIS
    matrix,cells,genes=load_counts(args.project_root,sets[GO]|sets[REACTOME],provenance)
    prior_manifest=INP/"source_provenance_sha256.tsv"
    if prior_manifest.exists():
        prior=pd.read_csv(prior_manifest,sep="\t").to_dict("records")
        known={r["source_path"] for r in provenance}
        provenance.extend(r for r in prior if r["source_path"] not in known)
    if not args.snapshots_only:
        pd.DataFrame(provenance).to_csv(INP/"source_provenance_sha256.tsv",sep="\t",index=False)
        pd.DataFrame([{"filename":name,"bytes":(INP/name).stat().st_size,"sha256":sha256(INP/name)} for name in ["union_phagocytosis_raw_counts.npz","scoring_cells.tsv","scoring_genes.tsv"]]).to_csv(INP/"count_cache_sha256.tsv",sep="\t",index=False)
    gindex={g:i for i,g in enumerate(genes)}
    variants={
      "GO_original":(GO,sets[GO]-AXIS),
      "GO_matched":(GO,(sets[GO]&universe)-AXIS),
      "Reactome_full":(REACTOME,sets[REACTOME]-AXIS),
      "Reactome_matched":(REACTOME,(sets[REACTOME]&universe)-AXIS),
    }
    member_rows=[]
    for program in [GO,REACTOME]:
        for gene in sorted(sets[program]):
            member_rows.append({"program":program,"gene_symbol":gene,"present_in_matrix":gene in gindex,"focal_axis_excluded":gene in AXIS,"in_original_GSEA_background":gene in set(ranked.gene_symbol),"in_matched_background":gene in universe,"used_in_full_score":gene in gindex and gene not in AXIS,"used_in_matched_score":gene in gindex and gene in universe})
    pd.DataFrame(member_rows).to_csv(RESULT/"gene_membership_audit.tsv",sep="\t",index=False)
    patient_rows=[]; gene_rows=[]; coverage=[]
    for patient, cm in cells.groupby("Patient_ID",sort=True):
        idx=cm.index.to_numpy(); norm=matrix[:,idx].toarray()*1e4/cm.library_umi.to_numpy()[None,:]
        x=np.log1p(norm); sd=x.std(axis=1,ddof=1); informative=np.isfinite(sd)&(sd>0)
        z=np.zeros_like(x);z[informative]=(x[informative]-x[informative].mean(axis=1,keepdims=True))/sd[informative,None]
        target=cm.cluster.eq(3).to_numpy(); comp=cm.cluster.isin([1,4,5]).to_numpy()
        eligible=target.sum()>=5 and comp.sum()>=20
        assert eligible
        delta=z[:,target].mean(axis=1)-z[:,comp].mean(axis=1)
        for v,(program,members) in variants.items():
            gi=np.array([gindex[g] for g in sorted(members&set(gindex)) if informative[gindex[g]]])
            score=z[gi].mean(axis=0)
            patient_rows.append({"variant":v,"program":program,"Patient_ID":patient,"target_n_cells":int(target.sum()),"comparator_n_cells":int(comp.sum()),"informative_gene_n":len(gi),"target_mean":float(score[target].mean()),"comparator_mean":float(score[comp].mean()),"delta_target_minus_comparator":float(delta[gi].mean())})
            coverage.append({"variant":v,"Patient_ID":patient,"informative_gene_n":len(gi)})
        for i,gene in enumerate(genes):
            gene_rows.append({"Patient_ID":patient,"gene_symbol":gene,"informative":bool(informative[i]),"delta_z":float(delta[i]),"target_mean_log1p_CP10K":float(x[i,target].mean()),"comparator_mean_log1p_CP10K":float(x[i,comp].mean()),"delta_log1p_CP10K":float(x[i,target].mean()-x[i,comp].mean())})
    paired=pd.DataFrame(patient_rows);paired.to_csv(RESULT/"patient_paired_scores.tsv",sep="\t",index=False)
    gene_effect=pd.DataFrame(gene_rows);gene_effect.to_csv(RESULT/"gene_patient_effects.tsv.gz",sep="\t",index=False)
    rows=[]; rng=np.random.default_rng(SEED)
    for variant,(program,members) in variants.items():
        d=paired[paired.variant.eq(variant)];vals=d.delta_target_minus_comparator.to_numpy()
        w=wilcoxon(vals,alternative="two-sided",zero_method="wilcox",correction=True,method="approx")
        boot=np.median(rng.choice(vals,size=(9999,len(vals)),replace=True),axis=1)
        rows.append({"variant":variant,"program":program,"n_patients":len(vals),"frozen_set_n":len(sets[program]),"gene_n":len(members&set(gindex)),"informative_gene_n_min":int(d.informative_gene_n.min()),"informative_gene_n_max":int(d.informative_gene_n.max()),"median_delta":np.median(vals),"mean_delta":np.mean(vals),"median_delta_bootstrap_CI_low":np.quantile(boot,.025),"median_delta_bootstrap_CI_high":np.quantile(boot,.975),"positive_patients":int(np.sum(vals>0)),"negative_patients":int(np.sum(vals<0)),"paired_Wilcoxon_W":w.statistic,"paired_Wilcoxon_P":w.pvalue})
    summary=pd.DataFrame(rows);summary["BH_q_four_posthoc_score_tests"]=bh(summary.paired_Wilcoxon_P)
    summary.to_csv(RESULT/"scoring_summary.tsv",sep="\t",index=False)
    orig=pd.read_csv(INP/"original_patient_differences.tsv",sep="\t");orig=orig[orig.program.eq(GO)]
    check=paired[paired.variant.eq("GO_original")].merge(orig,on="Patient_ID",suffixes=("_new","_old"))
    max_error=float(np.max(abs(check.delta_target_minus_comparator-check.delta_target_minus_other_recurrent)))
    assert max_error<1e-10, max_error
    audit={"status":"post_hoc_standalone_reviewer_sensitivity","n_myeloid_cells":len(cells),"n_patients":cells.Patient_ID.nunique(),"target_cluster":3,"comparator_clusters":[1,4,5],"within_patient_scaling_cells":"all myeloid cells including cluster 2","eligibility":"at least 5 target and 20 pooled comparator cells per patient","normalization":"log1p(CP10K), full-library denominator; gene-wise sample-SD Z scores within patient across all myeloid cells; equal mean across informative module genes; patient-state means","score_test":"two-sided asymptotic paired Wilcoxon signed-rank, continuity corrected","score_multiplicity":"BH across all 4 explicitly named post hoc score tests; does not replace original six-programme family","bootstrap":"9999 patient bootstrap resamples, percentile CI for median delta; seed2026110305","raw_DE_tested_n":len(de),"original_GSEA_ranked_n":len(ranked),"technical_removed_n":int(de.technical_excluded.sum()),"matched_GSEA_ranked_n":len(universe),"axis_removed_from_original_GSEA_background":sorted(set(ranked.gene_symbol)&AXIS),"original_GO_score_reproduction_max_abs_error":max_error,"versions":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"pandas":pd.__version__,"h5py":h5py.__version__}}
    # Reproduce the frozen 14-set analysis once before the matched-membership analysis.
    print("Reproducing original 14-set GSEA",flush=True)
    old=gsea(ranked,sets,20260801);old["analysis"]="original_reproduced"
    oldfile=pd.read_csv(INP/"original_gsea.tsv",sep="\t")
    match=old.merge(oldfile,on="program",suffixes=("_new","_old"))
    assert np.max(abs(match.NES_new-match.NES_old))<1e-10
    assert np.max(abs(match.nominal_p_new-match.nominal_p_old))<1e-10
    print("Running matched 14-set GSEA",flush=True)
    new=gsea(ranked[ranked.gene_symbol.isin(universe)].reset_index(drop=True),{k:v-AXIS for k,v in sets.items()},SEED);new["analysis"]="matched_axis_excluded_posthoc"
    gsea_all=pd.concat([old,new],ignore_index=True)
    gsea_all.to_csv(RESULT/"gsea_all14_original_and_matched.tsv",sep="\t",index=False)
    gsea_all[gsea_all.program.isin([GO,REACTOME])].to_csv(RESULT/"gsea_summary.tsv",sep="\t",index=False)
    gs=gene_effect.groupby("gene_symbol").agg(mean_patient_delta_z=("delta_z","mean"),median_patient_delta_z=("delta_z","median"),positive_patients=("delta_z",lambda v:int(np.sum(v>0))),mean_delta_log1p_CP10K=("delta_log1p_CP10K","mean")).reset_index()
    gs=gs.merge(de,on="gene_symbol",how="left")
    gs["ranked_GSEA_background"]=gs.gene_symbol.isin(set(ranked.gene_symbol))
    gs["focal_axis_excluded"]=gs.gene_symbol.isin(AXIS)
    for program,label in [(GO,"GO"),(REACTOME,"Reactome")]:
        gs[f"{label}_member"]=gs.gene_symbol.isin(sets[program])
        le=set(oldfile.loc[oldfile.program.eq(program),"leading_edge_genes"].iloc[0].split("|"))
        gs[f"{label}_original_GSEA_leading_edge"]=gs.gene_symbol.isin(le)
    gs.to_csv(RESULT/"gene_effects_and_GSEA_weights.tsv",sep="\t",index=False)
    group_rows=[]
    for program,label in [(GO,"GO"),(REACTOME,"Reactome")]:
        d=gs[gs[f"{label}_member"]&gs.ranked_GSEA_background&~gs.focal_axis_excluded].copy()
        audit[f"{label}_matched_gene_direction_summary"]={"gene_n":len(d),"positive_mean_delta_z":int(np.sum(d.mean_patient_delta_z>0)),"positive_pseudobulk_logFC":int(np.sum(d.logFC_target_minus_comparator>0)),"negative_pseudobulk_logFC":int(np.sum(d.logFC_target_minus_comparator<0))}
        for side,q in d.groupby(d.logFC_target_minus_comparator.gt(0).map({True:"target_up",False:"comparator_up"})):
            group_rows.append({"program":program,"pseudobulk_direction":side,"gene_n":len(q),"mean_patient_delta_z":q.mean_patient_delta_z.mean(),"sum_abs_GSEA_ranking_weight":float(np.sqrt(q.QLF).sum()),"median_log2FC":q.logFC_target_minus_comparator.median(),"genes":"|".join(q.gene_symbol)})
    pd.DataFrame(group_rows).to_csv(RESULT/"gene_direction_weight_summary.tsv",sep="\t",index=False)
    audit["GSEA_reproduction_max_NES_error"]=float(np.max(abs(match.NES_new-match.NES_old)))
    audit["new_GSEA_multiplicity"]="BH across fixed14 external sets on matched axis-excluded background; standalone post hoc sensitivity"
    (RESULT/"analysis_audit.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding="utf-8")
    print(summary.to_string(index=False),flush=True)
    print(gsea_all[gsea_all.program.isin([GO,REACTOME])].to_string(index=False),flush=True)

if __name__=="__main__":main()
