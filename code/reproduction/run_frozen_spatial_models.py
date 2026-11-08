"""Run preserved spatial source functions using explicit portable paths.

The AST adapter changes module-level path assignments only. Statistical function
bodies and fixed parameters remain those in the byte-identical frozen sources.
It does not derive a new model or alter the original correction families.
"""
from __future__ import annotations
import argparse, ast, json, types, hashlib
from pathlib import Path
import numpy as np
import pandas as pd

def file_sha256(path):
    h=hashlib.sha256()
    with path.open('rb')as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def load_source(path, overrides, injected=None, skip_base_import=False):
    tree=ast.parse(path.read_text(encoding='utf-8-sig'),filename=str(path))
    body=[]
    for node in tree.body:
        names=[t.id for t in node.targets if isinstance(t,ast.Name)] if isinstance(node,ast.Assign) else []
        if skip_base_import and any(n in {'SPEC','M'} for n in names):continue
        if skip_base_import and isinstance(node,ast.Assert) and 'SPEC.loader' in ast.unparse(node):continue
        if skip_base_import and isinstance(node,ast.Expr) and 'SPEC.loader.exec_module' in ast.unparse(node):continue
        if names and any(n in overrides for n in names):
            assert len(names)==1
            node.value=ast.Subscript(value=ast.Name(id='_PATH_OVERRIDES',ctx=ast.Load()),slice=ast.Constant(names[0]),ctx=ast.Load())
        body.append(node)
    tree.body=body;ast.fix_missing_locations(tree)
    module=types.ModuleType('frozen_'+path.stem)
    module.__dict__.update(__file__=str(path),__name__='frozen_'+path.stem,_PATH_OVERRIDES=overrides)
    module.__dict__.update(injected or {})
    exec(compile(tree,str(path),'exec'),module.__dict__)
    return module

def configure(repo,raw,out,definitions=None):
    defs=definitions or repo
    sources=repo/'code/frozen_sources';cov=repo/'resources/spatial_covariates'
    base_path=sources/'35_analyze_GSE270678_PLAT_spatial_chain.py'
    candidate_path=sources/'54_analyze_SPON2_ITGB2_spatial_candidate_context.py'
    base=load_source(base_path,{'PROJECT':repo,'PUBLIC':raw,'OUT':out/'primary_chain_helpers'})
    candidate=load_source(candidate_path,{
        'PROJECT':repo,'PUBLIC':raw,'BASE_SCRIPT':base_path,'OUT':out,
        'LOCK_FILE':cov/'analysis_locks/SPON2_ITGB2_SPATIAL_CANDIDATE_CONTEXT_LOCK_20260826.md',
        'SIGNATURE_FILE':defs/'resources/gene_sets/GSE183904_PLAT_ecosystem_bulk_signatures.tsv',
        'PROGRAM_FILE':defs/'resources/gene_sets/myeloid_receiver_programs_frozen.tsv',
        'MANIFEST_270678':cov/'GSE270678_spatial_sample_manifest.tsv',
        'MANIFEST_251950':cov/'GSE251950_spatial_sample_manifest.tsv',
    },{'M':base},True)
    candidate.DATASETS['GSE270678_main_ST'].update(raw=raw/'GSE270678/main_public_files',spot=cov/'GSE270678_main_ST_frozen_covariates.tsv.gz')
    candidate.DATASETS['GSE270678_ST_R_sensitivity'].update(raw=raw/'GSE270678/ST_R_public_files',spot=cov/'GSE270678_ST_R_frozen_covariates.tsv.gz')
    candidate.DATASETS['GSE251950_primary'].update(raw=raw/'GSE251950/samples',spot=cov/'GSE251950_frozen_covariates.tsv.gz')
    return candidate

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[2])
    p.add_argument('--raw-root',type=Path,required=True,help='Root containing GSE270678/ and GSE251950/, arranged per spatial_section_h5_map.tsv')
    p.add_argument('--definition-root',type=Path,help='Optional explicit root for existing signature/program resources; default repo root')
    p.add_argument('--output-dir',type=Path,required=True)
    action=p.add_mutually_exclusive_group(required=True)
    action.add_argument('--check-inputs',action='store_true')
    action.add_argument('--verify-first-sections',action='store_true',help='Recompute one existing section per dataset; not a new cohort analysis')
    action.add_argument('--run-all',action='store_true',help='Run preserved source54 main() for all original datasets, including original inference')
    p.add_argument('--expected-section-effects',type=Path,help='Frozen existing section-effect TSV for bounded verification')
    args=p.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    mod=configure(args.repo_root.resolve(),args.raw_root.resolve(),args.output_dir.resolve(),args.definition_root.resolve() if args.definition_root else None)
    raw_index=pd.read_csv(args.repo_root/'resources/spatial_covariates/GEO_raw_H5_file_index.tsv',sep='\t')
    expected_raw={str(Path(r.raw_h5_relative)):r for r in raw_index.itertuples(index=False)}
    for f in [mod.SIGNATURE_FILE,mod.PROGRAM_FILE,mod.LOCK_FILE,mod.MANIFEST_270678,mod.MANIFEST_251950]:
        if not f.is_file():raise FileNotFoundError(f)
    input_report=[]
    for name,cfg in mod.DATASETS.items():
        manifest=mod.dataset_manifest(name,cfg)
        assert len(manifest)==cfg['expected_sections'] and manifest.donor_id.nunique()==cfg['expected_donors']
        spots=pd.read_csv(cfg['spot'],sep='\t')
        assert not spots.duplicated(['sample_id' if cfg['kind']=='270678' else 'GSM','barcode']).any()
        for row in manifest.itertuples(index=False):
            h5=cfg['raw']/row.filtered_h5_file if cfg['kind']=='270678' else cfg['raw']/row.folder_id/(row.folder_id+'_filtered_feature_bc_matrix.h5')
            if not h5.is_file():raise FileNotFoundError(h5)
            rel=str(h5.relative_to(args.raw_root.resolve()))
            expected=expected_raw[rel]
            assert h5.stat().st_size==expected.raw_h5_bytes and file_sha256(h5)==expected.raw_h5_sha256,('H5 integrity mismatch',h5)
        input_report.append({'dataset':name,'manifest_sections':len(manifest),'individuals':int(manifest.donor_id.nunique()),'frozen_spots':len(spots),'all_raw_h5_present':True,'all_raw_h5_sizes_and_sha256_match':True})
    if args.check_inputs:
        (args.output_dir/'input_validation.json').write_text(json.dumps(input_report,indent=2),encoding='utf-8')
        print(json.dumps(input_report,indent=2));return
    if args.run_all:
        mod.main();return
    sender,receiver,programmes=mod.load_fixed_definitions()
    needed=set(sender+receiver+mod.M.MYELOID_GENES+sorted(mod.FOCAL_GENES))
    for genes in programmes.values():needed.update(genes)
    expected=pd.read_csv(args.expected_section_effects,sep='\t') if args.expected_section_effects else None
    verification=[];effects_all=[]
    for name,cfg in mod.DATASETS.items():
        row=next(mod.dataset_manifest(name,cfg).itertuples(index=False))
        spots=pd.read_csv(cfg['spot'],sep='\t')
        base,h5=mod.section_inputs(cfg,row,spots)
        norm=mod.M.read_selected_expression(h5,base.barcode.tolist(),needed,base.total_counts.to_numpy(float))
        effects,coverage=mod.analyze_section(base,norm,sender,receiver,programmes)
        effect=pd.DataFrame(effects);effect.insert(0,'section_id',row.section_id);effect.insert(0,'dataset',name)
        effects_all.append(effect)
        report={'dataset':name,'section_id':row.section_id,'spots':len(base),'effects':len(effects),'section_gate':bool(coverage['section_gate'])}
        if expected is not None:
            frozen=expected[(expected.dataset==name)&(expected.section_id==row.section_id)]
            merged=effect.merge(frozen,on=['dataset','section_id','endpoint','ring'],suffixes=('_recomputed','_frozen'),validate='one_to_one')
            assert len(merged)==len(effect)==len(frozen)
            a=merged.rho_recomputed.to_numpy(float);b=merged.rho_frozen.to_numpy(float)
            assert np.array_equal(np.isfinite(a),np.isfinite(b))
            err=float(np.max(np.abs(a[np.isfinite(a)]-b[np.isfinite(b)])))
            assert err<1e-10,(name,err)
            assert (merged.n_spots_recomputed==merged.n_spots_frozen).all()
            report.update(max_rho_absolute_difference=err,all_effects_and_denominators_match=True)
        verification.append(report);print(report,flush=True)
    pd.concat(effects_all,ignore_index=True).to_csv(args.output_dir/'bounded_recomputed_section_effects.tsv',sep='\t',index=False)
    (args.output_dir/'bounded_verification.json').write_text(json.dumps(verification,indent=2),encoding='utf-8')

if __name__=='__main__':main()
