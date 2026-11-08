"""Reconstruct original ECM and broad-lineage covariates for first sections.

Uses raw GEO H5 counts and the frozen QC spot/barcode selections. This is a
bounded verification of existing numeric values, not a new association model.
"""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd,h5py
from scipy.sparse import csc_matrix

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[2])
    p.add_argument('--raw-root',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
    cov=a.repo_root/'resources/spatial_covariates'
    fixed=pd.read_csv(a.repo_root/'resources/ecm_direction/frozen_ECM_highlow50.tsv',sep='\t')
    panels=pd.read_csv(cov/'spatial_model_broad_lineage_markers.tsv',sep='\t')
    sections=pd.read_csv(cov/'spatial_section_h5_map.tsv',sep='\t')
    cases=[('GSE270678_main_ST','sample_id','FROZEN_ECM_C2_HIGH50','FROZEN_ECM_C1_HIGH50'),('GSE270678_ST_R','sample_id','FROZEN_ECM_C2_HIGH50','FROZEN_ECM_C1_HIGH50'),('GSE251950','GSM','FROZEN_C2_HIGH50','FROZEN_C1_HIGH50')]
    out=[]
    for name,key,high,low in cases:
        data=pd.read_csv(cov/(name+'_frozen_covariates.tsv.gz'),sep='\t')
        row=sections[(sections.dataset==name)&sections.included_script54_manifest.astype(bool)].iloc[0]
        d=data[data[key]==row.section_id].copy()
        gs={high:fixed[fixed.ecm_group=='C2_HIGH50'].gene_symbol.tolist(),low:fixed[fixed.ecm_group=='C1_HIGH50'].gene_symbol.tolist()}
        for panel in ['MARKER_Fibroblast','MARKER_Epithelial']:gs[panel]=panels[panels.panel==panel].gene_symbol.tolist()
        with h5py.File(a.raw_root/row.raw_h5_relative,'r')as f:
            m=f['matrix'];symbols=np.asarray(m['features/name']).astype('U');bars=np.asarray(m['barcodes']).astype('U')
            counts=csc_matrix((m['data'][...],m['indices'][...],m['indptr'][...]),shape=tuple(m['shape'][...]))
        lookup={b:i for i,b in enumerate(bars)};j=[lookup[b]for b in d.barcode]
        total=np.asarray(counts[:,j].sum(axis=0)).ravel()
        original_total=d['total_counts' if 'total_counts'in d else'total_UMI'].to_numpy(float)
        assert np.array_equal(total,original_total)
        vectors={};coverage={}
        for panel,genes in gs.items():
            values=[];measured=0;informative=0
            for g in genes:
                index=np.flatnonzero(symbols==g)
                if not len(index):
                    if name=='GSE251950':continue
                    x=np.zeros(len(j))
                else:
                    measured+=1;x=np.asarray(counts[index,:][:,j].sum(axis=0)).ravel()
                x=np.log1p(1e4*x/np.maximum(total,1));sd=x.std(ddof=1)
                if np.isfinite(sd)and sd>0:z=(x-x.mean())/sd;informative+=1
                else:z=np.zeros_like(x)
                values.append(z)
            assert len(values)>=2
            vectors[panel]=np.mean(values,axis=0)
            coverage[panel]={'requested':len(genes),'measured':measured,'informative':informative,'averaged':len(values)}
        vectors['frozen_ECM_direction']=vectors[high]-vectors[low]
        err={k:float(np.max(np.abs(v-d[k].to_numpy(float))))for k,v in vectors.items()}
        assert max(err.values())<1e-10,(name,err)
        r={'dataset':name,'section_id':row.section_id,'spots':len(d),'max_absolute_error_by_covariate':err,'coverage':coverage,'all_totals_exact':True};out.append(r);print(r,flush=True)
    (a.output_dir/'covariate_reconstruction_verification.json').write_text(json.dumps(out,indent=2),encoding='utf-8')

if __name__=='__main__':main()
