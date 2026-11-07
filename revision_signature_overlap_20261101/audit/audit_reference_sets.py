from pathlib import Path
import sys,csv,json,re
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from symbols import harmonize
with (BASE.parent/'WJGOsubmit/11_Code_Release_Staging/gene_sets/GSE183904_primary_sender_receiver_signatures_30plus30.tsv').open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
A={r:{harmonize(x['gene_symbol'])[0] for x in rows if x['state_role']==r.lower()} for r in ['Sender','Receiver']}
with (BASE.parent/'analysis/processed/GSE183904_global_atlas/GSE183904_global_atlas_genes.tsv').open(encoding='utf-8-sig',newline='') as f:M={harmonize(x['gene_symbol'])[0] for x in csv.DictReader(f,delimiter='\t')}
sets=[]
for p in (BASE/'sources_caf').glob('*.json'):
    if 'metadata' in p.name:continue
    s=json.loads(p.read_text(encoding='utf-8'))
    if isinstance(s,dict) and 'gene_symbols' in s:sets.append({**s,'genes':s['gene_symbols'],'role':s['signature']})
sets+=json.loads((BASE/'sources_macrophage/verified_macrophage_sets.json').read_text(encoding='utf-8'))
summary=[]
for s in sets:
    raw=set(s['genes']);B={harmonize(g)[0] for g in raw};tested=B&M;k=len(A[s['role']]&tested)
    summary.append({'id':s['id'],'raw_n':len(raw),'harmonized_n':len(B),'in_M':len(tested),'k':k,'jaccard':k/(30+len(tested)-k),'recovery':k/30,'overlap':sorted(A[s['role']]&tested),'numeric_identifiers':[g for g in raw if re.fullmatch(r'[0-9.]+',str(g))],'harmonization_changes':[(g,harmonize(g)[0]) for g in sorted(raw) if g.upper()!=harmonize(g)[0]]})
(BASE/'audit/reference_set_independent_audit.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
for r in summary:print(json.dumps({k:v for k,v in r.items() if k!='harmonization_changes'},ensure_ascii=False))
