from pathlib import Path
import csv
from collections import defaultdict

BASE=Path(__file__).parent
HGNC=BASE/'sources_symbols/hgnc_complete_set.txt'
with HGNC.open(encoding='utf-8',newline='') as f:
    records=[r for r in csv.DictReader(f,delimiter='\t') if r['status']=='Approved']
APPROVED={r['symbol'].upper():r['symbol'].upper() for r in records}
PREV=defaultdict(set)
ALIAS=defaultdict(set)
for r in records:
    for value in r['prev_symbol'].split('|'):
        if value:PREV[value.upper()].add(r['symbol'].upper())
    for value in r['alias_symbol'].split('|'):
        if value:ALIAS[value.upper()].add(r['symbol'].upper())

def harmonize(g):
    value=str(g).strip().upper()
    if value in APPROVED:return value,'approved'
    previous=PREV.get(value,set())
    if len(previous)==1:return next(iter(previous)),'previous_symbol'
    if len(previous)>1:return value,'ambiguous_previous_retained'
    aliases=ALIAS.get(value,set())
    if len(aliases)==1:return next(iter(aliases)),'alias_symbol'
    if len(aliases)>1:return value,'ambiguous_alias_retained'
    return value,'unmapped_retained'

if __name__=='__main__':
    import json,sys
    from collections import Counter
    sys.stdout.reconfigure(encoding='utf-8')
    root=BASE.parent
    p=root/'analysis/processed/GSE183904_global_atlas/GSE183904_global_atlas_genes.tsv'
    with p.open(encoding='utf-8-sig',newline='') as f:genes=[r['gene_symbol'] for r in csv.DictReader(f,delimiter='\t')]
    mapped=[(g,*harmonize(g)) for g in genes]
    with (BASE/'audit/symbol_harmonization.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['original_symbol','harmonized_symbol','mapping_status']);w.writerows(mapped)
    print(json.dumps({'raw_unique':len(set(genes)),'harmonized_unique':len({r[1] for r in mapped}),'status_counts':dict(Counter(r[2] for r in mapped))},indent=2))
