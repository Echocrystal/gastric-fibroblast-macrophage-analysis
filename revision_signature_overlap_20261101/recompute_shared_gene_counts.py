"""Recompute the seven displayed overlap rows from frozen gene lists.

Standard library only. No clinical outcomes or fitted models are used here.
The archived symbol-map snapshot fixes the harmonisation used on 2026-11-01.
"""
from pathlib import Path
import argparse
import csv
import json

BASE = Path(__file__).resolve().parent

def read_tsv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=BASE/'recomputed')
    args = parser.parse_args()
    mapping = {r['original_symbol']:r['harmonized_symbol']
               for r in read_tsv(BASE/'audit/symbol_harmonization.tsv')}
    def harmonize(gene):
        # Fail on a new input: this is deliberately a frozen-map recount.
        return mapping[gene]
    primary_rows = read_tsv(BASE/'data/frozen_primary_signatures_30plus30.tsv')
    primary = {role:{harmonize(r['gene_symbol']) for r in primary_rows
                     if r['state_role']==role.lower()}
               for role in ['Sender','Receiver']}
    universe = {r['harmonized_symbol'] for r in
                read_tsv(BASE/'data/harmonized_measured_universe.tsv')}
    assert len(universe)==26523
    assert all(len(genes)==30 and genes<=universe for genes in primary.values())
    sets = json.loads((BASE/'frozen_final_seven_sets.json').read_text(encoding='utf-8'))
    output=[]
    for ref in sets:
        genes = {harmonize(g) for g in ref['source_symbols']}
        tested=genes & universe
        shared=sorted(primary[ref['role']] & tested)
        assert len(tested)==ref['expected_reference_n'], ref['id']
        assert shared==ref['expected_overlap_genes'], ref['id']
        output.append({'Our signature':ref['role'],
                       'Published reference set':ref['published_reference_set'],
                       'Gene list used':ref['gene_list_used'],
                       'Reference n*':len(tested),'Shared genes (n)':len(shared),
                       'Overlapping genes':', '.join(shared)})
    args.output_dir.mkdir(parents=True,exist_ok=True)
    with (args.output_dir/'Table6_signature_overlap_recomputed.tsv').open(
            'w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(output[0]),delimiter='\t')
        writer.writeheader();writer.writerows(output)
    (args.output_dir/'verification.json').write_text(json.dumps({
        'rows':len(output),'reference_counts':[r['Reference n*'] for r in output],
        'shared_counts':[r['Shared genes (n)'] for r in output],
        'all_frozen_expectations_match':True},indent=2),encoding='utf-8')
    print(json.dumps({'rows':len(output),'shared_counts':[
        r['Shared genes (n)'] for r in output],'output_dir':str(args.output_dir)}))

if __name__=='__main__':
    main()
