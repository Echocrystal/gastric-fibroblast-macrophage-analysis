from pathlib import Path
import csv, json, hashlib
from collections import defaultdict
from openpyxl import load_workbook

P=Path('I:/ECM-new/ECM_CFD_CXCL12_BCELL_GC_2026')
OUT=P/'review_signature_overlap_20261101/audit'
OUT.mkdir(exist_ok=True,parents=True)
primary=P/'WJGOsubmit/11_Code_Release_Staging/gene_sets/GSE183904_primary_sender_receiver_signatures_30plus30.tsv'
files=[primary,P/'WJGOsubmit/11_Code_Release_Staging/gene_sets/GSE183904_PLAT_ecosystem_bulk_signatures_frozen.tsv',P/'analysis/results/PLAT_axis/GSE183904/state_boundary_and_signatures/GSE183904_PLAT_ecosystem_bulk_signatures_frozen.tsv',P/'resources/gene_sets/GSE183904_PLAT_ecosystem_bulk_signatures.tsv',P/'WJGOsubmit/11_Code_Release_Staging/resources/gene_sets/GSE183904_PLAT_ecosystem_bulk_signatures.tsv']
def read_tsv(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
sets=defaultdict(list)
for row in read_tsv(primary): sets[row['signature_id']].append(row['gene_symbol'])
axis={'PLAT','ANXA2','ITGAM','ITGB2','LRP1','ILK'}
report={'primary_file':str(primary),'primary':{},'files':[]}
for sid,genes in sets.items(): report['primary'][sid]={'n':len(genes),'unique_n':len(set(genes)),'genes':genes,'focal_axis_genes_present':sorted(set(genes)&axis)}
for f in files:
    rows=read_tsv(f)
    check={'file':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'rows':len(rows),'signature_sizes':{},'primary_comparisons':{}}
    by=defaultdict(list)
    for row in rows:by[row['signature_id']].append((int(row['signature_rank']),row['gene_symbol']))
    for sid, items in by.items():
        genes=[g for _,g in sorted(items)]
        check['signature_sizes'][sid]=len(genes)
        if sid in sets:check['primary_comparisons'][sid]={'set_equal':set(genes)==set(sets[sid]),'ordered_equal':genes==sets[sid]}
    report['files'].append(check)
universe=P/'analysis/processed/GSE183904_global_atlas/GSE183904_global_atlas_genes.tsv'
M=[r['gene_symbol'] for r in read_tsv(universe)]
report['universe']={'file':str(universe),'column':'gene_symbol','rows':len(M),'unique_genes':len(set(M)),'sha256':hashlib.sha256(universe.read_bytes()).hexdigest(),'derivation_script':str(P/'scripts/45_export_GSE183904_primary_counts_for_python.R'),'derivation':'Rows of counts(sce) after secondary cell QC in the 26 primary tumour samples; identical gene order asserted for all samples; no additional expression threshold applied to the feature universe.','coverage':{sid:{'retained':len(set(genes)&set(M)),'total':len(genes),'missing':sorted(set(genes)-set(M))} for sid,genes in sets.items()}}
report['workbooks']=[]
for f in [P/'WJGOsubmit/00_Ready_to_Upload/Supplemental_Tables_1-10.xlsx',P/'MedComm-OncologySubmit/00_Ready_to_Upload/Supplementary_Tables_S1-S9.xlsx']:
    if not f.exists():continue
    wb=load_workbook(f,read_only=True,data_only=True)
    record={'file':str(f),'sheetnames':wb.sheetnames,'table6_sheets':[]}
    for ws in wb:
        if not ws.title.startswith('S6_'):continue
        nonempty=[]
        for i,row in enumerate(ws.iter_rows(values_only=True),1):
            vals=[str(v)[:180] for v in row if v is not None]
            if vals:nonempty.append({'row':i,'values':vals})
        record['table6_sheets'].append({'name':ws.title,'rows':ws.max_row,'columns':ws.max_column,'header':nonempty[0] if nonempty else {},'overlap_or_signature_rows':[{'row':r['row'],'values':r['values'][:12]} for r in nonempty if any(k in ' '.join(r['values']).lower() for k in ['overlap','signature','jaccard'])][:8]})
    gene_sheet='S6_FigureS6a_TCGA_gene_annotati'
    if gene_sheet in wb.sheetnames:
        vals=list(wb[gene_sheet].values)
        hdr=vals[0]
        anno=[dict(zip(hdr,row)) for row in vals[1:]]
        record['primary_gene_annotation_equal']={sid:set(genes)=={r['gene_symbol'] for r in anno if r['signature_id']==sid} for sid,genes in sets.items()}
    report['workbooks'].append(record)
    wb.close()
report['method_corrections']=[
'If A is restricted to M, Jaccard denominator is |A intersect M| + |B intersect M| - overlap; the simplified 30 + n - overlap denominator is valid only if all 30 primary signature genes are in M.',
'Recovery overlap/30 is an explicitly fixed original-signature denominator, not coverage-conditioned recovery. If any signature genes are absent in M, report signature coverage and optionally overlap/|A intersect M| separately.',
'A generic review framework is not a fixed signature. A non-retrieved list means not computable from retrieved evidence, not proof that no gene list exists.',
'Do not assert substantial overlap before computing it. A small Jaccard can coexist with appreciable recovery when the published set is large; both metrics should be interpreted with set sizes.',
'The actual frozen signatures exclude all six focal-axis genes but do not necessarily contain their naming markers (CXCL14, POSTN, C1QC, APOE). Comparing only naming markers would test a different set.',
'Full harmonized gene sets without expression-space restriction can be a primary set-definition comparison; measured-space restriction may be a separate sensitivity analysis if its M is precisely documented.'
]
report['selection_rule_from_code']={'script':str(P/'scripts/44_resolve_GSE183904_C1QC_boundary_and_freeze_bulk_signatures.R'),'sender_contrast':'fib_c1 vs fib_c4; CXCL14_POSTN_PLAT_vs_activated_matrix','receiver_contrast':'my_c3 vs my_c145; C1QC_APOE_vs_other_recurrent','eligible':'FDR <= 0.01; logFC >= 0.5; target patient detection >= 0.5; paired direction consistency >= 0.65; available in TCGA; no technical excluded genes; focal axis six genes excluded','rank':'ascending FDR, descending paired direction consistency, descending logFC; first 30 genes','note':'Gene availability in TCGA was used; TCGA outcomes were not used according to code and audit.'}
(OUT/'signature_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
md=['# Signature and overlap audit — 2026-11-01','',
'Primary signature source: `'+str(primary)+'`.','',
'The five frozen/list source files contain exactly the same two primary signatures in the same selection order. Each primary set has 30 unique genes; neither contains PLAT, ANXA2, ITGAM, ITGB2, LRP1 or ILK. The 120-row frozen files additionally contain a PLAT-inclusive sender and a merged C1QC receiver sensitivity signature; those must not be substituted for the primary sets.','',
'The GSE183904 feature universe contains **26,571 unique gene symbols**; both primary signatures have **30/30** represented. Source: `'+str(universe)+'`, column `gene_symbol`. The export script (45, lines 71–77) obtains these rows from the raw count assay and asserts identical gene order across the 26 primary tumour samples. This is a measured feature universe, without imposing a new expression threshold.','',
'For this explicitly defined M, Jaccard = k/(30+nB−k) and recovery = k/30 are valid. For any other M, coverage must be rechecked before using 30 in the Jaccard denominator.','',
'The WJGO and MedComm workbooks contain the same primary 60 genes in the S6 gene annotation panel. Their original S6 panels cover clinicopathology, biomarkers, TIDE, CAMERA, MCP-counter and drug sensitivity. No dedicated published CAF/macrophage signature overlap panel is present. The existing MCP-counter overlap-removal panel records **Monocytic lineage: PLA2G7, n=1**, which is not a CAF/TAM novelty analysis.','',
'## Exact selection rule found in script 44','',
'Sender: fib_c1 versus fib_c4. Receiver: my_c3 versus my_c145. Patient-blocked edgeR quasi-likelihood estimates are filtered by FDR ≤0.01, logFC ≥0.5, target patient detection ≥0.5, paired direction consistency ≥0.65, TCGA gene availability, and technical exclusions. The six focal-axis genes are removed. Eligible genes are ordered by ascending FDR, descending consistency and descending logFC, then the first 30 are frozen. TCGA gene availability is used; outcome-based selection is not used in the code.','',
'## Corrections to the supplied template','']
md += ['- '+s for s in report['method_corrections']]
(OUT/'SIGNATURE_AUDIT.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
print(json.dumps({'universe':report['universe'],'all_primary_versions_equal':all(all(x['set_equal'] and x['ordered_equal'] for x in f['primary_comparisons'].values()) for f in report['files']),'workbooks':[{'file':w['file'],'table6_sheets':[{'name':s['name'],'rows':s['rows'],'columns':s['columns']} for s in w['table6_sheets']]} for w in report['workbooks']]},ensure_ascii=False,indent=2))
