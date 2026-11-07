from pathlib import Path
import json
BASE=Path(__file__).parent
CAF=BASE/'sources_caf';MAC=BASE/'sources_macrophage'
refs={'Li2022':'[7]','Yang2025':'[8]','Lee2025':'[9]','Sahai2020':'[6]','Zhang2025':'[37]','Ma2025':'[38]','Becht2016':'[27]'}
def norm(s,folder):
    s=dict(s)
    s['role']=s.get('role',s.get('signature'))
    s['label']=s.get('label',s.get('published_reference_set'))
    s['genes']=s.get('genes',s.get('gene_symbols'))
    if s.get('source_file'):
        p=Path(s['source_file'])
        if not p.is_absolute():s['source_file']=str(folder/p)
    s['source_url']=s.get('source_url',s.get('retrieval_url',s.get('url','')))
    key=next((k for k in refs if s['id'].startswith(k)),None)
    if 'source' not in s:
        citation={'Becht2016':'Becht et al., Genome Biol 2016;17:218','Zhang2025':'Zhang et al., Cell Discov 2025;11:60','Ma2025':'Ma et al., Front Immunol 2025;16:1618591'}
        s['source']=citation[key]
    if key:s['source']+=' '+refs[key]
    s['limitation']=s.get('limitation',s.get('reference_note',s.get('reason','')))
    s['set_type']=s.get('set_type','Fixed generic lineage marker panel' if s['id'].startswith('Becht') else 'Author-reported characteristic marker panel')
    if s['id'].startswith('Sahai'):
        s['set_type']='Review framework; no single fixed signature';s['genes']=None
        s['source_location']='Full review; Table 1 concerns CAF-related clinical trials'
        s['extraction_rule']='No gene list was constructed by pooling markers mentioned in the review.'
    s['status']='NC' if s['genes'] is None else 'Computed'
    if s['id'].startswith('Becht'):s['source_location']='Official MCP-counter Signatures/genes.txt (commit pinned in extraction rule)'
    if s['id'].startswith('Zhang'):s['source_location']=s.get('source_location','Methods: Definitions of RTMs')
    if s['id']=='Zhang2025_C1QC_RTM20':s['source_location']='Supplementary Table S4, C1QC+ RTM row, PDF p35 (Page 35 of 41)'
    if s['id'].startswith('Ma'):s['source_location']='Supplementary Table S2; positive cluster markers; source Excel cell-type damage audited'
    return s
def caf(filename):return norm(json.loads((CAF/filename).read_text(encoding='utf-8')),CAF)
mac=[norm(s,MAC) for s in json.loads((MAC/'verified_macrophage_sets.json').read_text(encoding='utf-8'))]
byid={s['id']:s for s in mac}
sets=[caf('Li2022_eCAF_cluster3_deposited.json'),caf('Lee2025_CCL2_fibroblast425.json'),caf('Yang2025_ECM_hub10.json'),byid['Becht2016_fibroblasts'],caf('Sahai2020_not_computable.json')]
zhang_full=[s for s in mac if s['id']=='Zhang2025_C1QC_RTM20']
sets += zhang_full or [byid['Zhang2025_RTM_markers']]
sets += [s for s in mac if s['id'].startswith('Ma')]
sets += [caf('Lee2025_STAT3_macrophage425.json'),byid['Becht2016_monocytic']]
sets += [{'id':'Generic_TAM_unspecified','role':'Receiver','label':'Generic TAM signature (not specified in the proposed source list)','source':'No fixed reference supplied','set_type':'Undefined reference set','genes':None,'status':'NC','limitation':'No unique generic TAM signature is defined. A fixed publication and gene list must be specified before quantitative comparison. No zero or imputed value is assigned.','source_location':'User-supplied proposed panel','extraction_rule':'No invented gene list or pooled collection of macrophage markers.'}]
(BASE/'comparison_sets.json').write_text(json.dumps(sets,ensure_ascii=False,indent=2),encoding='utf-8')
print([(s['id'],None if s['genes'] is None else len(s['genes'])) for s in sets])
