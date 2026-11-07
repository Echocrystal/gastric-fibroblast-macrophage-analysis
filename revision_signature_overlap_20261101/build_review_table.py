"""Build a standalone review panel. Never edit the submission workbook."""
from pathlib import Path
import csv, json, hashlib, html, sys, shutil, textwrap
from collections import Counter
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, KeepTogether
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import A3, landscape
from symbols import harmonize

sys.stdout.reconfigure(encoding='utf-8')
BASE=Path(__file__).parent
ROOT=BASE.parent
PREFIX='Supplemental_Table_6_signature_overlap_REVIEW'
PRIMARY=ROOT/'WJGOsubmit/11_Code_Release_Staging/gene_sets/GSE183904_primary_sender_receiver_signatures_30plus30.tsv'
UNIVERSE=ROOT/'analysis/processed/GSE183904_global_atlas/GSE183904_global_atlas_genes.tsv'
def read_tsv(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
primary_rows=read_tsv(PRIMARY)
primary={role:{harmonize(r['gene_symbol'])[0] for r in primary_rows if r['state_role']==role.lower()} for role in ['Sender','Receiver']}
assert all(len(g)==30 for g in primary.values())
raw_m={r['gene_symbol'] for r in read_tsv(UNIVERSE)}
M={harmonize(g)[0] for g in raw_m}
assert all(g<=M for g in primary.values())
sets=json.loads((BASE/'comparison_sets.json').read_text(encoding='utf-8'))
all_rows=[];gene_rows=[];full_rows=[];checks=[]
for s in sets:
    A=primary[s['role']]
    result=dict(s)
    if s.get('genes') is None:
        result.update(original_n=None,tested_n=None,overlap_n=None,jaccard=None,recovery=None,overlap_genes=[],removed_genes=[])
    else:
        raw=list(dict.fromkeys(str(g).strip() for g in s['genes'] if str(g).strip()))
        mapping=[(g,*harmonize(g)) for g in raw]
        B={r[1] for r in mapping};Btest=B&M;Atest=A&M;intersection=Atest&Btest
        n=len(intersection);union=len(Atest|Btest)
        assert union>0 and n<=min(len(Atest),len(Btest))
        assert len(Atest)==30
        result.update(original_n=len(B),tested_n=len(Btest),overlap_n=n,jaccard=n/union,recovery=n/30,
                      overlap_genes=sorted(intersection),removed_genes=sorted(B-M),raw_unique_n=len(raw),union_n=union,
                      signature_measured_n=len(Atest))
        fulln=len(A&B)
        full_rows.append([s['id'],s['role'],s['label'],len(A),len(B),fulln,fulln/len(A|B),fulln/30,'; '.join(sorted(A&B))])
        for orig,mapped,status in mapping:
            gene_rows.append([s['id'],s['role'],orig,mapped,status,mapped in M,mapped in A,mapped in intersection])
        checks.append({'set_id':s['id'],'integer_bound':n<=min(30,len(Btest)),
                       'jaccard_identity':abs(result['jaccard']-n/(30+len(Btest)-n))<1e-14,
                       'recovery_identity':abs(result['recovery']-n/30)<1e-14,
                       'full_vs_measured_overlap_equal':n==fulln})
    all_rows.append(result)

source_note=(f'M is the GSE183904 count-assay feature universe: {len(raw_m):,} unique exported gene symbols '
             f'and {len(M):,} unique symbols after harmonisation. All 30 Sender and all 30 Receiver genes are represented in M. '
             'No new expression-detection threshold was imposed. The same M was used for every comparison; it is not an intersection across the external bulk or spatial platforms.')
footnote=('For each role, A is the frozen primary 30-gene signature, B is the explicitly identified reference set, A\u2032 = A intersect M, '
          "and B' = B intersect M. Reference set sizes are deduplicated after symbol harmonisation; both the full valid-symbol set size and the size retained in M are shown. "
          'Overlap k = |A\u2032 intersect B\u2032|; Jaccard = k / |A\u2032 union B\u2032|; recovery = k/30. '
          'Both primary signatures have 30/30 coverage, so Jaccard = k/(30 + |B\u2032| - k) in this table. '
          'Zero denotes a verified empty intersection; NC denotes a reference without a specified or retrieved fixed set, and is not zero. '
          'Metrics are descriptive; no enrichment P value, state-exclusivity test or signature reselection is performed. '
          'Marker panels, source-derived marker sets, ECM programmes and generic lineage panels have distinct scope and are labelled accordingly. '
          'A zero overlap with a short marker panel does not establish a new cell state. '
          'Ma et al. reuse GSE183904, the derivation cohort of our signatures, so those comparisons do not constitute independent-cohort validation. '
          'The Ma source workbook contains two unresolved Excel date-formatted identifiers per macrophage set; these entries were excluded without guessing their gene identities. '
          'Ma metrics describe only the recoverable-symbol subsets (626 and 805 valid symbols from 628 and 807 threshold-passing source entries).')
symbol_note=('Symbols were trimmed and uppercased, then matched to the downloaded HGNC complete-set snapshot (2026-11-01). '
             'Approved symbols took priority; previous symbols and aliases were replaced only when the corresponding mapping was unique. '
             'Ambiguous and unrecognised symbols were retained for exact matching; no speculative mapping was applied. '
             'The gene-level mapping log and a full-gene-space sensitivity comparison are included.')
footnote=footnote.replace('\u2032', "'")
methods_draft=('Signature overlap with published CAF and macrophage programmes was assessed after freezing the primary sender and receiver gene lists. '
               'We extracted explicitly defined reference signatures, author-reported marker sets and fixed marker panels, retaining their distinct labels and original scope. '
               +source_note+' '+symbol_note+' '+footnote)
results={'date':'2026-11-01','primary_file':str(PRIMARY),'universe_file':str(UNIVERSE),
         'universe_raw_n':len(raw_m),'universe_harmonized_n':len(M),'rows':all_rows,'checks':checks,'footnote':footnote,'source_note':source_note,'symbol_note':symbol_note}
(BASE/'overlap_results.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')

HEADERS=['Signature','Published reference set','Source (manuscript ref.)','Set class','Full reference set n (harmonised)','Reference set n (in M)',
         'Overlap genes n','Jaccard index','Recovery (k/30)','Overlap gene symbols','Status / scope']
def display_row(r):
    calc=r['overlap_n'] is not None
    return [r['role']+' (30 genes)',r['label'],r['source'],r['set_type'],
            r['original_n'] if calc else 'NC',r['tested_n'] if calc else 'NC',r['overlap_n'] if calc else 'NC',
            r['jaccard'] if calc else 'NC',r['recovery'] if calc else 'NC',
            '; '.join(r['overlap_genes']) if r['overlap_genes'] else ('None' if calc else 'Not applicable'),
            r.get('status','Computed')+'; '+r.get('limitation','')]
def text_row(r):
    row=display_row(r)
    if r['jaccard'] is not None:
        row[7]=f"{r['jaccard']:.3f}"
        row[8]=f"{r['overlap_n']}/30 = {r['recovery']:.3f} ({r['recovery']:.1%})"
    return row

w=openpyxl.Workbook();ws=w.active;ws.title='S6_overlap_review'
ws.append(['Supplemental Table 6 — proposed new panel (review copy)'])
ws.append(['Signature overlap with published CAF and macrophage gene sets'])
ws.append(['Standalone review file; existing submission documents and workbooks have not been merged or edited.'])
for row in range(1,4):ws.merge_cells(start_row=row,start_column=1,end_row=row,end_column=len(HEADERS))
ws.append(HEADERS)
for r in all_rows:ws.append(display_row(r))
ws.freeze_panes='E5';ws.auto_filter.ref=f'A4:K{ws.max_row}'
for i,r in enumerate(all_rows,5):
    if r['jaccard'] is not None:
        ws.cell(i,8).number_format='0.000'
        ws.cell(i,9).number_format=f'"{r["overlap_n"]}/30 = "0.000'
    import math
    rowtexts=[str(v) for v in display_row(r)]
    charwidths=[19,47,29,31,16,16,14,14,24,46,76]
    ws.row_dimensions[i].height=max(82,16*max(math.ceil(len(t)/(cw*.9)) for t,cw in zip(rowtexts,charwidths))+12)
data_last=ws.max_row
for text in [source_note,footnote,symbol_note]:
    ws.append([text]);ws.merge_cells(start_row=ws.max_row,start_column=1,end_row=ws.max_row,end_column=11);ws.row_dimensions[ws.max_row].height=65
ws.sheet_properties.pageSetUpPr.fitToPage=True
ws.page_setup.orientation='landscape';ws.page_setup.paperSize=ws.PAPERSIZE_A3;ws.page_setup.fitToWidth=1;ws.page_setup.fitToHeight=0
ws.print_title_rows='1:4';ws.print_options.horizontalCentered=True
widths=[19,47,29,31,16,16,14,14,24,46,76]
for i,width in enumerate(widths,1):ws.column_dimensions[get_column_letter(i)].width=width

def add_sheet(name,headers,rows,widths=None):
    sh=w.create_sheet(name);sh.append(headers)
    for row in rows:sh.append(row)
    sh.freeze_panes='A2';sh.auto_filter.ref=sh.dimensions
    for i in range(1,len(headers)+1):sh.column_dimensions[get_column_letter(i)].width=(widths[i-1] if widths else min(70,max(18,len(headers[i-1])+3)))
    return sh
add_sheet('Full_gene_space',[ 'Set ID','Signature','Reference set','Signature n','Reference n','Overlap n','Jaccard','Recovery (k/30)','Overlap genes'],full_rows,[30,18,62,15,15,15,15,20,70])
add_sheet('Reference_gene_lists',['Set ID','Signature','Source symbol','Harmonised symbol','Mapping status','In measured M','In primary signature','In tested overlap'],gene_rows,[34,16,24,24,30,18,20,20])
add_sheet('Primary_60_genes',['Signature','Selection rank','Frozen symbol','Harmonised symbol','In M','Signature ID'],
          [[r['state_role'],int(r['signature_rank']),r['gene_symbol'],harmonize(r['gene_symbol'])[0],harmonize(r['gene_symbol'])[0] in M,r['signature_id']] for r in primary_rows],
          [18,20,24,25,14,65])
source_rows=[]
for s in sets:
    source_rows.append([s['id'],s['role'],s['source'],s['label'],s.get('doi',''),s['set_type'],s.get('source_location',''),s.get('extraction_rule',''),s.get('source_url',''),s.get('source_file',''),s.get('limitation',''),s.get('selected_raw_entries',''),json.dumps(s.get('unresolved_identifiers',[]),ensure_ascii=False)])
add_sheet('Sources_and_rules',['Set ID','Role','Citation','Reference set','DOI','Set class','Source location','Extraction rule','Source URL','Local evidence','Limitations','Raw entries passing threshold','Unresolved source identifiers'],source_rows,[34,15,36,60,36,36,65,95,80,100,100,24,80])
notes=[['Universe',source_note],['Table footnote',footnote],['Symbols',symbol_note],['Draft supplementary method',methods_draft],
       ['审阅说明','这是独立新表。原Supplemental Table 6、正文和补充说明均未合并修改。数值来自冻结签名和可追溯文献列表。NC不等于零。'],
       ['结果解释','不得提前写 overlap substantially。不同集合大小、marker panel与完整状态signature、以及Ma2025与本研究共用GSE183904，都影响解释。零重叠不能证明新细胞状态。'],
       ['Selection audit','Patient-blocked edgeR QL; sender fib_c1 vs fib_c4; receiver my_c3 vs my_c145; FDR <=0.01, logFC >=0.5, target patient detection >=0.5, paired direction consistency >=0.65, TCGA gene availability and technical exclusions. Exclude PLAT, ANXA2, ITGAM, ITGB2, LRP1, ILK; sort ascending FDR, descending consistency, descending logFC, retain first 30. See audit/SIGNATURE_AUDIT.md.']]
add_sheet('Methods_and_notes',['Item','Text'],notes,[30,150])
orig=json.loads((BASE/'original_input_hashes.json').read_text(encoding='utf-8'))
integrity=[]
for path,v in orig.items():
    unchanged=digest(Path(path))==v['sha256'];assert unchanged
    integrity.append([path,v['sha256'],digest(Path(path)),unchanged])
add_sheet('Input_integrity',['Original input','SHA256 before','SHA256 after','Unchanged'],integrity,[110,70,70,20])
add_sheet('Arithmetic_QC',['Set ID','Integer bound','Jaccard identity','Recovery identity','Full vs measured overlap equal'],
          [[c[k] for k in ['set_id','integer_bound','jaccard_identity','recovery_identity','full_vs_measured_overlap_equal']] for c in checks],[36,20,22,22,35])

navy='17324D';teal='EAF6F3';blue='EDF3FA';amber='FFF2D7';thin=Side(style='hair',color='D7E0E8')
for sh in w:
    header=4 if sh==ws else 1
    sh.sheet_view.showGridLines=False
    for row in sh:
        for c in row:
            c.font=Font(name='Calibri',size=11,color='203448')
            c.alignment=Alignment(vertical='top',wrap_text=True)
            c.border=Border(bottom=thin)
    for c in sh[header]:c.fill=PatternFill('solid',fgColor=navy);c.font=Font(name='Calibri',size=11,bold=True,color='FFFFFF')
    sh.row_dimensions[header].height=35
    if sh!=ws:
        for ri in range(2,sh.max_row+1):
            if ri%2==0:
                for c in sh[ri]:c.fill=PatternFill('solid',fgColor='F3F6FA')
        if sh.title in ['Sources_and_rules','Methods_and_notes']:
            for ri in range(2,sh.max_row+1):sh.row_dimensions[ri].height=112
for i,r in enumerate(all_rows,5):
    for c in ws[i]:c.fill=PatternFill('solid',fgColor=amber if r['overlap_n'] is None else blue if r['role']=='Sender' else teal)
for ri in [1,2]:
    ws.cell(ri,1).font=Font(name='Calibri',bold=True,size=17 if ri==1 else 14,color=navy)
    ws.row_dimensions[ri].height=28
ws.row_dimensions[3].height=25
for sh in ['Full_gene_space']:
    for row in w[sh].iter_rows(min_row=2):row[6].number_format='0.000';row[7].number_format='0.000'
xlsx=BASE/(PREFIX+'.xlsx');w.save(xlsx)

def esc(v):return html.escape(str(v))
body=[]
for r in all_rows:
    cls='nc' if r['overlap_n'] is None else r['role'].lower()
    body.append('<tr class="'+cls+'">'+''.join('<td>'+esc(v)+'</td>' for v in text_row(r))+'</tr>')
source_html=''.join('<li><b>'+esc(s['id'])+'</b>: '+esc(s.get('source_location',''))+' '+esc(s.get('extraction_rule',''))+' '+
                    ('<a href="https://doi.org/'+esc(s['doi'])+'">DOI</a>' if s.get('doi') else '')+'</li>' for s in sets)
page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>Supplemental Table 6 — review panel</title><style>
body{font:15px/1.65 "Segoe UI","Microsoft YaHei",Arial,sans-serif;color:#183148;background:#f4f7fa;margin:0;padding:28px}
main{max-width:1800px;margin:auto;background:white;border-radius:12px;padding:30px;box-shadow:0 2px 12px #17324d12}
h1{font-size:25px;margin:0 0 5px}h2{font-size:18px;margin:25px 0 10px}.tag{color:#846000;background:#fff2d7;padding:8px 12px;border-radius:5px}
.scroll{overflow:auto}table{border-collapse:collapse;width:100%;min-width:1450px;margin:15px 0}th{background:#17324d;color:white;text-align:left;position:sticky;top:0}
td,th{border:1px solid #d8e1e8;padding:10px;vertical-align:top}td{font-size:13px}.sender{background:#edf3fa}.receiver{background:#eaf6f3}.nc{background:#fff5df}
td:nth-child(n+5):nth-child(-n+9){white-space:nowrap;text-align:center}th:nth-child(2){min-width:190px}th:last-child{min-width:260px}
.note{font-size:13px;color:#40586d}.callout{border-left:4px solid #da9f26;background:#fffaed;padding:12px 16px}li{margin:8px 0}a{color:#126d9a}
@media print{body{padding:0;background:white}main{padding:0;box-shadow:none}.scroll{overflow:visible}table{min-width:0}td,th{font-size:8px;padding:4px}}
</style><main><h1>Supplemental Table 6 · proposed new panel</h1><p>Signature overlap with published CAF and macrophage gene sets</p>
<p class="tag">待审阅独立新表 · 2026-11-01 · 尚未并入正式补充表</p>
<div class="callout">数值按真实基因列表计算。NC 表示没有已定义或已取得的固定集合，不能当作 0。短 marker panel 的零重叠不能证明状态新颖；Ma 2025 与本研究共用 GSE183904，不能作为独立队列验证。</div>
<div class="scroll"><table><thead><tr>'''+''.join('<th>'+esc(h)+'</th>' for h in HEADERS)+'</tr></thead><tbody>'+''.join(body)+'''</tbody></table></div>
<h2>Table footnote</h2><p class="note">'''+esc(source_note)+'</p><p class="note">'+esc(footnote)+'</p><p class="note">'+esc(symbol_note)+'''</p>
<h2>可核对的来源与提取规则</h2><ol>'''+source_html+'''</ol><h2>与审稿回复的关系</h2>
<p>“overlap substantially”需要由具体集合的实测结果支持，不能作为预设结论。完整基因列表、逐基因映射、全基因空间敏感性结果、计算核查及正式文件未改动的校验记录，均在同名 Excel 文件中。</p></main></html>'''
(BASE/(PREFIX+'.html')).write_text(page,encoding='utf-8')

styles=getSampleStyleSheet();styles['Normal'].fontSize=9;styles['Normal'].leading=12
styles['Title'].fontSize=18;styles['Title'].leading=22
pdf=SimpleDocTemplate(str(BASE/(PREFIX+'.pdf')),pagesize=landscape(A3),rightMargin=28,leftMargin=28,topMargin=25,bottomMargin=25)
story=[Paragraph('Supplemental Table 6. Proposed new panel',styles['Title']),
       Paragraph('Signature overlap with published CAF and macrophage gene sets',styles['Heading2']),
       Paragraph('REVIEW COPY | 1 November 2026 | Not merged into the submission workbook',styles['Normal']),Spacer(1,10)]
pdf_headers=['Signature','Published reference set','Source','Set class','Full valid n','n in M','Overlap n','Jaccard','Recovery k/30','Overlap genes','Status / scope']
pstyle=styles['Normal']
pdfrows=[[Paragraph(esc(v),pstyle) for v in pdf_headers]]
pdfrows += [[Paragraph(esc(v),pstyle) for v in text_row(r)] for r in all_rows]
widths_pdf=[62,145,102,118,46,46,48,52,85,120,285]
assert sum(widths_pdf) <= landscape(A3)[0]-56
t=Table(pdfrows,colWidths=widths_pdf,repeatRows=1,hAlign='LEFT')
commands=[('BACKGROUND',(0,0),(-1,0),colors.HexColor('#DCE6F0')),('VALIGN',(0,0),(-1,-1),'TOP'),('GRID',(0,0),(-1,-1),.3,colors.HexColor('#CCD6DF')),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6)]
for i,r in enumerate(all_rows,1):commands.append(('BACKGROUND',(0,i),(-1,i),colors.HexColor('#FFF2D7' if r['overlap_n'] is None else '#EDF3FA' if r['role']=='Sender' else '#EAF6F3')))
t.setStyle(TableStyle(commands));story.append(t);story.append(Spacer(1,10))
for para in [source_note,footnote,symbol_note]:story.extend([Paragraph(esc(para),styles['Normal']),Spacer(1,7)])
story.append(PageBreak());story.append(Paragraph('Reference sets and extraction rules',styles['Heading1']))
for s in sets:
    text='<b>'+esc(s['id'])+' - '+esc(s['source'])+'</b><br/>'+esc(s['label'])+'<br/>'+esc(s.get('source_location',''))+'<br/>'+esc(s.get('extraction_rule',''))+'<br/>Scope: '+esc(s.get('limitation',''))
    if s.get('doi'):text+='<br/>DOI: '+esc(s['doi'])
    story.extend([Paragraph(text,styles['Normal']),Spacer(1,12)])
pdf.build(story)

data=BASE/'data';data.mkdir(exist_ok=True)
shutil.copyfile(PRIMARY,data/'frozen_primary_signatures_30plus30.tsv')
with (data/'harmonized_measured_universe.tsv').open('w',encoding='utf-8',newline='') as f:
    wr=csv.writer(f,delimiter='\t');wr.writerow(['harmonized_symbol']);wr.writerows([[g] for g in sorted(M)])
with (BASE/'review_panel.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    wr=csv.writer(f,delimiter='\t');wr.writerow(HEADERS);wr.writerows([text_row(r) for r in all_rows])
used={str(p):digest(p) for p in [PRIMARY,UNIVERSE,BASE/'comparison_sets.json',BASE/'sources_symbols/hgnc_complete_set.txt']}
for s in sets:
    p=Path(s.get('source_file',''))
    if p.is_file():used[str(p)]=digest(p)
(BASE/'source_hashes.json').write_text(json.dumps(used,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'xlsx':str(xlsx),'rows':len(all_rows),'computed':sum(r['overlap_n'] is not None for r in all_rows),'universe':len(M),'arithmetic_checks':len(checks),'original_inputs_unchanged':len(integrity),'summary':[[r['id'],r['tested_n'],r['overlap_n'],r['jaccard'],r['recovery']] for r in all_rows]},ensure_ascii=False,indent=2))
