import fs from 'node:fs/promises';
import { FileBlob, SpreadsheetFile } from '@oai/artifact-tool';
import { fileURLToPath } from 'node:url';
import path from 'node:path';
import { spawnSync } from 'node:child_process';

const base = path.dirname(fileURLToPath(import.meta.url));
const target = path.join(base, 'Signature_overlap_SIMPLE.xlsx');
const wb = await SpreadsheetFile.importXlsx(await FileBlob.load(target));
if (process.argv[2] !== 'edit') {
  console.log((await wb.inspect({kind:'sheet', include:'id,name', maxChars:1500})).ndjson);
  console.log((await wb.inspect({kind:'table', range:'Overlap_table!A1:G3', include:'values,formulas', maxChars:2500, tableMaxRows:3, tableMaxCols:7})).ndjson);
  console.log(wb.help('range.delete', {include:'index,examples,notes', maxChars:3000}).ndjson);
  const preview = await wb.render({sheetName:'Overlap_table',range:'A1:G4',scale:1,format:'png'});
  await fs.writeFile(path.join(base,'audit','overlap_before_shared_genes.png'),new Uint8Array(await preview.arrayBuffer()));
  process.exit(0);
}

const source = JSON.parse(await fs.readFile(path.join(base,'overlap_results.json'),'utf8'));
const byId = new Map(source.rows.map(r=>[r.id,r]));
const ids = [
 'Li2022_eCAF_cluster3_deposited','Lee2025_CCL2_fibroblast425','Zhang2025_C1QC_RTM20',
 'Ma2025_Macro_C1QC_author_threshold','Ma2025_Macro_SPP1_author_threshold',
 'Lee2025_STAT3_macrophage425','Becht2016_monocytic'
];
const s = wb.worksheets.getItem('Overlap_table');
const original = s.getRange('A2:D8').values;
const headers = ['Our signature','Published reference set','How the reference genes were obtained','Reference genes (n)','Shared genes (n)','Overlapping genes'];
const rows = ids.map((id,i)=>{
 const r=byId.get(id);
 if (new Set(r.overlap_genes).size!==r.overlap_n) throw new Error('Gene count mismatch: '+id);
 return [...original[i],r.overlap_n,r.overlap_genes.join(', ')];
});
s.getRange('A9:G11').unmerge();
s.getRange('G1:G11').clear({applyTo:'all'});
s.getRange('A1:F8').values = [headers,...rows];
s.getRange('A1:F8').format.wrapText = true;
s.getRange('A1:F8').format.verticalAlignment = 'top';
s.getRange('D2:E8').setNumberFormat('0');
for (const [col,width] of [['A',16],['B',42],['C',65],['D',18],['E',18],['F',92]]) {
 s.getRange(`${col}1:${col}11`).format.columnWidth = width;
}
s.getRange('A1:F1').format.rowHeight=34;
for (let i=0;i<rows.length;i++) {
 s.getRange(`A${i+2}:F${i+2}`).format.rowHeight = rows[i][4]===30 ? 92 : rows[i][4]>=13 ? 76 : 58;
}
const notes = [
 'Reference gene counts use harmonised symbols restricted to the GSE183904 measured gene space. All 30 genes in each primary signature are covered. Shared genes (n) is the number of genes common to our primary signature and the indicated reference set; all shared symbols are listed in the last column.',
 'Ma et al. share GSE183904 with our derivation analysis. Each Ma set excludes two unresolved gene identifiers corrupted in the source workbook; the reported numbers describe the recoverable-symbol subsets. These comparisons do not constitute independent-cohort validation.',
 'MCP-counter monocytic lineage is a generic lineage panel. Li et al. provide deposited positive cluster markers. Additional ECM and generic fibroblast panel comparisons are retained in the separate source record.'
];
for (let i=0;i<notes.length;i++) {
 const range=s.getRange(`A${i+9}:F${i+9}`);
 range.merge();
 s.getRange(`A${i+9}`).values=[[notes[i]]];
 range.format.wrapText=true;
 range.format.rowHeight=36;
}
const rationale='CXCL14+POSTN+ PLAT-rich fibroblast cluster 1 vs activated-matrix fibroblast cluster 4. Marker review flagged clusters 2, 3 and 5 for other-lineage contamination; cluster 6 was represented by only one patient. Cluster 4 supported within-patient analysis in 16 patients with at least five cells per state. Cluster 4 had uneven donor representation; enrichment is relative to this specified comparator.';
const selection=wb.worksheets.getItem('Our_marker_selection');
selection.getRange('B3').values=[[rationale]];
selection.getRange('A3:B3').format.rowHeight=84;
wb.recalculate();
console.log((await wb.inspect({kind:'table',range:'Overlap_table!A1:F8',include:'values',tableMaxRows:8,tableMaxCols:6,tableMaxCellChars:800,maxChars:7000})).ndjson);
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!',options:{useRegex:true,maxResults:10},maxChars:1000})).ndjson);
const preview=await wb.render({sheetName:'Overlap_table',range:'A1:F8',scale:1,format:'png'});
await fs.writeFile(path.join(base,'audit','overlap_shared_genes_final.png'),new Uint8Array(await preview.arrayBuffer()));
const out=await SpreadsheetFile.exportXlsx(wb);
await out.save(target);
// The import/export API does not expose sheet AutoFilters. Restore these native
// view controls directly from the pre-edit snapshot after the authored export.
const repair = spawnSync('C:/Users/ASUS/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe', ['-c', `
import json,sys,zipfile,os
from pathlib import Path
from xml.etree import ElementTree as ET
b=Path(sys.argv[1]); p=b/'Signature_overlap_SIMPLE.xlsx'
before=json.loads((b/'audit/shared_genes_edit_before.json').read_text(encoding='utf-8'))
ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
ET.register_namespace('',ns)
ET.register_namespace('r','http://schemas.openxmlformats.org/officeDocument/2006/relationships')
temp=p.with_suffix('.filters.tmp')
with zipfile.ZipFile(p) as original,zipfile.ZipFile(temp,'w',zipfile.ZIP_DEFLATED) as revised:
    book=ET.fromstring(original.read('xl/workbook.xml'))
    rels=ET.fromstring(original.read('xl/_rels/workbook.xml.rels'))
    targets={x.attrib['Id']:x.attrib['Target'] for x in rels}
    filters={}
    for sheet in book.find('{'+ns+'}sheets'):
        name=sheet.attrib['name']; key=sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
        raw=targets[key]; file=raw.lstrip('/') if raw.startswith('/') else 'xl/'+raw
        ref=before.get(name,{}).get('auto_filter')
        if name=='Overlap_table':ref='A1:F8'
        if ref:filters[file]=ref
    for item in original.infolist():
        data=original.read(item.filename)
        if item.filename in filters:
            root=ET.fromstring(data)
            old=root.find('{'+ns+'}autoFilter')
            if old is not None:root.remove(old)
            el=ET.Element('{'+ns+'}autoFilter',{'ref':filters[item.filename]})
            before_tags={'sortState','dataConsolidate','customSheetViews','mergeCells','phoneticPr','conditionalFormatting','dataValidations','hyperlinks','printOptions','pageMargins','pageSetup','headerFooter','rowBreaks','colBreaks','customProperties','cellWatches','ignoredErrors','smartTags','drawing','legacyDrawing','legacyDrawingHF','picture','oleObjects','controls','webPublishItems','tableParts','extLst'}
            at=next((i for i,c in enumerate(root) if c.tag.rsplit('}',1)[-1] in before_tags),len(root))
            root.insert(at,el)
            data=ET.tostring(root,encoding='utf-8',xml_declaration=True)
        revised.writestr(item,data)
os.replace(temp,p)
`, base],{encoding:'utf8'});
if (repair.status!==0) throw new Error(repair.stderr || 'Native filter preservation failed');
console.log('Updated existing workbook: 7 comparisons, integer shared-gene counts and complete gene names.');
