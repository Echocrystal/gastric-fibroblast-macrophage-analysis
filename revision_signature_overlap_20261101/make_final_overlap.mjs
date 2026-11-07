import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const base=path.dirname(fileURLToPath(import.meta.url));
const data=JSON.parse(await fs.readFile(path.join(base,'final_package_content.json'),'utf8'));
const outDir=data.output_dir;
const wb=Workbook.create();
const s=wb.worksheets.add('Signature_overlap');
s.showGridLines=false;
s.getRange('A1:F8').values=[data.headers,...data.rows];
s.getRange('A1:F8').format.font={name:'Arial',size:11};
s.getRange('A1:F8').format.wrapText=true;
s.getRange('A1:F8').format.verticalAlignment='top';
s.getRange('A1:F1').format={fill:'#17324D',font:{name:'Arial',size:11,bold:true,color:'#FFFFFF'},wrapText:true,verticalAlignment:'center',horizontalAlignment:'center',rowHeight:34};
s.getRange('D2:E8').setNumberFormat('0');
s.getRange('D2:E8').format.horizontalAlignment='right';
s.getRange('F2:F8').format.font={name:'Arial',size:11,italic:true};
for (const [col,width] of [['A',16],['B',44],['C',60],['D',18],['E',18],['F',90]]){
 s.getRange(`${col}1:${col}12`).format.columnWidth=width;
}
for(let i=0;i<data.rows.length;i++){
 const k=data.rows[i][4];
 s.getRange(`A${i+2}:F${i+2}`).format.rowHeight=k===30?72:k>=13?58:46;
}
for(let i=0;i<data.notes.length;i++){
 const row=i+10;
 const range=s.getRange(`A${row}:F${row}`);
 range.merge();
 s.getRange(`A${row}`).values=[[data.notes[i]]];
 range.format={font:{name:'Arial',size:10},wrapText:true,verticalAlignment:'top',rowHeight:34};
}
s.freezePanes.freezeRows(1);
wb.recalculate();
console.log((await wb.inspect({kind:'table',range:'Signature_overlap!A1:F8',include:'values',tableMaxRows:8,tableMaxCols:6,tableMaxCellChars:100,maxChars:2300})).ndjson);
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!',options:{useRegex:true,maxResults:10},maxChars:600})).ndjson);
await fs.mkdir(outDir,{recursive:true});
const preview=await wb.render({sheetName:'Signature_overlap',range:'A1:F12',scale:1,format:'png'});
await fs.writeFile(path.join(base,'audit','final_delivery_overlap.png'),new Uint8Array(await preview.arrayBuffer()));
const output=await SpreadsheetFile.exportXlsx(wb);
await output.save(path.join(outDir,'Supplemental_Table6_overlap.xlsx'));
console.log('Saved final single-sheet Excel table.');
