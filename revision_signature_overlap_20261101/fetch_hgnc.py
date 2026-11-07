from pathlib import Path
import urllib.request, hashlib, json, sys
sys.stdout.reconfigure(encoding='utf-8')
out=Path(__file__).parent/'sources_symbols'
out.mkdir(exist_ok=True)
url='https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'
p=out/'hgnc_complete_set.txt'
try:
    if not p.exists():
        with urllib.request.urlopen(url,timeout=45) as r:p.write_bytes(r.read())
    text=p.read_text(encoding='utf-8')
    assert text.startswith('hgnc_id\tsymbol\t'),text[:100]
    (out/'source.json').write_text(json.dumps({'url':url,'accessed_date':'2026-11-01','sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'file':str(p),'rows':len(text.splitlines())-1},indent=2),encoding='utf-8')
    print('HGNC saved',p.stat().st_size,'bytes',len(text.splitlines())-1,'records')
except Exception as e:print(type(e).__name__,str(e))
