#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
python - <<'PY'
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import json
root=Path.cwd()
completed=[]
for p in (root/'results').glob('*/protocol.json'):
    if json.loads(p.read_text()).get('status')=='completed_full_inference':
        completed.append(p.parent)
if not completed:
    raise SystemExit('Nenhuma inferência completa encontrada em results/.')
latest=max(completed,key=lambda p:p.stat().st_mtime_ns)
out=root/f'ADMET_Resultados_Local_{latest.name}.zip'
if out.exists():
    raise SystemExit(f'Arquivo já existe; preservado: {out.name}')
with ZipFile(out,'w',ZIP_DEFLATED) as z:
    for p in sorted(latest.rglob('*')):
        if p.is_file(): z.write(p,p.relative_to(root))
    log=root/'logs'/f'{latest.name}.log'
    if log.is_file(): z.write(log,log.relative_to(root))
    for name in ['README_PT.md','checksums.json','run_admet.py','requirements.txt']:
        z.write(root/name,'executed_package/'+name)
print(f'Envie este ZIP para a revisão do artigo: {out}')
PY
