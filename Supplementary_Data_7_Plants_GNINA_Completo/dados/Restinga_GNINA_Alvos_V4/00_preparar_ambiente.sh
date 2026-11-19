#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
task_python="${PYTHON_BIN:-python3}"
"$task_python" -c 'import sys; assert sys.version_info >= (3,12), "Use Python 3.12 ou mais recente"'
if [[ ! -d .venv ]]; then "$task_python" -m venv .venv; fi
.venv/bin/python -m pip install -r requirements-analysis.lock
if [[ "${INSTALL_REPORT:-0}" == "1" ]]; then .venv/bin/python -m pip install -r requirements-report.lock; fi
mkdir -p tools
if [[ -z "${GNINA_BIN:-}" ]]; then
.venv/bin/python - <<'PY'
import hashlib,json,os,urllib.request
from pathlib import Path
meta=json.loads(Path('gnina_binary_manifest.json').read_text());p=Path('tools/gnina-1.3.3')
if not p.exists():
 print('Obtendo binário oficial GNINA, aproximadamente 2 GB.',flush=True)
 with urllib.request.urlopen(meta['url'],timeout=90) as r,p.with_suffix('.part').open('wb') as f:
  while b:=r.read(8*1024*1024):f.write(b)
 os.replace(p.with_suffix('.part'),p)
h=hashlib.sha256()
with p.open('rb') as f:
 while b:=f.read(8*1024*1024):h.update(b)
assert h.hexdigest()==meta['sha256'],'Executável diferente do utilizado nesta rodada'
p.chmod(0o755);print('GNINA conferido.')
PY
fi
printf '%s\n' 'Ambiente preparado. Consulte README.md para preservar/reutilizar os resultados entregues.'
