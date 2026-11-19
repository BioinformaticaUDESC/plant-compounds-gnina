#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
command -v conda >/dev/null || { echo 'Conda não encontrado. Abra um terminal com conda disponível.' >&2; exit 1; }
if conda env list --json | python3 -c 'import json,sys; from pathlib import Path; sys.exit(0 if any(Path(p).name=="admet_plantas_v1" for p in json.load(sys.stdin)["envs"]) else 1)'; then
    echo 'O ambiente admet_plantas_v1 já existe. Não foi alterado. Ative-o e confira as versões antes de executar.' >&2
    exit 1
fi
# Ambiente novo; não altera ambientes de docking/MD existentes.
conda create -y -n admet_plantas_v1 python=3.12 pip
conda run --no-capture-output -n admet_plantas_v1 python -m pip install 'torch==2.10.0' --index-url https://download.pytorch.org/whl/cpu
conda run --no-capture-output -n admet_plantas_v1 python -m pip install -r requirements.txt
conda run --no-capture-output -n admet_plantas_v1 python -m pip check
echo 'Instalação concluída. Execute: conda activate admet_plantas_v1'
echo 'Depois: bash run_local.sh'
