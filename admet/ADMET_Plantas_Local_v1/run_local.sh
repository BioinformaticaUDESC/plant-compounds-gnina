#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export CUDA_VISIBLE_DEVICES=''
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
run_id="$(date -u +%Y%m%dT%H%M%SZ)_$$"
mkdir -p logs results
python run_admet.py --mode full --threads 4 --out "results/$run_id" 2>&1 | tee "logs/$run_id.log"
echo "Resultados: results/$run_id"
