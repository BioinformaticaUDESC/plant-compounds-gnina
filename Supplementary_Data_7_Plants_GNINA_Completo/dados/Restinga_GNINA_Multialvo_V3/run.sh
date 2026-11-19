#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
task_stage="${1:-validate}"
if [[ $# -gt 0 ]]; then shift; fi
if [[ ! -x .venv/bin/python ]]; then printf '%s\n' 'Execute: bash 00_preparar_ambiente.sh'; exit 1; fi
task_binary="${GNINA_BIN:-$PWD/tools/gnina-1.3.3}"
case "$task_stage" in
 validate|screen|all) exec .venv/bin/python -u run_v3.py "$task_stage" --gnina "$task_binary" "$@" ;;
 audit) .venv/bin/python -u audit_v3.py "$@"; .venv/bin/python synthesize_v3.py ;;
 report) exec .venv/bin/python build_report.py ;;
 *) printf '%s\n' 'Uso: bash run.sh validate|screen|all|audit|report [--resume]'; exit 2 ;;
esac
