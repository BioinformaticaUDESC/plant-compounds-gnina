"""Copy fixed inputs and code into a new study folder without computed results."""
import argparse,shutil,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--destination',required=True,type=Path);a=ap.parse_args();dst=a.destination.expanduser().resolve()
 if dst.exists():raise FileExistsError('Destination already exists; preserve original study: '+str(dst))
 dst.mkdir(parents=True)
 for p in ROOT.iterdir():
  if p.is_file() and p.suffix in ['.py','.csv','.json','.lock','.freeze','.sh','.md'] and p.name not in ['pipeline_complete.json','delivery_checks.json','final_code_sha256.json','method_checks.json','execution_environment.json']:shutil.copy2(p,dst/p.name)
 for name in ['prepared','raw','baseline_v3','previous_delivery']:
  if (ROOT/name).exists():shutil.copytree(ROOT/name,dst/name,ignore=shutil.ignore_patterns('__pycache__'))
 (dst/'targets').mkdir()
 for source in sorted((ROOT/'targets').iterdir()):
  if not source.is_dir() or not (source/'config.json').exists():continue
  target=dst/'targets'/source.name;target.mkdir()
  for name in ['config.json','preparation_audit.json','receptor_H.pdb','native_reference.sdf','redock_input.sdf','native_crystal_heavy.pdb']:shutil.copy2(source/name,target/name)
 for name in ['reports','results','execution_logs']:(dst/name).mkdir()
 shutil.copy2(ROOT/'reports/input_identity_qc.csv',dst/'reports/input_identity_qc.csv')
 (dst/'clean_reproduction_origin.json').write_text(json.dumps({'created_from':ROOT.name,'new_V4_computed_results_copied':False,'previous_V3_results_retained_for_provenance':True,'verified_inputs':59,'target_configs':8,'scope':'Fresh rerun of fixed V4_A/B, preserving previous V3 baseline provenance'},indent=2))
 print(dst,flush=True)
if __name__=='__main__':main()
