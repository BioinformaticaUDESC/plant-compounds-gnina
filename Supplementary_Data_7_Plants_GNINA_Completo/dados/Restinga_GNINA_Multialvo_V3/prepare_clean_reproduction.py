"""Create a fresh portable study containing only the verified primary input library."""
import argparse,csv,json,shutil,subprocess,sys
from pathlib import Path
from rdkit import Chem
import gnina_core as core
ROOT=Path(__file__).resolve().parent
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--destination',required=True);args=ap.parse_args();dst=Path(args.destination).expanduser().resolve()
    if dst.exists():raise RuntimeError('Destination exists; choose a new directory')
    rows=core.read_csv(ROOT/'primary_docking_inputs.csv');ss={s['structure_id']:s for s in core.read_csv(ROOT/'structures.csv')}
    for row in rows:assert core.chemical_key(core.load(ROOT/row['input_path'])[0])==core.chemical_key(Chem.MolFromSmiles(ss[row['structure_id']]['isomeric_smiles']))
    dst.mkdir(parents=True);(dst/'prepared/ligands').mkdir(parents=True);(dst/'reports').mkdir();(dst/'results').mkdir()
    for name in ['targets','inputs','raw']:shutil.copytree(ROOT/name,dst/name)
    for p in ROOT.glob('*.py'):shutil.copy2(p,dst/p.name)
    for p in ROOT.glob('*.sh'):shutil.copy2(p,dst/p.name)
    for name in ['structures.csv','physicochemical_descriptors.csv','literature_species_links.csv','species_coverage.csv','literature_sources.json','target_rationale.json','gnina_binary_manifest.json','requirements-analysis.lock','requirements-report.lock']:
        shutil.copy2(ROOT/name,dst/name)
    for row in rows:shutil.copy2(ROOT/row['input_path'],dst/row['input_path'])
    for name in ['docking_inputs.csv','primary_docking_inputs.csv','all_attempted_docking_inputs.csv']:core.write_csv(dst/name,rows)
    with (dst/'correction_docking_inputs.csv').open('w',newline='') as f:csv.DictWriter(f,fieldnames=list(rows[0])).writeheader()
    for name in ['all_inputs.sdf','primary_inputs.sdf']:shutil.copy2(ROOT/'prepared/primary_inputs.sdf',dst/'prepared'/name)
    qc=[]
    for row in rows:
        expected=core.chemical_key(core.load(ROOT/row['input_path'])[0]);qc.append({'input_id':row['input_id'],'structure_id':row['structure_id'],'compound':row['compound'],'expected_chemical_key':expected,'observed_3D_input_chemical_key':expected,'status':'PASS','input_sha256':row['sha256'],'primary_analysis_eligible':True})
    core.write_csv(dst/'reports/input_identity_qc.csv',qc)
    (dst/'README_CLEAN.md').write_text('Fresh reproduction with only the verified primary library. No historical quarantined or diagnostic input is included in docking manifests. Run validate then screen, audit and synthesis using the local pinned GNINA and Python environment.\n')
    (dst/'source_primary_manifest_sha256.json').write_text(json.dumps({'source_manifest_sha256':core.sha(ROOT/'primary_docking_inputs.csv'),'inputs':len(rows),'compounds':len(set(r['compound'] for r in rows)),'source_study':'GNINA_V3_curated_primary_library'},indent=2));print(dst)

if __name__=='__main__':main()
