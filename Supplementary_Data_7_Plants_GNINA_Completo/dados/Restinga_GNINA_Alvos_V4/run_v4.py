"""V4: Validate six new targets and screen the clean curated library in three seeds."""
import argparse,concurrent.futures,json,statistics,subprocess
from pathlib import Path
import numpy as np
from rdkit import Chem,RDLogger
import gnina_core as core
ROOT=Path(__file__).resolve().parent
def validate(binary,cfg,resume):
    target=cfg['target_id'];dst=ROOT/'results'/target/'redocking';ref=core.load(ROOT/cfg['native_reference'])[0];probe=Chem.Mol(ref)
    for i,p in enumerate(probe.GetConformer().GetPositions()):probe.GetConformer().SetAtomPosition(i,(float(p[0])+10,float(p[1]),float(p[2])))
    assert abs(core.rmsd(ref,ref))<1e-6 and abs(core.rmsd(ref,probe)-10)<1e-5
    rows=[];runs=[]
    def one(seed):return seed,*core.run(binary,cfg,ROOT/cfg['redock_input'],dst/f'{target}_native_seed{seed}.sdf',seed,resume)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        for seed,poses,elapsed in ex.map(one,[101,202,303,404,505]):
            ds=[]
            for rank,m in enumerate(poses,1):
                d=core.rmsd(ref,m);ds.append(d);rows.append({'target_id':target,'seed':seed,'energy_rank':rank,'heavy_atom_fixed_frame_symmetry_RMSD_A':d,'minimizedAffinity_kcal_mol':float(m.GetProp('minimizedAffinity')),'CNNscore_fast':float(m.GetProp('CNNscore'))})
            runs.append({'target_id':target,'seed':seed,'top1_RMSD_A':ds[0],'best_sampled_RMSD_A':min(ds),'top1_within_2A':ds[0]<=2,'elapsed_s':elapsed});core.write_csv(dst/'poses.csv',rows);core.write_csv(dst/'runs.csv',runs)
    passed=sum(r['top1_within_2A'] for r in runs)
    validation={'target_id':target,'runs':5,'top1_passes':passed,'required_passes':4,'median_top1_RMSD_A':statistics.median(r['top1_RMSD_A'] for r in runs),'prespecified_rule':'Energy_top1_heavy_atom_symmetry_fixed_frame_RMSD_le_2A_in_at_least_4_of_5_seeds','geometry_status':'PASS' if passed>=4 else 'FAIL','screening_eligible':passed>=4,'limits':cfg.get('preparation_note',cfg.get('biological_scope','')),'input_chemical_key':core.chemical_key(core.load(ROOT/cfg['redock_input'])[0]),'reference_chemical_key':core.chemical_key(ref)}
    (dst.parent/'validation.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation),flush=True)
def screen(binary,cfg,resume):
    target=cfg['target_id'];dst=ROOT/'results'/target;validation=json.loads((dst/'validation.json').read_text())
    if not validation['screening_eligible']:print('SKIP failed redocking '+target,flush=True);return
    inputs=core.read_csv(ROOT/'docking_inputs.csv');metas={r['input_id']:r for r in inputs};expected={r['input_id']:core.chemical_key(core.load(ROOT/r['input_path'])[0]) for r in inputs}
    rows=[]
    def one(seed):return seed,*core.run(binary,cfg,ROOT/'prepared/all_inputs.sdf',dst/'docking'/f'{target}_batch_seed{seed}.sdf',seed,resume)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        for seed,poses,elapsed in ex.map(one,[42,142,242]):
            groups={}
            for m in poses:groups.setdefault(m.GetProp('_Name'),[]).append(m)
            assert set(groups)==set(expected),'Missing/extra batch compound groups'
            for iid,mm in groups.items():
                assert all(core.chemical_key(m)==expected[iid] for m in mm),'Chemical identity mismatch '+iid
                energies=[float(m.GetProp('minimizedAffinity')) for m in mm];assert all(a<=b+1e-4 for a,b in zip(energies,energies[1:]))
                cnnrank={idx:rank for rank,idx in enumerate(sorted(range(len(mm)),key=lambda i:float(mm[i].GetProp('CNNscore')),reverse=True),1)}
                for rank,m in enumerate(mm,1):rows.append({'target_id':target,'input_id':iid,'structure_id':metas[iid]['structure_id'],'compound':metas[iid]['compound'],'panel':metas[iid]['panel'],'seed':seed,'energy_rank':rank,'CNNscore_rank_within_energy_exported_pool':cnnrank[rank-1],'minimizedAffinity_kcal_mol':float(m.GetProp('minimizedAffinity')),'CNNscore_fast':float(m.GetProp('CNNscore')),'CNNaffinity_fast':float(m.GetProp('CNNaffinity')),'sample_presence':'not_confirmed','native_redocking_geometry_status':validation['geometry_status']})
            core.write_csv(dst/'docking_poses.csv',rows);print(f'SCREEN_SEED_COMPLETE {target} seed={seed} compounds={len(groups)} poses={len(poses)}',flush=True)
    (dst/'screening_complete.json').write_text(json.dumps({'target_id':target,'inputs_per_seed':len(inputs),'seeds':[42,142,242],'exported_poses':len(rows),'status':'COMPLETE_EXPLORATORY'},indent=2))
def main():
    RDLogger.DisableLog('rdApp.warning');RDLogger.DisableLog('rdApp.error')
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['validate','screen','all']);ap.add_argument('--gnina',required=True);ap.add_argument('--resume',action='store_true');ap.add_argument('--targets',nargs='+');ap.add_argument('--historical-input-reproduction',action='store_true');args=ap.parse_args();binary=str(Path(args.gnina).resolve())
    manifest=json.loads((ROOT/'gnina_binary_manifest.json').read_text());assert core.sha(binary)==manifest['sha256']
    if args.stage in ['all','screen']:
        ss={r['structure_id']:r for r in core.read_csv(ROOT/'structures.csv')};mismatches=[]
        for row in core.read_csv(ROOT/'docking_inputs.csv'):
            expected=core.chemical_key(Chem.MolFromSmiles(ss[row['structure_id']]['isomeric_smiles']))
            observed=core.load(ROOT/row['input_path'])[0];declared=Chem.MolFromSmiles(ss[row['structure_id']]['isomeric_smiles'])
            if core.chemical_key(observed)!=expected or core.microstate_graph_key(observed)!=core.microstate_graph_key(declared):mismatches.append(row['input_id'])
        if mismatches and not args.historical_input_reproduction:raise RuntimeError('Quarantined input identity mismatch. Use the clean curated reproduction folder; historical replay requires explicit --historical-input-reproduction. Inputs: '+','.join(mismatches))
    round_ids=set(json.loads((ROOT/'selection_plan.json').read_text())['target_ids'])
    assert len(round_ids)==6 and all((ROOT/'targets'/tid/'config.json').exists() for tid in round_ids),'Expected six initial round target preparations'
    cfgs=[json.loads(p.read_text()) for p in sorted((ROOT/'targets').glob('*/config.json'))];cfgs=[c for c in cfgs if c['target_id'] in round_ids and (not args.targets or c['target_id'] in args.targets)]
    assert cfgs,'No targets selected in this round'
    protocol={'GNINA_version':subprocess.run([binary,'--version'],env=core.env(),capture_output=True,text=True,check=True).stdout.strip(),'binary_sha256':manifest['sha256'],'cnn_scoring':'rescore','cnn_model':'fast_builtin_pinned_GNINA_1.3.3','primary_pose_order':'Energy','exhaustiveness':16,'num_modes':20,'cpu':4,'gpu':False,'docking_seeds':[42,142,242],'redocking_seeds':[101,202,303,404,505],'redocking_rule':'top1_fixed_frame_symmetry_heavy_atom_RMSD<=2A_in_at_least_4/5_seeds','library_manifest_sha256':core.sha(ROOT/'docking_inputs.csv'),'batch_policy':'GNINA internally resets reported random seed for each ligand; same ordered complete library in each batch','comparison_limits':'Scores not compared as affinity across targets; fast-CNN scores not numerically comparable with previous builtin-default ensemble','targets':{c['target_id']:{'config':c,'receptor_sha256':core.sha(ROOT/c['receptor']),'reference_sha256':core.sha(ROOT/c['native_reference']),'redock_input_sha256':core.sha(ROOT/c['redock_input'])} for c in cfgs}}
    pp=ROOT/'protocol.json'
    if pp.exists():
        old=json.loads(pp.read_text());assert old==protocol,'Protocol mismatch; use separate study folder'
    else:pp.write_text(json.dumps(protocol,indent=2))
    for c in cfgs:
        if args.stage in ['all','validate']:validate(binary,c,args.resume)
    if args.stage in ['all','screen']:
        for c in cfgs:screen(binary,c,args.resume)

if __name__=='__main__':main()
