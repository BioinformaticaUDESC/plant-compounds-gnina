"""Run docking first, then independent redocking, and report fixed-frame RMSD."""
import argparse,csv,hashlib,json,os,shutil,site,subprocess,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

ROOT=Path(__file__).resolve().parent
def env_for_gnina(cpu):
    env=os.environ.copy()
    libs=[str(p) for s in site.getsitepackages() for p in (Path(s)/'nvidia').glob('*/lib')]
    env['LD_LIBRARY_PATH']=':'.join(libs+[env.get('LD_LIBRARY_PATH','')])
    for name in ['OMP_NUM_THREADS','MKL_NUM_THREADS']:env[name]=str(cpu)
    env['OPENBLAS_NUM_THREADS']='1'
    return env

def save_csv(path,rows):
    if not rows:raise RuntimeError('Empty result table')
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def load_poses(path):
    mols=list(Chem.SDMolSupplier(str(path),removeHs=False))
    if not mols or any(m is None for m in mols):raise RuntimeError(f'Invalid SDF poses: {path}')
    return mols

def score(m,key):
    if not m.HasProp(key):return ''
    try:return float(m.GetProp(key))
    except ValueError:return m.GetProp(key)

def rmsd(reference,probe):
    a=Chem.RemoveHs(Chem.Mol(reference));b=Chem.RemoveHs(Chem.Mol(probe))
    if Chem.MolToSmiles(a,isomericSmiles=True)!=Chem.MolToSmiles(b,isomericSmiles=True):
        raise RuntimeError('Reference and pose chemical graphs differ; RMSD not computed')
    return float(rdMolAlign.CalcRMS(b,a,maxMatches=1000000,symmetrizeConjugatedTerminalGroups=True))

def run_one(binary,config,ligand,out,seed,args,env):
    out.parent.mkdir(parents=True,exist_ok=True)
    if out.exists() and not args.overwrite:raise RuntimeError(f'Refusing to overwrite {out}')
    cmd=[binary,'-r',str(ROOT/config['receptor']),'-l',str(ligand)]
    for axis,val in zip('xyz',config['center']):cmd += [f'--center_{axis}',str(val)]
    for axis,val in zip('xyz',config['size']):cmd += [f'--size_{axis}',str(val)]
    cmd += ['--cnn_scoring','rescore','--exhaustiveness',str(args.exhaustiveness),'--num_modes',str(args.modes),'--seed',str(seed),'--cpu',str(args.cpu),'--no_gpu','-o',str(out)]
    out.with_suffix('.command.json').write_text(json.dumps({'command':cmd,'start_utc':datetime.now(timezone.utc).isoformat()},indent=2))
    print(f'START {out.stem}',flush=True);start=time.monotonic()
    with out.with_suffix('.log').open('w') as log:
        r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env=env,timeout=args.timeout)
    if r.returncode!=0:raise RuntimeError(f'GNINA exited {r.returncode}; see {out.with_suffix(".log")}')
    poses=load_poses(out);elapsed=round(time.monotonic()-start,3)
    print(f'COMPLETE {out.stem} poses={len(poses)} seconds={elapsed}',flush=True)
    return poses,elapsed

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('stage',choices=['preflight','dock','redock','all'])
    ap.add_argument('--gnina',default=os.environ.get('GNINA_BIN') or shutil.which('gnina'))
    ap.add_argument('--cpu',type=int,default=4)
    ap.add_argument('--exhaustiveness',type=int,default=16)
    ap.add_argument('--modes',type=int,default=20)
    ap.add_argument('--timeout',type=int,default=900)
    ap.add_argument('--overwrite',action='store_true')
    ap.add_argument('--redock-seeds',nargs='+',type=int,default=[101,202,303,404,505])
    args=ap.parse_args()
    if not args.gnina:ap.error('Specify --gnina /absolute/path/to/gnina')
    if args.cpu<1 or args.exhaustiveness<1 or args.modes<1:ap.error('Positive execution parameters required')
    binary=str(Path(args.gnina).expanduser().resolve())
    cfg=json.loads((ROOT/'config.json').read_text());env=env_for_gnina(args.cpu)
    results=ROOT/'results';results.mkdir(exist_ok=True)
    version=subprocess.run([binary,'--version'],capture_output=True,text=True,env=env,check=True).stdout.strip()
    fingerprint=hashlib.sha256((ROOT/cfg['receptor']).read_bytes()).hexdigest()
    protocol={'gnina_version':version,'receptor_sha256':fingerprint,'exhaustiveness':args.exhaustiveness,'num_modes':args.modes,'cnn_scoring':'rescore','cnn_model':'builtin_default_of_pinned_binary','pose_order':'GNINA default CNNscore','cpu':args.cpu,'gpu':False,'box':cfg,'cofactor_model':'deposited COH retained; general empirical/CNN metal representation, no custom electronic model','production_status':'exploratory_no_GCMS_assignment'}
    reference=load_poses(ROOT/'prepared/native_reference.sdf')[0]
    # A translated ligand must remain displaced. An aligning RMSD implementation fails this check.
    shifted=Chem.Mol(reference)
    for i,pos in enumerate(shifted.GetConformer().GetPositions()):shifted.GetConformer().SetAtomPosition(i,(float(pos[0])+10,float(pos[1]),float(pos[2])))
    if abs(rmsd(reference,reference))>1e-6 or abs(rmsd(reference,shifted)-10)>1e-5:raise RuntimeError('Fixed-frame RMSD verification failed')
    if args.stage=='preflight':
        cmd=[binary,'-r',str(ROOT/cfg['receptor']),'-l',str(ROOT/'prepared/native_reference.sdf'),'--score_only','--cnn_scoring','rescore','--no_gpu']
        with (results/'preflight.log').open('w') as log:
            r=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout)
        if r.returncode:raise RuntimeError('GNINA preflight failed; inspect results/preflight.log')
        print('PREFLIGHT_OK',flush=True);return
    if args.stage in ['dock','all']:
        if args.overwrite:
            (results/'validation.json').unlink(missing_ok=True)
        (results/'protocol.json').write_text(json.dumps(protocol,indent=2))
        rows=[]
        for path in sorted((ROOT/'prepared/ligands').glob('*.sdf')):
            poses,elapsed=run_one(binary,cfg,path,results/'docking'/(path.stem+'_seed42.sdf'),42,args,env)
            for rank,m in enumerate(poses,1):rows.append({'ligand_id':path.stem,'seed':42,'pose_rank':rank,'minimizedAffinity_kcal_mol':score(m,'minimizedAffinity'),'CNNscore':score(m,'CNNscore'),'CNNaffinity':score(m,'CNNaffinity'),'elapsed_s':elapsed,'sample_presence':'not_confirmed','validation_status':'pending_redocking'})
            save_csv(results/'docking_poses.csv',rows)
        (results/'docking_stage.json').write_text(json.dumps({'completed_ligands':len(set(r['ligand_id'] for r in rows)),'completed_poses':len(rows),'status':'complete_exploratory_pending_redocking'},indent=2))
    if args.stage in ['redock','all']:
        if not (results/'docking_stage.json').exists():raise RuntimeError('Docking stage must finish before redocking in this workflow')
        previous=json.loads((results/'protocol.json').read_text())
        if previous!=protocol:raise RuntimeError('Docking and redocking protocols differ; cannot validate prior docking')
        rows=[];runs=[]
        for seed in args.redock_seeds:
            poses,elapsed=run_one(binary,cfg,ROOT/'prepared/redock_input.sdf',results/'redocking'/f'ID8_seed{seed}.sdf',seed,args,env)
            distances=[]
            for rank,m in enumerate(poses,1):
                d=rmsd(reference,m);distances.append(d)
                rows.append({'seed':seed,'pose_rank':rank,'heavy_atom_fixed_frame_symmetry_RMSD_A':d,'minimizedAffinity_kcal_mol':score(m,'minimizedAffinity'),'CNNscore':score(m,'CNNscore'),'CNNaffinity':score(m,'CNNaffinity')})
            runs.append({'seed':seed,'top1_RMSD_A':distances[0],'best_sampled_RMSD_A':min(distances),'top1_within_2A':distances[0]<=2.0,'elapsed_s':elapsed})
            save_csv(results/'redocking_poses.csv',rows);save_csv(results/'redocking_runs.csv',runs)
        passed=sum(r['top1_within_2A'] for r in runs)
        threshold=max(1,int(np.ceil(.8*len(runs))))
        validation={'runs':len(runs),'top1_passes':passed,'required_top1_passes':threshold,'predefined_rule':'top1 heavy-atom symmetry-aware fixed-frame RMSD <=2 Å in >=80% of seeds','top1_median_RMSD_A':float(np.median([r['top1_RMSD_A'] for r in runs])),'status':'GEOMETRY_CRITERION_MET_PILOT_ONLY' if passed>=threshold else 'GEOMETRY_CRITERION_NOT_MET','limitations':['literature candidate ligands not confirmed in samples','single rigid receptor preparation','waters and glycans excluded','deposited cobalt porphyrin approximated by standard docking atom types','redocking does not validate binding affinity selectivity safety or efficacy']}
        (results/'validation.json').write_text(json.dumps(validation,indent=2))
        print(json.dumps(validation,indent=2),flush=True)

if __name__=='__main__':main()
