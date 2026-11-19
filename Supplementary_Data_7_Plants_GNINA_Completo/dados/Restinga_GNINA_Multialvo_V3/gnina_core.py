"""Pinned GNINA execution and symmetry-aware, unaligned pose comparisons."""
import csv,hashlib,json,os,site,subprocess,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from rdkit import Chem,RDLogger
from rdkit.Chem import rdMolAlign
ROOT=Path(__file__).resolve().parent
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def write_csv(p,rows):
    assert rows
    with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def read_csv(p):return list(csv.DictReader(p.open()))
def load(p):
    mols=list(Chem.SDMolSupplier(str(p),removeHs=False));assert mols and all(m is not None for m in mols),'Invalid SDF '+str(p)
    return mols
def chemical_key(m):
    # InChI normalizes redundant stereotags sometimes generated for symmetric cage/ring atoms.
    return Chem.MolToInchiKey(Chem.RemoveHs(Chem.Mol(m)))
def rmsd(ref,probe):
    a=Chem.RemoveHs(Chem.Mol(ref));b=Chem.RemoveHs(Chem.Mol(probe));assert chemical_key(a)==chemical_key(b),'Chemical identity/microstate mismatch'
    return float(rdMolAlign.CalcRMS(b,a,maxMatches=1000000,symmetrizeConjugatedTerminalGroups=True))
def env(cpu=4):
    e=os.environ.copy();libs=[str(p) for s in site.getsitepackages() for p in (Path(s)/'nvidia').glob('*/lib')]
    e['LD_LIBRARY_PATH']=':'.join(libs+[e.get('LD_LIBRARY_PATH','')]);e['OMP_NUM_THREADS']=str(cpu);e['MKL_NUM_THREADS']=str(cpu);e['OPENBLAS_NUM_THREADS']='1';return e
def portable_fingerprint(fp):
    result=dict(fp);cmd=list(result['command']);cmd[0]='<PINNED_GNINA_BINARY>'
    for opt,label in [('-r','<RECEPTOR_WITH_RECORDED_HASH>'),('-l','<LIGAND_WITH_RECORDED_HASH>'),('-o','<OUTPUT>')]:
        cmd[cmd.index(opt)+1]=label
    result['command']=cmd;return result
def run(binary,cfg,ligand,out,seed,resume=False,exhaust=16,modes=20):
    out.parent.mkdir(parents=True,exist_ok=True);done=out.with_suffix('.complete.json')
    cmd=[binary,'-r',str(ROOT/cfg['receptor']),'-l',str(ligand)]
    for ax,v in zip('xyz',cfg['center']):cmd += ['--center_'+ax,str(v)]
    for ax,v in zip('xyz',cfg['size']):cmd += ['--size_'+ax,str(v)]
    cmd += ['--cnn_scoring','rescore','--cnn','fast','--pose_sort_order','Energy','--exhaustiveness',str(exhaust),'--num_modes',str(modes),'--seed',str(seed),'--cpu','4','--no_gpu','-o',str(out)]
    fingerprint={'receptor_sha256':sha(ROOT/cfg['receptor']),'ligand_sha256':sha(ligand),'command':cmd}
    if out.exists():
        assert resume and done.exists(),'Refuse overwrite/incomplete output '+str(out)
        d=json.loads(done.read_text());assert portable_fingerprint(d['fingerprint'])==portable_fingerprint(fingerprint) and sha(out)==d['output_sha256'];return load(out),d['elapsed_s']
    out.with_suffix('.command.json').write_text(json.dumps({'command':cmd,'start_utc':datetime.now(timezone.utc).isoformat(),'fingerprint':fingerprint},indent=2))
    print('START '+out.stem,flush=True);start=time.monotonic()
    with out.with_suffix('.log').open('w') as log:r=subprocess.run(cmd,env=env(),stdout=log,stderr=subprocess.STDOUT,timeout=10800)
    assert r.returncode==0,'GNINA failed '+str(out.with_suffix('.log'))
    poses=load(out);elapsed=round(time.monotonic()-start,3)
    done.write_text(json.dumps({'fingerprint':fingerprint,'output_sha256':sha(out),'elapsed_s':elapsed},indent=2));print(f'COMPLETE {out.stem} poses={len(poses)} seconds={elapsed}',flush=True);return poses,elapsed
