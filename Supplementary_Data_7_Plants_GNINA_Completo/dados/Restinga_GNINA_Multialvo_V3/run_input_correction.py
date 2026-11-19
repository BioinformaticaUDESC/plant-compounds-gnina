"""Run the chemically verified replacement as a separately labeled diagnostic."""
import argparse,concurrent.futures,json
from pathlib import Path
from rdkit import RDLogger
import gnina_core as core
ROOT=Path(__file__).resolve().parent
def main():
    RDLogger.DisableLog('rdApp.warning');RDLogger.DisableLog('rdApp.error')
    ap=argparse.ArgumentParser();ap.add_argument('--gnina',required=True);ap.add_argument('--resume',action='store_true');ap.add_argument('--targets',nargs='+');ap.add_argument('--workers',type=int,choices=[1,2],default=2);args=ap.parse_args();binary=str(Path(args.gnina).resolve())
    protocol=json.loads((ROOT/'protocol.json').read_text());assert core.sha(binary)==protocol['binary_sha256'];inp=core.read_csv(ROOT/'correction_docking_inputs.csv')[0];assert core.sha(ROOT/inp['input_path'])==inp['sha256'];expected=core.chemical_key(core.load(ROOT/inp['input_path'])[0])
    plan=[]
    for cfgpath in sorted((ROOT/'targets').glob('*/config.json')):
        cfg=json.loads(cfgpath.read_text());tid=cfg['target_id'];dst=ROOT/'results'/tid;v=json.loads((dst/'validation.json').read_text())
        if args.targets and tid not in args.targets:continue
        if not v['screening_eligible']:continue
        assert (dst/'screening_complete.json').exists(),'Run initial screen before correction diagnostic'
        for seed in [42,142,242]:plan.append((cfg,seed))
    corrected_protocol={'input_curation_sha256':core.sha(ROOT/'reports/input_curation.json'),'corrected_input_sha256':inp['sha256'],'original_protocol_sha256':core.sha(ROOT/'protocol.json'),'purpose':'separate verified-identity diagnostic; inclusion/exclusion from primary determined by identity/diversity curation before corrected-input docking scores','seeds':[42,142,242]};(ROOT/'correction_protocol.json').write_text(json.dumps(corrected_protocol,indent=2))
    rows={}
    def one(job):
        cfg,seed=job;tid=cfg['target_id'];out=ROOT/'results'/tid/'docking'/f'{tid}_correction_seed{seed}.sdf'
        poses,elapsed=core.run(binary,cfg,ROOT/inp['input_path'],out,seed,args.resume);assert all(core.chemical_key(m)==expected for m in poses)
        cnnrank={i:r for r,i in enumerate(sorted(range(len(poses)),key=lambda i:float(poses[i].GetProp('CNNscore')),reverse=True),1)}
        rr=[{'target_id':tid,'input_id':inp['input_id'],'structure_id':inp['structure_id'],'compound':inp['compound'],'panel':inp['panel'],'seed':seed,'energy_rank':i+1,'CNNscore_rank_within_energy_exported_pool':cnnrank[i],'minimizedAffinity_kcal_mol':float(m.GetProp('minimizedAffinity')),'CNNscore_fast':float(m.GetProp('CNNscore')),'CNNaffinity_fast':float(m.GetProp('CNNaffinity')),'sample_presence':'not_confirmed','native_redocking_geometry_status':'PASS'} for i,m in enumerate(poses)]
        return tid,rr
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        for tid,rr in ex.map(one,plan):rows.setdefault(tid,[]).extend(rr)
    for tid,rr in rows.items():
        dst=ROOT/'results'/tid;core.write_csv(dst/'correction_docking_poses.csv',rr);(dst/'correction_complete.json').write_text(json.dumps({'corrected_input':inp['input_id'],'input_sha256':inp['sha256'],'seeds':[42,142,242],'exported_poses':len(rr),'protocol_sha256':core.sha(ROOT/'correction_protocol.json')},indent=2))
    print(json.dumps({tid:len(rr) for tid,rr in rows.items()}),flush=True)

if __name__=='__main__':main()
