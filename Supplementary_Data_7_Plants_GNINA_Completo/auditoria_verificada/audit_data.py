from pathlib import Path
from collections import defaultdict,Counter
import json,csv,zipfile,hashlib,statistics,itertools,re,time,argparse,shutil
import numpy as np
from rdkit import Chem,RDLogger
from rdkit.Chem import rdMolAlign,rdMolDescriptors

ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser(description='Read-only verification of archived GNINA data; never runs GNINA or modifies source results.')
ap.add_argument('--data-root',type=Path,default=ROOT/'source_final/unpacked/Restinga_GNINA_Alvos_V4')
ap.add_argument('--v3-root',type=Path,default=ROOT/'figure_work/V3_source/Restinga_GNINA_Multialvo_V3')
ap.add_argument('--out',type=Path,default=Path(__file__).resolve().parent)
args=ap.parse_args();OUT=args.out.resolve();B=args.data_root.resolve();V3=args.v3_root.resolve();OUT.mkdir(parents=True,exist_ok=True)
fixture=Path(__file__).resolve().parent/'manuscript_tables.json'
if not (OUT/'manuscript_tables.json').exists():shutil.copy2(fixture,OUT/'manuscript_tables.json')
Z=zipfile.ZipFile(B/'previous_delivery/Restinga_GNINA_Multialvo_V3_40_Compostos_25-10-2026.zip')
RDLogger.DisableLog('rdApp.warning');RDLogger.DisableLog('rdApp.error')
checks=[];issues=[]
def check(label,a,b,tol=None):
    ok=abs(float(a)-float(b))<=tol if tol is not None else a==b
    r={'check':label,'status':'PASS' if ok else 'DIVERGENCE','recomputed':a,'archived_or_manuscript':b};checks.append(r)
    if not ok:issues.append(r);print('DIVERGENCE',label,a,b,flush=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return list(csv.DictReader(open(p)))
def write(name,rs):
    if rs:
        with open(OUT/name,'w') as f:w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
def load(p):
    ms=list(Chem.SDMolSupplier(str(p),removeHs=False));assert ms and all(m is not None for m in ms),str(p);return ms
def key(m):return Chem.MolToInchiKey(m)
def graph(m):return Chem.MolToSmiles(m,isomericSmiles=False,canonical=True)
def rms(a,b):return float(rdMolAlign.CalcRMS(b,a,maxMatches=1000000,symmetrizeConjugatedTerminalGroups=True))
def archived_json(src,path):
    return json.loads(Z.read('Restinga_GNINA_Multialvo_V3/'+str(path))) if src==V3 else json.load(open(src/path))
def archived_bytes(src,path):return Z.read('Restinga_GNINA_Multialvo_V3/'+str(path)) if src==V3 else (src/path).read_bytes()

# Every source-package checksum, without changing the archive.
hashes=0
for line in (B/'SHA256SUMS.txt').read_text().splitlines():
    if not line.strip():continue
    digest,name=line.split(None,1);name=name.lstrip('* ');p=B/name
    check('archive SHA256 '+name,sha(p) if p.is_file() else 'MISSING',digest);hashes+=1
print('SOURCE_CHECKSUMS',hashes,flush=True)
inputs={r['input_id']:r for r in read(B/'primary_docking_inputs.csv')};structures={r['structure_id']:r for r in read(B/'structures.csv')};inm={};alltop=read(B/'reports/all_top1_with_baseline.csv');topmap={(r['target_id'],r['input_id'],int(r['seed'])):r for r in alltop}
tables=json.load(open(OUT/'manuscript_tables.json'));name_to_tid={};vals={r['target_id']:r for r in read(B/'reports/all_targets_validation.csv')}
check('primary input count',len(inputs),59);check('declared structure count',len(structures),51);check('named compound count',len({r['compound'] for r in structures.values()}),40)
for iid,r in inputs.items():
    p=B/r['input_path'];check('input SHA256 '+iid,sha(p),r['sha256']);m=Chem.RemoveHs(load(p)[0]);inm[iid]=m;s=structures[r['structure_id']];check('input InChIKey '+iid,key(m),s['structure_inchikey']);check('input formula '+iid,rdMolDescriptors.CalcMolFormula(m),s['formula'])
native=[];candidate=[];pose_counts=[];stability=[];execution=[];configs={};pt=Chem.GetPeriodicTable();NOFLAG='no_flag_from_limited_geometry_audit';topm={};metalset={'Fe','Co','Zn','Mg','Ca','Mn','Ni','Cu'}
for tid in sorted(vals):
    ts=time.monotonic();src=V3 if tid in ['COX2_5IKR','ACHE_4EY7','CB2_5ZTY','ALOX5_6N2W'] else B
    cfg=json.load(open(src/'targets'/tid/'config.json'));configs[tid]=cfg;v=vals[tid]
    # All 60 native runs; all 1,200 poses and independent starting-state identity.
    ref=Chem.RemoveHs(load(src/cfg['native_reference'])[0]);ind=Chem.RemoveHs(load(src/cfg['redock_input'])[0]);check(tid+' independent native input identity',key(ind),key(ref));check(tid+' native input microstate',graph(ind),graph(ref))
    xyz=ref.GetConformer().GetPositions()
    if tid!='COX2_5IKR':
        check(tid+' box center',np.allclose((xyz.min(0)+xyz.max(0))/2,cfg['center'],atol=1e-8),True)
        check(tid+' box size',np.allclose(np.maximum(22,np.ptp(xyz,axis=0)+12),cfg['size'],atol=1e-8),True)
    for seed in [101,202,303,404,505]:
        rel=Path('results')/tid/'redocking'/f'{tid}_native_seed{seed}.sdf';p=src/rel;ms=load(p);check(tid+f' native modes seed{seed}',len(ms),20)
        for rank,raw in enumerate(ms,1):
            m=Chem.RemoveHs(raw);check(tid+f' native identity {seed}/{rank}',key(m),key(ref));check(tid+f' native microstate {seed}/{rank}',graph(m),graph(ref))
            native.append({'target_id':tid,'seed':seed,'rank':rank,'RMSD_A':rms(ref,m),'Energy_kcal_mol':float(raw.GetProp('minimizedAffinity'))})
        meta=archived_json(src,rel.with_suffix('.complete.json'));fp=meta['fingerprint'];check(str(rel)+' output SHA256',sha(p),meta['output_sha256']);check(str(rel)+' receptor SHA256',sha(src/cfg['receptor']),fp['receptor_sha256']);check(str(rel)+' ligand SHA256',sha(src/cfg['redock_input']),fp['ligand_sha256'])
        execution.append({'target_id':tid,'stage':'native','seed':seed,'command':fp['command'],'output_sha256':meta['output_sha256']})
    nt=[r for r in native if r['target_id']==tid and r['rank']==1];passes=sum(r['RMSD_A']<=2 for r in nt);median=statistics.median(r['RMSD_A'] for r in nt)
    check(tid+' five-seed pass count',passes,int(v['top1_passes']));check(tid+' median top1 RMSD',median,float(v['median_top1_RMSD_A']),1e-8);check(tid+' eligibility','PASS' if passes>=4 else 'FAIL',v['geometry_status'])
    sdfpaths=sorted((src/'results'/tid/'docking').glob('*.sdf'))
    if passes<4:
        check(tid+' failed control has no primary screen',len(sdfpaths),0);print('TARGET_AUDITED',tid,'native only',round(time.monotonic()-ts,2),flush=True);continue
    rec=[];rr=[];metals=[]
    for l in (src/cfg['receptor']).read_text().splitlines():
        if not l.startswith(('ATOM  ','HETATM')):continue
        el=l[76:78].strip().title()
        if el=='H':continue
        rec.append([float(l[30:38]),float(l[38:46]),float(l[46:54])]);rr.append(pt.GetRvdw(pt.GetAtomicNumber(el)));metals.append(el in metalset)
    rec=np.array(rec);rr=np.array(rr);nonmetal=~np.array(metals);center=np.array(cfg['center']);half=np.array(cfg['size'])/2
    count=0
    for p in sdfpaths:
        seed=int(p.stem.rsplit('seed',1)[1]);groups=defaultdict(list)
        for m in load(p):groups[m.GetProp('_Name')].append(m)
        if '_batch_seed' in p.name:check(tid+f' primary coverage seed{seed}',sorted(set(groups)&set(inputs)),sorted(inputs))
        rel=p.relative_to(src);meta=archived_json(src,rel.with_suffix('.complete.json'));check(str(rel)+' output SHA256',sha(p),meta['output_sha256']);fp=meta['fingerprint'];check(str(rel)+' receptor SHA256',sha(src/cfg['receptor']),fp['receptor_sha256'])
        ligand_rel=fp['command'][fp['command'].index('-l')+1].split('Restinga_GNINA_Multialvo_V3/' if src==V3 else 'Restinga_GNINA_Alvos_V4/')[-1]
        if Path(ligand_rel).is_absolute():
            # Source paths can name the original unpacked root; retain the portable suffix.
            for marker in ['prepared/','targets/']:
                if marker in ligand_rel:ligand_rel=marker+ligand_rel.split(marker,1)[1];break
        check(str(rel)+' batch ligand SHA256',hashlib.sha256(archived_bytes(src,Path(ligand_rel))).hexdigest(),fp['ligand_sha256'])
        execution.append({'target_id':tid,'stage':'candidate','seed':seed,'command':fp['command'],'output_sha256':meta['output_sha256']})
        for iid,ms in groups.items():
            if iid not in inputs:continue
            check(tid+f' {iid} modes seed{seed}',len(ms),20);energies=[float(m.GetProp('minimizedAffinity')) for m in ms];check(tid+f' {iid} Energy ordering {seed}',energies,sorted(energies));count+=len(ms)
            inp=inm[iid];ix=inp.GetConformer().GetPositions()
            for rank,raw in enumerate(ms,1):
                m=Chem.RemoveHs(raw);check(tid+f' {iid} stereo {seed}/{rank}',key(m),key(inp));check(tid+f' {iid} graph {seed}/{rank}',graph(m),graph(inp));x=m.GetConformer().GetPositions();check(tid+f' {iid} finite coordinates {seed}/{rank}',bool(np.isfinite(x).all()),True)
                if rank!=1:continue
                stored=topmap[(tid,iid,seed)]
                for prop,col in [('minimizedAffinity','minimizedAffinity_kcal_mol'),('CNNscore','CNNscore_fast'),('CNNaffinity','CNNaffinity_fast')]:check(tid+f' {iid} {prop} seed{seed}',float(raw.GetProp(prop)),float(stored[col]),1e-9)
                d=np.linalg.norm(x[:,None,:]-rec[None,:,:],axis=-1);lr=np.array([pt.GetRvdw(a.GetAtomicNum()) for a in m.GetAtoms()]);ratio=d/(lr[:,None]+rr[None,:]);close=int((d[:,nonmetal]<2).sum());over=int((ratio[:,nonmetal]<.7).sum());outside=float(np.maximum(np.abs(x-center)-half,0).max())
                mappings=inp.GetSubstructMatches(m,useChirality=False,uniquify=False,maxMatches=10000);assert mappings
                delta=min(max(abs(np.linalg.norm(x[b.GetBeginAtomIdx()]-x[b.GetEndAtomIdx()])-np.linalg.norm(ix[mp[b.GetBeginAtomIdx()]]-ix[mp[b.GetEndAtomIdx()]])) for b in m.GetBonds()) for mp in mappings)
                sd=np.linalg.norm(x[:,None,:]-x[None,:,:],axis=-1);td=Chem.GetDistanceMatrix(m);intra=int(((sd/(lr[:,None]+lr[None,:])<.7)&np.triu(td>3,k=1)).sum());flags=[]
                for condition,label in [(close,'nonmetal_heavy_atom_contacts_below_2A'),(over,'nonmetal_vdw_ratio_below_0.70'),(energies[0]>0,'positive_empirical_energy'),(outside>.1,'outside_box_above_0.1A'),(delta>.15,'bond_length_change_above_0.15A'),(intra,'severe_intraligand_nonbonded_overlap')]:
                    if condition:flags.append(label)
                q=';'.join(flags) or NOFLAG;check(tid+f' {iid} flags seed{seed}',q,stored['quality_flags'])
                for actual,col in [(close,'nonmetal_heavy_pairs_below_2A'),(over,'nonmetal_vdw_pairs_ratio_below_0.70'),(outside,'maximum_outside_box_A'),(delta,'maximum_mapped_bond_length_change_A'),(intra,'severe_intraligand_nonbonded_overlaps')]:check(tid+f' {iid} {col} seed{seed}',actual,float(stored[col]),1e-8)
                topm[(tid,iid,seed)]=m
                candidate.append({'target_id':tid,'input_id':iid,'compound':inputs[iid]['compound'],'structure_id':inputs[iid]['structure_id'],'panel':inputs[iid]['panel'],'seed':seed,'Energy':energies[0],'CNNscore':float(raw.GetProp('CNNscore')),'flags':q})
    check(tid+' primary candidate poses',count,3540);pose_counts.append({'target_id':tid,'primary_poses':count})
    print('TARGET_AUDITED',tid,'candidate poses',count,'seconds',round(time.monotonic()-ts,2),flush=True)

# Exact command parameters for all archived candidate and native runs.
for r in execution:
    cmd=r['command'];prefix=r['target_id']+' '+r['stage']+' '+str(r['seed'])
    for arg,expected in [('--cnn_scoring','rescore'),('--cnn','fast'),('--pose_sort_order','Energy'),('--exhaustiveness','16'),('--num_modes','20'),('--cpu','4'),('--seed',str(r['seed']))]:check(prefix+' '+arg,cmd[cmd.index(arg)+1] if arg in cmd else 'MISSING',expected)
    check(prefix+' CPU only','--no_gpu' in cmd,True)
    cfg=configs[r['target_id']]
    for axis,val in zip('xyz',cfg['center']):check(prefix+' center_'+axis,float(cmd[cmd.index('--center_'+axis)+1]),val,1e-10)
    for axis,val in zip('xyz',cfg['size']):check(prefix+' size_'+axis,float(cmd[cmd.index('--size_'+axis)+1]),val,1e-10)
stored_stability={(r['target_id'],r['input_id']):r for r in read(B/'reports/all_seed_stability_with_baseline.csv')}
for tid in sorted({r['target_id'] for r in candidate}):
    for iid in inputs:
        pairs=[rms(topm[(tid,iid,a)],topm[(tid,iid,b)]) for a,b in itertools.combinations([42,142,242],2)];mx=max(pairs);rec=stored_stability[(tid,iid)]
        check(tid+' '+iid+' maximum pairwise seed RMSD',mx,float(rec['maximum_top1_pairwise_fixed_frame_RMSD_between_seeds_A']),1e-8)
        stability.append({'target_id':tid,'input_id':iid,'max_pairwise_RMSD_A':mx,'all_seeds_without_flags':all(r['flags']==NOFLAG for r in candidate if r['target_id']==tid and r['input_id']==iid)})
stored_compounds={(r['target_id'],r['compound']):r for r in read(B/'reports/all_compound_summaries_with_baseline.csv')};groups=defaultdict(list);compound=[]
for r in candidate:groups[(r['target_id'],r['compound'])].append(r)
for (tid,name),rs in sorted(groups.items()):
    es=[r['Energy'] for r in rs];stored=stored_compounds[(tid,name)];md=statistics.median(es)
    for actual,col in [(md,'Energy_median_kcal_mol'),(min(es),'Energy_min_kcal_mol'),(max(es),'Energy_max_kcal_mol'),(len(rs),'top1_observations'),(sum(r['flags']!=NOFLAG for r in rs),'observations_with_geometry_flags')]:check(tid+' '+name+' '+col,actual,float(stored[col]),1e-9)
    compound.append({'target_id':tid,'compound':name,'panel':rs[0]['panel'],'observations':len(rs),'median_Energy':md,'minimum_Energy':min(es),'maximum_Energy':max(es),'flagged':sum(r['flags']!=NOFLAG for r in rs)})
check('all primary candidate poses',sum(r['primary_poses'] for r in pose_counts),24780);check('all top1 candidate observations',len(candidate),1239);check('all native control poses',len(native),1200);check('all top1 flagged observations',sum(r['flags']!=NOFLAG for r in candidate),33)
names={'COX-2':'COX2_5IKR','AChE':'ACHE_4EY7','CB2':'CB2_5ZTY','5-LOX':'ALOX5_6N2W','MAO-A':'MAOA_2Z5X','MAO-B':'MAOB_2V61','sEH':'SEH_3ANS','KEAP1':'KEAP1_4XMB','EGFR':'EGFR_1M17','GyrB':'GYRB_1KZN','PARP1':'PARP1_7KK4','TyrS':'TYRS_1JIL'}
def canonical(s):
    for a,b in [('α-','alpha-'),('β-','beta-'),('γ-','gamma-'),('δ-','delta-')]:s=s.replace(a,b)
    return s.lower()
compound_map={(r['target_id'],canonical(r['compound'])):r for r in compound}
for row in tables[1][1:]:
    tid=names[row[0]];rs=[r for r in native if r['target_id']==tid and r['rank']==1];check('Table 2 '+row[0]+' PDB',configs[tid]['pdb'],row[1]);check('Table 2 '+row[0]+' successes',str(sum(r['RMSD_A']<=2 for r in rs))+'/5',row[3]);check('Table 2 '+row[0]+' median',f'{statistics.median(r["RMSD_A"] for r in rs):.3f}',row[4]);check('Table 2 '+row[0]+' screened','Yes' if sum(r['RMSD_A']<=2 for r in rs)>=4 else 'No',row[5])
for row in tables[4][1:]:
    for protein,value in zip(tables[4][0][1:],row[1:]):check('Table S1 '+row[0]+' '+protein,f'{compound_map[(names[protein],canonical(row[0]))]["median_Energy"]:.2f}',value.replace('−','-'))
for row in tables[5][1:]:
    name=row[0].split('\n')[0];cfg=configs[names[name]]
    for field,values in [(3,cfg['center']),(4,cfg['size'])]:check('Table S2 '+name+' '+str(field),', '.join(f'{v:.4f}' for v in values),row[field])
for row in tables[7][1:]:
    tid=names[row[0]];rs={r['seed']:r['RMSD_A'] for r in native if r['target_id']==tid and r['rank']==1}
    for seed,value in zip([101,202,303,404,505],row[1:6]):check('Table S4 '+row[0]+' seed'+str(seed),f'{rs[seed]:.3f}',value)
for row in tables[8][1:]+tables[9][1:]:
    rs=[r for r in structures.values() if canonical(r['compound'])==canonical(row[0])];irs=[r for r in inputs.values() if canonical(r['compound'])==canonical(row[0])]
    check('Table S5 '+row[0]+' CID',set(r['pubchem_cid'] for r in rs),{row[1]});check('Table S5 '+row[0]+' charge',set(r['formal_charge'] for r in rs),{row[2]});check('Table S5 '+row[0]+' structures',len(rs),int(row[3]));check('Table S5 '+row[0]+' inputs',len(irs),int(row[4]));check('Table S5 '+row[0]+' observations',len(irs)*3,int(row[5]))
for row in tables[3][1:]:
    tid=names[row[0]];rs=[r for r in candidate if r['target_id']==tid];sr=[r for r in stability if r['target_id']==tid]
    check('Table 4 '+row[0]+' flagged',f'{sum(r["flags"]!=NOFLAG for r in rs)} / {len(rs)}',row[1])
    if row[2]!='Not summarized':check('Table 4 '+row[0]+' clear inputs',f'{sum(r["all_seeds_without_flags"] for r in sr)} / 59',row[2]);check('Table 4 '+row[0]+' repeatable inputs',f'{sum(r["max_pairwise_RMSD_A"]<=2 for r in sr)} / 59',row[3])
write('recomputed_native_poses.csv',native);write('recomputed_candidate_top1.csv',candidate);write('recomputed_compound_summaries.csv',compound);write('recomputed_seed_stability.csv',stability)
# Sets are serialized explicitly for an auditable check log.
(OUT/'data_checks.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2,default=lambda x:sorted(x) if isinstance(x,set) else str(x)))
(OUT/'data_divergences.json').write_text(json.dumps(issues,ensure_ascii=False,indent=2,default=lambda x:sorted(x) if isinstance(x,set) else str(x)))
(OUT/'verified_execution_records.json').write_text(json.dumps(execution,indent=2))
summary={'checks':len(checks),'divergences':len(issues),'archive_hashes':hashes,'native_poses':len(native),'primary_candidate_poses':sum(r['primary_poses'] for r in pose_counts),'candidate_top1':len(candidate),'flagged_top1':sum(r['flags']!=NOFLAG for r in candidate),'compound_target_summaries':len(compound),'seed_stability_entries':len(stability),'executed_commands_verified':len(execution),'new_docking_simulations':0}
(OUT/'data_audit_summary.json').write_text(json.dumps(summary,indent=2));print('AUDIT_COMPLETE',json.dumps(summary),flush=True)
