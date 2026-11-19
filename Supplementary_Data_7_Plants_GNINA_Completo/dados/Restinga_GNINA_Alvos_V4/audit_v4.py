"""Full pose QC, seed variation and residue proximities for completed target batches."""
import argparse,csv,json,statistics,itertools
from pathlib import Path
from collections import defaultdict
import numpy as np
from rdkit import Chem,RDLogger
import gnina_core as core
ROOT=Path(__file__).resolve().parent
NOFLAG='no_flag_from_limited_geometry_audit'
METALS={'Fe','Co','Zn','Mg','Ca','Mn','Ni','Cu'}
def graph_order_signature(m):
    return (tuple((a.GetAtomicNum(),a.GetFormalCharge()) for a in m.GetAtoms()),tuple((b.GetBeginAtomIdx(),b.GetEndAtomIdx(),str(b.GetBondType())) for b in m.GetBonds()))
def main():
    RDLogger.DisableLog('rdApp.warning');RDLogger.DisableLog('rdApp.error')
    ap=argparse.ArgumentParser();ap.add_argument('--targets',nargs='+');ap.add_argument('--reuse-audited',action='store_true');args=ap.parse_args()
    inputs=core.read_csv(ROOT/'primary_docking_inputs.csv');attempted=core.read_csv(ROOT/'all_attempted_docking_inputs.csv');initial_ids={r['input_id'] for r in core.read_csv(ROOT/'docking_inputs.csv')};correction_ids={r['input_id'] for r in core.read_csv(ROOT/'correction_docking_inputs.csv')}
    meta={r['input_id']:r for r in attempted};inm={r['input_id']:Chem.RemoveHs(core.load(ROOT/r['input_path'])[0]) for r in attempted};identity_qc={r['input_id']:r for r in core.read_csv(ROOT/'reports/input_identity_qc.csv')}
    primary_ids={r['input_id'] for r in inputs};assert all(identity_qc[i]['primary_analysis_eligible']=='True' for i in primary_ids)
    pt=Chem.GetPeriodicTable();all_top=[];all_struct=[];all_valid=[];contact_summary=[]
    for cfgpath in sorted((ROOT/'targets').glob('*/config.json')):
        cfg=json.loads(cfgpath.read_text());tid=cfg['target_id']
        if args.targets and tid not in args.targets:continue
        dst=ROOT/'results'/tid;validation=json.loads((dst/'validation.json').read_text());all_valid.append(validation)
        if not validation['screening_eligible']:continue
        marker=dst/'screening_complete.json'
        if not marker.exists():raise RuntimeError('Incomplete target '+tid)
        if correction_ids:assert (dst/'correction_complete.json').exists(),'Corrected input must finish for '+tid
        audit_sources={'core_script_sha256':core.sha(Path(core.__file__)),'audit_script_sha256':core.sha(Path(__file__)),'protocol_sha256':core.sha(ROOT/cfg.get('protocol_file','protocol.json')),'library_sha256':core.sha(ROOT/'primary_docking_inputs.csv'),'input_identity_qc_sha256':core.sha(ROOT/'reports/input_identity_qc.csv'),'receptor_sha256':core.sha(ROOT/cfg['receptor']),'validation_sha256':core.sha(dst/'validation.json'),'raw_pose_sha256':{p.name:core.sha(p) for p in sorted((dst/'docking').glob('*.sdf'))}}
        cached=dst/'audit_complete.json'
        if args.reuse_audited and cached.exists():
            previous=json.loads(cached.read_text());assert previous['source_fingerprint']==audit_sources,'Stale audit: rerun without --reuse-audited'
            for name,digest in previous['audited_output_sha256'].items():assert core.sha(dst/name)==digest,'Changed audit output '+name
            all_top.extend(core.read_csv(dst/'top1_by_input_seed.csv'));all_struct.extend(core.read_csv(dst/'structure_summary.csv'));contact_summary.extend(core.read_csv(dst/'residue_contact_recurrence.csv'));print('AUDIT_REUSED '+tid,flush=True);continue
        print('AUDIT_START '+tid,flush=True)
        rec=[];radii=[];elements=[];labels=[]
        for l in (ROOT/cfg['receptor']).read_text().splitlines():
            if not l.startswith(('ATOM  ','HETATM')):continue
            el=l[76:78].strip().title()
            if el=='H':continue
            anum=pt.GetAtomicNumber(el);assert anum>0
            rec.append([float(l[30:38]),float(l[38:46]),float(l[46:54])]);radii.append(pt.GetRvdw(anum));elements.append(el);labels.append((l[21],l[17:20].strip(),l[22:27].strip(),l[12:16].strip()))
        rec=np.array(rec);rr=np.array(radii);is_metal=np.array([e in METALS for e in elements]);center=np.array(cfg['center']);half=np.array(cfg['size'])/2
        # All atoms relevant to the box plus a 6 Å buffer; nearest distances still use full receptor.
        pose_rows=core.read_csv(dst/'docking_poses.csv')+(core.read_csv(dst/'correction_docking_poses.csv') if correction_ids else []);rows_by={(int(r['seed']),r['input_id'],int(r['energy_rank'])):r for r in pose_rows}
        audit=[];tops=[];all_attempted_tops=[];proximities=[];topm={};map_cache={}
        for sdf in sorted((dst/'docking').glob('*.sdf')):
            seed=int(sdf.stem.rsplit('seed',1)[1]);groups=defaultdict(list)
            for raw in core.load(sdf):groups[raw.GetProp('_Name')].append(raw)
            assert set(groups)==(initial_ids if '_batch_seed' in sdf.name else correction_ids)
            for iid,mm in groups.items():
                inp=inm[iid];ix=inp.GetConformer().GetPositions();inpkey=core.chemical_key(inp)
                for rank,raw in enumerate(mm,1):
                    m=Chem.RemoveHs(raw);assert core.chemical_key(m)==inpkey,'Chemical identity differs '+iid
                    assert core.microstate_graph_key(m)==core.microstate_graph_key(inp),'Charge/bond/H-state differs '+iid
                    x=m.GetConformer().GetPositions();assert np.isfinite(x).all()
                    d=np.linalg.norm(x[:,None,:]-rec[None,:,:],axis=-1);lr=np.array([pt.GetRvdw(a.GetAtomicNum()) for a in m.GetAtoms()]);ratio=d/(lr[:,None]+rr[None,:])
                    # Metal coordination is recorded separately, not classified by organic vdW radii.
                    close=int((d[:,~is_metal]<2).sum());overlap=int((ratio[:,~is_metal]<.70).sum());outside=float(np.maximum(np.abs(x-center)-half,0).max())
                    sig=(iid,graph_order_signature(m))
                    if sig not in map_cache:
                        mappings=inp.GetSubstructMatches(m,useChirality=False,uniquify=False,maxMatches=10000);assert mappings,'No atom mapping'
                        map_cache[sig]=mappings
                    delta=min(max(abs(np.linalg.norm(x[b.GetBeginAtomIdx()]-x[b.GetEndAtomIdx()])-np.linalg.norm(ix[mp[b.GetBeginAtomIdx()]]-ix[mp[b.GetEndAtomIdx()]])) for b in m.GetBonds()) for mp in map_cache[sig])
                    # Ignore graph-near atoms when checking severe intramolecular overlap.
                    selfd=np.linalg.norm(x[:,None,:]-x[None,:,:],axis=-1);topdist=Chem.GetDistanceMatrix(m);mask=np.triu(topdist>3,k=1);selfratio=selfd/(lr[:,None]+lr[None,:]);intra=int(((selfratio<.70)&mask).sum())
                    row=rows_by[(seed,iid,rank)];energy=float(row['minimizedAffinity_kcal_mol']);flags=[];eligible=iid in primary_ids
                    quarantined=identity_qc[iid]['status'].startswith('QUARANTINE')
                    if quarantined:flags.append('QUARANTINED_INPUT_STEREOCHEMICAL_IDENTITY_MISMATCH')
                    if close:flags.append('nonmetal_heavy_atom_contacts_below_2A')
                    if overlap:flags.append('nonmetal_vdw_ratio_below_0.70')
                    if energy>0:flags.append('positive_empirical_energy')
                    if outside>.1:flags.append('outside_box_above_0.1A')
                    if delta>.15:flags.append('bond_length_change_above_0.15A')
                    if intra:flags.append('severe_intraligand_nonbonded_overlap')
                    q=dict(row);q.update({'primary_analysis_eligible':eligible,'input_identity_status':identity_qc[iid]['status'],'primary_exclusion_reason':('stereochemical_identity_mismatch' if quarantined else ('corrected_but_not_distinct_diagnostic_conformer' if not eligible else '')),'minimum_receptor_heavy_atom_distance_A':float(d.min()),'minimum_nonmetal_heavy_atom_distance_A':float(d[:,~is_metal].min()),'minimum_metal_distance_A':float(d[:,is_metal].min()) if is_metal.any() else '', 'nonmetal_heavy_pairs_below_2A':close,'nonmetal_vdw_pairs_ratio_below_0.70':overlap,'maximum_outside_box_A':outside,'maximum_mapped_bond_length_change_A':float(delta),'severe_intraligand_nonbonded_overlaps':intra,'quality_flags':';'.join(flags) or NOFLAG});audit.append(q)
                    if rank==1:
                        all_attempted_tops.append(q)
                        if not eligible:continue
                        tops.append(q);topm[(iid,seed)]=raw;resdist={}
                        for ai,ri in np.argwhere(d<=4):
                            key=labels[ri][:3]
                            if key not in resdist or d[ai,ri]<resdist[key][0]:resdist[key]=(float(d[ai,ri]),int(ai),labels[ri][3])
                        for key,(distance,ai,rname) in resdist.items():proximities.append({'target_id':tid,'input_id':iid,'structure_id':meta[iid]['structure_id'],'compound':meta[iid]['compound'],'seed':seed,'receptor_chain':key[0],'residue_name':key[1],'residue_id':key[2],'minimum_distance_A':distance,'ligand_atom_index_0based':ai,'ligand_atom_element':m.GetAtomWithIdx(ai).GetSymbol(),'receptor_atom':rname,'interpretation':'geometric_proximity_le_4A_only; no_Hbond_mechanism_or_affinity_assignment'})
            print(f'AUDIT_SEED_COMPLETE {tid} {seed} poses={len(audit)}',flush=True)
        core.write_csv(dst/'pose_quality_audit.csv',audit);core.write_csv(dst/'primary_docking_poses.csv',[r for r in audit if r['primary_analysis_eligible']]);core.write_csv(dst/'top1_all_attempted_inputs.csv',all_attempted_tops);core.write_csv(dst/'top1_by_input_seed.csv',tops);core.write_csv(dst/'top1_residue_proximities.csv',proximities)
        expected_seeds={42,142,242};grouped=defaultdict(list)
        for row in tops:grouped[row['input_id']].append(row)
        seedsummary=[]
        for iid,rrr in grouped.items():
            assert {int(r['seed']) for r in rrr}==expected_seeds
            ee=[float(r['minimizedAffinity_kcal_mol']) for r in rrr];cc=[float(r['CNNscore_fast']) for r in rrr];pair=[core.rmsd(topm[(iid,a)],topm[(iid,b)]) for a,b in itertools.combinations(sorted(expected_seeds),2)]
            seedsummary.append({'target_id':tid,'input_id':iid,'structure_id':meta[iid]['structure_id'],'compound':meta[iid]['compound'],'panel':meta[iid]['panel'],'seeds':3,'top1_energy_median_kcal_mol':statistics.median(ee),'top1_energy_min_kcal_mol':min(ee),'top1_energy_max_kcal_mol':max(ee),'top1_energy_range_kcal_mol':max(ee)-min(ee),'top1_CNNscore_fast_median':statistics.median(cc),'seeds_without_limited_geometry_flags':sum(r['quality_flags']==NOFLAG for r in rrr),'maximum_top1_pairwise_fixed_frame_RMSD_between_seeds_A':max(pair),'median_top1_pairwise_fixed_frame_RMSD_between_seeds_A':statistics.median(pair),'native_redocking_geometry_status':validation['geometry_status'],'sample_presence':'not_confirmed','interpretation':'repeated_search_variation_not_measurement_error_or_experimental_affinity'})
        core.write_csv(dst/'seed_stability_by_input.csv',seedsummary)
        sb=defaultdict(list)
        for r in tops:sb[r['structure_id']].append(r)
        structsummary=[]
        for sid,rrr in sb.items():
            ee=[float(r['minimizedAffinity_kcal_mol']) for r in rrr]
            structsummary.append({'target_id':tid,'structure_id':sid,'compound':rrr[0]['compound'],'panel':rrr[0]['panel'],'input_conformers':len(set(r['input_id'] for r in rrr)),'top1_observations':len(rrr),'top1_energy_median_kcal_mol':statistics.median(ee),'top1_energy_min_kcal_mol':min(ee),'top1_energy_max_kcal_mol':max(ee),'observations_without_limited_geometry_flags':sum(r['quality_flags']==NOFLAG for r in rrr),'sample_presence':'not_confirmed','note':'conformer_and_seed_sensitivity; no_target_or_species_potency_comparison'})
        core.write_csv(dst/'structure_summary.csv',structsummary)
        # Residue recurrence is calculated separately for each unchanged input conformation.
        cg=defaultdict(list)
        for r in proximities:cg[(r['input_id'],r['receptor_chain'],r['residue_name'],r['residue_id'])].append(r)
        target_contacts=[]
        for key,rrr in cg.items():target_contacts.append({'target_id':tid,'input_id':key[0],'structure_id':meta[key[0]]['structure_id'],'compound':meta[key[0]]['compound'],'receptor_chain':key[1],'residue_name':key[2],'residue_id':key[3],'seeds_with_proximity_le_4A':len(set(r['seed'] for r in rrr)),'total_seeds':3,'minimum_distance_A':min(r['minimum_distance_A'] for r in rrr),'interpretation':'contact_recurrence_not_binding_probability'})
        core.write_csv(dst/'residue_contact_recurrence.csv',target_contacts);contact_summary.extend(target_contacts)
        counts={'target_id':tid,'input_conformers':len(inputs),'attempted_inputs_including_quarantine':len(attempted),'stereochemical_structures':len(set(r['structure_id'] for r in inputs)),'compounds':len(set(r['compound'] for r in inputs)),'screening_seeds':3,'docking_poses':sum(r['primary_analysis_eligible'] for r in audit),'all_exported_docking_poses':len(audit),'quarantined_docking_poses':sum(r['input_identity_status'].startswith('QUARANTINE') for r in audit),'diagnostic_correction_poses_outside_primary':sum(not r['primary_analysis_eligible'] and not r['input_identity_status'].startswith('QUARANTINE') for r in audit),'top1_observations':len(tops),'top1_with_flags':sum(r['quality_flags']!=NOFLAG for r in tops),'top1_without_flags':sum(r['quality_flags']==NOFLAG for r in tops),'input_conformers_all_3_seeds_without_flags':sum(r['seeds_without_limited_geometry_flags']==3 for r in seedsummary),'input_conformers_pairwise_RMSD_le_2A':sum(r['maximum_top1_pairwise_fixed_frame_RMSD_between_seeds_A']<=2 for r in seedsummary),'audit_limits':'Operational geometry checks; absence of flags does not establish pose validity or biological activity'}
        (dst/'audit_counts.json').write_text(json.dumps(counts,indent=2));all_top.extend(tops);all_struct.extend(structsummary)
        output_names=['pose_quality_audit.csv','primary_docking_poses.csv','top1_all_attempted_inputs.csv','top1_by_input_seed.csv','top1_residue_proximities.csv','seed_stability_by_input.csv','structure_summary.csv','residue_contact_recurrence.csv','audit_counts.json']
        cached.write_text(json.dumps({'source_fingerprint':audit_sources,'audited_output_sha256':{name:core.sha(dst/name) for name in output_names}},indent=2));print(json.dumps(counts),flush=True)
    if not args.targets:
        core.write_csv(ROOT/'reports/all_top1_by_input_seed.csv',all_top);core.write_csv(ROOT/'reports/all_structure_summaries.csv',all_struct);core.write_csv(ROOT/'reports/residue_contact_recurrence.csv',contact_summary)
        (ROOT/'reports/redocking_summary.json').write_text(json.dumps(all_valid,indent=2))
        acs=[json.loads(p.read_text()) for p in (ROOT/'results').glob('*/audit_counts.json')]
        (ROOT/'reports/study_counts.json').write_text(json.dumps({'compounds':len(set(r['compound'] for r in inputs)),'structures':len(set(r['structure_id'] for r in inputs)),'input_conformers':len(inputs),'attempted_inputs_including_quarantine':len(attempted),'targets_tested':len(all_valid),'targets_geometry_pass':sum(v['geometry_status']=='PASS' for v in all_valid),'targets_geometry_fail':sum(v['geometry_status']=='FAIL' for v in all_valid),'docking_top1_observations':len(all_top),'docking_poses':sum(a['docking_poses'] for a in acs),'all_exported_docking_poses':sum(a['all_exported_docking_poses'] for a in acs),'quarantined_docking_poses':sum(a['quarantined_docking_poses'] for a in acs),'diagnostic_correction_poses_outside_primary':sum(a['diagnostic_correction_poses_outside_primary'] for a in acs),'native_redocking_poses':sum(len(core.read_csv(p)) for p in (ROOT/'results').glob('*/redocking/poses.csv'))},indent=2))

if __name__=='__main__':main()
