"""Quarantine inherited stereo mismatch and prepare a verified replacement conformer."""
import json,math
from pathlib import Path
from rdkit import Chem,RDLogger
from rdkit.Chem import AllChem,rdMolAlign
import gnina_core as core
ROOT=Path(__file__).resolve().parent
def roundtrip(m,confId=-1):return Chem.MolFromMolBlock(Chem.MolToMolBlock(m,confId=confId),removeHs=False)
def main():
    RDLogger.DisableLog('rdApp.warning');RDLogger.DisableLog('rdApp.error')
    original=core.read_csv(ROOT/'docking_inputs.csv');ss={s['structure_id']:s for s in core.read_csv(ROOT/'structures.csv')};qc=[]
    for row in original:
        m=core.load(ROOT/row['input_path'])[0];expected=core.chemical_key(Chem.MolFromSmiles(ss[row['structure_id']]['isomeric_smiles']));observed=core.chemical_key(m)
        qc.append({'input_id':row['input_id'],'structure_id':row['structure_id'],'compound':row['compound'],'expected_chemical_key':expected,'observed_3D_input_chemical_key':observed,'status':'PASS' if expected==observed else 'QUARANTINE_STEREOCHEMICAL_IDENTITY_MISMATCH','input_sha256':core.sha(ROOT/row['input_path']),'primary_analysis_eligible':expected==observed})
    bad=[q for q in qc if not q['primary_analysis_eligible']];assert len(bad)==1 and bad[0]['input_id']=='caryophyllene_oxide_conf03',bad
    badrow=next(r for r in original if r['input_id']==bad[0]['input_id']);structure=ss[badrow['structure_id']];wanted=Chem.MolFromSmiles(structure['isomeric_smiles']);expected=core.chemical_key(wanted)
    m=Chem.AddHs(wanted);params=AllChem.ETKDGv3();params.randomSeed=2026102503;params.numThreads=1;params.pruneRmsThresh=.20;params.enforceChirality=True
    ids=list(AllChem.EmbedMultipleConfs(m,numConfs=160,params=params));assert ids
    values=AllChem.MMFFOptimizeMoleculeConfs(m,numThreads=1,maxIters=2000)
    valid=[];rejected=[]
    for ci,(status,energy) in zip(ids,values):
        if status!=0 or not math.isfinite(energy):continue
        rt=roundtrip(m,ci)
        if rt is None or core.chemical_key(rt)!=expected:rejected.append(ci);continue
        valid.append((float(energy),ci))
    assert valid;valid.sort()
    existing=[Chem.RemoveHs(core.load(ROOT/r['input_path'])[0]) for r in original if r['structure_id']==badrow['structure_id'] and r['input_id']!=badrow['input_id']]
    choices=[]
    for energy,ci in valid:
        if energy>valid[0][0]+8:continue
        probe=Chem.RemoveHs(roundtrip(m,ci));distances=[float(rdMolAlign.GetBestRMS(Chem.Mol(probe),ref,maxMatches=100000)) for ref in existing]
        choices.append((min(distances),energy,ci,distances))
    distinct=[x for x in choices if x[0]>.5]
    # If no genuinely distinct third conformer exists, retain two primary conformers;
    # keep a chemically correct replacement only as a separate methodological diagnostic.
    use_in_primary=bool(distinct)
    separation,energy,ci,distances=min(distinct if distinct else choices,key=lambda x:x[1])
    corrected=roundtrip(m,ci);iid='caryophyllene_oxide_conf03_corrected';corrected.SetProp('_Name',iid);corrected.SetProp('structure_id',badrow['structure_id']);corrected.SetProp('sample_presence','not_confirmed');corrected.SetProp('preparation_note','replacement_of_inherited_stereochemical_mismatch; 3D_roundtrip_identity_verified')
    path=ROOT/'prepared/ligands'/(iid+'.sdf');w=Chem.SDWriter(str(path));w.write(corrected);w.close();assert core.chemical_key(core.load(path)[0])==expected
    row=dict(badrow);row.update({'input_id':iid,'input_path':str(path.relative_to(ROOT)),'MMFF_energy_kcal_mol':energy,'delta_MMFF_kcal_mol':energy-valid[0][0],'source_conformer_id':ci,'sha256':core.sha(path)})
    core.write_csv(ROOT/'correction_docking_inputs.csv',[row]);core.write_csv(ROOT/'all_attempted_docking_inputs.csv',original+[row]);primary=[r for r in original if r['input_id']!=badrow['input_id']]
    if use_in_primary:primary.append(row)
    core.write_csv(ROOT/'primary_docking_inputs.csv',primary)
    qc.append({'input_id':iid,'structure_id':badrow['structure_id'],'compound':badrow['compound'],'expected_chemical_key':expected,'observed_3D_input_chemical_key':core.chemical_key(core.load(path)[0]),'status':'PASS_CORRECTED_INPUT' if use_in_primary else 'PASS_CHEMICAL_IDENTITY_DIAGNOSTIC_NOT_DISTINCT','input_sha256':core.sha(path),'primary_analysis_eligible':use_in_primary});core.write_csv(ROOT/'reports/input_identity_qc.csv',qc)
    w=Chem.SDWriter(str(ROOT/'prepared/primary_inputs.sdf'))
    for r in primary:w.write(core.load(ROOT/r['input_path'])[0])
    w.close()
    report={'original_input_count':len(original),'quarantined_inputs':bad,'corrected_input_id':iid,'corrected_input_in_primary_analysis':use_in_primary,'final_primary_input_count':len(primary),'final_primary_structure_count':len(set(r['structure_id'] for r in primary)),'corrected_chemical_key':expected,'MMFF_converged_stereo_valid_candidates':len(valid),'MMFF_converged_stereo_invalid_candidates_rejected':len(rejected),'corrected_MMFF_energy_kcal_mol':energy,'aligned_heavy_atom_RMSDs_to_retained_conformers_A':distances,'selection':'lowest_MMFF_energy_verified_conformer_within_8kcal_of_verified_min; primary inclusion additionally requires RMSD above 0.5A to both retained inputs; otherwise diagnostic only','interpretation':'Post-screen identity QC identified an inherited mismatch. Curation used chemical identity and conformer-diversity criteria independently of docking-score favorability; original outputs retained, mismatched input excluded from primary analysis and chemically correct but redundant replacement kept as a separate diagnostic.'}
    (ROOT/'reports/input_curation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
