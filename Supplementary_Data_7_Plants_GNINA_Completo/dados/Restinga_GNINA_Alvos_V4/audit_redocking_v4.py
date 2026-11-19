"""Secondary geometry audit of every native control; does not redefine the RMSD gate."""
from pathlib import Path
import json,csv
import numpy as np
from rdkit import Chem,RDLogger
import gnina_core as core
ROOT=Path(__file__).resolve().parent
METALS={'Fe','Co','Zn','Mg','Ca','Mn','Ni','Cu','Na','K'}
def main():
 RDLogger.DisableLog('rdApp.warning');RDLogger.DisableLog('rdApp.error');pt=Chem.GetPeriodicTable();allrows=[];summaries=[];checks=[]
 for cp in sorted((ROOT/'targets').glob('*/config.json')):
  cfg=json.loads(cp.read_text());tid=cfg['target_id'];dst=ROOT/'results'/tid;assert (dst/'validation.json').exists()
  ref=Chem.RemoveHs(core.load(ROOT/cfg['native_reference'])[0]);inp=Chem.RemoveHs(core.load(ROOT/cfg['redock_input'])[0]);expected=core.chemical_key(ref);assert expected==core.chemical_key(inp);assert core.microstate_graph_key(ref)==core.microstate_graph_key(inp)
  shifted=Chem.Mol(ref)
  for j,p in enumerate(shifted.GetConformer().GetPositions()):shifted.GetConformer().SetAtomPosition(j,(float(p[0])+10,float(p[1]),float(p[2])))
  selfr=core.rmsd(ref,ref);tr=core.rmsd(ref,shifted);assert abs(selfr)<1e-6 and abs(tr-10)<1e-5
  # Reorder all atoms: the same coordinates must still yield symmetry-aware0A.
  reordered=Chem.RenumberAtoms(ref,list(reversed(range(ref.GetNumAtoms()))));rr=core.rmsd(ref,reordered);assert rr<1e-6
  checks.append({'target_id':tid,'reference_self_RMSD_A':selfr,'reference_translated_10A_RMSD_A':tr,'atom_reordered_reference_RMSD_A':rr,'native_identity_key':expected,'native_charge_bond_H_graph_key':core.microstate_graph_key(ref),'independent_input_same_charge_bond_H_graph':core.microstate_graph_key(ref)==core.microstate_graph_key(inp),'status':'PASS'})
  rx=[];rad=[];metal=[]
  for l in (ROOT/cfg['receptor']).read_text().splitlines():
   if not l.startswith(('ATOM  ','HETATM')):continue
   el=l[76:78].strip().title()
   if el=='H':continue
   n=pt.GetAtomicNumber(el);rx.append([float(l[30:38]),float(l[38:46]),float(l[46:54])]);rad.append(pt.GetRvdw(n));metal.append(el in METALS)
  rx=np.array(rx);rad=np.array(rad);metal=np.array(metal);centre=np.array(cfg['center']);half=np.array(cfg['size'])/2
  target=[]
  for sp in sorted((dst/'redocking').glob('*.sdf')):
   seed=int(sp.stem.rsplit('seed',1)[1])
   for rank,raw in enumerate(core.load(sp),1):
    m=Chem.RemoveHs(raw);assert core.chemical_key(m)==expected;xx=m.GetConformer().GetPositions();assert np.isfinite(xx).all();lr=np.array([pt.GetRvdw(a.GetAtomicNum()) for a in m.GetAtoms()]);d=np.linalg.norm(xx[:,None]-rx[None,:],axis=-1);vd=d/(lr[:,None]+rad[None,:]);outside=float(np.maximum(np.abs(xx-centre)-half,0).max())
    selfd=np.linalg.norm(xx[:,None]-xx[None,:],axis=-1);td=Chem.GetDistanceMatrix(m);intra=int((((selfd/(lr[:,None]+lr[None,:]))<.70)&np.triu(td>3,k=1)).sum())
    flags=[];close=int((d[:,~metal]<2).sum());overlap=int((vd[:,~metal]<.7).sum())
    if close:flags.append('nonmetal_contacts_below2A')
    if overlap:flags.append('nonmetal_vdw_ratio_below0.70')
    if outside>.1:flags.append('outside_box_above0.1A')
    if intra:flags.append('severe_intraligand_nonbonded_overlap')
    row={'target_id':tid,'seed':seed,'energy_rank':rank,'RMSD_A':core.rmsd(ref,m),'Energy_kcal_mol':float(raw.GetProp('minimizedAffinity')),'minimum_nonmetal_distance_A':float(d[:,~metal].min()),'close_nonmetal_pairs':close,'vdw_overlap_pairs':overlap,'maximum_outside_box_A':outside,'intraligand_overlaps':intra,'flags':';'.join(flags) or 'no_flag_from_limited_geometry_audit','role':'secondary_geometry_diagnostic; does_not_redefine_primary_native_RMSD_gate'};target.append(row);allrows.append(row)
  core.write_csv(dst/'redocking/native_pose_quality_audit.csv',target);top=[x for x in target if x['energy_rank']==1];summaries.append({'target_id':tid,'exported_control_poses':len(target),'top1_observations':len(top),'top1_with_secondary_geometry_flags':sum(x['flags']!='no_flag_from_limited_geometry_audit' for x in top),'primary_gate_unchanged':True})
 core.write_csv(ROOT/'reports/native_control_geometry_audit.csv',allrows);core.write_csv(ROOT/'reports/native_control_geometry_counts.csv',summaries);(ROOT/'method_checks.json').write_text(json.dumps(checks,indent=2));print(json.dumps(summaries,indent=2),flush=True)
if __name__=='__main__':main()
