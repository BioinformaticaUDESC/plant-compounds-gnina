"""Prepare six prospectively selected targets without moving deposited heavy atoms."""
import io,json,random,hashlib
from pathlib import Path
import numpy as np
from rdkit import Chem,RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from pdbfixer import PDBFixer
from openmm import Platform,unit
from openmm.app import Modeller,PDBFile
ROOT=Path(__file__).resolve().parent
SPECS=[
 ('MAOA_2Z5X','2Z5X','HRM','A','harmine_pyridine_protonated',406,'Monoamine oxidase A; deposited monomer, oxidized covalent FAD retained; no membrane or redox model'),
 ('MAOB_2V61','2V61','C18','A','ccd_neutral',397,'Monoamine oxidase B; C18 coumarin-related reference with exact deposited CCD connectivity retained, not reconstructed from common name; chain A of dimer and oxidized covalent FAD; no membrane or redox model'),
 ('SEH_3ANS','3ANS','S38','A','ccd_neutral',None,'Human soluble epoxide hydrolase C-terminal catalytic domain; neutral S38 reference; no full-length dynamics'),
 ('KEAP1_4XMB','4XMB','41P','A','ccd_neutral',None,'Human wild-type Kelch domain; native neutral diamide41P; noncovalent NRF2-binding-pocket hypothesis only, no electrophilic Cys151 mechanism; no missing standard sidechain heavy atoms'),
 ('EGFR_1M17','1M17','AQ4','A','ccd_neutral',None,'Human EGFR kinase domain; erlotinib-bound active-like state, no resistance mutations or cellular ATP competition modeled'),
 ('GYRB_1KZN','1KZN','CBN','A','ccd_neutral',None,'Escherichia coli GyrB 24kDa ATPase-domain pocket; deposited neutral clorobiocin state retained as explicit approximation; not whole-cell antimicrobial efficacy')]
def xyz(l):return np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
def key(m):return Chem.MolToInchiKey(Chem.RemoveHs(Chem.Mol(m)))
def dominant(lines):
 d={}
 for l in lines:
  k=(l[21],l[22:27],l[12:16]);priority=(float(l[54:60]),l[16]==' ',l[16]=='A')
  if k not in d or priority>d[k][0]:d[k]=(priority,l[:16]+' '+l[17:])
 return [v[1] for v in d.values()]
def sdf(m,p):
 w=Chem.SDWriter(str(p));w.write(m);w.close()
def native_from_pdb(lines,lig):
 native=Chem.MolFromPDBBlock('\n'.join(lines+['END'])+'\n',removeHs=True,proximityBonding=True)
 template=Chem.RemoveHs(Chem.SDMolSupplier(str(ROOT/'raw'/f'{lig}_ideal.sdf'),removeHs=False)[0])
 assert native is not None and native.GetNumAtoms()==template.GetNumAtoms()
 native=AllChem.AssignBondOrdersFromTemplate(template,native)
 Chem.AssignStereochemistryFrom3D(native,replaceExistingTags=True)
 assert key(native)==key(template),'CCD-versus-crystal chemistry mismatch '+lig
 return native,template
def independent(m,seed,dst):
 mol=Chem.AddHs(Chem.Mol(m));mol.RemoveAllConformers();params=AllChem.ETKDGv3();params.randomSeed=seed;params.numThreads=1;params.enforceChirality=True
 ids=list(AllChem.EmbedMultipleConfs(mol,numConfs=64,params=params));assert ids
 vals=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=2000)
 valid=[]
 for cid,(code,e) in zip(ids,vals):
  if code!=0 or not np.isfinite(e):continue
  trial=Chem.Mol(mol);trial.RemoveAllConformers();trial.AddConformer(Chem.Conformer(mol.GetConformer(cid)),assignId=True)
  back=Chem.MolFromMolBlock(Chem.MolToMolBlock(trial),removeHs=False)
  if back is not None and key(back)==key(m):valid.append((float(e),int(cid),trial))
 assert valid,'No converged stereo-valid independent input'
 e,cid,trial=min(valid,key=lambda t:t[0]);trial.SetProp('_Name',dst.name+'_native_independent_input');sdf(trial,dst/'redock_input.sdf')
 assert key(Chem.SDMolSupplier(str(dst/'redock_input.sdf'),removeHs=False)[0])==key(m)
 return {'attempted_conformers':len(ids),'converged_stereo_valid':len(valid),'MMFF_energy_kcal_mol':e,'selected_conformer_id':cid,'seed':seed,'method':'ETKDGv3 then MMFF2000; native coordinates removed before embedding; chemical key checked after SDF roundtrip'}
def fad_pdb(lines,protein,cys):
 fl=dominant([l for l in lines if l.startswith('HETATM') and l[17:20]=='FAD' and l[21]=='A'])
 fad,template=native_from_pdb(fl,'FAD');original_key=key(fad)
 phosphate=[]
 for atom in fad.GetAtoms():
  if atom.GetAtomicNum()==8 and atom.GetDegree()==1 and atom.GetBonds()[0].GetBondType()==Chem.BondType.SINGLE and atom.GetNeighbors()[0].GetAtomicNum()==15:
   atom.SetFormalCharge(-1);atom.SetNumExplicitHs(0);atom.SetNoImplicit(True);phosphate.append(atom.GetPDBResidueInfo().GetName().strip())
 assert len(phosphate)==2;Chem.SanitizeMol(fad)
 sg=next(l for l in protein if int(l[22:26])==cys and l[12:16].strip()=='SG');cb=next(l for l in protein if int(l[22:26])==cys and l[12:16].strip()=='CB')
 ci=next(a.GetIdx() for a in fad.GetAtoms() if a.GetPDBResidueInfo().GetName().strip()=='C8M')
 # Temporary cysteine S/CB anchors give linked C8M two H, not an isolated methyl.
 rw=Chem.RWMol(fad);si=rw.AddAtom(Chem.Atom(16));bi=rw.AddAtom(Chem.Atom(6));rw.AddBond(ci,si,Chem.BondType.SINGLE);rw.AddBond(si,bi,Chem.BondType.SINGLE)
 linked=rw.GetMol();conf=linked.GetConformer();conf.SetAtomPosition(si,tuple(xyz(sg)));conf.SetAtomPosition(bi,tuple(xyz(cb)));Chem.SanitizeMol(linked)
 h=Chem.AddHs(linked,addCoords=True,addResidueInfo=True)
 hlines=Chem.MolToPDBBlock(h).splitlines();kept=[]
 for atom,line in zip(h.GetAtoms(),[l for l in hlines if l.startswith(('ATOM  ','HETATM'))]):
  if atom.GetIdx() in [si,bi]:continue
  if atom.GetAtomicNum()==1 and atom.GetNeighbors()[0].GetIdx() in [si,bi]:continue
  line=line[:17]+'FAD'+' '+ 'A'+f'{600:4d}'+line[26:];kept.append(line)
 assert len([l for l in kept if l[76:78].strip()!='H'])==fad.GetNumAtoms()
 return kept,{'CCD_neutral_key':original_key,'modeled_FAD_charge':-2,'deprotonated_phosphate_oxygen_names':phosphate,'covalent_Cys_residue':cys,'covalent_FAD_atom':'C8M','deposited_SG_C8M_distance_A':float(np.linalg.norm(xyz(sg)-fad.GetConformer().GetAtomPosition(ci))),'hydrogen_policy':'Cys CYX removes SG-H; explicit FAD hydrogens from CCD connectivity and temporary cysteine anchors; oxidized flavin; phosphate dianion approximation','electronic_policy':'fixed docking cofactor; no redox, QM, metal, or custom FAD scoring model'}
def prepare(spec,index):
 tid,pdb,lig,chain,state,cys,note=spec;dst=ROOT/'targets'/tid;dst.mkdir(exist_ok=True)
 if (dst/'config.json').exists():print('ALREADY_PREPARED',tid,flush=True);return
 lines=(ROOT/'raw'/f'{pdb}.pdb').read_text().splitlines()
 ll=dominant([l for l in lines if l.startswith('HETATM') and l[17:20]==lig and l[21]==chain and l[76:78].strip()!='H'])
 assert len({l[22:27] for l in ll})==1
 native,template=native_from_pdb(ll,lig);ccdkey=key(native)
 if state=='harmine_pyridine_protonated':
  aa=[a for a in native.GetAtoms() if a.GetAtomicNum()==7 and a.GetIsAromatic() and a.GetTotalNumHs()==0];assert len(aa)==1
  aa[0].SetFormalCharge(1);aa[0].SetNumExplicitHs(1);aa[0].SetNoImplicit(True)
 elif state=='safinamide_secondary_amine_protonated':
  aa=[a for a in native.GetAtoms() if a.GetAtomicNum()==7 and a.GetDegree()==2 and all(b.GetBondType()==Chem.BondType.SINGLE for b in a.GetBonds())];assert len(aa)==1
  aa[0].SetFormalCharge(1);aa[0].SetNumExplicitHs(2);aa[0].SetNoImplicit(True)
 elif state=='amino_acid_zwitterion':
  matches=native.GetSubstructMatches(Chem.MolFromSmarts('[CX3](=O)[OX2H1]'));assert len(matches)==1
  a=native.GetAtomWithIdx(matches[0][2]);a.SetFormalCharge(-1);a.SetNumExplicitHs(0);a.SetNoImplicit(True)
  amines=[a for a in native.GetAtoms() if a.GetAtomicNum()==7 and a.GetDegree()==1 and a.GetNeighbors()[0].GetHybridization()==Chem.HybridizationType.SP3];assert len(amines)==1
  a=amines[0];a.SetFormalCharge(1);a.SetNumExplicitHs(3);a.SetNoImplicit(True)
 elif state=='dicarboxylate_minus2':
  matches=native.GetSubstructMatches(Chem.MolFromSmarts('[CX3](=O)[OX2H1]'));assert len(matches)==2
  for match in matches:
   a=native.GetAtomWithIdx(match[2]);a.SetFormalCharge(-1);a.SetNumExplicitHs(0);a.SetNoImplicit(True)
 Chem.SanitizeMol(native)
 assert key(rdMolStandardize.Uncharger().uncharge(Chem.Mol(native)))==key(rdMolStandardize.Uncharger().uncharge(Chem.Mol(template)))
 native.SetProp('_Name',tid+'_native_reference');native.SetProp('microstate_hypothesis',state);sdf(native,dst/'native_reference.sdf')
 assert key(Chem.SDMolSupplier(str(dst/'native_reference.sdf'))[0])==key(native)
 (dst/'native_crystal_heavy.pdb').write_text('\n'.join(ll+['END'])+'\n')
 embed=independent(native,2026102500+index,dst);nc=native.GetConformer().GetPositions()
 protein=dominant([l for l in lines if l.startswith('ATOM  ') and l[21]==chain and l[76:78].strip()!='H'])
 fixer=PDBFixer(pdbfile=io.StringIO('\n'.join(protein+['END'])+'\n'));fixer.findMissingResidues();fixer.missingResidues={};fixer.findNonstandardResidues();assert not fixer.nonstandardResidues
 fixer.findMissingAtoms();missing=[]
 for r,atoms in fixer.missingAtoms.items():
  x=np.array([fixer.positions[a.index].value_in_unit(unit.angstrom) for a in r.atoms()]);d=float(np.linalg.norm(x[:,None]-nc[None,:],axis=-1).min());missing.append({'residue':r.id,'name':r.name,'atoms':[a.name for a in atoms],'distance_existing_atoms_to_native_A':d})
  assert d>=8,'Incomplete standard sidechain within8A: '+tid+' '+r.id
 terminals=[{'residue':r.id,'atoms':list(atoms)} for r,atoms in fixer.missingTerminals.items()]
 fixer.addMissingAtoms(seed=2026102500+index)
 variants=['CYX' if cys is not None and int(r.id)==cys and r.name=='CYS' else None for r in fixer.topology.residues()]
 random.seed(2026102500+index);mod=Modeller(fixer.topology,fixer.positions);chosen=mod.addHydrogens(pH=7.4,variants=variants,platform=Platform.getPlatformByName('CPU'))
 out=io.StringIO();PDBFile.writeFile(mod.topology,mod.positions,out,keepIds=True);ph=[l for l in out.getvalue().splitlines() if l.startswith(('ATOM  ','TER   '))]
 before={(l[21],l[22:27],l[12:16].strip()):xyz(l) for l in protein};after={(l[21],l[22:27],l[12:16].strip()):xyz(l) for l in ph if l.startswith('ATOM  ')}
 maxdelta=max(float(np.linalg.norm(x-after[k])) for k,x in before.items());assert maxdelta<.002
 cof=[];cofaudit={};con=[]
 if cys is not None:
  cof,cofaudit=fad_pdb(lines,protein,cys)
  serial=max(int(l[6:11]) for l in ph if l.startswith('ATOM  '))+10
  cof=[l[:6]+f'{serial+i:5d}'+l[11:] for i,l in enumerate(cof)]
  sg=next(l for l in ph if l.startswith('ATOM  ') and int(l[22:26])==cys and l[12:16].strip()=='SG')
  cm=next(l for l in cof if l[12:16].strip()=='C8M');a=int(sg[6:11]);b=int(cm[6:11]);con=[f'CONECT{a:5d}{b:5d}',f'CONECT{b:5d}{a:5d}']
  assert not any(l[12:16].strip()=='HG' and int(l[22:26])==cys for l in ph if l.startswith('ATOM  '))
 (dst/'receptor_H.pdb').write_text('\n'.join(ph+cof+con+['END'])+'\n')
 modeled={int(l[22:26]):[] for l in protein}
 for l in protein:modeled[int(l[22:26])].append(xyz(l))
 missing_res=[]
 for l in lines:
  if not l.startswith('REMARK 465') or len(l)<27 or l[19:20]!=chain:continue
  try:rid=int(l[21:26])
  except ValueError:continue
  flank=[v for v in [max((i for i in modeled if i<rid),default=None),min((i for i in modeled if i>rid),default=None)] if v is not None]
  distance=min(float(np.linalg.norm(np.array(modeled[v])[:,None]-nc[None,:],axis=-1).min()) for v in flank)
  missing_res.append({'residue':rid,'residue_author_number_with_insertion_code':l[21:27].strip(),'name':l[15:18].strip(),'nearest_existing_flank_distance_to_native_A':distance})
 other=[l for l in lines if l.startswith('ATOM  ') and l[21]!=chain and l[76:78].strip()!='H']
 otherdistance=float(np.linalg.norm(np.array([xyz(l) for l in other])[:,None]-nc[None,:],axis=-1).min()) if other else None
 center=(nc.max(0)+nc.min(0))/2;size=np.maximum(nc.max(0)-nc.min(0)+12,22)
 audit={'target_id':tid,'pdb':pdb,'chain':chain,'native_CCD':lig,'CCD_key_before_microstate':ccdkey,'native_key':key(native),'native_formal_charge':Chem.GetFormalCharge(native),'native_microstate_hypothesis':state,'native_heavy_atoms':native.GetNumAtoms(),'independent_embedding':embed,'missing_standard_heavy_atoms_added':missing,'terminal_atoms_added':terminals,'missing_residues_not_rebuilt':missing_res,'max_deposited_heavy_atom_displacement_A':maxdelta,'cofactor':cofaudit,'retained_cofactor_heavy_atoms':sum(l[76:78].strip()!='H' for l in cof),'selected_hydrogen_variants':chosen,'all_waters_removed':True,'water_sensitivity_tested':False,'minimum_other_ASU_chain_distance_to_native_A':otherdistance,'symmetry_mates_modeled':False,'sequence_differences':[l for l in lines if l.startswith('SEQADV') and l[16]==chain],'note':note}
 (dst/'preparation_audit.json').write_text(json.dumps(audit,indent=2))
 cfg={'target_id':tid,'pdb':pdb,'chain':chain,'ligand':lig,'receptor':str((dst/'receptor_H.pdb').relative_to(ROOT)),'native_reference':str((dst/'native_reference.sdf').relative_to(ROOT)),'redock_input':str((dst/'redock_input.sdf').relative_to(ROOT)),'center':center.tolist(),'size':size.tolist(),'pH_model':7.4,'preparation_note':note,'near_site_missing_residue_warning':any(r['nearest_existing_flank_distance_to_native_A']<8 for r in missing_res)}
 (dst/'config.json').write_text(json.dumps(cfg,indent=2));print('PREPARED',tid,'nativecharge',audit['native_formal_charge'],'missing',missing,'near_site_missing_loop',cfg['near_site_missing_residue_warning'],flush=True)
def main():
 RDLogger.DisableLog('rdApp.warning')
 for i,s in enumerate(SPECS):prepare(s,i)
if __name__=='__main__':main()
