"""Prepare deposited receptor frames and independent native-ligand inputs."""
import io,json,hashlib,random
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pdbfixer import PDBFixer
from openmm import Platform,unit
from openmm.app import Modeller,PDBFile

ROOT=Path(__file__).resolve().parent
SPECS=[('ACHE_4EY7','4EY7','E20','B',[], 'donepezil tertiary amine protonated +1; reference stereochemistry from deposited coordinates'),('ALOX5_6N2W','6N2W','30Z','B',['FE2'],'NDGA neutral deposited reference; engineered stable construct and disordered active-site capping region'),('CB2_5ZTY','5ZTY','9JU','A',[],'AM10257 neutral antagonist; engineered inactive CB2-T4 lysozyme fusion; membrane excluded')]
def xyz(l):return np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
def write_sdf(m,p):
    w=Chem.SDWriter(str(p));w.write(m);w.close()
def dominant(lines):
    chosen={}
    for l in lines:
        key=(l[21],l[22:27],l[12:16])
        priority=(float(l[54:60]),l[16]==' ',l[16]=='A')
        if key not in chosen or priority>chosen[key][0]:chosen[key]=(priority,l[:16]+' '+l[17:])
    return [v[1] for v in chosen.values()]
def independent_embed(m,seed):
    mol=Chem.AddHs(Chem.Mol(m));mol.RemoveAllConformers();p=AllChem.ETKDGv3();p.randomSeed=seed;p.numThreads=1
    ids=list(AllChem.EmbedMultipleConfs(mol,numConfs=48,params=p));assert ids
    vals=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=2000)
    valid=sorted((e,i) for i,(code,e) in zip(ids,vals) if code==0 and np.isfinite(e));assert valid
    conf=Chem.Conformer(mol.GetConformer(valid[0][1]));mol.RemoveAllConformers();mol.AddConformer(conf,assignId=True)
    return mol
def main():
    random.seed(20261025)
    for index,(tid,pdb,lig,chain,keep,note) in enumerate(SPECS):
        dst=ROOT/'targets'/tid;dst.mkdir(exist_ok=True)
        if (dst/'config.json').exists():raise RuntimeError('Refuse overwrite '+tid)
        lines=(ROOT/'raw'/(pdb+'.pdb')).read_text().splitlines()
        ligandlines=dominant([l for l in lines if l.startswith('HETATM') and l[17:20]==lig and l[21]==chain and l[76:78].strip()!='H'])
        (dst/'native_crystal_heavy.pdb').write_text('\n'.join(ligandlines+['END'])+'\n')
        native=Chem.MolFromPDBBlock('\n'.join(ligandlines+['END'])+'\n',removeHs=True,proximityBonding=True)
        template=Chem.RemoveHs(Chem.SDMolSupplier(str(ROOT/'raw'/(lig+'_ideal.sdf')),removeHs=False)[0]);native=AllChem.AssignBondOrdersFromTemplate(template,native)
        if tid=='ACHE_4EY7':
            ns=[a for a in native.GetAtoms() if a.GetAtomicNum()==7];assert len(ns)==1 and ns[0].GetDegree()==3
            ns[0].SetFormalCharge(1);ns[0].SetNumExplicitHs(1);ns[0].SetNoImplicit(True);Chem.SanitizeMol(native)
        Chem.AssignStereochemistryFrom3D(native,replaceExistingTags=True)
        native.SetProp('_Name',tid+'_native_reference');native.SetProp('microstate_note',note);write_sdf(native,dst/'native_reference.sdf')
        redock=independent_embed(native,923100+index);redock.SetProp('_Name',tid+'_native_independent_input');write_sdf(redock,dst/'redock_input.sdf')
        nativecoords=native.GetConformer().GetPositions()
        protein=dominant([l for l in lines if l.startswith('ATOM  ') and l[21]==chain and l[76:78].strip()!='H'])
        fixer=PDBFixer(pdbfile=io.StringIO('\n'.join(protein+['END'])+'\n'));fixer.findMissingResidues();fixer.missingResidues={};fixer.findNonstandardResidues();assert not fixer.nonstandardResidues
        fixer.findMissingAtoms();missing=[]
        for residue,atoms in fixer.missingAtoms.items():
            x=np.array([fixer.positions[a.index].value_in_unit(unit.angstrom) for a in residue.atoms()]);d=float(np.linalg.norm(x[:,None,:]-nativecoords[None,:,:],axis=-1).min())
            missing.append({'residue':residue.id,'name':residue.name,'atoms':[a.name for a in atoms],'minimum_existing_residue_distance_to_native_A':d})
            if d<8:raise RuntimeError('Missing standard heavy atoms near reference; review '+tid)
        terminals=[{'residue':r.id,'name':r.name,'atoms':list(atoms)} for r,atoms in fixer.missingTerminals.items()]
        fixer.addMissingAtoms(seed=20261025+index)
        cof=dominant([l for l in lines if l.startswith('HETATM') and l[21]==chain and l[17:20] in keep])
        metals=np.array([xyz(l) for l in cof if l[76:78].strip().upper() in ['FE','CO','ZN']])
        variants=[];metal_his=[]
        for r in fixer.topology.residues():
            variant=None
            if r.name=='HIS' and len(metals):
                for a in r.atoms():
                    if a.name not in ['ND1','NE2']:continue
                    x=np.array(fixer.positions[a.index].value_in_unit(unit.angstrom))
                    if np.linalg.norm(metals-x,axis=1).min()<3:
                        variant='HID' if a.name=='NE2' else 'HIE';metal_his.append({'residue':r.id,'coordinating_atom':a.name,'variant':variant})
            variants.append(variant)
        modeller=Modeller(fixer.topology,fixer.positions);selected=modeller.addHydrogens(pH=7.4,variants=variants,platform=Platform.getPlatformByName('CPU'))
        s=io.StringIO();PDBFile.writeFile(modeller.topology,modeller.positions,s,keepIds=True)
        ph=[l for l in s.getvalue().splitlines() if l.startswith(('ATOM  ','TER   '))]
        before={(l[21],l[22:27],l[12:16].strip()):xyz(l) for l in protein};after={(l[21],l[22:27],l[12:16].strip()):xyz(l) for l in ph if l.startswith('ATOM  ')}
        deviations=[float(np.linalg.norm(x-after[k])) for k,x in before.items()];assert max(deviations)<.002
        nextserial=max(int(l[6:11]) for l in ph if l.startswith('ATOM  '))+10
        cof=[l[:6]+f'{nextserial+i:5d}'+l[11:] for i,l in enumerate(cof)]
        receptor=dst/'receptor_H.pdb';receptor.write_text('\n'.join(ph+cof+['END'])+'\n')
        missing_res=[]
        modeled={int(l[22:26]):[] for l in protein}
        for l in protein:modeled[int(l[22:26])].append(xyz(l))
        for l in lines:
            if not l.startswith('REMARK 465') or len(l)<27 or l[19:20]!=chain:continue
            try:rid=int(l[21:26])
            except ValueError:continue
            flanks=[v for v in [max((i for i in modeled if i<rid),default=None),min((i for i in modeled if i>rid),default=None)] if v is not None]
            nearest=min(float(np.linalg.norm(np.array(modeled[v])[:,None,:]-nativecoords[None,:,:],axis=-1).min()) for v in flanks)
            missing_res.append({'residue':rid,'name':l[15:18].strip(),'nearest_flanking_residue_distance_to_native_A':nearest})
        center=(nativecoords.max(axis=0)+nativecoords.min(axis=0))/2;size=np.maximum(nativecoords.max(axis=0)-nativecoords.min(axis=0)+12,22)
        audit={'pdb':pdb,'chain':chain,'native_ligand':lig,'native_heavy_atoms':native.GetNumAtoms(),'native_smiles':Chem.MolToSmiles(native,isomericSmiles=True),'formal_charge':Chem.GetFormalCharge(native),'note':note,'missing_standard_heavy_atoms_added':missing,'terminal_atoms_added':terminals,'missing_residues_not_rebuilt_from_REMARK465':missing_res,'max_deposited_heavy_atom_displacement_A':max(deviations),'retained_cofactor_atoms':len(cof),'metal_coordinating_histidines':metal_his,'selected_hydrogen_variants':selected,'protein_residues_modeled':sum(1 for _ in modeller.topology.residues()),'excluded_heterogens':sorted(set(l[17:20].strip() for l in lines if l.startswith('HETATM') and l[17:20] not in keep)),'water_policy':'all removed; sensitivity not tested','construct_policy':'deposited modeled chain retained including fusion if present; no missing-loop modeling or heavy-atom relaxation','cofactor_policy':'deposited metal retained with standard docking atom types; no custom coordination/redox/electronic model'}
        (dst/'preparation_audit.json').write_text(json.dumps(audit,indent=2))
        cfg={'target_id':tid,'pdb':pdb,'chain':chain,'ligand':lig,'receptor':str(receptor.relative_to(ROOT)),'native_reference':str((dst/'native_reference.sdf').relative_to(ROOT)),'redock_input':str((dst/'redock_input.sdf').relative_to(ROOT)),'center':center.tolist(),'size':size.tolist(),'pH_model':7.4,'preparation_note':note,'near_site_missing_residue_warning':any(x['nearest_flanking_residue_distance_to_native_A']<8 for x in missing_res)}
        (dst/'config.json').write_text(json.dumps(cfg,indent=2));print(json.dumps({'target':tid,'native_atoms':native.GetNumAtoms(),'missing_heavy':missing,'near_site_missing_residues':cfg['near_site_missing_residue_warning'],'center':cfg['center'],'size':cfg['size']},indent=2),flush=True)

if __name__=='__main__':main()
