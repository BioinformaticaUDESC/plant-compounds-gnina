"""Extend the previous structures without altering their input coordinates."""
import csv,hashlib,json,math
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem,Descriptors,rdMolDescriptors
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers,StereoEnumerationOptions

ROOT=Path(__file__).resolve().parent
NEW={
 'eucalyptol':('L1','Lippia alba','chemotype_specific_essential_oils'),
 'sabinene':('L1','Lippia alba','chemotype_specific_essential_oils'),
 '(E)-beta-ocimene':('L1','Lippia alba','chemotype_specific_essential_oils'),
 'gamma-terpinene':('L1','Lippia alba','chemotype_specific_essential_oils'),
 'alpha-terpineol':('L1','Lippia alba','chemotype_specific_essential_oils'),
 'citronellol':('L1','Lippia alba','chemotype_specific_essential_oils'),
 'beta-elemene':('L1','Lippia alba','chemotype_specific_essential_oils'),
 'gamma-muurolene':('L1','Lippia alba','chemotype_specific_essential_oils'),
 'elemol':('L1','Lippia alba','chemotype_specific_essential_oils'),
 'delta-cadinene':('P2','Psidium guajava','leaf_essential_oil'),
 'uvaol':('E1','Eugenia catharinae','ethanolic_leaf_extract_NONVOLATILE'),
 'erythrodiol':('E1','Eugenia catharinae','ethanolic_leaf_extract_NONVOLATILE'),
 'rotundic acid':('E1','Eugenia catharinae','ethanolic_leaf_extract_NONVOLATILE'),
 'beta-sitosterol':('E1','Eugenia catharinae','ethanolic_leaf_extract_NONVOLATILE'),
 'stigmasterol':('E1','Eugenia catharinae','ethanolic_leaf_extract_NONVOLATILE')}
def read(p):return list(csv.DictReader(p.open()))
def write(p,rows):
    with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    if (ROOT/'structures.csv').exists():raise RuntimeError('Refuse overwrite library')
    ss=read(ROOT/'baseline_structures.csv');inputs=read(ROOT/'baseline_docking_inputs.csv');links=read(ROOT/'baseline_literature_species_links.csv');sources=json.loads((ROOT/'baseline_literature_sources.json').read_text())
    sources['sources']['P2']={'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC11763297/','title':'Chemical Characterization and Antimicrobial Activity of Essential Oils and Nanoemulsions of Eugenia uniflora and Psidium guajava','matrix':'leaf_essential_oil','species':'Psidium guajava','compounds':['delta-cadinene']}
    for ni,(name,(sid,species,matrix)) in enumerate(NEW.items()):
        p=json.loads((ROOT/'inputs'/(name.replace(' ','_')+'.json')).read_text())['PropertyTable']['Properties'][0]
        mol=Chem.MolFromSmiles(p.get('SMILES') or p['IsomericSMILES']);assert rdMolDescriptors.CalcMolFormula(mol)==p['MolecularFormula']
        # The carboxylic acid in rotundic acid is represented as carboxylate for this pH 7.4 hypothesis.
        if name=='rotundic acid':
            matches=mol.GetSubstructMatches(Chem.MolFromSmarts('[CX3](=O)[OX2H1]'));assert len(matches)==1
            o=mol.GetAtomWithIdx(matches[0][2]);o.SetFormalCharge(-1);o.SetNumExplicitHs(0);o.SetNoImplicit(True);Chem.SanitizeMol(mol)
        isomers=list(EnumerateStereoisomers(mol,options=StereoEnumerationOptions(onlyUnassigned=True,unique=True,tryEmbedding=True,maxIsomers=32)));assert 0<len(isomers)<=8
        panel='nonvolatile_extract_comparator' if sid=='E1' else 'volatile_literature_panel'
        for ii,isomer in enumerate(isomers,1):
            ident=name.replace(' ','_').replace('(','').replace(')','')+(f'_isomer{ii}' if len(isomers)>1 else '')
            m=Chem.AddHs(Chem.Mol(isomer));m.RemoveAllConformers();params=AllChem.ETKDGv3();params.randomSeed=20271025+ni*100+ii;params.numThreads=1;params.pruneRmsThresh=.25
            ids=list(AllChem.EmbedMultipleConfs(m,numConfs=48,params=params));assert ids and AllChem.MMFFHasAllMoleculeParams(m)
            vals=AllChem.MMFFOptimizeMoleculeConfs(m,numThreads=1,maxIters=2000);valid=sorted((float(e),ci) for ci,(c,e) in zip(ids,vals) if c==0 and math.isfinite(e));assert valid
            e,ci=valid[0];conf=Chem.Conformer(m.GetConformer(ci));m.RemoveAllConformers();m.AddConformer(conf,assignId=True);inp=ident+'_conf01';m.SetProp('_Name',inp);m.SetProp('sample_presence','not_confirmed');m.SetProp('structure_id',ident);m.SetProp('panel',panel)
            path=ROOT/'prepared/ligands'/(inp+'.sdf');w=Chem.SDWriter(str(path));w.write(m);w.close()
            row=dict(ss[0]);row.update({'structure_id':ident,'compound':name,'panel':panel,'pubchem_cid':str(p['CID']),'pubchem_reference_inchikey':p['InChIKey'],'structure_inchikey':Chem.MolToInchiKey(isomer),'formula':rdMolDescriptors.CalcMolFormula(isomer),'isomeric_smiles':Chem.MolToSmiles(isomer,isomericSmiles=True),'formal_charge':str(Chem.GetFormalCharge(isomer)),'enumerated_isomers_from_reference':str(len(isomers)),'microstate_policy':'carboxylate_hypothesis_at_pH7.4' if name=='rotundic acid' else 'PubChem_microstate_retained; no_pKa_enumeration','source_ids':sid,'embedded_conformers':str(len(ids)),'MMFF_converged':str(len(valid)),'selected_input_conformers':'1'});ss.append(row)
            inputs.append({'input_id':inp,'structure_id':ident,'compound':name,'panel':panel,'input_path':str(path.relative_to(ROOT)),'conformer_rank':1,'MMFF_energy_kcal_mol':e,'delta_MMFF_kcal_mol':0,'macrocycle':False,'source_conformer_id':ci,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
            links.append({'structure_id':ident,'compound':name,'literature_species':species,'reported_matrix':matrix,'source_id':sid,'source_url':sources['sources'][sid]['url'],'association':'reported_in_literature_not_confirmed_in_project_samples','stereochemical_assignment':'database_reference_or_enumeration_not_sample_evidence'})
        if name not in sources['sources'][sid]['compounds']:sources['sources'][sid]['compounds'].append(name)
    # Add independently documented shared occurrences without repeating chemical structures/docking inputs.
    shared={'L1':['beta-caryophyllene','limonene','myrcene','geraniol','alpha-pinene','p-cymene'],'P2':['eucalyptol','alpha-terpineol','alpha-pinene','p-cymene','limonene','beta-caryophyllene','beta-selinene','caryophyllene oxide']}
    for sid,compounds in shared.items():
        src=sources['sources'][sid]
        for name in compounds:
            if name not in src['compounds']:src['compounds'].append(name)
            for s in ss:
                if s['compound']!=name:continue
                if any(l['structure_id']==s['structure_id'] and l['source_id']==sid for l in links):continue
                links.append({'structure_id':s['structure_id'],'compound':name,'literature_species':src['species'],'reported_matrix':src['matrix'],'source_id':sid,'source_url':src['url'],'association':'reported_in_literature_not_confirmed_in_project_samples','stereochemical_assignment':'database_reference_or_enumeration_not_sample_evidence'})
                s['source_ids']=';'.join(sorted(set(s['source_ids'].split(';'))|{sid}))
    write(ROOT/'structures.csv',ss);write(ROOT/'docking_inputs.csv',inputs);write(ROOT/'literature_species_links.csv',links)
    cov=[]
    for species in ['Varronia curassavica','Campomanesia littoralis','Myrcia palustris','Eugenia catharinae','Schinus terebinthifolia','Psidium guajava','Lippia alba']:
        ll=[l for l in links if l['literature_species']==species];names=sorted(set(l['compound'] for l in ll));cov.append({'species':species,'literature_compounds':';'.join(names),'compound_count':len(names),'matrices':';'.join(sorted(set(l['reported_matrix'] for l in ll))) or 'pending_GCMS_and_taxonomic_voucher','sample_presence':'not_confirmed','note':'extract_comparators_only; volatile_composition_pending' if species=='Eugenia catharinae' else ('no_chemical_assignment; no_generic_C_xanthocarpa_substitution' if not names else 'literature_selection_not_exhaustive')})
    write(ROOT/'species_coverage.csv',cov);sources['retrieved_date']='2026-10-25';(ROOT/'literature_sources.json').write_text(json.dumps(sources,indent=2,ensure_ascii=False))
    descriptors=[]
    for s in ss:
        m=Chem.MolFromSmiles(s['isomeric_smiles']);descriptors.append({'structure_id':s['structure_id'],'compound':s['compound'],'panel':s['panel'],'MW':Descriptors.MolWt(m),'cLogP_RDKit':Descriptors.MolLogP(m),'TPSA_A2':Descriptors.TPSA(m),'HBD':rdMolDescriptors.CalcNumHBD(m),'HBA':rdMolDescriptors.CalcNumHBA(m),'rotatable_bonds':rdMolDescriptors.CalcNumRotatableBonds(m),'formal_charge':Chem.GetFormalCharge(m),'interpretation':'calculated_descriptors_only_not_experimental_ADMET_or_toxicology'})
    write(ROOT/'physicochemical_descriptors.csv',descriptors)
    w=Chem.SDWriter(str(ROOT/'prepared/all_inputs.sdf'))
    for row in inputs:
        m=Chem.SDMolSupplier(str(ROOT/row['input_path']),removeHs=False)[0];assert m.GetProp('_Name')==row['input_id'];w.write(m)
    w.close()
    print(json.dumps({'compounds':len(set(s['compound'] for s in ss)),'structures':len(ss),'input_conformers':len(inputs),'species_counts':{r['species']:r['compound_count'] for r in cov}},indent=2),flush=True)

if __name__=='__main__':main()
