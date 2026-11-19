#!/usr/bin/env python3
"""ADMET complementar: identidade química, descritores e inferência local auditável."""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.metadata as md
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n')

def version(name):
    try: return md.version(name)
    except md.PackageNotFoundError: return 'not_installed'

def validate_bundle():
    expected = json.loads((ROOT/'checksums.json').read_text())
    for name, digest in expected.items():
        p = ROOT/name
        if not p.is_file() or sha(p) != digest:
            raise ValueError(f'Arquivo ausente ou alterado: {name}. Recupere o pacote original.')

def prepare_inputs():
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    structures = pd.read_csv(ROOT/'data/original/structures.csv', keep_default_na=False)
    docking = pd.read_csv(ROOT/'data/original/primary_docking_inputs.csv', keep_default_na=False)
    if (len(structures), structures.compound.nunique(), len(docking)) != (51, 40, 59):
        raise ValueError('A biblioteca não corresponde a 51 estruturas / 40 compostos / 59 SDF.')
    if structures.structure_id.duplicated().any() or docking.input_id.duplicated().any():
        raise ValueError('Identificadores duplicados.')
    if set(docking.structure_id) != set(structures.structure_id):
        raise ValueError('Mapeamento structures.csv / primary_docking_inputs.csv incompleto.')
    canonical = []
    for row in structures.itertuples():
        mol = Chem.MolFromSmiles(row.isomeric_smiles)
        if mol is None or len(Chem.GetMolFrags(mol)) != 1:
            raise ValueError(f'SMILES inválido ou multicomponente: {row.structure_id}')
        if Chem.MolToInchiKey(mol) != row.structure_inchikey:
            raise ValueError(f'InChIKey divergente: {row.structure_id}')
        if Chem.GetFormalCharge(mol) != int(row.formal_charge):
            raise ValueError(f'Carga divergente: {row.structure_id}')
        if rdMolDescriptors.CalcMolFormula(mol) != row.formula:
            raise ValueError(f'Fórmula divergente: {row.structure_id}')
        canonical.append(Chem.MolToSmiles(mol, isomericSmiles=True))
    structures['model_smiles'] = canonical
    checks = []
    mapping = structures.set_index('structure_id')
    for row in docking.itertuples():
        p = ROOT/'data/sdf'/Path(row.input_path).name
        if sha(p) != row.sha256:
            raise ValueError(f'SHA256 divergente: {p.name}')
        mols = list(Chem.SDMolSupplier(str(p), removeHs=True))
        if len(mols) != 1 or mols[0] is None:
            raise ValueError(f'SDF ilegível: {p.name}')
        ik = Chem.MolToInchiKey(mols[0])
        if ik != mapping.loc[row.structure_id, 'structure_inchikey']:
            raise ValueError(f'Identidade SDF/CSV divergente: {p.name}')
        checks.append({'input_id': row.input_id, 'structure_id': row.structure_id,
                       'sha256': row.sha256, 'sdf_inchikey': ik, 'status': 'PASS'})
    return structures, docking, pd.DataFrame(checks)

def descriptors(structures):
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors, FilterCatalog
    catalogs = {}
    for name in ['PAINS', 'BRENK', 'NIH']:
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(getattr(FilterCatalog.FilterCatalogParams.FilterCatalogs, name))
        catalogs[name] = FilterCatalog.FilterCatalog(params)
    rows = []
    for row in structures.itertuples():
        m = Chem.MolFromSmiles(row.model_smiles)
        mw, logp = Descriptors.MolWt(m), Crippen.MolLogP(m)
        hba, hbd = Lipinski.NumHAcceptors(m), Lipinski.NumHDonors(m)
        result = {'structure_id': row.structure_id, 'compound': row.compound, 'panel': row.panel,
                  'formula': rdMolDescriptors.CalcMolFormula(m), 'formal_charge': Chem.GetFormalCharge(m),
                  'MW_Da': mw, 'cLogP_RDKit': logp, 'HBA': hba, 'HBD': hbd,
                  'TPSA_A2': rdMolDescriptors.CalcTPSA(m),
                  'rotatable_bonds_RDKit': Lipinski.NumRotatableBonds(m),
                  'heavy_atoms': m.GetNumHeavyAtoms(), 'fraction_Csp3': rdMolDescriptors.CalcFractionCSP3(m),
                  'QED': QED.qed(m),
                  'Lipinski_violations': sum([mw > 500, logp > 5, hba > 10, hbd > 5])}
        for name, cat in catalogs.items():
            hits = [match.GetDescription() for match in cat.GetMatches(m)]
            result[name+'_alert_count'] = len(hits)
            result[name+'_alerts'] = '; '.join(hits)
        rows.append(result)
    return pd.DataFrame(rows)

def sensitivity_inputs(structures):
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    uncharger = rdMolStandardize.Uncharger()
    rows = []
    for row in structures.itertuples():
        parent = uncharger.uncharge(Chem.MolFromSmiles(row.model_smiles))
        smiles = Chem.MolToSmiles(parent, isomericSmiles=True)
        if smiles != row.model_smiles:
            rows.append({'structure_id': row.structure_id, 'compound': row.compound, 'panel': row.panel,
                         'primary_smiles': row.model_smiles, 'neutral_sensitivity_smiles': smiles,
                         'primary_formal_charge': row.formal_charge,
                         'sensitivity_formal_charge': Chem.GetFormalCharge(parent)})
    return pd.DataFrame(rows, columns=['structure_id','compound','panel','primary_smiles',
        'neutral_sensitivity_smiles','primary_formal_charge','sensitivity_formal_charge'])

def predict(structures, sensitivity, out, torch_threads):
    import numpy as np
    import pandas as pd
    import torch
    import admet_ai
    from admet_ai import ADMETModel
    from admet_ai.constants import DEFAULT_ADMET_PATH, DEFAULT_MODELS_DIR
    if version('admet-ai') != '2.0.1':
        raise ValueError('Este protocolo requer admet-ai==2.0.1. Use o ambiente indicado.')
    torch.set_num_threads(torch_threads)
    torch.manual_seed(2026)
    # CPU deliberadamente: painel pequeno, sem alterar drivers/CUDA da máquina.
    model = ADMETModel(include_physchem=False, drugbank_path=None, num_workers=0)
    if model.device != 'cpu':
        raise RuntimeError('O protocolo requer CPU; CUDA_VISIBLE_DEVICES não foi aplicado.')
    tasks = [task for group in model.task_lists for task in group]
    if len(tasks) != len(set(tasks)) or len(tasks) != 41:
        raise ValueError(f'Esperados 41 endpoints únicos; encontrados {len(tasks)}.')
    metadata = pd.read_csv(DEFAULT_ADMET_PATH)
    endpoint_meta = metadata[metadata.id.isin(tasks)].copy()
    if set(endpoint_meta.id) != set(tasks):
        raise ValueError('Metadados dos endpoints incompletos.')
    metadata.to_csv(out/'vendor_endpoint_metadata.csv', index=False)
    # Não usar percentis DrugBank como probabilidade, AD ou incerteza.
    unique_smiles = list(dict.fromkeys(structures.model_smiles.tolist() +
                        sensitivity.neutral_sensitivity_smiles.tolist()))
    raw = model.predict(smiles=unique_smiles)
    if raw.index.duplicated().any() or set(raw.index) != set(unique_smiles):
        raise ValueError('Saída do modelo perdeu ou duplicou identidades moleculares.')
    if set(raw.columns) != set(tasks) or not np.isfinite(raw.to_numpy()).all():
        raise ValueError('Saída incompleta, não finita ou com endpoints inesperados.')
    for task in endpoint_meta[endpoint_meta.task_type=='classification'].id:
        if not raw[task].between(-1e-6, 1+1e-6).all():
            raise ValueError(f'Score de classificação fora de [0,1]: {task}')
    raw.index.name = 'model_smiles'
    raw.to_csv(out/'raw_model_predictions.csv')
    prefix = structures[['structure_id','compound','panel','model_smiles']].copy()
    primary = prefix.merge(raw.reset_index(), on='model_smiles', how='left', validate='many_to_one')
    if len(primary) != 51:
        raise ValueError('O resultado primário deve ter 51 estruturas.')
    primary.to_csv(out/'admet_by_structure.csv', index=False)
    long = primary.melt(id_vars=prefix.columns.tolist(), value_vars=tasks,
                        var_name='id', value_name='predicted_value')
    long = long.merge(endpoint_meta, on='id', how='left', validate='many_to_one')
    long['evidence_type'] = 'model_prediction_not_experimental'
    long.to_csv(out/'admet_long_with_units.csv', index=False)
    summary = []
    for (compound, panel), group in primary.groupby(['compound','panel'], sort=False):
        for task in tasks:
            summary.append({'compound': compound, 'panel': panel, 'endpoint': task,
                'n_declared_structures': len(group), 'min_prediction': group[task].min(),
                'max_prediction': group[task].max(),
                'range_interpretation': 'representation_range_not_confidence_interval'})
    pd.DataFrame(summary).to_csv(out/'admet_by_compound_ranges.csv', index=False)
    changed = []
    for row in sensitivity.itertuples():
        for task in tasks:
            a, b = raw.loc[row.primary_smiles,task], raw.loc[row.neutral_sensitivity_smiles,task]
            changed.append({'structure_id': row.structure_id, 'compound': row.compound,
                'endpoint': task, 'primary_value': a, 'neutral_sensitivity_value': b,
                'delta_neutral_minus_primary': b-a})
    pd.DataFrame(changed).to_csv(out/'charge_sensitivity_predictions.csv', index=False)
    weights = {str(p.relative_to(DEFAULT_MODELS_DIR)): sha(p)
               for p in sorted(Path(DEFAULT_MODELS_DIR).rglob('*.pt'))}
    if not weights: raise ValueError('Nenhum peso de modelo encontrado para registro.')
    write_json(out/'model_weights_sha256.json', weights)
    shutil.copy2(DEFAULT_ADMET_PATH, out/'vendor_endpoint_metadata_original.csv')
    manifest = {'engine': 'ADMET-AI', 'version': version('admet-ai'), 'device': model.device,
                'num_workers': 0, 'torch_threads': torch_threads, 'tasks': tasks,
                'n_endpoints':len(tasks), 'n_unique_model_smiles':len(unique_smiles),
                'n_primary_structures':len(primary), 'n_charge_sensitivity_rows':len(sensitivity),
                'include_physchem':False, 'drugbank_reference':None,
                'package_path':str(Path(admet_ai.__file__).parent)}
    return manifest

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--mode', choices=['full','descriptors'], default='full')
    ap.add_argument('--threads', type=int, default=4)
    args = ap.parse_args()
    if args.threads < 1: ap.error('--threads deve ser positivo')
    os.environ['CUDA_VISIBLE_DEVICES'] = ''
    os.environ.setdefault('OMP_NUM_THREADS', str(args.threads))
    os.environ.setdefault('MKL_NUM_THREADS', str(args.threads))
    os.environ.setdefault('OPENBLAS_NUM_THREADS', str(args.threads))
    if args.out.exists() and any(args.out.iterdir()):
        ap.error('A pasta de saída não está vazia; escolha um novo caminho.')
    args.out.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    protocol = {'status':'running','started_utc':started,'mode':args.mode,
        'structures_source':'original docking library, no silent neutralization',
        'analysis_unit':'declared chemical structure; conformers not replicates',
        'scope':'post hoc exploratory ADMET; volatile and extract panels kept separate',
        'charge_policy':'primary = original microstate; neutralized changed molecules = sensitivity only',
        'applicability_domain':'not established for this panel; no validated per-molecule AD filter',
        'classification_policy':'continuous positive-class scores; no binary safe/toxic labels',
        'regression_policy':'raw model scale and vendor units retained; no back-transform/clipping',
        'uncertainty_policy':'structure ranges are not uncertainty estimates or confidence intervals',
        'versions':{x:version(x) for x in ['admet-ai','rdkit','torch','chemprop','lightning','numpy','pandas','scipy']},
        'python':sys.version,'platform':platform.platform(),'command':sys.argv,
        'package_file_sha256':sha(Path(__file__))}
    try:
        validate_bundle()
        structures, docking, checks = prepare_inputs()
        structures.to_csv(args.out/'validated_structures.csv', index=False)
        docking.to_csv(args.out/'docking_input_mapping.csv', index=False)
        checks.to_csv(args.out/'input_identity_audit.csv', index=False)
        descriptors(structures).to_csv(args.out/'physchem_by_structure.csv', index=False)
        sensitivity = sensitivity_inputs(structures)
        sensitivity.to_csv(args.out/'charge_sensitivity_inputs.csv', index=False)
        if args.mode == 'full':
            protocol['model'] = predict(structures, sensitivity, args.out, args.threads)
            protocol['status'] = 'completed_full_inference'
        else:
            protocol['status'] = 'completed_descriptors_only_no_ADMET_predictions'
        protocol['counts'] = {'compounds':structures.compound.nunique(),
                              'structures':len(structures),'sdf_entries':len(docking)}
        shutil.copy2(ROOT/'checksums.json', args.out/'input_checksums.json')
        proc = subprocess.run([sys.executable,'-m','pip','freeze'], capture_output=True, text=True)
        (args.out/'environment_freeze.txt').write_text(proc.stdout)
        protocol['pip_freeze_returncode'] = proc.returncode
    except Exception as exc:
        protocol['status']='failed'; protocol['error']=repr(exc)
        raise
    finally:
        protocol['finished_utc']=datetime.now(timezone.utc).isoformat()
        write_json(args.out/'protocol.json', protocol)
    outputs = {p.name:sha(p) for p in sorted(args.out.iterdir()) if p.is_file()}
    write_json(args.out/'output_sha256.json', outputs)
    print(f'Concluído: {args.out.resolve()} | {protocol["status"]}')

if __name__ == '__main__':
    main()
