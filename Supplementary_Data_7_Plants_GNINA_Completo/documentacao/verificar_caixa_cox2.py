#!/usr/bin/env python3
"""Confere a caixa da COX-2 usando somente coordenadas e a biblioteca padrão.

Esta conta não executa docking e não estima o efeito de uma caixa maior.
"""
from pathlib import Path
import argparse,json

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--v3-root',type=Path,default=Path(__file__).resolve().parents[1]/'dados/Restinga_GNINA_Multialvo_V3')
    args=parser.parse_args()
    p=args.v3_root/'targets/COX2_5IKR'
    cfg=json.loads((p/'config.json').read_text())
    lines=(p/'native_reference.sdf').read_text().splitlines()
    count=int(lines[3][:3])
    xyz=[tuple(float(line[j:j+10]) for j in [0,10,20]) for line in lines[4:4+count] if line[31:34].strip()!='H']
    lo=[min(v[j] for v in xyz) for j in range(3)]
    hi=[max(v[j] for v in xyz) for j in range(3)]
    extent=[b-a for a,b in zip(lo,hi)]
    center=[(a+b)/2 for a,b in zip(lo,hi)]
    size=[max(22,x+12) for x in extent]
    equal=all(abs(a-b)<1e-8 for a,b in zip(center,cfg['center'])) and size==cfg['size']
    result={'native_heavy_atom_extent_xyz_A':extent,'native_extent_plus_12_A':[x+12 for x in extent],
            'rule_center_A':center,'inherited_center_A':cfg['center'],'rule_size_A':size,'inherited_size_A':cfg['size'],
            'identical':equal,'new_docking_runs':0}
    print(json.dumps(result,ensure_ascii=False,indent=2))
    if not equal:raise SystemExit(1)

if __name__=='__main__':main()
