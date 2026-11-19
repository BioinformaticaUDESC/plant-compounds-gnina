"""Prepare additional targets chosen after round A controls, before any B outcome."""
import json
from pathlib import Path
from prepare_targets_v4 import prepare
ROOT=Path(__file__).resolve().parent
SPECS_B=[
 ('PARP1_7KK4','7KK4','09L','A','ccd_neutral',None,'Human wild-type PARP1 catalytic-domain monomer, olaparib neutral reference; no DNA-binding/trapping, catalytic inhibition, synthetic lethality or tumor-cell selectivity demonstrated'),
 ('TYRS_1JIL','1JIL','485','A','amino_acid_zwitterion',None,'S.aureus tyrosyl-tRNA synthetase catalytic construction; ligandSB284485 amino-acid zwitterion; modeled323/420residues, single protomer, no full tRNA-recognition/dimer dynamics; not whole-cell antibacterial efficacy')]
def main():
 for i,spec in enumerate(SPECS_B,100):
  prepare(spec,i)
  p=ROOT/'targets'/spec[0]/'config.json';cfg=json.loads(p.read_text());cfg['protocol_file']='protocol_round_b.json';cfg['selection_round']='V4_B_complementary_exploratory';p.write_text(json.dumps(cfg,indent=2))
if __name__=='__main__':main()
