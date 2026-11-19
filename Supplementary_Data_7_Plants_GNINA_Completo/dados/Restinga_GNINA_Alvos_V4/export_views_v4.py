"""Export actual top1 native poses and portable PyMOL instructions without ligand fitting."""
from pathlib import Path
import json
from rdkit import Chem
import gnina_core as c
ROOT=Path(__file__).resolve().parent
def main():
 for cp in sorted((ROOT/'targets').glob('*/config.json')):
  cfg=json.loads(cp.read_text());tid=cfg['target_id'];pose=c.load(ROOT/'results'/tid/'redocking'/f'{tid}_native_seed101.sdf')[0];p=cp.parent/'native_top1_seed101.sdf';w=Chem.SDWriter(str(p));w.write(pose);w.close()
  text=f'''# Execute from project root: @targets/{tid}/view_redocking.pml
# Actual coordinates; no align, fit, or super command.
reinitialize
load targets/{tid}/receptor_H.pdb, receptor
load targets/{tid}/native_reference.sdf, crystal
load targets/{tid}/native_top1_seed101.sdf, redocking
hide everything
remove hydro
show cartoon, receptor
show sticks, crystal or redocking
show sticks, receptor within 4 of crystal
color gray80, receptor
color marine, crystal
color tv_orange, redocking
set stick_radius, 0.16
set cartoon_transparency, 0.55
bg_color white
orient crystal
zoom crystal, 8
set ray_opaque_background, on
'''
  (cp.parent/'view_redocking.pml').write_text(text)
 print('Portable PyMOL scripts and actual seed101 Energy-top1 native SDFs exported',flush=True)
if __name__=='__main__':main()
