from pathlib import Path
import json,math
from PIL import Image,ImageDraw,ImageFont
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parent
M=json.loads((ROOT/'figure_manifest.json').read_text())
regular=font_manager.findfont('DejaVu Sans');bold=font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans',weight='bold'))
def ft(s,b=False):return ImageFont.truetype(bold if b else regular,s)
def cname(s):
 for a,b in [('alpha-','α-'),('beta-','β-'),('gamma-','γ-'),('delta-','δ-')]:s=s.replace(a,b)
 return s
def compose(part,filename,start_letter='A'):
 rows=math.ceil(len(part)/2);pw,head,ih,gap=1200,155,900,24
 W=2*pw+3*gap;H=rows*(head+ih+gap)+175
 im=Image.new('RGB',(W,H),'white');d=ImageDraw.Draw(im)
 for i,e in enumerate(part):
  x=gap+(i%2)*(pw+gap);y=gap+(i//2)*(head+ih+gap)
  letter=chr(ord(start_letter)+i)
  d.text((x+15,y),f'{letter}  {e["name"]}  |  {e["pdb"]}',font=ft(44,True),fill='black')
  if e['kind']=='native':
   d.text((x+15,y+58),f'RMSD {e["RMSD_A"]:.3f} Å  |  seed 101',font=ft(36),fill='#34495D')
   color='#246961' if e['system_status']=='PASS' else '#9C4D37'
   d.text((x+15,y+103),f'System control: {e["native_passes"]}/5  |  {e["system_status"]}',font=ft(33),fill=color)
  else:
   d.text((x+15,y+58),cname(e['compound']),font=ft(39,True),fill='#864A13')
   d.text((x+15,y+107),f'Pose Energy {e["pose_Energy_kcal_mol"]:.2f} kcal/mol  |  seed 42',font=ft(32),fill='#34495D')
  pic=Image.open(ROOT/e['panel_file']).convert('RGB');assert pic.size==(1200,900)
  im.paste(pic,(x,y+head))
 y=H-120;kind=part[0]['kind'];keycolor='#AB388C' if kind=='native' else '#E07A14'
 for x,color,text in [(35,'#147D87','Experimental ligand carbons'),(865,keycolor,'Redocked carbons' if kind=='native' else 'Candidate carbons'),(1630,'#BCC9D1','Receptor context')]:
  d.rectangle((x,y+8,x+32,y+40),fill=color);d.text((x+48,y),text,font=ft(33),fill='#222222')
 d.text((35,y+62),'Original coordinates  |  No ligand fitting  |  Oxygen red and nitrogen blue  |  Gray surface is not electrostatic potential',font=ft(30),fill='#44515C')
 im.save(ROOT/'figures'/filename,dpi=(300,300))
 return {'file':'figures/'+filename,'entries':[{'target_id':e['target_id'],'kind':e['kind'],'panel':chr(ord(start_letter)+i)} for i,e in enumerate(part)]}
entries=M['native_and_candidate_entries'];passed=[e for e in entries if e['kind']=='native' and e['system_status']=='PASS'];candidates=[e for e in entries if e['kind']=='candidate'];failed=[e for e in entries if e['kind']=='native' and e['system_status']=='FAIL']
out=[]
out.append(compose(passed[:4],'Figure_3A_Native_Redocking.png'))
out.append(compose(passed[4:],'Figure_3B_Native_Redocking.png','E'))
out.append(compose(candidates[:4],'Figure_4A_Candidate_Poses.png'))
out.append(compose(candidates[4:],'Figure_4B_Candidate_Poses.png','E'))
out.append(compose(failed[:3],'Figure_S1_Failed_Controls.png'))
out.append(compose(failed[3:],'Figure_S2_Failed_Controls.png','D'))
(ROOT/'pymol_figure_layouts.json').write_text(json.dumps(out,indent=2))
print('Composed',len(out),'PyMOL plates')
