from pathlib import Path
import json,io,textwrap
from collections import Counter
from PIL import Image,ImageDraw,ImageFont
from rdkit import Chem,RDLogger
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from matplotlib import font_manager
RDLogger.DisableLog('rdApp.warning')
ROOT=Path(__file__).resolve().parent
M=json.loads((ROOT/'figure_manifest.json').read_text())
structures=sorted(M['chemical_structures'],key=lambda r:(r['panel']!='volatile_literature_panel',r['compound'],r['structure_id']))
regular=font_manager.findfont('DejaVu Sans');bold=font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans',weight='bold'))
def font(size,b=False):return ImageFont.truetype(bold if b else regular,size)
def name(s):
 for a,b in [('alpha-','α-'),('beta-','β-'),('gamma-','γ-'),('delta-','δ-')]:s=s.replace(a,b)
 return s
def draw_mol(row,w,h,svg=False):
 mol=Chem.MolFromSmiles(row['isomeric_smiles']);assert mol is not None
 key=Chem.MolToInchiKey(mol);assert key==row['structure_inchikey']
 Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
 rdDepictor.SetPreferCoordGen(row['compound']=='eucalyptol')
 rdDepictor.Compute2DCoords(mol)
 assert Chem.MolToInchiKey(mol)==key
 dr=rdMolDraw2D.MolDraw2DSVG(w,h) if svg else rdMolDraw2D.MolDraw2DCairo(w,h)
 op=dr.drawOptions();op.bondLineWidth=2.7;op.minFontSize=20;op.maxFontSize=36;op.padding=.10;op.addStereoAnnotation=True
 dr.DrawMolecule(mol);dr.FinishDrawing()
 return dr.GetDrawingText()
dep=[]
for idx,row in enumerate(structures,1):
 row=dict(row);row['figure_id']=f'C{idx:02d}'
 stem=row['figure_id']+'_'+row['structure_id']
 (ROOT/'molecules'/(stem+'.svg')).write_text(draw_mol(row,950,600,True))
 (ROOT/'molecules'/(stem+'.png')).write_bytes(draw_mol(row,950,600,False))
 row['png']='molecules/'+stem+'.png';row['svg']='molecules/'+stem+'.svg';dep.append(row)
def grid(part,path,cols=3,cellw=900,cellh=650):
 rows=(len(part)+cols-1)//cols;canvas=Image.new('RGB',(cols*cellw,rows*cellh+75),'white');d=ImageDraw.Draw(canvas)
 for i,r in enumerate(part):
  x=(i%cols)*cellw;y=(i//cols)*cellh
  title=r['figure_id']+'  '+name(r['compound'])
  d.text((x+26,y+16),title,font=font(31,True),fill='black')
  if len([z for z in dep if z['compound']==r['compound']])>1:
   subtype=r['structure_id'].split('_isomer')[-1];d.text((x+26,y+58),'Declared stereochemical structure '+subtype,font=font(24),fill='#56616B')
  else:d.text((x+26,y+58),'Declared library structure',font=font(24),fill='#56616B')
  mol=Image.open(ROOT/r['png']).convert('RGB');mol.thumbnail((cellw-30,cellh-110),Image.Resampling.LANCZOS)
  canvas.paste(mol,(x+(cellw-mol.width)//2,y+95+(cellh-110-mol.height)//2))
 d.text((26,rows*cellh+20),'Chemical structures from the curated library  |  Occurrence in project samples is unconfirmed',font=font(25),fill='#45515D')
 canvas.save(ROOT/'figures'/path,dpi=(300,300))
vol=[r for r in dep if r['panel']=='volatile_literature_panel'];non=[r for r in dep if r['panel']!='volatile_literature_panel']
assert len(vol)==43 and len(non)==8
pages=[]
for i in range(0,len(vol),12):
 fn=f'Figure_S{3+i//12}_Volatile_Structures.png';grid(vol[i:i+12],fn);pages.append({'file':'figures/'+fn,'ids':[r['figure_id'] for r in vol[i:i+12]],'panel':'volatile'})
grid(non,'Figure_S7_Nonvolatile_Structures.png',cols=2,cellw=1250,cellh=610)
pages.append({'file':'figures/Figure_S7_Nonvolatile_Structures.png','ids':[r['figure_id'] for r in non],'panel':'nonvolatile'})
high=['delta-cadinene','alpha-bulnesene','caryophyllene oxide','alpha-cadinol','(E)-nerolidol','beta-caryophyllene','alpha-humulene','quercetin','myricitrin']
chosen=[next(r for r in dep if r['compound']==c) for c in high]
grid(chosen,'Figure_5_Highlighted_Molecules.png')
(ROOT/'molecule_depictions.json').write_text(json.dumps({'structures':dep,'supplementary_pages':pages,'main_figure_ids':[r['figure_id'] for r in chosen]},indent=2))
print(json.dumps({'depicted_structures':len(dep),'volatile':len(vol),'nonvolatile':len(non),'supplementary_pages':len(pages)}))
