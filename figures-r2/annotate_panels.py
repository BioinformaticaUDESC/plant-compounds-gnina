"""Regenerate selected PyMOL panels from archived sessions, using verified pose coordinates."""
import csv, hashlib, json, sys, shutil
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
import pymol2
from PIL import Image,ImageDraw,ImageFont
from itertools import permutations
from matplotlib import font_manager
root=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
m=json.loads((root/'figure_manifest.json').read_text());entries=m['native_and_candidate_entries']
entries=[e for e in entries if e['kind']=='candidate' or e['system_status']=='PASS']
assert len(entries)==14
rows=[];audit=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
with pymol2.PyMOL() as pm:
 c=pm.cmd
 for k in ('native','candidate'):
  for i,e in enumerate(x for x in entries if x['kind']==k):
   panel=str(3 if k=='native' else 7)+chr(65+i)
   c.reinitialize()
   base=root/'sessions'/f"{e['target_id']}_{k}.pse"
   c.load(str(base));camera=c.get_view()
   sels=['receptor','native','pose']
   coords={sel:c.get_coords(sel+' and not hydro').copy() for sel in sels}
   for sel,key in zip(sels,['receptor_file','reference_file','pose_file']):
    c.load(str(root/e[key]),'check')
    z=c.get_coords('check and not hydro')
    assert z.shape==coords[sel].shape and np.allclose(z,coords[sel],atol=0.0001,rtol=0),(panel,sel)
    c.delete('check')
   c.set_view(camera)
   atoms=c.get_model('receptor and polymer.protein and not hydro').atom
   nearest=cdist([a.coord for a in atoms],coords['pose']).min(axis=1)
   bins={}
   for a,d in zip(atoms,nearest):
    key=(a.chain,a.resn,a.resi)
    if d<=4 and (key not in bins or d<bins[key][0]):bins[key]=(float(d),a)
   ranked=sorted(bins.items(),key=lambda x:(x[1][0],x[0]))
   assert len(ranked)>=3
   picks=ranked[:3]
   for j,((chain,resn,resi),(d,a)) in enumerate(ranked):
    rows.append(dict(panel=panel,pdb=e['pdb'],chain=chain,residue=resn,number=resi,min_distance_A=f'{d:.6f}',labelled=int(j<3)))
   assert c.get_view()==camera
   assert all(np.array_equal(c.get_coords(sel+' and not hydro'),coords[sel])for sel in sels)
   # Project each exact atom through the archived camera. A marker-only ray
   # render provides the pixel coordinates without altering the published plate.
   c.disable('all')
   points=[]
   for j,(_,(_,a)) in enumerate(picks):
    name='anchor_'+str(j)
    c.pseudoatom(name,pos=a.coord)
    c.show('spheres',name);c.set('sphere_scale',.065,name);c.enable(name)
   c.bg_color('white');c.set('ray_opaque_background',1)
   c.set('ambient',1);c.set('specular',0);c.set('depth_cue',0)
   c.set('ray_trace_mode',0);c.set('max_threads',4)
   c.viewport(600,450);c.set_view(camera)
   for j in range(3):
    for k2 in range(3):
     (c.enable if j==k2 else c.disable)('anchor_'+str(k2))
    c.color('red','anchor_'+str(j))
    raster=out/'marker.png'
    c.png(str(raster),width=600,height=450,ray=1,quiet=1)
    pix=np.asarray(Image.open(raster).convert('RGB'))
    iy,ix=np.where((pix[:,:,0]>100)&(pix[:,:,1]<80)&(pix[:,:,2]<80))
    assert len(ix)>0,(panel,j)
    points.append((float(ix.mean()*2+1),float(iy.mean()*2+1)))
   raster.unlink(missing_ok=True)
   original=root/e['panel_file']
   assert sha(original)==e['panel_sha256']
   im=Image.open(original).convert('RGB')
   assert im.size==(1200,900)
   draw=ImageDraw.Draw(im)
   font=ImageFont.truetype(font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans',weight='bold')),37)
   def btw(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
   def crosses(a,b,c,d):return btw(a,b,c)*btw(a,b,d)<0 and btw(c,d,a)*btw(c,d,b)<0
   slots=[(134,130),(1066,130),(134,770),(1066,770),(130,450),(1070,450)]
   def cost(plan):
    return sum(np.linalg.norm(np.array(a)-b) for a,b in zip(points,plan))+4000*sum(crosses(points[i],plan[i],points[j],plan[j])for i in range(3)for j in range(i))
   plan=min(permutations(slots,3),key=cost)
   figs=[]
   for ((chain,resn,resi),_),pt,box in zip(picks,points,plan):
    lab=f'{resn}{resi} ({chain})' if chain else f'{resn}{resi}'
    draw.line([box,pt],fill='white',width=10)
    draw.line([box,pt],fill='#18282d',width=4)
    x,y=pt;draw.ellipse((x-7,y-7,x+7,y+7),fill='#18282d',outline='white',width=2)
    figs.append((lab,box,pt))
   for lab,(x,y),pt in figs:
    w=draw.textlength(lab,font=font)/2+11
    draw.rounded_rectangle((x-w,y-26,x+w,y+26),radius=9,fill='white',outline='#18282d',width=3)
    draw.text((x,y),lab,font=font,anchor='mm',fill='#111d20')
   png=out/e['panel_file'];png.parent.mkdir(exist_ok=True)
   im.save(png,dpi=(300,300))
   audit.append(dict(panel=panel,pdb=e['pdb'],seed=e['seed'],source_session_sha256=sha(base),source_panel_sha256=sha(original),labelled_panel_sha256=sha(png),pose_sha256=sha(root/e['pose_file']),receptor_sha256=sha(root/e['receptor_file']),reference_sha256=sha(root/e['reference_file']),labels=[dict(name=lab,text_pixel=pos,atom_pixel=pt) for lab,pos,pt in figs],coordinate_check='archive files match session objects; camera unchanged'))
   print(panel,len(ranked),[x[0] for x in picks],flush=True)
with (out/'residue_proximities.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
(out/'figure_label_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
m['native_and_candidate_entries']=entries
(out/'figure_manifest.json').write_text(json.dumps(m,indent=2)+'\n')
