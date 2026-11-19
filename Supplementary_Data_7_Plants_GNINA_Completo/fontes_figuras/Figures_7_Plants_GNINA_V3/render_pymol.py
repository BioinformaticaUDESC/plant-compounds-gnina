from pathlib import Path
import json,sys,hashlib
import pymol2

ROOT=Path(__file__).resolve().parent
manifest=json.loads((ROOT/'figure_manifest.json').read_text())
entries=manifest['native_and_candidate_entries']
if len(sys.argv)>1: entries=[e for e in entries if e['target_id']==sys.argv[1] and (len(sys.argv)<3 or e['kind']==sys.argv[2])]

with pymol2.PyMOL() as pm:
 c=pm.cmd;manifest['PyMOL_version']=list(c.get_version())
 for e in entries:
  c.reinitialize()
  c.set('max_threads',4);c.set('retain_order',1)
  c.load(str(ROOT/e['receptor_file']),'receptor')
  c.load(str(ROOT/e['reference_file']),'native')
  c.load(str(ROOT/e['pose_file']),'pose')
  c.hide('everything','all');c.remove('hydro')
  c.bg_color('white');c.set('orthoscopic',1);c.set('ray_opaque_background',1)
  c.set('antialias',2);c.set('ray_shadows',0);c.set('ambient',.55);c.set('specular',.15)
  c.set('depth_cue',0);c.set('ray_trace_fog',0);c.set('ray_trace_mode',1)
  c.set('cartoon_fancy_helices',1);c.set('cartoon_smooth_loops',1)
  c.set_color('protein_gray',[.75,.79,.82]);c.set_color('pocket_gray',[.57,.62,.67])
  c.set_color('native_teal',[.08,.49,.53]);c.set_color('redock_magenta',[.67,.22,.55]);c.set_color('candidate_amber',[.88,.48,.08])
  c.color('protein_gray','receptor');c.show('cartoon','receptor');c.set('cartoon_transparency',.60,'receptor')
  c.select('pocket','byres (receptor within 4.0 of (native or pose))')
  c.select('shell','byres (receptor within 6.0 of (native or pose))')
  c.show('sticks','pocket and not name N+C+O');c.set('stick_radius',.12,'pocket');c.color('pocket_gray','pocket')
  c.show('surface','shell');c.set('transparency',.78,'receptor');c.set('surface_quality',1);c.set('surface_color','protein_gray','receptor')
  c.show('sticks','native or pose');c.set('stick_radius',.22,'native or pose')
  c.color('native_teal','native');c.color('redock_magenta' if e['kind']=='native' else 'candidate_amber','pose')
  for sel in ['native','pose','pocket']:
   for elem,col in [('N','blue'),('O','red'),('S','yellow'),('F','green'),('Cl','green'),('P','orange')]:c.color(col,f'{sel} and elem {elem}')
  c.orient('native or pose');c.zoom('native or pose',buffer=6.0,complete=1)
  c.clip('slab',70)
  c.viewport(1200,900)
  c.save(str(ROOT/'sessions'/f'{e["target_id"]}_{e["kind"]}.pse'))
  c.png(str(ROOT/e['panel_file']),width=1200,height=900,dpi=300,ray=1)
  e['panel_sha256']=hashlib.sha256((ROOT/e['panel_file']).read_bytes()).hexdigest()
  print('RENDERED',e['target_id'],e['kind'],flush=True)
(ROOT/'figure_manifest.json').write_text(json.dumps(manifest,indent=2))
