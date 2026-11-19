"""Aggregate completed V4 outcomes, keeping previous delivery and new runs distinct."""
import json,csv,statistics,collections,hashlib
from pathlib import Path
import numpy as np
from rdkit import Chem
import gnina_core as core
ROOT=Path(__file__).resolve().parent
MARKERS=['alpha-humulene','beta-caryophyllene','alpha-guaiene','alpha-bulnesene','beta-selinene','linalool','geranial','carvone','eucalyptol','elemol','quercetin','rotundic_acid']
NOFLAG='no_flag_from_limited_geometry_audit'
def rows(p):return list(csv.DictReader(p.open()))
def write(p,data):
 if data:core.write_csv(p,data)
def median(v):return statistics.median(float(x) for x in v)
def main():
 reports=ROOT/'reports';rationale=json.loads((ROOT/'target_rationale.json').read_text());names={x['target_id']:x for x in rationale['targets']}
 controls=[];new=[];runrows=[];sensitivity=[];parents=[];contacts=[];ranks=[];all_flags=[]
 for cp in sorted((ROOT/'targets').glob('*/config.json')):
  tid=cp.parent.name;d=ROOT/'results'/tid;v=json.loads((d/'validation.json').read_text());n=names[tid]
  controls.append({'target_id':tid,'origin':('NEW_V4_B_COMPLEMENTARY' if json.loads(cp.read_text()).get('selection_round') else 'NEW_V4_A_INITIAL'),'name_pt':n['name_pt'],'pdb':n['pdb'],'top1_passes':v['top1_passes'],'runs':v['runs'],'median_top1_RMSD_A':v['median_top1_RMSD_A'],'geometry_status':v['geometry_status'],'screening_eligible':v['screening_eligible']});new.append(v)
  rr=rows(d/'redocking/poses.csv');grouped=collections.defaultdict(list)
  for row in rr:grouped[row['seed']].append(row)
  for seed,pp in grouped.items():
   a=min(pp,key=lambda x:int(x['energy_rank']));b=max(pp,key=lambda x:float(x['CNNscore_fast']));best=min(pp,key=lambda x:float(x['heavy_atom_fixed_frame_symmetry_RMSD_A']))
   sensitivity.append({'target_id':tid,'seed':seed,'Energy_top1_RMSD_A':a['heavy_atom_fixed_frame_symmetry_RMSD_A'],'best_sampled_RMSD_A':best['heavy_atom_fixed_frame_symmetry_RMSD_A'],'CNN_top1_within_Energy_exported_pool_RMSD_A':b['heavy_atom_fixed_frame_symmetry_RMSD_A'],'Energy_top1_CNN_ranked_pose_Energy_rank':b['energy_rank'],'interpretation':'secondary_reordering_same_exported_pose_pool; not_independent_search_or_primary_gate; not_prospective_validation'})
  runrows.extend(rows(d/'redocking/runs.csv'))
  if not v['screening_eligible']:continue
  assert (d/'audit_complete.json').exists(),'Audit incomplete '+tid
  tops=rows(d/'top1_by_input_seed.csv');groups=collections.defaultdict(list)
  for row in tops:groups[(row['compound'],row['panel'])].append(row)
  for (compound,panel),pp in groups.items():
   ee=[float(x['minimizedAffinity_kcal_mol']) for x in pp];cc=[float(x['CNNscore_fast']) for x in pp];free=[x for x in pp if x['quality_flags']==NOFLAG]
   ss=[]
   for seed in [42,142,242]:ss.append(median(x['minimizedAffinity_kcal_mol'] for x in pp if int(x['seed'])==seed))
   parents.append({'target_id':tid,'compound':compound,'panel':panel,'stereochemical_structures':len({x['structure_id'] for x in pp}),'conformer_inputs':len({x['input_id'] for x in pp}),'seed_conformer_observations':len(pp),'Energy_median_kcal_mol':statistics.median(ee),'Energy_min_kcal_mol':min(ee),'Energy_max_kcal_mol':max(ee),'median_of_seedwise_medians_kcal_mol':statistics.median(ss),'CNNscore_fast_median':statistics.median(cc),'observations_without_limited_geometry_flags':len(free),'observations_with_flags':len(pp)-len(free),'Energy_median_excluding_flagged_diagnostic_only':median(x['minimizedAffinity_kcal_mol'] for x in free) if free else '', 'sample_presence':'not_confirmed','interpretation':'descriptive_sensitivity_over_stereo_conformers_and_search_seeds; not_statistical_replicates_or_potency'})
  counts=collections.Counter(flag for x in tops for flag in x['quality_flags'].split(';') if flag!=NOFLAG)
  for flag,count in counts.items():all_flags.append({'target_id':tid,'flag':flag,'top1_observations':count,'total_top1':len(tops)})
  contacts.extend(rows(d/'residue_contact_recurrence.csv'))
  poses=rows(d/'docking_poses.csv');pg=collections.defaultdict(list)
  for row in poses:pg[(row['input_id'],row['seed'])].append(row)
  for (iid,seed),pp in pg.items():
   a=min(pp,key=lambda x:int(x['energy_rank']));b=max(pp,key=lambda x:float(x['CNNscore_fast']))
   ranks.append({'target_id':tid,'input_id':iid,'compound':a['compound'],'seed':seed,'Energy_top1_kcal_mol':a['minimizedAffinity_kcal_mol'],'Energy_top1_CNNscore_fast':a['CNNscore_fast'],'CNN_best_within_exported_pool_Energy_rank':b['energy_rank'],'CNN_best_within_exported_pool_score':b['CNNscore_fast'],'interpretation':'restricted_exported_pool_secondary_diagnostic; no_independent_search_or_activity_probability'})
 baseline=[]
 for p in sorted((ROOT/'baseline_v3/results').glob('*/validation.json')):
  v=json.loads(p.read_text());baseline.append(v);controls.append({'target_id':v['target_id'],'origin':'PREVIOUS_V3_NOT_RERUN','name_pt':v['target_id'].split('_')[0],'pdb':v['target_id'].split('_')[1],'top1_passes':v['top1_passes'],'runs':v['runs'],'median_top1_RMSD_A':v['median_top1_RMSD_A'],'geometry_status':v['geometry_status'],'screening_eligible':v['screening_eligible']})
 write(reports/'all_targets_validation.csv',controls);write(reports/'new_redocking_runs.csv',runrows);write(reports/'redocking_secondary_rank_diagnostic.csv',sensitivity);write(reports/'new_compound_summaries.csv',parents);write(reports/'new_residue_contact_recurrence.csv',contacts);write(reports/'new_pose_ranking_comparison.csv',ranks);write(reports/'new_top1_flag_counts.csv',all_flags)
 current=[json.loads(p.read_text()) for p in (ROOT/'results').glob('*/audit_counts.json')];old=[json.loads(p.read_text()) for p in (ROOT/'baseline_v3/results').glob('*/audit_counts.json')]
 total={'compounds':40,'stereochemical_structures':51,'verified_conformer_inputs':59,'new_targets_tested':len(new),'new_targets_primary_redocking_pass':sum(x['screening_eligible'] for x in new),'new_targets_primary_redocking_fail':sum(not x['screening_eligible'] for x in new),'new_primary_screening_poses':sum(x['docking_poses'] for x in current),'new_primary_top1_observations':sum(x['top1_observations'] for x in current),'new_redocking_poses':sum(len(rows(p)) for p in (ROOT/'results').glob('*/redocking/poses.csv')),'new_top1_observations_with_flags':sum(x['top1_with_flags'] for x in current),'previous_targets_tested':len(baseline),'previous_targets_pass':sum(x['screening_eligible'] for x in baseline),'previous_primary_screening_poses':sum(x['docking_poses'] for x in old),'previous_redocking_poses':400,'all_target_systems_tested':len(new)+len(baseline),'all_target_systems_pass':sum(x['screening_eligible'] for x in new+baseline),'all_primary_screening_poses':sum(x['docking_poses'] for x in current+old),'all_redocking_poses':400+sum(len(rows(p)) for p in (ROOT/'results').glob('*/redocking/poses.csv')),'previous_results_recomputed':False,'computation_mode':'CPU;GNINA1.3.3;CNNfast;Energyprimary;exhaustiveness16;20modes;3screeningseeds;5controlseeds'}
 (reports/'study_counts.json').write_text(json.dumps(total,indent=2))
 lines=['# Ampliação de alvos GNINA V4','',f"Biblioteca curada: 40 compostos, 51 estruturas, 59 entradas. Novos alvos testados: {len(new)}; aprovados: {total['new_targets_primary_redocking_pass']}; reprovados: {total['new_targets_primary_redocking_fail']}.",'','|Alvo|Origem|Controle|RMSD mediano (A)|Triagem|','|---|---|---|---:|---|']
 for v in controls:lines.append(f"|{v['target_id']}|{v['origin']}|{v['top1_passes']}/5|{float(v['median_top1_RMSD_A']):.3f}|{'sim' if v['screening_eligible'] else 'não'}|")
 lines+=['',f"Novas poses de candidatos: {total['new_primary_screening_poses']}; novas poses de redocking: {total['new_redocking_poses']}. Acumulado com V3: {total['all_primary_screening_poses']} poses de candidatos e{total['all_redocking_poses']} de controles. Resultados anteriores não recalculados.",'','Aprovação refere-se à reprodução geométrica do ligante nativo no sistema e protocolo definidos. Não comprova atividade dos candidatos. Reprovação não demonstra inatividade biológica. Escores entre alvos não são escala comum de afinidade.','', 'Todas as 40 substâncias permanecem candidatas de literatura; ausência de GC-MS do projeto. Dados químicos de Campomanesia littoralis e composição completa do óleo de Eugenia catharinae continuam pendentes. Sem previsão validada de toxicologia/ADMET, efeito antidepressivo, antitumoral ou ação antimicrobiana celular.','', 'KEAP1 avalia sítio não covalente deNRF2, não modificação de Cys151 ou capacidade antioxidante direta. GyrB é uma construção com segmentos ausentes e exige cautela estrutural. Os controles falhos permanecem disponíveis e não geram rankings principais de candidatos.']
 (reports/'RESULTADOS_V4.md').write_text('\n'.join(lines)+'\n');print(json.dumps(total,indent=2),flush=True)
if __name__=='__main__':
 main()
 import combine_baseline_v4
 combine_baseline_v4.main()
