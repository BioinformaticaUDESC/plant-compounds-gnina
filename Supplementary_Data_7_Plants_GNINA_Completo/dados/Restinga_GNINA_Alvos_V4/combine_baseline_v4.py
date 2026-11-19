"""Harmonized long tables: prior curated V3 and new V4 results, never pooled across targets."""
import csv,json,statistics
from pathlib import Path
from collections import defaultdict
import gnina_core as core
ROOT=Path(__file__).resolve().parent
NOFLAG='no_flag_from_limited_geometry_audit'
def read(p):return list(csv.DictReader(p.open()))
def main():
 compounds=[];alltop=[];seeds=[]
 sources=[(ROOT/'baseline_v3/results','PREVIOUS_V3_NOT_RERUN'),(ROOT/'results','NEW_V4')]
 for parent,origin in sources:
  for p in sorted(parent.glob('*/top1_by_input_seed.csv')):
   tid=p.parent.name;rr=read(p);assert len(rr)==177;groups=defaultdict(list)
   for row in rr:
    row['result_origin']=origin;alltop.append(row);groups[(row['compound'],row['panel'])].append(row)
   for (compound,panel),gg in groups.items():
    ee=[float(x['minimizedAffinity_kcal_mol']) for x in gg]
    compounds.append({'target_id':tid,'result_origin':origin,'compound':compound,'panel':panel,'input_conformers':len({x['input_id'] for x in gg}),'stereochemical_structures':len({x['structure_id'] for x in gg}),'top1_observations':len(gg),'Energy_median_kcal_mol':statistics.median(ee),'Energy_min_kcal_mol':min(ee),'Energy_max_kcal_mol':max(ee),'CNNscore_fast_median':statistics.median(float(x['CNNscore_fast']) for x in gg),'observations_with_geometry_flags':sum(x['quality_flags']!=NOFLAG for x in gg),'sample_presence':'not_confirmed','interpretation':'descriptive_sensitivity_within_this_target; not_cross_target_affinity_or_selectivity; same_curated59inputs'})
   for row in read(p.parent/'seed_stability_by_input.csv'):row['result_origin']=origin;seeds.append(row)
 def write_union(p,rr):
  keys=list(dict.fromkeys(k for row in rr for k in row));core.write_csv(p,[{k:row.get(k,'') for k in keys} for row in rr])
 write_union(ROOT/'reports/all_compound_summaries_with_baseline.csv',compounds);write_union(ROOT/'reports/all_top1_with_baseline.csv',alltop);write_union(ROOT/'reports/all_seed_stability_with_baseline.csv',seeds)
 counts=json.loads((ROOT/'reports/study_counts.json').read_text());assert len(compounds)==40*counts['all_target_systems_pass'];assert len(alltop)==177*counts['all_target_systems_pass']
 print('COMBINED',len(compounds),'compound-target summaries;',len(alltop),'top1 observations;prior results explicitly labeled',flush=True)
if __name__=='__main__':main()
