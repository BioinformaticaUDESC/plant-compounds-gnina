"""Reproduce the revision's post hoc sensitivity analyses from rounded Table S1a medians."""
import argparse,csv
from pathlib import Path
import numpy as np
from scipy.stats import kendalltau
from scipy.optimize import least_squares,linprog
parser=argparse.ArgumentParser()
parser.add_argument('--out',type=Path,required=True)
out=parser.parse_args().out
out.mkdir(parents=True,exist_ok=True)
with Path(__file__).with_name('table_s1a_volatile_medians.csv').open(newline='') as f: rows=list(csv.reader(f))
targets=rows[0][2:]
names=[r[0] for r in rows[1:]]
x=np.array([int(r[1]) for r in rows[1:]])
y=np.array([[float(v) for v in r[2:]] for r in rows[1:]])
X=np.column_stack([np.ones(32),x])
res=y-X@np.linalg.lstsq(X,y,rcond=None)[0]
idx=[0,1,5,6]
pairs=[(i,j)for j in idx for i in range(32) if res[i,j]-res[:,j].min()<=.34]
def gaps(r):return r-r.min(axis=0)
assert set(pairs)=={(i,j)for j in idx for i in range(32)if gaps(res)[i,j]<=.34}
rng=np.random.default_rng(4725052); hits=np.zeros(len(pairs),int)
for _ in range(2000):
 sample=rng.integers(0,32,32); g=gaps(y-X@np.linalg.lstsq(X[sample],y[sample],rcond=None)[0]); hits += [g[i,j]<=.34+1e-12 for i,j in pairs]
freq=hits/20
sweep=[]
for w in sorted({.34,*[i/100 for i in range(10,101,5)]}):
 c=(gaps(res)<=w+1e-12).sum(axis=0);sweep.append([f'{w:.2f}',*map(str,c),str(sum(c[idx]))])
persistent=sum(gaps(res)[i,j]<=.10 for i,j in pairs)
alt=[]
for j in idx:
 z=y[:,j];b=np.linalg.lstsq(X,z,rcond=None)[0]
 h=least_squares(lambda k:X@k-z,b,loss='huber',f_scale=.5).x
 # Minimize the sum of absolute residuals.
 l=linprog(np.r_[0,0,np.ones(32)],A_ub=np.block([[X,-np.eye(32)],[-X,-np.eye(32)]]),b_ub=np.r_[z,-z],bounds=[(None,None)]*2+[(0,None)]*32,method='highs')
 assert l.success
 ranks=[z,z/x,z-X@h,z-X@l.x[:2]]
 alt.append([targets[j],*[f'{kendalltau(res[:,j],v).statistic:.3f}' for v in ranks]])
assert persistent==4 and sum(freq>=90)==10 and len(pairs)==11
def write_csv(name,header,rows):
    with (out/name).open('w',newline='') as f: csv.writer(f).writerows([header,*rows])
write_csv('s15_window_sweep.csv',['window',*targets,'four_targets'],sweep)
write_csv('s15_bootstrap.csv',['target','compound','selected_pct'],[[targets[j],names[i],f'{v:.2f}']for (i,j),v in zip(pairs,freq)])
write_csv('s16_rank_agreement.csv',['target','raw_or_z','energy_per_heavy_atom','huber_residual','median_residual'],alt)
