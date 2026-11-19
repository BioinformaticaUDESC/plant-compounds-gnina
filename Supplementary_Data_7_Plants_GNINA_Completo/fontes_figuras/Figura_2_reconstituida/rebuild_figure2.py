#!/usr/bin/env python3
"""Reconstituição de 30/10/2026 da Figura 2; não é o código histórico original.

Recalcula os ranks com empates médios a partir das medianas de duas casas
decimais efetivamente apresentadas na Tabela S1. Compara os 224 valores à
matriz arquivada e verificada antes de produzir PNG, SVG e PDF. A aparência
é uma reconstrução; não se afirma identidade de pixels com a figura antiga.
"""
from pathlib import Path
import argparse, csv, hashlib, json, platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
TARGETS=[('COX2_5IKR','COX-2'),('ACHE_4EY7','AChE'),('CB2_5ZTY','CB2'),
         ('SEH_3ANS','sEH'),('KEAP1_4XMB','KEAP1'),('PARP1_7KK4','PARP1'),('TYRS_1JIL','TyrS')]

def read_csv(path):
    with path.open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))

def average_ranks(values):
    ordered=sorted(range(len(values)),key=values.__getitem__)
    result=[None]*len(values)
    i=0
    while i<len(ordered):
        j=i+1
        while j<len(ordered) and values[ordered[j]]==values[ordered[i]]:j+=1
        rank=((i+1)+j)/2
        for k in ordered[i:j]:result[k]=rank
        i=j
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=HERE/'figure2_reported_volatile_medians.csv')
    parser.add_argument('--expected',type=Path,default=HERE/'verified_figure2_rank_matrix.csv')
    parser.add_argument('--output-dir',type=Path,default=HERE/'reconstructed')
    args=parser.parse_args()
    data=read_csv(args.input);expected={r['compound']:r for r in read_csv(args.expected)}
    if len(data)!=32 or len(expected)!=32:raise ValueError('São necessários 32 voláteis únicos.')
    ranks=np.column_stack([average_ranks([float(r[tid]) for r in data]) for tid,_ in TARGETS])
    discrepancies=[]
    for i,row in enumerate(data):
        for j,(tid,_) in enumerate(TARGETS):
            if ranks[i,j]!=float(expected[row['compound']][tid]):
                discrepancies.append({'compound':row['compound'],'target':tid,'rank':float(ranks[i,j]),'expected':expected[row['compound']][tid]})
    if discrepancies:raise ValueError(json.dumps(discrepancies,indent=2))
    args.output_dir.mkdir(parents=True,exist_ok=True)
    with (args.output_dir/'figure2_reconstructed_ranks.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['compound']+[tid for tid,_ in TARGETS])
        for row,values in zip(data,ranks):w.writerow([row['compound']]+values.tolist())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(8.8,11.4))
    im=ax.imshow(ranks,aspect='auto',interpolation='nearest',cmap='Blues_r',vmin=1,vmax=32)
    ax.set_xticks(range(7),[label for _,label in TARGETS]);ax.xaxis.tick_top()
    ax.set_yticks(range(32),[r['display_name'] for r in data])
    ax.tick_params(axis='both',which='both',length=0)
    ax.set_xticks(np.arange(-.5,7,1),minor=True);ax.set_yticks(np.arange(-.5,32,1),minor=True)
    ax.grid(which='minor',color='white',linewidth=.45)
    for i in range(32):
        for j in range(7):
            if ranks[i,j]<=3:ax.text(j,i,f'{ranks[i,j]:g}',ha='center',va='center',color='white',fontweight='bold',fontsize=8.8)
    for side in ['top','right']:ax.spines[side].set_visible(False)
    bar=fig.colorbar(im,ax=ax,fraction=.042,pad=.025)
    bar.set_label('Within-protein rank among 32 volatiles\n1 = most negative median Energy')
    bar.set_ticks([1,5,10,15,20,25,30,32])
    fig.subplots_adjust(left=.27,right=.89,bottom=.035,top=.96)
    for ext in ['png','svg','pdf']:
        fig.savefig(args.output_dir/f'Figure_2_Rank_Map_Reconstructed.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)
    metadata={'provenance':'new reconstruction dated 2026-10-30; original plotting source not recovered',
              'scope':'32 volatile compounds x 7 eligible targets; within-target ranks only',
              'median_precision':'two decimal places as presented in manuscript Table S1',
              'tie_rule':'average rank; lower, more negative Energy first',
              'rank_comparisons':224,'rank_discrepancies':0,'pixel_identity_to_original_claimed':False,
              'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),
              'expected_sha256':hashlib.sha256(args.expected.read_bytes()).hexdigest(),
              'python':platform.python_version(),'numpy':np.__version__,'matplotlib':matplotlib.__version__}
    (args.output_dir/'rank_reconstruction_check.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(metadata,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
