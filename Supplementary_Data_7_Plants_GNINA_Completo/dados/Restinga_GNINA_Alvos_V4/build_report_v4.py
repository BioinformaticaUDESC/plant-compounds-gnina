"""Create a source-grounded report only after all simulations and audits finish."""
import csv,json,statistics,math,html,re
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image,KeepTogether
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
ROOT=Path(__file__).resolve().parent
REPORTS=ROOT/'reports';FIG=REPORTS/'figures';FIG.mkdir(exist_ok=True)
BLUE='#23465B';TEAL='#34796E';BROWN='#A76348';GRAY='#68747B';LIGHT='#EAF0F2'
W,H=A4;PW=W-92
def rows(p):return list(csv.DictReader(p.open()))
def esc(t):return html.escape(str(t))
def num(v,nd=2):return f'{float(v):.{nd}f}'.replace('.',',')
def count(v):return format(int(v),',').replace(',','.')
def savefig(fig,name):
 fig.savefig(FIG/(name+'.png'),dpi=240,bbox_inches='tight',facecolor='white');fig.savefig(FIG/(name+'.svg'),bbox_inches='tight',facecolor='white');plt.close(fig)
def plots(targets):
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','axes.labelcolor':GRAY,'text.color':'#17252C'})
 fig,ax=plt.subplots(figsize=(8.2,4.15));labels=[]
 for i,t in enumerate(targets):
  d=rows(ROOT/'results'/t['target_id']/'redocking/runs.csv');v=json.loads((ROOT/'results'/t['target_id']/'validation.json').read_text());xx=np.linspace(-.13,.13,len(d))+i
  yy=[float(x['top1_RMSD_A']) for x in d];ax.scatter(xx,yy,s=28,color=[TEAL if y<=2 else BROWN for y in yy],zorder=3)
  ax.plot([i-.23,i+.23],[v['median_top1_RMSD_A']]*2,lw=2,color=BLUE)
  labels.append(t['target_id'].split('_')[0]+'\n'+t['pdb']+f"\n{v['top1_passes']}/5")
 ax.axhline(2,color=GRAY,ls='--',lw=1);ax.set_xticks(range(len(targets)),labels);ax.set_ylabel('RMSD pesado, referencial fixo (Å)');ax.set_title('Redocking: primeira pose ordenada por Energy',loc='left');ax.grid(axis='y',alpha=.15);fig.tight_layout();savefig(fig,'01_redocking_v4')
 try:
  from rdkit import Chem
  import gnina_core as core
 except ImportError:return
 nrows=math.ceil(len(targets)/3);fig,axes=plt.subplots(nrows,3,figsize=(9,2.7*nrows));export=[]
 for ax,t in zip(axes.flat,targets):
  tid=t['target_id'];ref=Chem.RemoveHs(core.load(ROOT/'targets'/tid/'native_reference.sdf')[0]);probe=Chem.RemoveHs(core.load(ROOT/'results'/tid/'redocking'/f'{tid}_native_seed101.sdf')[0]);a=ref.GetConformer().GetPositions();b=probe.GetConformer().GetPositions()
  # One common display projection; no molecular fitting or coordinate changes.
  allx=np.vstack([a,b]);origin=allx.mean(0);_,_,vh=np.linalg.svd(allx-origin,full_matrices=False);basis=vh[:2].T
  for m,x,color,label in [(ref,a,BLUE,'Cristal'),(probe,b,BROWN,'Redocking')]:
   xy=(x-origin)@basis
   for bond in m.GetBonds():
    u,v=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx();ax.plot(xy[[u,v],0],xy[[u,v],1],color=color,lw=1.2,alpha=.85)
   ax.scatter(xy[:,0],xy[:,1],s=7,color=color,label=label)
   for j,p in enumerate(x):export.append({'target_id':tid,'model':label,'atom_index_0based':j,'element':m.GetAtomWithIdx(j).GetSymbol(),'x_A':float(p[0]),'y_A':float(p[1]),'z_A':float(p[2]),'display_x_A':float(xy[j,0]),'display_y_A':float(xy[j,1])})
  d=core.rmsd(ref,probe);ax.set_title(tid.split('_')[0]+f' | semente 101 | {d:.2f} Å',fontsize=9);ax.axis('equal');ax.axis('off')
 for ax in list(axes.flat)[len(targets):]:ax.axis('off')
 handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=2,frameon=False);fig.tight_layout(rect=(0,.05,1,1));savefig(fig,'02_native_coordinates_v4');core.write_csv(REPORTS/'native_overlay_coordinates.csv',export)
def main():
 data=json.loads((ROOT/'target_rationale.json').read_text());targets=data['targets'];counts=json.loads((REPORTS/'study_counts.json').read_text());controls=rows(REPORTS/'all_targets_validation.csv');parents=rows(REPORTS/'new_compound_summaries.csv') if (REPORTS/'new_compound_summaries.csv').exists() else []
 for t in targets:
  v=json.loads((ROOT/'results'/t['target_id']/'validation.json').read_text())
  if v['screening_eligible']:assert (ROOT/'results'/t['target_id']/'audit_complete.json').exists()
 plots(targets)
 fontroot=Path(matplotlib.get_data_path())/'fonts/ttf'
 for name,file in [('D','DejaVuSans.ttf'),('DB','DejaVuSans-Bold.ttf'),('DI','DejaVuSans-Oblique.ttf')]:pdfmetrics.registerFont(TTFont(name,str(fontroot/file)))
 pdfmetrics.registerFontFamily('D',normal='D',bold='DB',italic='DI',boldItalic='DB')
 styles={
 'title':ParagraphStyle('title',fontName='DB',fontSize=21,leading=25,textColor=colors.HexColor('#101E26'),spaceAfter=14),
 'h':ParagraphStyle('h',fontName='DB',fontSize=12,leading=16,textColor=colors.HexColor('#101E26'),spaceBefore=12,spaceAfter=8),
 'p':ParagraphStyle('p',fontName='D',fontSize=9.1,leading=13.3,spaceAfter=8),
 'small':ParagraphStyle('small',fontName='D',fontSize=7.6,leading=10.7,spaceAfter=6),
 'cell':ParagraphStyle('cell',fontName='D',fontSize=7.3,leading=10.2),
 'head':ParagraphStyle('head',fontName='DB',fontSize=7.3,leading=10.2,textColor=colors.white),
 'lead':ParagraphStyle('lead',fontName='D',fontSize=11,leading=16,textColor=colors.HexColor(BLUE),spaceAfter=12)}
 def P(s,style='p'):
  # Typography only: preserve identifiers such as 2Z5X, KEAP1, and SHA256.
  s=re.sub(r'(\d)(?=(?:compostos|estruturas|entradas|poses|sistemas|alvos|primeiras|sementes|tentativas|conformações|CPUs|inputs|kcal|Å|anos|resíduos)\b)',r'\1 ',s)
  for a,b in [('emCPU','em CPU'),('semGPU','sem GPU'),('comV3','com V3'),('resultadosV3','resultados V3'),('bibliotecaV3','biblioteca V3'),('daV3','da V3'),('aoTCC','ao TCC'),('doTCC','do TCC'),('deGC-MS','de GC-MS'),('deSAG','de SAG'),('deMAO','de MAO'),('EmMAO','Em MAO'),('emMAO','em MAO'),('oCCD','o CCD'),('químicaCCD','química CCD'),('conectividadeCCD','conectividade CCD'),('porEnergy','por Energy'),('porSHA256','por SHA256'),('emCPU','em CPU'),('modeloCNN','modelo CNN'),('ensembleCNN','ensemble CNN'),('sítioNRF2','sítio NRF2'),('deNRF2','de NRF2'),('comNRF2','com NRF2'),('deCys151','de Cys151'),('a3,02','a 3,02'),('ácido16','ácido 16'),('a4Å','a 4 Å'),('em3D','em 3D'),('SDF3D','SDF 3D'),('emkcal','em kcal'),('noHOME','no HOME'),('pasta00','pasta 00'),('script00','script 00'),('semGPU','sem GPU'),('asTesla','as Tesla'),('nasTesla','nas Tesla'),('deC18','de C18'),('controleC18','controle C18'),('SAG/2V5Z','SAG/2V5Z'),('emMAO','em MAO'),('aoCCD','ao CCD'),('doSAG','do SAG'),('deGyrB','de GyrB'),('deE.coli','de E. coli'),('deEugenia','de Eugenia'),('deCampomanesia','de Campomanesia'),('de40','de 40'),('de59','de 59'),('as40','as 40'),('todos40','todos 40'),('nosCSVs','nos CSVs'),('uma semente101','uma semente 101'),('semente101','semente 101'),('em pelo menos4/5','em pelo menos 4/5'),('em pelo menos4','em pelo menos 4'),('até2','até 2'),('semH','sem H'),('doFAD','do FAD'),('oFAD','o FAD'),('de2','de 2'),('de5','de 5'),('de20','de 20'),('de24','de 24'),('de2026','de 2026'),('emCPU','em CPU'),('V3em','V3 em'),('a14','a 14'),('NRF2-binding','NRF2-binding')]:s=s.replace(a,b)
  s=re.sub(r'\b(com|de|em|por|para|até|a|e|os|as|todos|todas|mantém|somam|foram|inclui|contém|são|sementes?|RMSD|pH|GNINA)(?=\d)',r'\1 ',s)
  s=re.sub(r'(\d)(?=(?:satisfizeram|aprovados|reprovados|novos|resíduos|kDa)\b)',r'\1 ',s)
  for a,b in [('eTyrS','e TyrS'),('deS.aureus','de S. aureus'),('etal.','et al.'),('RMSD4','RMSD 4')]:s=s.replace(a,b)
  s=re.sub(r'\b(mediana|marca|cerca|translação|alternativa|selecionada|ausentes|Exhaustiveness|Sementes|Os|Todos)(?=\d)',r'\1 ',s,flags=re.I)
  for term in ['Energy','CNN','CCD','SDF','ADMET','RMSD','PyMOL','GC-MS','SHA256','PNG','SVG','HOME','V3','V100','TCC','IA','AChE','SB284485','C18','Arg380','R354D','CZ','vdW']:
   s=re.sub(r'(?<=[a-zÀ-ÿ])(?='+re.escape(term)+r')',' ',s)
  for a,b in [('Energytop1','Energy top1'),('CNNfast','CNN fast'),('resolverCampomanesia','resolver Campomanesia'),(',323/420',', 323/420')]:s=s.replace(a,b)
  s=re.sub(r'(\d)(?=processos\b)',r'\1 ',s)
  s=re.sub(r'(?<=[A-Za-zÀ-ÿ]),(?=[A-Za-zÀ-ÿ0-9])',', ',s)
  s=re.sub(r'([:;])(?=[A-Za-zÀ-ÿ0-9])',r'\1 ',s)
  return Paragraph(s,styles[style])
 story=[]
 def add(s,style='p'):story.append(P(s,style))
 def table(headers,body,widths):
  vv=[[P(esc(h),'head') for h in headers]]+[[P(esc(v),'cell') for v in row] for row in body]
  q=Table(vv,colWidths=widths,repeatRows=1,hAlign='LEFT');q.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor(BLUE)),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#F2F5F6')]),('LINEBELOW',(0,0),(-1,0),.5,colors.HexColor(BLUE))]));story.append(q);story.append(Spacer(1,8))
 def page():story.append(PageBreak())
 def image(name,height=None):
  p=FIG/(name+'.png')
  if p.exists():
   from PIL import Image as PILImage
   iw,ih=PILImage.open(p).size;story.append(Image(str(p),width=PW,height=PW*ih/iw if height is None else height));story.append(Spacer(1,8))
 add('Restinga | GNINA V4','lead');add('Ampliação de alvos<br/>e controles de redocking','title')
 add('Estudo exploratório computacional | 25 de outubro de2026 | Execução emCPU','small')
 add(f"<b>{counts['new_targets_tested']} novos sistemas testados; {counts['new_targets_primary_redocking_pass']} aprovados para triagem.</b> A biblioteca mantém40compostos de literatura,51estruturas estereoquímicas e59entradas conformacionais verificadas.",'lead')
 table(['Sistema novo','PDB','Redocking','RMSD mediano','Triagem principal'],[[t['name_pt'],t['pdb'],str(next(c for c in controls if c['target_id']==t['target_id'])['top1_passes'])+'/5',num(next(c for c in controls if c['target_id']==t['target_id'])['median_top1_RMSD_A'])+' Å','Sim' if next(c for c in controls if c['target_id']==t['target_id'])['screening_eligible']=='True' else 'Não'] for t in targets],[145,46,60,83,PW-334])
 add(f"<b>Novas simulações:</b> {count(counts['new_primary_screening_poses'])}poses de candidatos e{counts['new_redocking_poses']}poses de redocking. <b>Acumulado comV3:</b> {counts['all_target_systems_tested']}sistemas testados,{counts['all_target_systems_pass']}aprovados e{count(counts['all_primary_screening_poses'])}poses principais de candidatos. Os resultadosV3 foram reaproveitados com procedência registrada; não foram recalculados.")
 add('O redocking testa a reprodução geométrica de um ligante experimental. Sua aprovação não comprova atividade dos compostos vegetais. Sua reprovação impede o ranking principal naquele sistema e não demonstra inatividade biológica.')
 add('Não foram recebidos resultados deGC-MS ou de composição dos extratos deste projeto. Os40compostos são candidatos bibliográficos. Esta entrega não contém validação experimental de eficácia, segurança ou toxicologia.','small')
 page();add('O que cada alvo permite investigar','title')
 table(['Alvo e organismo','Hipótese exploratória','Limite central'],[[t['name_pt']+'\n'+t['species'],t['axis_pt'],t['caution_pt']] for t in targets],[122,142,PW-264])
 add('Seleção estrutural e cronologia','h')
 add('KEAP1/4L7B foi preterida por conterR354D. A alternativa6HWS não foi simulada: faltamCZ/NH1/NH2 daArg380, a3,02Å do ligante. Foi selecionada4XMB, humana, sem mutações declaradas e sem cadeias laterais incompletas. A decisão não dependeu dos escores dos candidatos.')
 add('EmMAO-B, oCCD deSAG/2V5Z apresenta conectividade de imina apesar do sinônimo safinamida. O controle escolhido foiC18/2V61, mantendo exatamente sua conectividadeCCD depositada, sem reconstrução pelo nome comum. As verificações químicas e a seleção estão documentadas.')
 add('Os seis sistemas iniciais foram definidos antes dos respectivos resultados. PARP1 eTyrS foram acrescentados após conhecer esses seis controles, antes de seus próprios redockings. É uma expansão exploratória documentada, sem substituir ou ocultar reprovações.','small')
 page();add('Protocolo e rastreabilidade','title')
 table(['Etapa','Definição'],[
 ['Programa','GNINA1.3.3 oficial, binário fixado porSHA256; semGPU. CNNfast com rescore; ordenação principalEnergy.'],
 ['Busca','Exhaustiveness16; até20poses por entrada;4CPUs por processo, até2processos simultâneos. Receptor rígido; busca localizada no ligante experimental.'],
 ['Controle nativo','Cinco sementes101/202/303/404/505. Aprovação: primeira poseEnergy comRMSD pesado e simetria até2Å em pelo menos4/5sementes. Melhor pose amostrada é diagnóstico separado.'],
 ['RMSD','Referencial fixo do receptor, sem alinhamento isolado. Testes com molécula idêntica e translação10Å. Não usar sobreposição ajustada do ligante para validar redocking.'],
 ['Entrada independente','QuímicaCCD e coordenadas cristalográficas para referência;64tentativasETKDGv3/MMFF com coordenadas nativas removidas para entrada. Identidade conferida após leituraSDF3D.'],
 ['Microestados','Harmina+1 eSB284485 zwitteriônico (amina+1/carboxilato-1) como hipóteses explícitas; demais controles em estadoCCD neutro. FAD oxidado com fosfatos dianiónicos emMAO-A/B;Cys-SG semH indevido. Sem cálculo populacional de pKa.'],
 ['Biblioteca','59entradas curadas, com todos40compostos e51estruturas. EntradaV3em quarentena e diagnóstico redundante excluídos. Nenhuma presença confirmada nas amostras.'],
 ['Candidatos','Sementes42/142/242 apenas nos sistemas aprovados. Variação da busca sobre a mesma proteína, sem réplicas biológicas ou intervalos de confiança.'],
 ['Auditoria','Identidade estereoquímica e grafo de cargas/ligações/hidrogênios; coordenadas; contatos<2Å; razão de distância/raiosvdW<0,70; energia positiva; fora da caixa>0,1Å; alteração de ligação mapeada>0,15Å; sobreposição intramolecular grave.']
 ],[105,PW-105])
 add('O modeloCNNfast é igual ao daV3. Escores de proteínas distintas não são escala comum de afinidade ou seletividade. Comparações com as primeiras entregas que usaram outro ensembleCNN não são equivalentes. Não se auditou sobreposição com dados de treinamento do CNN; estes controles internos não medem generalização independente. O pacote registra comandos, sementes, hashes e saídas brutas.','small')
 add('Situação acumulada, incluindo os controles anteriores','h')
 table(['Alvo anterior','Redocking','RMSD mediano','Origem'],[[c['target_id'],c['top1_passes']+'/5',num(c['median_top1_RMSD_A'])+' Å','V3; sem recálculo'] for c in controls if c['origin'].startswith('PREVIOUS')],[125,80,100,PW-305])
 page();add('Redocking: todos os controles novos','title');image('01_redocking_v4')
 add('Cada ponto representa a primeira pose de uma semente, ordenada porEnergy: verde até2Å, marrom acima de2Å. A linha azul curta indica a mediana; a linha tracejada marca2Å. O resultado principal não pode ser substituído pela melhor pose encontrada em qualquer posição do ranking.','small')
 table(['Alvo','Passes/5','Mediana top1','Melhor amostrada: mediana','Estado'],[[t['target_id'].split('_')[0],str(json.loads((ROOT/'results'/t['target_id']/'validation.json').read_text())['top1_passes']),num(json.loads((ROOT/'results'/t['target_id']/'validation.json').read_text())['median_top1_RMSD_A'])+' Å',num(statistics.median(float(x['best_sampled_RMSD_A']) for x in rows(ROOT/'results'/t['target_id']/'redocking/runs.csv')))+' Å',json.loads((ROOT/'results'/t['target_id']/'validation.json').read_text())['geometry_status']] for t in targets],[70,65,90,160,PW-385])
 add('Uma pose próxima ao cristal amostrada em posição inferior indica que a busca encontrou uma geometria compatível, mas a ordenação principal não a priorizou. O diagnóstico de reordenaçãoCNN usa apenas o conjunto exportado porEnergy e não constitui outra busca nem aprovação alternativa.','small')
 page();add('Coordenadas experimentais e redocking','title');image('02_native_coordinates_v4')
 add('Azul: ligante cristalográfico. Marrom: primeira poseEnergy da semente101. Visualização por projeção comum em duas dimensões das coordenadas originais; não houve ajuste ou alinhamento isolado dos ligantes. OsRMSDs foram calculados em3D. A figura mostra uma semente fixa, escolhida para apresentação consistente, e não a melhor das cinco.','small')
 add('O arquivo native_overlay_coordinates.csv permite reproduzir os painéis. Cada alvo também inclui um scriptPyMOL para inspecionar referência e poses no mesmo referencial. A figura é derivada das coordenadas efetivas, sem imagens geradas porIA.')
 for t in targets:
  page();tid=t['target_id'];dst=ROOT/'results'/tid;v=json.loads((dst/'validation.json').read_text());prep=json.loads((ROOT/'targets'/tid/'preparation_audit.json').read_text());cfg=json.loads((ROOT/'targets'/tid/'config.json').read_text())
  add(t['name_pt'],'title');add(t['pdb']+' | '+t['species']+' | '+t['axis_pt'],'lead')
  add('<b>Referência:</b> '+esc(t['ligand'])+'. Resolução:'+num(t['resolution_A'])+'Å. Cadeia:'+cfg['chain']+'. '+esc(t['caution_pt']))
  add(f"<b>Redocking: {v['top1_passes']}/5; mediana{num(v['median_top1_RMSD_A'])}Å. {'Aprovado para triagem exploratória.' if v['screening_eligible'] else 'Reprovado: candidatos não triados neste sistema.'}</b>")
  if tid=='TYRS_1JIL':add('A semente101 apresentouRMSD4,36Å; as demais ficaram próximas de0,62Å. A aprovação4/5 inclui uma falha e não equivale a acerto em todas as repetições.','small')
  runs=rows(dst/'redocking/runs.csv')
  if v['screening_eligible']:
   add('RMSD top1 por semente: '+ '; '.join(x['seed']+': '+num(x['top1_RMSD_A'])+' Å' for x in runs)+'.','small')
  else:
   table(['Semente','RMSD top1','Melhor amostrada','Top1 até2Å'],[[x['seed'],num(x['top1_RMSD_A'])+' Å',num(x['best_sampled_RMSD_A'])+' Å','Sim' if x['top1_within_2A']=='True' else 'Não'] for x in runs],[78,115,160,PW-353])
  if v['screening_eligible']:
   ac=json.loads((dst/'audit_counts.json').read_text());pp=[p for p in parents if p['target_id']==tid]
   add(f"Triagem:{count(ac['docking_poses'])}poses,{ac['top1_observations']}primeiras poses. Alertas geométricos em{ac['top1_with_flags']}primeiras poses;{ac['input_conformers_all_3_seeds_without_flags']}/59entradas sem alertas nas três sementes. RMSD máximo entre sementes até2Å em{ac['input_conformers_pairwise_RMSD_le_2A']}/59entradas.",'small')
   # Fixed markers inherited from V3, not selected by current favorable results.
   markers=['alpha-humulene','beta-caryophyllene','alpha-guaiene','alpha-bulnesene','beta-selinene','linalool','geranial','carvone','eucalyptol','elemol','quercetin','rotundic_acid']
   chosen=[]
   for marker in markers:
    matches=[p for p in pp if p['compound']==marker or p['compound'].replace(' ','_').replace('-','_')==marker.replace('-','_')]
    if matches:chosen.append(matches[0])
   add('Marcadores fixos da rodada anterior','h')
   table(['Composto','Energy mediana','Intervalo observado','Top1 sem alerta'],[[x['compound'],num(x['Energy_median_kcal_mol']),num(x['Energy_min_kcal_mol'])+' a '+num(x['Energy_max_kcal_mol']),x['observations_without_limited_geometry_flags']+'/'+x['seed_conformer_observations']] for x in chosen],[150,100,140,PW-390])
   add('Energy emkcal/mol é escore empírico. Mediana e intervalo abrangem sementes, estereoisômeros e conformações disponíveis; não são estimativas de afinidade experimental, potência ou incerteza estatística. A coluna de alertas não foi usada para ocultar resultados desfavoráveis. Todos40compostos constam nosCSVs; all_compound_summaries_with_baseline.csv também reúne os resultados anteriores, identificados por origem.','small')
  else:
   sampled_near=sum(float(x['best_sampled_RMSD_A'])<=2 for x in runs)
   add(f'No conjunto exportado, {sampled_near}/5 sementes amostraram pelo menos uma pose a até2Å do cristal. Esse diagnóstico não altera a reprovação da primeira pose porEnergy.')
   add('O controle nativo não satisfez o critério definido antes da busca. Não foram geradas poses ou ranking principal de candidatos. A falha pode envolver busca, ordenação, estado conformacional ou simplificações da preparação; os dados não isolam uma causa única.')
   add('Não se conclui que os compostos não interagem com o alvo. Uma tentativa futura de aprimorar esse protocolo deve permanecer separada e repetir o controle antes de avaliar candidatos. Mudar o critério após ver estes resultados não validaria esta rodada.')
  if prep['missing_standard_heavy_atoms_added']:add('Cadeias laterais completadas: '+ '; '.join(x['name']+x['residue']+' (distância dos átomos existentes ao ligante: '+num(x['distance_existing_atoms_to_native_A'],1)+' Å)' for x in prep['missing_standard_heavy_atoms_added'])+'.','small')
  if cfg['near_site_missing_residue_warning']:add('<b>Limitação estrutural:</b> segmentos não resolvidos próximos do ligante; sem reconstrução de alças. A aprovação geométrica, se houver, é restrita a esta construção incompleta.','small')
  add('Fonte estrutural:'+esc(t['url'])+' | DOI:'+esc(t['structural_doi']),'small')
 page();add('Como interpretar e priorizar os próximos passos','title')
 add('O ganho desta ampliação é a cobertura de mecanismos e a documentação dos sistemas que suportam uma triagem geométrica sob parâmetros definidos. Nenhum resultado permite afirmar qual espécie é mais potente, qual óleo será eficaz ou qual substância é segura.')
 table(['Tema','Conclusão permitida','Evidência ainda necessária'],[
 ['Anti-inflamatório','Hipóteses de interação com alvos aprovados, consideradas as construções utilizadas.','Composição confirmada; ensaio enzimático e celular compatível com o alvo; controles e exposição.'],
 ['Neurotransmissão','AChE teve controle aprovado na V3; MAO-A/B reprovaram nesta ampliação.','Inibição funcional, seletividade, acesso cerebral, efeitos comportamentais e segurança; docking não demonstra antidepressão.'],
 ['Antioxidante','KEAP1 explora uma interação regulatória não covalente.','Competição comNRF2, ativação de genes e ensaios celulares apropriados;DPPH/ABTS não confirmam esse mecanismo.'],
 ['Antitumoral','EGFR reprovado não gera candidatos principais; PARP1 aprovado sustenta apenas hipótese de interação no domínio catalítico.','Alvo funcional, viabilidade e seletividade tumoral, contexto de sinalização e relação concentração-resposta.'],
 ['Antibacteriano','GyrB representa um domínio deE.coli eTyrS uma construção catalítica deS.aureus. Somente protocolos aprovados recebem candidatos.','CIM/atividade celular, permeabilidade e efluxo; a espécie bacteriana e a concentração importam.'],
 ['Toxicologia','Descritores físico-químicos calculados daV3 podem contextualizar propriedades moleculares.','Não houve previsãoADMET validada nem ensaio toxicológico; estabelecer exposição, matriz e margem de segurança.']
 ],[90,184,PW-274])
 add('Antes de escolher candidatos para oTCC e a dissertação, integrarGC-MS/índices de retenção/áreas e composição dos extratos, confirmar vouchers e resolverCampomanesia littoralis e o óleo deEugenia catharinae. Comparadores de extratos não devem ser apresentados como constituintes voláteis comprovados.')
 add('Para comunicação ao professor','h')
 add(f"A análise foi ampliada de quatro para{counts['all_target_systems_tested']}sistemas estruturais, mantendo a biblioteca curada de40compostos. Nesta ampliação,{counts['new_targets_tested']}novos alvos passaram por controles de redocking e{counts['new_targets_primary_redocking_pass']}satisfizeram o critério geométrico, recebendo a triagem em três sementes. Os controles reprovados foram preservados e excluídos da análise principal de candidatos. Os resultados orientam hipóteses para validação e dependem da confirmação química das amostras; não demonstram eficácia ou segurança.")
 page();add('Fontes e reprodução','title')
 refs=[('GNINA1.3.3: versão executada','https://github.com/gnina/gnina/releases/tag/v1.3.3'),('GNINA1.3: McNutt et al.(2025)','https://doi.org/10.1186/s13321-025-00973-x'),('AutoDock Vina: Trott e Olson(2010)','https://doi.org/10.1002/jcc.21334'),('ETKDG: Wang et al.(2020)','https://doi.org/10.1021/acs.jcim.0c00025'),('MAO-A:Son etal.(2008)','https://doi.org/10.1073/pnas.0710626105'),('MAO-B:Binda etal.(2007)','https://doi.org/10.1021/jm070677y'),('sEH:Tanaka etal.(2011)','https://doi.org/10.1021/jm101273e'),('sEH/inflamação:Schmelzer etal.(2005)','https://doi.org/10.1073/pnas.0503279102'),('KEAP1:Jain etal.(2015)','https://doi.org/10.1016/j.ejmech.2015.08.049'),('EGFR:Stamos etal.(2002)','https://doi.org/10.1074/jbc.M207135200'),('GyrB:Lafitte etal.(2002)','https://doi.org/10.1021/bi0159837'),('PARP1:Ryan etal.(2021)','https://doi.org/10.1074/jbc.RA120.016573'),('TyrS:Qiu etal.(2001)','https://doi.org/10.1110/ps.18001'),('KEAP1/4L7B:fonte estrutural preterida','https://doi.org/10.1002/cmdc.201300525'),('IdentidadeCCD doSAG:fonte verificada na seleção','https://www.rcsb.org/ligand/SAG')]
 for label,url in refs:add('<b>'+esc(label)+'</b> | '+esc(url),'small')
 add('As fontes de composição vegetal são as mesmas da bibliotecaV3 e estão emliterature_sources.json e literature_species_links.csv. Elas justificam candidatos de literatura, não identificação nas amostras deste projeto. Os artigos estruturais sustentam seus próprios controles; não demonstram atividade dos40compostos contra esses alvos.','small')
 add('Conteúdo reprodutível','h')
 add('O pacote inclui todas as posesSDF novas, inclusive controles reprovados; receptores; referências; inputs; comandos e logs; parâmetros e sementes; identidades químicas; auditorias; tabelas completas por composto/estrutura/entrada; figurasPNG/SVG; scripts e inventárioSHA256. As síntesesV3 são identificadas como anteriores e não recalculadas.')
 add('Extraia a pasta no seuHOME. O script00_preparar_ambiente.sh cria um ambiente local e obtém o mesmo binário oficial. Não instala driver, não altera kernel e não exigeGPU. Reproduza em uma pasta nova, preservando a entrega. Use os arquivos preparados para manter a definição do sistema; protocolos diferentes requerem outra pasta.','small')
 add('Proveniência e limites','h')
 add('Rodada executada neste ambiente emCPU, não na estação local ou nasTeslaV100 do usuário. O conjunto curado exclui a entrada histórica divergente de óxido de cariofileno e o diagnóstico corrigido redundante. Não foram realizadas dinâmica molecular, cálculo de energia livre, ensaio de afinidade ou avaliação toxicológica.','small')
 def footer(can,doc):
  can.saveState();can.setStrokeColor(colors.HexColor('#D2DDE2'));can.line(46,39,W-46,39);can.setFont('D',7);can.setFillColor(colors.HexColor(GRAY));can.drawString(46,26,'Restinga | GNINA V4 | resultados exploratórios');can.drawRightString(W-46,26,str(doc.page));can.restoreState()
 path=REPORTS/'Relatorio_GNINA_Expansao_Alvos_V4_25-10-2026.pdf';doc=SimpleDocTemplate(str(path),pagesize=A4,rightMargin=46,leftMargin=46,topMargin=43,bottomMargin=52,title='Restinga GNINA V4 - Expansão de alvos',author='Juan Philippe Teixeira | análise computacional')
 doc.build(story,onFirstPage=footer,onLaterPages=footer);print(path,flush=True)
if __name__=='__main__':main()
