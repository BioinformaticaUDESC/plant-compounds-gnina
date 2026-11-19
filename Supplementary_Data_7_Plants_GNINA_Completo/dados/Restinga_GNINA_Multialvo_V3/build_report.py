"""Create a source-backed PDF and scientific figures from the completed V3 audit."""
import csv,json,math,statistics,html
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image,KeepTogether
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.pagesizes import A4

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'reports';FIG=OUT/'figures';FIG.mkdir(exist_ok=True)
NAMES={'ACHE_4EY7':'Acetilcolinesterase','ALOX5_6N2W':'5-lipoxigenase','CB2_5ZTY':'CB2','COX2_5IKR':'COX-2'}
ORDER=['COX2_5IKR','ACHE_4EY7','CB2_5ZTY','ALOX5_6N2W']
COLORS=['#224867','#3B7778','#806554','#9C4C4C']
def read(p):return list(csv.DictReader(p.open()))
def pretty(s):
    return s.replace('alpha-','α-').replace('beta-','β-').replace('gamma-','γ-').replace('delta-','δ-').replace('_',' ')
def figsave(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=240,bbox_inches='tight',facecolor='white');fig.savefig(FIG/(name+'.svg'),bbox_inches='tight');plt.close(fig)
def build_figures(native_only=False):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    validations={v['target_id']:v for v in json.loads((OUT/'redocking_summary.json').read_text())}
    fig,ax=plt.subplots(figsize=(8.1,3.6))
    for i,tid in enumerate(ORDER):
        rows=read(ROOT/'results'/tid/'redocking/runs.csv');ys=[float(r['top1_RMSD_A']) for r in rows]
        ax.scatter(np.linspace(i-.16,i+.16,5),ys,c=COLORS[i],s=38,zorder=3);ax.plot([i-.23,i+.23],[statistics.median(ys)]*2,c=COLORS[i],lw=2)
        ax.text(i,max(ys)+.20,f"{validations[tid]['top1_passes']}/5",ha='center',fontsize=10,color=COLORS[i])
    ax.axhline(2,c='#777777',ls='--',lw=1,label='Limite predefinido: 2 Å');ax.set_xticks(range(4),[NAMES[t]+'\n'+t.split('_')[1] for t in ORDER]);ax.set_ylabel('RMSD da primeira pose (Å)');ax.set_ylim(0,4.2);ax.grid(axis='y',alpha=.18);ax.legend(loc='upper left',frameon=False);fig.tight_layout();figsave(fig,'01_redocking')
    overlays={r['target_id']:r for r in json.loads((OUT/'native_overlay_coordinates.json').read_text())}
    fig,axes=plt.subplots(2,2,figsize=(8.1,6.0))
    for ax,tid in zip(axes.flat,ORDER):
        d=overlays[tid];ref=np.array(d['reference']['xyz']);pr=np.array(d['redock']['xyz']);origin=ref.mean(axis=0);_,_,vt=np.linalg.svd(ref-origin,full_matrices=False);basis=vt[:2].T
        for mol,color,style in [(d['reference'],'#224867','-'),(d['redock'],'#AF6160','--')]:
            x=(np.array(mol['xyz'])-origin)@basis
            for a,b in mol['bonds']:ax.plot(x[[a,b],0],x[[a,b],1],style,c=color,lw=1.8,alpha=.85)
        ax.set_title(f"{NAMES[tid]} | RMSD {d['top1_RMSD_A']:.2f} Å",fontsize=10);ax.set_aspect('equal');ax.set_xlabel('Projeção comum (Å)');ax.set_ylabel('Projeção comum (Å)');ax.grid(alpha=.12)
    fig.legend(handles=[Line2D([0],[0],color='#224867',lw=2,label='Referência cristalográfica'),Line2D([0],[0],color='#AF6160',ls='--',lw=2,label='Redocking, semente 101')],loc='lower center',ncol=2,frameon=False);fig.tight_layout(rect=[0,.05,1,1]);figsave(fig,'02_native_overlays')
    if native_only:return
    quality=read(OUT/'pose_order_quality_comparison.csv');active=[t for t in ORDER if validations[t]['screening_eligible']];qm={r['target_id']:r for r in quality}
    fig,ax=plt.subplots(figsize=(8.1,3.4));x=np.arange(len(active));width=.32
    ax.bar(x-width/2,[int(qm[t]['Energy_top1_with_flags']) for t in active],width,color='#3B7778',label='Primeira pose por Energy')
    ax.bar(x+width/2,[int(qm[t]['best_CNNscore_in_exported_Energy_pool_with_flags']) for t in active],width,color='#A67F69',label='Maior CNNscore no conjunto exportado')
    ax.set_xticks(x,[NAMES[t] for t in active]);ax.set_ylabel('Primeiras poses com alerta (n)');ax.grid(axis='y',alpha=.16);ax.set_axisbelow(True);ax.legend(frameon=False,fontsize=8,loc='upper right');fig.tight_layout();figsave(fig,'03_pose_quality')
    desc=read(ROOT/'physicochemical_descriptors.csv');one={}
    for r in desc:one.setdefault(r['compound'],r)
    fig,ax=plt.subplots(figsize=(8.1,3.8))
    for panel,color,label in [('volatile_literature_panel','#3B7778','Voláteis da literatura'),('nonvolatile_extract_comparator','#806554','Comparadores de extratos')]:
        rows=[r for r in one.values() if r['panel']==panel];ax.scatter([float(r['MW']) for r in rows],[float(r['cLogP_RDKit']) for r in rows],c=color,s=45,alpha=.8,label=label,edgecolors='white',linewidth=.5)
    for name in ['myricitrin','quercetin','rotundic acid','beta-sitosterol','linalool']:
        r=one[name];ax.annotate(pretty(name),(float(r['MW']),float(r['cLogP_RDKit'])),xytext=(4,5),textcoords='offset points',fontsize=7)
    ax.set_xlabel('Massa molecular calculada (g/mol)');ax.set_ylabel('cLogP calculado (RDKit)');ax.grid(alpha=.16);ax.legend(frameon=False,fontsize=8);fig.tight_layout();figsave(fig,'04_descriptors')

def main():
    build_figures()
    fonts=Path(matplotlib.get_data_path())/'fonts/ttf'
    for name,file in [('DV','DejaVuSans.ttf'),('DV-Bold','DejaVuSans-Bold.ttf'),('DV-Italic','DejaVuSans-Oblique.ttf')]:pdfmetrics.registerFont(TTFont(name,str(fonts/file)))
    pdfmetrics.registerFontFamily('DV',normal='DV',bold='DV-Bold',italic='DV-Italic',boldItalic='DV-Bold')
    styles=getSampleStyleSheet();styles.add(ParagraphStyle(name='BodyV3',fontName='DV',fontSize=9.2,leading=13,spaceAfter=7,textColor=colors.HexColor('#20242A')));styles.add(ParagraphStyle(name='TitleV3',fontName='DV-Bold',fontSize=22,leading=28,spaceAfter=14,textColor=colors.black));styles.add(ParagraphStyle(name='H1V3',fontName='DV-Bold',fontSize=15,leading=20,spaceBefore=3,spaceAfter=12,textColor=colors.black));styles.add(ParagraphStyle(name='H2V3',fontName='DV-Bold',fontSize=10.5,leading=15,spaceAfter=8,textColor=colors.black));styles.add(ParagraphStyle(name='SmallV3',fontName='DV',fontSize=7.7,leading=10.3,spaceAfter=6,textColor=colors.HexColor('#414A54')));styles.add(ParagraphStyle(name='CellV3',fontName='DV',fontSize=7.8,leading=10.5,spaceAfter=0));styles.add(ParagraphStyle(name='HeadCellV3',fontName='DV-Bold',fontSize=7.8,leading=10.5,textColor=colors.white,spaceAfter=0))
    story=[];W=A4[0]-88
    def p(text,style='BodyV3'):story.append(Paragraph(text,styles[style]))
    def h(text):p(text,'H1V3')
    def small(text):p(text,'SmallV3')
    def page():story.append(PageBreak())
    def table(headers,rows,widths):
        data=[[Paragraph(html.escape(str(x)),styles['HeadCellV3']) for x in headers]]+[[Paragraph(html.escape(str(x)),styles['CellV3']) for x in row] for row in rows]
        t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#224867')),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#224867')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#F0F3F5')]),('LINEBELOW',(0,-1),(-1,-1),.5,colors.HexColor('#D4DBDF'))]));story.append(t);story.append(Spacer(1,10))
    def picture(name,height):story.append(Image(str(FIG/(name+'.png')),width=W,height=height,kind='proportional'));story.append(Spacer(1,6))
    counts=json.loads((OUT/'study_counts.json').read_text());vs=json.loads((OUT/'redocking_summary.json').read_text());vm={v['target_id']:v for v in vs};active=[t for t in ORDER if vm[t]['screening_eligible']];coverage=read(ROOT/'species_coverage.csv');parents=read(OUT/'parent_compound_sensitivity.csv');markers=read(OUT/'preselected_markers.csv');contacts=read(OUT/'residue_contact_recurrence.csv');rationale=json.loads((ROOT/'target_rationale.json').read_text())
    p('Docking e redocking<br/>GNINA multialvo','TitleV3');p('Biblioteca de espécies da Restinga | Relatório técnico V3','H2V3');small('25 de outubro de 2026 | Estudo exploratório com compostos da literatura')
    p(f"A ampliação avaliou <b>{counts['compounds']} compostos</b>, representados por <b>{counts['structures']} estruturas químicas</b> e <b>{counts['input_conformers']} entradas conformacionais</b>. Foram testados quatro alvos por redocking; três atenderam ao critério geométrico e receberam a triagem em três sementes independentes.")
    table(['Alvo / estrutura','Redocking top1 ≤2 Å','RMSD mediano','Decisão'],[[NAMES[t]+' / '+t.split('_')[1],f"{vm[t]['top1_passes']}/5",f"{vm[t]['median_top1_RMSD_A']:.3f} Å",'Triagem executada' if vm[t]['screening_eligible'] else 'Excluído da triagem'] for t in ORDER],[W*.35,W*.22,W*.20,W*.23])
    p(f"A análise principal contém <b>{counts['docking_poses']:,}</b> poses de docking e <b>{counts['native_redocking_poses']}</b> de redocking. O critério predefinido exigiu RMSD da primeira pose até 2 Å em pelo menos quatro das cinco sementes, sem alinhar o ligante ao cristal.")
    p('O resultado mais importante do controle foi a reprovação de 5-LOX neste protocolo: houve amostragem de poses próximas ao cristal, mas elas não ocuparam a primeira posição por energia. Esse alvo não foi usado para interpretar os candidatos.')
    p('<b>Limite central:</b> não foi fornecida a GC-MS do projeto. As substâncias e suas relações com espécies são bibliográficas. Os resultados não demonstram composição das amostras, atividade farmacológica de óleos/hidrolatos ou segurança.')
    if counts['quarantined_docking_poses']:
        small('Curadoria de identidade: uma conformação herdada de óxido de cariofileno foi colocada em quarentena por divergência estereoquímica 3D. A geometria corrigida foi testada separadamente; como reproduziu uma conformação já retida, não recebeu peso adicional na análise principal. Permaneceram 59 entradas principais quimicamente verificadas.')
    else:small('Esta reprodução limpa utilizou somente a biblioteca principal quimicamente verificada. Entradas históricas em quarentena e diagnósticas não foram simuladas nesta execução.')
    small('Os arquivos CSV, SDF, protocolos, logs, fontes e scripts acompanham o pacote técnico. Não se trata de validação experimental nem de estudo clínico.')
    page();h('Cobertura botânica e matrizes')
    table(['Espécie','Compostos','Matriz e alcance'],[[r['species'],r['compound_count'],'Extratos não voláteis; óleo pendente' if r['species']=='Eugenia catharinae' else ('Sem atribuição química; GC-MS/voucher pendentes' if r['species']=='Campomanesia littoralis' else 'Constituintes voláteis selecionados da literatura')] for r in coverage],[W*.36,W*.13,W*.51])
    p('A soma por espécie inclui substâncias compartilhadas. Um mesmo composto não foi contado como uma nova molécula apenas por ocorrer em outra planta. A seleção não é um inventário exaustivo de nenhum óleo, extrato ou hidrolato.')
    structures=read(ROOT/'structures.csv');unique={}
    for s in structures:unique.setdefault(s['compound'],s)
    for panel,title in [('volatile_literature_panel','Voláteis selecionados (32 compostos)'),('nonvolatile_extract_comparator','Comparadores de extratos (8 compostos)')]:
        p(title,'H2V3');p('; '.join(html.escape(pretty(n)) for n,s in sorted(unique.items()) if s['panel']==panel)+'.','SmallV3')
    p('Para <i>E. catharinae</i>, o artigo de óleo de 2023 foi localizado, inclusive um fragmento indexado de sua tabela; a composição completa não foi acessada. O painel executado continua identificado como comparadores de extratos. Para <i>C. littoralis</i>, não houve lista nominal suficientemente confirmada para atribuição; a química genérica de <i>C. xanthocarpa</i> não foi transferida automaticamente.')
    small('Fontes e matrizes: literature_sources.json e literature_species_links.csv. Os nomes, CIDs, SMILES, InChIKeys e estados de carga estão em structures.csv.')
    page();h('Preparação e protocolo')
    table(['Componente','Procedimento'],[
        ['Software e pontuação','GNINA 1.3.3; CPU; CNN fast como rescoring; ordenação principal por Energy; exhaustiveness 16; até 20 poses por entrada.'],
        ['Repetições','Docking: sementes 42, 142 e 242 em toda a biblioteca dos alvos aprovados. Redocking: 101, 202, 303, 404 e 505.'],
        ['Identidade dos controles','Estruturas do RCSB/CCD; referências químicas conferidas por InChIKey, considerando a protonação. Conformação inicial de redocking gerada independentemente do cristal.'],
        ['Ligantes','Entradas anteriores preservadas; novas estruturas com ETKDGv3 e MMFF, até 48 conformações iniciais. Posições estereoquímicas indefinidas enumeradas. Anéis grandes herdados mantêm até três conformações selecionadas.'],
        ['Microestados','Estados do banco mantidos; ácido rotúndico como carboxilato. Não foi feita previsão de pKa ou distribuição de microestados.'],
        ['Receptores','Rígidos, sem relaxação dos átomos pesados depositados ou reconstrução de alças; hidrogênios no modelo pH 7,4; águas/glicanos/solventes excluídos.'],
        ['RMSD','Átomos pesados, simetria, mesmo referencial do receptor; nenhum alinhamento do ligante. Controles de RMSD zero e translação de 10 Å verificados.']],[W*.26,W*.74])
    p('O modelo CNN fast reduz o custo do rescoring da biblioteca ampliada. Seus valores não devem ser comparados numericamente aos do ensemble padrão da entrega anterior. A CNN atua após a busca neste protocolo; a energia empírica foi escolhida como ordenação principal antes dos resultados desta rodada.')
    p('As repetições medem a variação da busca sob um receptor e uma preparação fixos. Elas não representam réplicas experimentais, diferentes conformações da proteína ou incerteza completa da afinidade.')
    if counts['quarantined_docking_poses']:small('A identidade foi conferida após escrever e reler a geometria 3D. A entrada herdada divergente foi excluída; a substituta correta ficou como diagnóstico por não constituir uma terceira conformação suficientemente distinta. Saídas originais, curadoria e diagnóstico acompanham o pacote.')
    else:small('A biblioteca usada na reprodução limpa é a biblioteca principal curada, com identidades verificadas após escrever e reler as geometrias 3D.')
    small('Rastreabilidade: protocol.json, native_identity_checks.json, auditorias de preparação por alvo e registros dos comandos/saídas. O executável é fixado por SHA256.')
    page();h('Redocking: critério e resultado')
    picture('01_redocking',235)
    small('Figura 1. RMSD tridimensional da primeira pose de cada semente, sem alinhamento. Pontos: cinco sementes; traço: mediana. A linha pontilhada indica o limiar predefinido de 2 Å. O alvo precisava passar em pelo menos 4/5 sementes.')
    rows=[]
    for t in ORDER:
        runs=read(ROOT/'results'/t/'redocking/runs.csv');first=[float(r['top1_RMSD_A']) for r in runs];best=[float(r['best_sampled_RMSD_A']) for r in runs]
        rows.append([NAMES[t],f'{min(first):.2f} - {max(first):.2f}',f'{min(best):.2f} - {max(best):.2f}',vm[t]['geometry_status']])
    table(['Alvo','Faixa RMSD top1 (Å)','Melhor pose amostrada (Å)','Critério'],rows,[W*.26,W*.26,W*.29,W*.19])
    p('Na 5-LOX, poses abaixo de 2 Å apareceram entre as soluções exportadas, mas nenhuma foi a primeira por energia. Portanto, a falha observada envolve a priorização das poses sob este protocolo; não autoriza concluir que o alvo ou qualquer composto seja biologicamente inativo. O resultado reprovado foi preservado e não foi convertido em aprovação escolhendo retrospectivamente outra pose.')
    p('A aprovação dos demais controles demonstra recuperação geométrica dos respectivos ligantes cristalográficos. Ela não valida automaticamente a pose, a afinidade ou o mecanismo dos demais esqueletos químicos.')
    page();h('Sobreposição das coordenadas nativas')
    picture('02_native_overlays',375)
    small('Figura 2. Ligações dos ligantes de referência e da primeira pose de redocking da semente 101. As coordenadas vêm dos arquivos estruturais e de saída do GNINA. Cada par usa a mesma projeção linear em dois eixos do cristal; não houve ajuste entre ligantes. A projeção pode ocultar diferenças fora do plano; os RMSDs indicados foram calculados em 3D.')
    p('No PyMOL, os scripts view_redocking.pml permitem examinar cada controle junto ao receptor no referencial original. A inspeção química visual completa de todas as poses da biblioteca permanece necessária.')
    page();h('Auditoria geométrica e repetibilidade')
    qr=[]
    for t in active:
        d=json.loads((ROOT/'results'/t/'audit_counts.json').read_text());qr.append([NAMES[t],f"{d['top1_with_flags']}/{d['top1_observations']}",f"{d['input_conformers_all_3_seeds_without_flags']}/{d['input_conformers']}",f"{d['input_conformers_pairwise_RMSD_le_2A']}/{d['input_conformers']}"])
    table(['Alvo','Top1 com alerta','Entradas sem alerta nas 3 sementes','Entradas com RMSD entre sementes ≤2 Å'],qr,[W*.23,W*.20,W*.28,W*.29])
    picture('03_pose_quality',210)
    small('Figura 3. Alertas das primeiras poses por energia e das poses com maior CNNscore dentro do conjunto exportado por energia. A segunda série é uma comparação restrita ao mesmo conjunto de poses; não corresponde a uma nova triagem com ordenação global por CNNscore.')
    p('Foram verificados contatos não metálicos abaixo de 2 Å, distância menor que 70% da soma de raios de van der Waals, energia positiva, saída da caixa acima de 0,1 Å, alteração de ligação acima de 0,15 Å e sobreposição intramolecular grave entre átomos separados por mais de três ligações. Contatos com metais foram registrados separadamente.')
    p('O número de entradas com RMSD entre sementes até 2 Å é descritivo. Nenhum desses indicadores comprova afinidade, efeito farmacológico ou validade completa da pose. As proximidades de resíduos até 4 Å também são geométricas, sem classificação automática de ligação de hidrogênio.')
    for tid in active:
        page();h('Sensibilidade dos marcadores: '+NAMES[tid])
        source=rationale['targets'][tid];p(source['scope'])
        mm=[r for r in markers if r['target_id']==tid];rows=[]
        for r in mm:
            rows.append([pretty(r['compound']),r['top1_observations'],f"{float(r['pooled_top1_energy_median_kcal_mol']):.2f}",f"{float(r['pooled_top1_energy_min_kcal_mol']):.2f} a {float(r['pooled_top1_energy_max_kcal_mol']):.2f}",f"{r['top1_observations_without_limited_flags']}/{r['top1_observations']}"])
        table(['Marcador','N','Mediana energia','Mín. a máx. (kcal/mol)','Sem alerta'],rows,[W*.30,W*.08,W*.20,W*.24,W*.18])
        small('Os 12 marcadores foram definidos antes de examinar seus resultados. Cada linha reúne primeiras poses de sementes, isômeros e conformações do mesmo composto. N é o número de observações agrupadas. A faixa mostra sensibilidade dessas escolhas; não é intervalo de confiança, potência ou afinidade experimental.')
        p('Proximidades recorrentes nos marcadores selecionados','H2V3')
        cr=[]
        for name in ['beta-caryophyllene','linalool','quercetin']:
            candidates=sorted(set(r['input_id'] for r in contacts if r['target_id']==tid and r['compound']==name))
            if not candidates:continue
            iid=candidates[0];rr=[r for r in contacts if r['target_id']==tid and r['input_id']==iid and int(r['seeds_with_proximity_le_4A'])==3];rr=sorted(rr,key=lambda r:float(r['minimum_distance_A']))[:5]
            cr.append([pretty(name),' / '.join(r['residue_name']+r['residue_id']+'('+r['receptor_chain']+')' for r in rr) or 'Nenhum resíduo em 3/3 sementes'])
        table(['Marcador','Resíduos a ≤4 Å em 3/3 sementes'],cr,[W*.30,W*.70])
        small('Primeira entrada em ordem alfabética de cada marcador, escolhida sem usar o score. Recorrência de proximidade não é probabilidade de ligação. Os resíduos usam a numeração e a cadeia do PDB preparado, sem conversão para UniProt. A tabela completa acompanha o pacote.')
        p('<b>Limites específicos:</b> '+html.escape(source['limitations']))
    page();h('Descritores físico-químicos calculados')
    picture('04_descriptors',245)
    small('Figura 4. Uma representação por composto para evitar duplicar isômeros com os mesmos descritores. Massa e cLogP foram calculados por RDKit a partir das estruturas preparadas. Os pontos distinguem voláteis bibliográficos de comparadores de extratos; não representam biodisponibilidade, permeação cerebral ou toxicidade medidas.')
    desc=read(ROOT/'physicochemical_descriptors.csv');one={}
    for r in desc:one.setdefault(r['compound'],r)
    rows=[]
    for panel,title in [('volatile_literature_panel','Voláteis (32)'),('nonvolatile_extract_comparator','Extratos (8)')]:
        rr=[r for r in one.values() if r['panel']==panel]
        fmt=lambda k:f"{min(float(r[k]) for r in rr):.1f} - {max(float(r[k]) for r in rr):.1f}"
        rows.append([title,fmt('MW'),fmt('cLogP_RDKit'),fmt('TPSA_A2')])
    table(['Painel','Massa (g/mol)','cLogP','TPSA (Å²)'],rows,[W*.28,W*.25,W*.20,W*.27])
    p('A separação por matriz importa: esteróis, triterpenos e flavonoides descritos em extratos não foram apresentados como constituintes do óleo ou do hidrolato. Sem a composição experimental e sem quantificação, não é possível converter o docking das moléculas em resposta de uma mistura.')
    p('Os microestados foram definidos por uma política simples de preparação. A ionização de fenólicos e outras formas tautoméricas não foi explorada. Assim, até a comparação entre moléculas isoladas depende de uma hipótese química que precisa ser refinada antes de priorizar candidatos.')
    page();h('Interpretação para o projeto')
    p('O que esta rodada entrega','H2V3')
    p('Uma biblioteca química com fontes e matrizes rastreáveis, três alvos com controle geométrico aprovado, um controle reprovado preservado, repetição integral da busca em três sementes, auditoria das poses e tabelas de contatos e descritores. Esses dados ajudam a formular e revisar hipóteses estruturais para as moléculas isoladas.')
    p('O que os resultados ainda não estabelecem','H2V3')
    p('Não há demonstração de atividade anti-inflamatória, antimicrobiana, anticancerígena, antioxidante ou antidepressiva das amostras. Scores de proteínas diferentes não são uma escala comum de afinidade. A ligação prevista a CB2 não distingue agonismo de antagonismo, e a aproximação a AChE não equivale a efeito clínico ou antidepressivo.')
    p('As três sementes variam a busca sob uma estrutura rígida; não exploram a dinâmica completa do receptor. O redocking avalia apenas a recuperação do ligante de referência. Ausência de alerta geométrico não prova ligação, seletividade ou segurança, e os descritores calculados não substituem toxicologia.')
    p('Dados necessários para a próxima decisão','H2V3')
    table(['Dado','Finalidade'],[
        ['GC-MS por espécie e matriz','Confirmar quais moléculas pertencem efetivamente a cada óleo/hidrolato, com índices de retenção, identificação e áreas relativas.'],
        ['Voucher e composição pendente','Resolver a atribuição de C. littoralis e obter a tabela completa do óleo de E. catharinae.'],
        ['Revisão estrutural e microestados','Inspecionar poses, conformações, protonação e ambiente do sítio antes de selecionar hipóteses.'],
        ['Validação experimental específica','Confrontar as hipóteses com ensaios adequados ao alvo e à matriz, incluindo seletividade e toxicidade.']],[W*.33,W*.67])
    p('Texto breve para discussão com o orientador','H2V3')
    p(f"A triagem foi ampliada para 40 compostos de literatura e quatro sistemas de referência. COX-2, acetilcolinesterase e CB2 atenderam ao critério de redocking em cinco sementes; 5-LOX não atendeu e ficou fora da triagem. Nos três alvos aprovados, toda a biblioteca foi analisada em três sementes, com controle geométrico das poses e registro das fontes. Os resultados sustentam hipóteses estruturais, ainda dependentes de confirmação da composição das amostras e de validação experimental.")
    page();h('Referências e rastreabilidade')
    p('Fontes da biblioteca','H2V3');sources=json.loads((ROOT/'literature_sources.json').read_text())['sources']
    for sid,s in sources.items():small(f"<b>[{sid}]</b> {html.escape(s['title'])}. <link href=\"{html.escape(s['url'],quote=True)}\" color=\"#224867\">{html.escape(s['url'])}</link>")
    p('Estruturas e justificativas dos alvos','H2V3')
    for tid in ORDER:
        s=rationale['targets'][tid];small(f"<b>{NAMES[tid]}:</b> {html.escape(s['structure_url'])}; base biológica: {html.escape(s['evidence_url'])}.")
    small('Estruturas adicionais: AChE/donepezila, DOI 10.1021/jm300871x; CB2/AM10257, DOI 10.1016/j.cell.2018.12.011; 5-LOX/NDGA, DOI 10.1038/s41589-020-0544-7.')
    small('O artigo de 5-LOX tem nota editorial sobre a estrutura companheira ligada a AKBA (6NCF), que não foi utilizada aqui. O limite operacional desta rodada para 6N2W foi a reprovação do redocking, registrada integralmente.')
    p('Métodos e arquivos reproduzíveis','H2V3')
    small('GNINA: https://github.com/gnina/gnina. RDKit: https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html. O pacote inclui entradas, saídas, logs, hashes, código, protocolo, tabelas completas e figuras PNG/SVG. Um inventário SHA256 permite conferir a integridade da entrega.')
    small('As referências fundamentam a seleção de constituintes e hipóteses; as simulações e estatísticas apresentadas neste relatório são resultados computacionais desta execução. Nenhuma identidade química foi confirmada nas amostras do projeto.')
    pdf=OUT/'Relatorio_GNINA_Multialvo_V3_25-10-2026.pdf'
    def decorate(c,doc):
        c.saveState();c.setStrokeColor(colors.HexColor('#B7C4CD'));c.setLineWidth(.5);c.line(44,A4[1]-36,A4[0]-44,A4[1]-36);c.setFont('DV',7);c.setFillColor(colors.HexColor('#51606A'));c.drawString(44,A4[1]-28,'RESTINGA | GNINA MULTIALVO V3');c.drawString(44,27,'Relatório exploratório | 25/10/2026');c.drawRightString(A4[0]-44,27,f'Página {doc.page}');c.restoreState()
    doc=SimpleDocTemplate(str(pdf),pagesize=A4,leftMargin=44,rightMargin=44,topMargin=52,bottomMargin=44,title='Docking e redocking GNINA multialvo V3',author='Projeto Restinga',pageCompression=1);doc.build(story,onFirstPage=decorate,onLaterPages=decorate);print(pdf.resolve())

if __name__=='__main__':main()
