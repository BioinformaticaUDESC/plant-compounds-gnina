"""Build reviewable summaries, native overlays and PyMOL files from audited results."""
import json,statistics
from pathlib import Path
from collections import defaultdict
from rdkit import Chem,RDLogger
import gnina_core as core
ROOT=Path(__file__).resolve().parent
ANCHORS=['alpha-humulene','beta-caryophyllene','alpha-guaiene','alpha-bulnesene','beta-selinene','linalool','geranial','carvone','eucalyptol','elemol','quercetin','rotundic acid']
def main():
    RDLogger.DisableLog('rdApp.warning');RDLogger.DisableLog('rdApp.error')
    counts=json.loads((ROOT/'reports/study_counts.json').read_text());valid=json.loads((ROOT/'reports/redocking_summary.json').read_text());tops=core.read_csv(ROOT/'reports/all_top1_by_input_seed.csv');primary_ids={r['input_id'] for r in core.read_csv(ROOT/'primary_docking_inputs.csv')}
    parents=defaultdict(list)
    for row in tops:parents[(row['target_id'],row['compound'])].append(row)
    pooled=[]
    for (tid,name),rows in sorted(parents.items()):
        ee=[float(r['minimizedAffinity_kcal_mol']) for r in rows];pooled.append({'target_id':tid,'compound':name,'panel':rows[0]['panel'],'structures':len(set(r['structure_id'] for r in rows)),'input_conformers':len(set(r['input_id'] for r in rows)),'top1_observations':len(rows),'pooled_top1_energy_median_kcal_mol':statistics.median(ee),'pooled_top1_energy_min_kcal_mol':min(ee),'pooled_top1_energy_max_kcal_mol':max(ee),'top1_observations_without_limited_flags':sum(r['quality_flags']=='no_flag_from_limited_geometry_audit' for r in rows),'interpretation':'pooled_isomer_conformer_seed_sensitivity_not_confidence_interval_or_potency'})
    core.write_csv(ROOT/'reports/parent_compound_sensitivity.csv',pooled);core.write_csv(ROOT/'reports/preselected_markers.csv',[r for r in pooled if r['compound'] in ANCHORS])
    native_data=[];quality=[]
    for v in valid:
        tid=v['target_id'];cfg=json.loads((ROOT/'targets'/tid/'config.json').read_text());dst=ROOT/'results'/tid
        ref=Chem.RemoveHs(core.load(ROOT/cfg['native_reference'])[0]);red=Chem.RemoveHs(core.load(dst/'redocking'/f'{tid}_native_seed101.sdf')[0])
        def moldata(m):return {'xyz':m.GetConformer().GetPositions().tolist(),'bonds':[[b.GetBeginAtomIdx(),b.GetEndAtomIdx()] for b in m.GetBonds()]}
        native_data.append({'target_id':tid,'reference':moldata(ref),'redock':moldata(red),'seed':101,'top1_RMSD_A':core.rmsd(ref,red),'status':v['geometry_status']})
        pml=['reinitialize',f'load {cfg["receptor"]}, receptor',f'load {cfg["native_reference"]}, reference',f'load results/{tid}/redocking/{tid}_native_seed101.sdf, redocking','hide everything','show cartoon, receptor','color gray80, receptor','show sticks, reference or redocking','color cyan, reference','color magenta, redocking','frame 1','select site, byres (receptor within 5 of reference)','show sticks, site','zoom reference, 8','bg_color white','# No alignment or fitting: common deposited receptor frame.']
        (ROOT/'targets'/tid/'view_redocking.pml').write_text('\n'.join(pml)+'\n')
        if not v['screening_eligible']:continue
        audit=[r for r in core.read_csv(dst/'pose_quality_audit.csv') if r['primary_analysis_eligible']=='True'];cnn=[r for r in audit if int(r['CNNscore_rank_within_energy_exported_pool'])==1];energy=[r for r in audit if int(r['energy_rank'])==1]
        quality.append({'target_id':tid,'top1_observations':len(energy),'Energy_top1_with_flags':sum(r['quality_flags']!='no_flag_from_limited_geometry_audit' for r in energy),'best_CNNscore_in_exported_Energy_pool_with_flags':sum(r['quality_flags']!='no_flag_from_limited_geometry_audit' for r in cnn),'note':'CNN comparison restricted to Energy-exported pool; not a second independent sampling run'})
        molecules=defaultdict(list)
        for path in sorted((dst/'docking').glob('*seed42.sdf')):
            for m in core.load(path):
                if m.GetProp('_Name') in primary_ids:molecules[m.GetProp('_Name')].append(m)
        w=Chem.SDWriter(str(dst/'top1_seed42.sdf'))
        for iid,mm in molecules.items():w.write(mm[0])
        w.close()
        if 'beta-caryophyllene_conf01' in molecules:
            w=Chem.SDWriter(str(dst/'beta_caryophyllene_conf01_seed42_top1.sdf'));w.write(molecules['beta-caryophyllene_conf01'][0]);w.close()
    (ROOT/'reports/native_overlay_coordinates.json').write_text(json.dumps(native_data,indent=2));core.write_csv(ROOT/'reports/pose_order_quality_comparison.csv',quality)
    lines=['# Relatório técnico - GNINA multialvo V3','',f"Biblioteca: {counts['compounds']} compostos de literatura, {counts['structures']} estruturas e {counts['input_conformers']} entradas conformacionais. Quatro alvos avaliados por redocking e três aprovados para a triagem. Resultados exploratórios de moléculas isoladas; nenhum composto confirmado nas amostras.",'',
           '## Redocking','', '| Alvo | Critério top1 <=2 A | RMSD mediano (A) | Situação |','|---|---:|---:|---|']
    for v in valid:lines.append(f"| {v['target_id']} | {v['top1_passes']}/5 | {v['median_top1_RMSD_A']:.3f} | {v['geometry_status']} |")
    lines += ['',f"A análise principal contém {counts['docking_poses']} poses de docking e {counts['native_redocking_poses']} de redocking. Foram exportadas {counts['all_exported_docking_poses']} poses de docking no total: {counts['quarantined_docking_poses']} ficaram em quarentena por divergência estereoquímica e {counts['diagnostic_correction_poses_outside_primary']} pertencem a uma correção metodológica mantida separada. A triagem completa dos alvos aprovados usou as sementes 42, 142 e 242, com até 20 poses por entrada e semente.",
              ('Uma entrada herdada do óxido de cariofileno apresentou geometria 3D com identidade estereoquímica diferente da estrutura declarada. Foi excluída da análise principal. A substituta correta ficou como diagnóstico separado por não ser uma conformação suficientemente distinta das já retidas. A biblioteca principal tem 59 entradas verificadas.' if counts['quarantined_docking_poses'] else 'Reprodução limpa com somente as entradas principais verificadas; não incluiu entradas históricas em quarentena ou diagnósticas.'),
              'A 5-LOX apresentou poses amostradas abaixo de 2 A, mas a primeira pose ordenada por energia falhou nas cinco sementes. O alvo foi excluído da triagem de candidatos neste protocolo, sem escolher retrospectivamente uma pose para aprovar o controle. Isso não demonstra ausência de atividade biológica.','',
              '## Cobertura botânica','', '| Espécie | Compostos selecionados | Observação |','|---|---:|---|']
    for r in core.read_csv(ROOT/'species_coverage.csv'):lines.append(f"| {r['species']} | {r['compound_count']} | {r['note']} |")
    lines += ['', 'Não houve acesso à GC-MS do projeto. E. catharinae está representada por comparadores de extratos; fragmentos indexados de sua tabela de óleo foram registrados como pistas bibliográficas e não como composição completa confirmada. C. littoralis permanece sem atribuição química.','',
              '## Método','', 'GNINA 1.3.3, CPU; CNN fast como rescoring, Energy como ordenação principal, exhaustiveness 16, até 20 poses; mesmo protocolo de pontuação por alvo no redocking e na triagem. O CNN fast foi escolhido para a triagem ampliada e seus valores não devem ser comparados numericamente aos do ensemble padrão da entrega anterior.',
              'Redocking independente da conformação cristalográfica, cinco sementes, RMSD de átomos pesados com simetria no referencial fixo da proteína, sem alinhar o ligante. Regra definida antes dos resultados: primeira pose ate 2 A em pelo menos 4/5 sementes. Identidades químicas dos controles comparadas com o CCD do RCSB, considerando a protonação de donepezila.',
              'Receptores rígidos; coordenadas pesadas depositadas preservadas; hidrogênios no modelo pH 7,4; águas e glicanos excluídos. CB2 mantém a construção mutante com fusão T4-lisozima, sem membrana, em estado ligado a antagonista. 5-LOX mantém ferro com aproximação empírica e a desordem da região que recobre o sítio, sem reconstrução.',
              'Estruturas antigas preservadas e novas estruturas preparadas com ETKDGv3/MMFF; estereoquímica indefinida enumerada. Microestados do banco mantidos, com a exceção documentada do carboxilato de ácido rotúndico. Não foi calculada a distribuição por pKa, inclusive para fenólicos.',
              'Auditoria: conservação de identidade química, coordenadas finitas, contatos não metálicos <2 A, distância inferior a 70% da soma dos raios de van der Waals, energia empírica positiva, saída da caixa >0,1 A, variação de ligação >0,15 A e sobreposição intramolecular grave entre átomos separados por mais de três ligações. Contatos com metais registrados separadamente. Ausência de alerta é uma verificação limitada, não prova de ligação.',
              'Comprimentos de ligação comparados com correspondência de átomos; RMSD entre sementes calculado sem alinhamento. As proximidades de resíduos ate 4 A são geométricas e não foram automaticamente rotuladas como ligações de hidrogênio ou mecanismos.','',
              '## O que se pode concluir','', 'A ampliação oferece hipóteses estruturais rastreáveis e permite avaliar dependência de sementes, conformações e ordenação. A falha de um controle ficou registrada e impediu a triagem naquele alvo. Não há ranking de potência entre plantas ou proteínas, nem estimativa de atividade de misturas.',
              'Os scores não comprovam atividade anti-inflamatória, antimicrobiana, anticancerígena, antioxidante ou antidepressiva. A ligação prevista a CB2 não distingue agonismo de antagonismo. Os descritores físico-químicos calculados não são ensaios nem previsões validadas de toxicidade, biodisponibilidade ou segurança.',
              'Próximos passos: integrar GC-MS/índices de retenção por matriz e espécie, confirmar vouchers, resolver composição pendente, revisar as poses e confrontar as hipóteses com ensaios específicos.','',
              '## Arquivos','', 'O PDF reúne a síntese visual e os métodos. As tabelas completas ficam em reports/ e results/<alvo>/. Identidades e fontes estão em structures.csv e literature_species_links.csv. target_rationale.json documenta as justificativas e ressalvas estruturais. Os scripts PyMOL de cada alvo usam o mesmo referencial, sem alinhamento.','']
    (ROOT/'reports/RESULTADOS_V3.md').write_text('\n'.join(lines))
    print(json.dumps(counts,indent=2))

if __name__=='__main__':main()
