# Relatório técnico - GNINA multialvo V3

Biblioteca: 40 compostos de literatura, 51 estruturas e 59 entradas conformacionais. Quatro alvos avaliados por redocking e três aprovados para a triagem. Resultados exploratórios de moléculas isoladas; nenhum composto confirmado nas amostras.

## Redocking

| Alvo | Critério top1 <=2 A | RMSD mediano (A) | Situação |
|---|---:|---:|---|
| ACHE_4EY7 | 5/5 | 1.070 | PASS |
| ALOX5_6N2W | 0/5 | 3.549 | FAIL |
| CB2_5ZTY | 5/5 | 0.537 | PASS |
| COX2_5IKR | 5/5 | 0.481 | PASS |

A análise principal contém 10620 poses de docking e 400 de redocking. Foram exportadas 10980 poses de docking no total: 180 ficaram em quarentena por divergência estereoquímica e 180 pertencem a uma correção metodológica mantida separada. A triagem completa dos alvos aprovados usou as sementes 42, 142 e 242, com até 20 poses por entrada e semente.
Uma entrada herdada do óxido de cariofileno apresentou geometria 3D com identidade estereoquímica diferente da estrutura declarada. Foi excluída da análise principal. A substituta correta ficou como diagnóstico separado por não ser uma conformação suficientemente distinta das já retidas. A biblioteca principal tem 59 entradas verificadas.
A 5-LOX apresentou poses amostradas abaixo de 2 A, mas a primeira pose ordenada por energia falhou nas cinco sementes. O alvo foi excluído da triagem de candidatos neste protocolo, sem escolher retrospectivamente uma pose para aprovar o controle. Isso não demonstra ausência de atividade biológica.

## Cobertura botânica

| Espécie | Compostos selecionados | Observação |
|---|---:|---|
| Varronia curassavica | 2 | literature_selection_not_exhaustive |
| Campomanesia littoralis | 0 | no_chemical_assignment; no_generic_C_xanthocarpa_substitution |
| Myrcia palustris | 3 | literature_selection_not_exhaustive |
| Eugenia catharinae | 8 | extract_comparators_only; volatile_composition_pending |
| Schinus terebinthifolia | 8 | literature_selection_not_exhaustive |
| Psidium guajava | 11 | literature_selection_not_exhaustive |
| Lippia alba | 19 | literature_selection_not_exhaustive |

Não houve acesso à GC-MS do projeto. E. catharinae está representada por comparadores de extratos; fragmentos indexados de sua tabela de óleo foram registrados como pistas bibliográficas e não como composição completa confirmada. C. littoralis permanece sem atribuição química.

## Método

GNINA 1.3.3, CPU; CNN fast como rescoring, Energy como ordenação principal, exhaustiveness 16, até 20 poses; mesmo protocolo de pontuação por alvo no redocking e na triagem. O CNN fast foi escolhido para a triagem ampliada e seus valores não devem ser comparados numericamente aos do ensemble padrão da entrega anterior.
Redocking independente da conformação cristalográfica, cinco sementes, RMSD de átomos pesados com simetria no referencial fixo da proteína, sem alinhar o ligante. Regra definida antes dos resultados: primeira pose ate 2 A em pelo menos 4/5 sementes. Identidades químicas dos controles comparadas com o CCD do RCSB, considerando a protonação de donepezila.
Receptores rígidos; coordenadas pesadas depositadas preservadas; hidrogênios no modelo pH 7,4; águas e glicanos excluídos. CB2 mantém a construção mutante com fusão T4-lisozima, sem membrana, em estado ligado a antagonista. 5-LOX mantém ferro com aproximação empírica e a desordem da região que recobre o sítio, sem reconstrução.
Estruturas antigas preservadas e novas estruturas preparadas com ETKDGv3/MMFF; estereoquímica indefinida enumerada. Microestados do banco mantidos, com a exceção documentada do carboxilato de ácido rotúndico. Não foi calculada a distribuição por pKa, inclusive para fenólicos.
Auditoria: conservação de identidade química, coordenadas finitas, contatos não metálicos <2 A, distância inferior a 70% da soma dos raios de van der Waals, energia empírica positiva, saída da caixa >0,1 A, variação de ligação >0,15 A e sobreposição intramolecular grave entre átomos separados por mais de três ligações. Contatos com metais registrados separadamente. Ausência de alerta é uma verificação limitada, não prova de ligação.
Comprimentos de ligação comparados com correspondência de átomos; RMSD entre sementes calculado sem alinhamento. As proximidades de resíduos ate 4 A são geométricas e não foram automaticamente rotuladas como ligações de hidrogênio ou mecanismos.

## O que se pode concluir

A ampliação oferece hipóteses estruturais rastreáveis e permite avaliar dependência de sementes, conformações e ordenação. A falha de um controle ficou registrada e impediu a triagem naquele alvo. Não há ranking de potência entre plantas ou proteínas, nem estimativa de atividade de misturas.
Os scores não comprovam atividade anti-inflamatória, antimicrobiana, anticancerígena, antioxidante ou antidepressiva. A ligação prevista a CB2 não distingue agonismo de antagonismo. Os descritores físico-químicos calculados não são ensaios nem previsões validadas de toxicidade, biodisponibilidade ou segurança.
Próximos passos: integrar GC-MS/índices de retenção por matriz e espécie, confirmar vouchers, resolver composição pendente, revisar as poses e confrontar as hipóteses com ensaios específicos.

## Arquivos

O PDF reúne a síntese visual e os métodos. As tabelas completas ficam em reports/ e results/<alvo>/. Identidades e fontes estão em structures.csv e literature_species_links.csv. target_rationale.json documenta as justificativas e ressalvas estruturais. Os scripts PyMOL de cada alvo usam o mesmo referencial, sem alinhamento.
