# Restinga: ampliação de alvos com GNINA V4

A etapa A amplia seis sistemas: MAO-A/2Z5X, MAO-B/2V61, epóxido-hidrolase solúvel/3ANS, KEAP1/4XMB, EGFR/1M17 e GyrB de E. coli/1KZN. A etapa B acrescenta PARP1/7KK4 e tirosil-tRNA sintetase de S. aureus/1JIL. A etapa B foi definida depois de conhecer os seis controles iniciais e antes de seus próprios resultados: expansão exploratória documentada, sem substituir reprovações.

COX-2, acetilcolinesterase, CB2 e 5-LOX permanecem identificados como resultados anteriores V3, sem recálculo. Consulte o PDF e `reports/RESULTADOS_V4.md` para os resultados finais, incluindo todos os controles reprovados.

## Biblioteca e limites

São 40 compostos de literatura, 51 estruturas estereoquímicas e 59 entradas conformacionais verificadas. Todos os arquivos de entrada correspondem aos hashes registrados. A conformação histórica com divergência estereoquímica e o diagnóstico corrigido redundante de óxido de cariofileno não integram a biblioteca V4.

Nenhuma análise química das amostras do projeto foi fornecida. As relações entre substâncias, plantas e matrizes são bibliográficas. Campomanesia littoralis permanece sem atribuição química nominal suficientemente confirmada; Eugenia catharinae está representada por comparadores de extratos, com a composição completa do óleo pendente.

Docking não demonstra mecanismo, inibição, agonismo, seletividade, eficácia, ADMET, toxicologia ou segurança. KEAP1 não equivale a um ensaio antioxidante, MAO não demonstra antidepressão, PARP1/EGFR não demonstram citotoxicidade seletiva e alvos bacterianos não demonstram CIM ou efeito bactericida. Escores de proteínas diferentes não são uma escala comum de afinidade.

## Métodos definidos antes de cada etapa

GNINA 1.3.3 oficial, fixado por SHA256, em CPU. CNN `fast` com `rescore`; ordenação principal `Energy`; exhaustiveness 16; até 20 poses por entrada; quatro CPUs por processo e no máximo dois processos GNINA simultâneos.

Redocking: sementes 101, 202, 303, 404 e 505. Aprovação exige RMSD pesado com simetria de até 2 Å na primeira pose por Energy em pelo menos quatro das cinco sementes. RMSD no referencial fixo do receptor, sem alinhamento isolado do ligante. Testes com referência idêntica, translação de 10 Å e átomos reordenados conferem a implementação. Melhor pose amostrada e reordenação por CNN no conjunto exportado são diagnósticos secundários, sem alterar a aprovação principal.

Triagem apenas nos aprovados: todas as 59 entradas, sementes 42, 142 e 242, mesmos parâmetros. São repetições da busca, não réplicas biológicas. Medianas e intervalos descrevem sensibilidade a busca, estereoquímica e conformação, sem inferência estatística de potência. O modelo CNN é igual ao da V3; não comparar seus números com entregas iniciais que usaram outro ensemble.

As caixas são definidas pelas coordenadas dos ligantes experimentais. Entradas independentes de redocking: ETKDGv3/MMFF, 64 tentativas, sem preservar coordenadas nativas e com identidade conferida após leitura SDF 3D. Hipóteses: harmina +1; SB284485 zwitteriônico (amina +1/carboxilato -1, carga líquida zero); demais controles em estado CCD neutro. Sem previsão populacional de pKa ou microestados.

## Preparação estrutural

Coordenadas pesadas preservadas; hidrogênios em modelo pH 7,4. Sem relaxação pesada ou reconstrução de alças. Cadeias laterais incompletas só completadas quando seus átomos existentes estavam a pelo menos 8 Å do ligante. Distâncias e resíduos constam nas auditorias.

FAD oxidado retido em MAO-A/B, fosfatos dianiónicos. Conexão Cys-SG/C8M preservada, sem H artificial em SG; coordenadas pesadas do FAD também preservadas. Sem modelo eletrônico/redox ou membrana. Todas as águas removidas; sensibilidade a essa decisão não testada.

GyrB conserva segmentos ausentes 82A-D e 101-115, com flancos a cerca de 5,9 Å do ligante. TyrS é construção catalítica com 323/420 resíduos modelados; o parceiro da montagem biológica fica a pelo menos 9,91 Å do ligante. Não se modelaram tRNA ou dinâmica do dímero. Numeração do PDB preparado, incluindo inserções, sem conversão automática para UniProt.

Antes das simulações, KEAP1/4L7B foi preterida por R354D e 6HWS por faltarem CZ/NH1/NH2 da Arg380 a 3,02 Å do ligante. Foi escolhida 4XMB, humana, sem mutações declaradas ou cadeias laterais incompletas. SAG/2V5Z tem conectividade de imina no CCD apesar do sinônimo safinamida; foi escolhido C18/2V61, mantendo sua química depositada. Preparações preliminares excluídas não contam como controles simulados.

## Auditoria e arquivos

Verificações: identidade, coordenadas finitas, contatos não metálicos abaixo de 2 Å, razão distância/raios vdW abaixo de 0,70, energia positiva, fora da caixa acima de 0,1 Å, alteração de ligação mapeada acima de 0,15 Å e sobreposição intramolecular grave. Metais separados. Ausência de alerta não prova pose correta; proximidade até 4 Å não demonstra ligação de hidrogênio ou mecanismo.

- `reports/`: PDF, síntese, figuras PNG/SVG e tabelas completas por composto, estrutura, entrada, semente e contatos.
- `targets/`: receptores, referências, inputs independentes, caixas e auditorias. Scripts PyMOL usam coordenadas reais, sem alinhamento isolado.
- `results/`: todas as poses novas, comandos, logs, hashes, controles reprovados e auditorias.
- `baseline_v3/`: sínteses anteriores identificadas, com procedência e integridade verificadas.
- `previous_delivery/`: pacote V3 original completo, inclusive sua quarentena documentada.
- `prepared/`, `structures.csv`, `docking_inputs.csv`: biblioteca curada utilizada.
- `raw/`: fontes RCSB/CCD, inclusive alternativas estruturais descartadas.
- `selection_plan.json` / `protocol.json`: etapa A; `selection_plan_round_b.json` / `protocol_round_b.json`: etapa B.
- `target_rationale.json`: hipóteses, limites e fontes primárias. Fontes vegetais em `literature_sources.json` / `literature_species_links.csv`.
- `SHA256SUMS.txt`: inventário da entrega.

## Reprodução no HOME

Extraia no HOME. Binário e ambiente virtual não estão incluídos; a instalação baixa vários GB, sem instalar drivers, alterar kernel ou exigir sudo/GPU.

```bash
cd "$HOME/Restinga_GNINA_Alvos_V4"
INSTALL_REPORT=1 bash 00_preparar_ambiente.sh
```

Para nova busca, crie uma pasta limpa, preservando a entrega. O script recusa destino existente e não copia poses calculadas da V4:

```bash
.venv/bin/python prepare_clean_reproduction_v4.py \
  --destination "$HOME/Restinga_GNINA_V4_Reproducao_Limpa"
```

Na pasta nova, instale o ambiente ou use o Python e o executável já preparados. Execute em sequência:

```bash
.venv/bin/python run_v4.py all --gnina "$HOME/Restinga_GNINA_Alvos_V4/tools/gnina-1.3.3"
.venv/bin/python run_v4_round_b.py all --gnina "$HOME/Restinga_GNINA_Alvos_V4/tools/gnina-1.3.3"
.venv/bin/python audit_v4.py
.venv/bin/python audit_redocking_v4.py
.venv/bin/python synthesize_v4.py
.venv/bin/python export_views_v4.py
.venv/bin/python build_report_v4.py
```

`--resume` reaproveita somente saídas concluídas com hashes e parâmetros compatíveis. Para recriar tabelas a partir das poses, use `audit_v4.py --reuse-audited`, depois `synthesize_v4.py` e `build_report_v4.py`. Os arquivos preparados mantêm a definição do sistema; nova preparação de hidrogênios pode variar numericamente. Nunca apague a entrega para recalcular.

Antes de priorizar candidatos, integrar os resultados químicos, confirmar vouchers/matrizes e confrontar hipóteses com ensaios específicos. As simulações foram executadas neste ambiente em CPU, não nas Tesla V100 ou estação local do usuário.
