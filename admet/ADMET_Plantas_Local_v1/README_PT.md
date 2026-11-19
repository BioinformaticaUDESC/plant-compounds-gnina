# ADMET complementar — plantas, execução local

Pacote para Ubuntu/Linux com Conda. A análise é exploratória e complementar ao docking do artigo. Instalação e inferência serão feitas na sua máquina. O pacote não implica que alguma substância seja segura, eficaz ou adequada para uso humano.

## Rodar

Extraia o ZIP na pasta Downloads. No terminal:

```bash
cd ~/Downloads/ADMET_Plantas_Local_v1
bash install_cpu.sh
conda activate admet_plantas_v1
bash run_local.sh
bash pack_results.sh
```

Se o navegador salvou a pasta em outro lugar, ajuste apenas o `cd`. O instalador precisa de internet e cria um ambiente separado chamado `admet_plantas_v1`. Se esse nome já existir, ele para sem substituir o ambiente. A inferência usa CPU, quatro threads e zero workers; não depende das V100 nem modifica drivers/CUDA. O pacote de modelos acompanha a distribuição oficial do ADMET-AI. Não são enviados SMILES a um servidor de predição.

Cada execução cria uma pasta nova em `results/` e um log em `logs/`. Não são sobrescritos resultados anteriores. Se houver erro, o processo termina e o `protocol.json` registra `failed`; arquivos parciais não constituem uma execução concluída.

`pack_results.sh` compacta a última inferência completa em `ADMET_Resultados_Local_*.zip`. Esse é o arquivo a enviar para incorporar a análise, as tabelas suplementares e os comentários ao artigo. O manuscrito ainda não foi alterado por este módulo.

Para validar somente identidade e calcular descritores, sem inferência ADMET:

```bash
python run_admet.py --mode descriptors --out results/somente_descritores
```

Esse modo NÃO produz resultados de absorção, distribuição, metabolismo, excreção ou toxicidade.

## Biblioteca preservada

Os dados vieram de `Supplementary_Data_7_Plants_GNINA_Final.zip`, pasta original `Restinga_GNINA_Alvos_V4`. O pacote inclui os CSV originais e os 59 SDF com seus hashes originais. Os arquivos históricos de associação com espécies são preservados como procedência e não foram novamente revisados neste módulo; qualquer associação ao texto deve seguir o manuscrito atual.

São 40 compostos, 51 estruturas declaradas e 59 entradas conformacionais de docking. A unidade primária aqui é cada estrutura química, não cada conformação ou semente. As 59 entradas não são 59 moléculas independentes. As 51 estruturas mantêm seus IDs, estereoquímica e microestados. O programa verifica SMILES, fórmula, carga, InChIKey e SHA256 de todos os SDF antes da análise. Identidades divergentes interrompem o processo.

Os 32 compostos voláteis e os oito comparadores não voláteis de extratos mantêm o campo `panel` separado em todas as tabelas. A análise não comprova presença em amostras, propriedades de um óleo ou efeitos da mistura. Não foram acrescentados compostos por inferência taxonômica.

O ácido rotúndico está com carga −1 na biblioteca de docking. O resultado primário mantém essa representação. Uma análise de sensibilidade adicional usa a forma neutralizada por RDKit Uncharger apenas quando a neutralização muda o SMILES. Essa forma não substitui a entrada primária e não representa enumeração de pKa, pH fisiológico ou distribuição de microestados. O esperado para esta biblioteca é uma estrutura adicional na sensibilidade.

## Métodos implementados

RDKit 2025.9.5: massa molecular média, cLogP, TPSA, doadores e aceptores de hidrogênio, ligações rotacionáveis, átomos pesados, fração Csp3, QED, violações de Lipinski (MW >500, cLogP >5, HBA >10, HBD >5) e alertas PAINS/BRENK/NIH. São descritores e filtros estruturais. Alertas não demonstram toxicidade, e Lipinski/QED não foram usados para excluir os produtos naturais nem para mudar a seleção do docking.

ADMET-AI **2.0.1**, Chemprop **2.2.2**, CPU: os 41 endpoints treinados são obtidos com `ADMETModel(include_physchem=False, drugbank_path=None, num_workers=0)`. O programa conserva as saídas contínuas, confere alinhamento de cada SMILES e arquiva hashes dos pesos e versões. As dependências centrais são fixadas; o ambiente resolvido completo é registrado em `environment_freeze.txt`. O ambiente de ADMET é independente do ambiente original de docking.

O modelo v2 foi retreinado e não reproduz numericamente o modelo v1 descrito no artigo de 2024 nem o servidor baseado em v1. Cite o artigo metodológico e identifique explicitamente a versão 2.0.1.

O painel cobre absorção intestinal, biodisponibilidade, solubilidade, permeabilidade, P-gp, BBB, ligação a proteínas, volume de distribuição, CYPs, meia-vida, clearance e endpoints toxicológicos como hERG, AMES e DILI. O arquivo oficial de endpoints contém definições, espécies, unidades, tamanhos e métricas de referência. Metadados e unidades são conservados na escala do fornecedor; não é feita conversão automática para doses ou concentrações. As métricas do fornecedor não são uma avaliação independente neste painel. Alguns endpoints têm desempenho limitado; precisam de interpretação individual.

Scores de classificação indicam a classe positiva de cada conjunto de treinamento; não são probabilidades clínicas calibradas. Valores altos têm significados diferentes entre endpoints. Não há regra global de 0,5 nem selo “seguro/tóxico”. Predições de regressão não são truncadas em limites físicos, para não esconder saídas problemáticas. O LD50 permanece na escala logarítmica do modelo; não se converte em dose humana. Ausência de alertas não demonstra segurança.

Não foi estabelecido domínio de aplicabilidade validado para essas estruturas, nem estimada incerteza individual. Os intervalos mínimo–máximo entre estruturas declaradas são amplitudes de representação, não intervalos de confiança. O programa não calcula um escore global ADMET, não combina endpoints com energias de docking e não promove ou exclui compostos por um corte escolhido após os resultados. Esta análise nova deve ser descrita como complementar pós hoc, mantendo o histórico da seleção original.

## Arquivos de saída

| Arquivo | Conteúdo |
|---|---|
| `protocol.json` | Estado da execução, versões, opções, contagens e limitações |
| `validated_structures.csv` | 51 estruturas com SMILES canônicos usados pelo modelo |
| `input_identity_audit.csv` | Auditoria de todas as 59 entradas SDF |
| `docking_input_mapping.csv` | Correspondência com os IDs do docking |
| `physchem_by_structure.csv` | Descritores e alertas estruturais |
| `admet_by_structure.csv` | 51 estruturas × 41 endpoints preditos |
| `admet_long_with_units.csv` | Predições com definições, unidades e espécie do endpoint |
| `admet_by_compound_ranges.csv` | Amplitudes por composto e painel, sem média entre estereoisômeros |
| `charge_sensitivity_inputs.csv` | Representações neutralizadas que diferem das primárias |
| `charge_sensitivity_predictions.csv` | Comparação forma original vs. neutralizada |
| `raw_model_predictions.csv` | Saída bruta para os SMILES únicos, sem percentis DrugBank |
| `vendor_endpoint_metadata*.csv` | Dicionário original do fornecedor |
| `model_weights_sha256.json` | Identidade dos pesos efetivamente usados |
| `environment_freeze.txt` | Ambiente resolvido da execução |
| `input_checksums.json`, `output_sha256.json` | Rastreabilidade de entradas e saídas |

Somente `protocol.json` com `completed_full_inference` permite declarar que a inferência ADMET foi concluída. O modo descritores registra estado distinto. Resultados de validação que acompanhem a entrega são identificados em `VALIDACAO.md` e não substituem a execução na sua máquina.

## Referências

- Swanson K. et al. ADMET-AI: a machine learning ADMET platform for evaluation of large-scale chemical libraries. Bioinformatics. 2024;40:btae416. https://doi.org/10.1093/bioinformatics/btae416
- Implementação usada: https://github.com/swansonk14/admet_ai/tree/v_2.0.1
- Alterações do v2: https://github.com/swansonk14/admet_ai/releases/tag/v_2.0.0
- RDKit: https://www.rdkit.org/docs/

Os metadados do ADMET-AI são reproduzidos com sua licença MIT em `reference/ADMET_AI_LICENSE.txt`. Os scripts deste pacote são novos e não são scripts históricos do docking.
