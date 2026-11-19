# GNINA multialvo V3 - Restinga

Estudo computacional exploratório realizado em 25/10/2026, em CPU. Leia `reports/Relatorio_GNINA_Multialvo_V3_25-10-2026.pdf` e `reports/RESULTADOS_V3.md` para a síntese. Não foi recebida a GC-MS do projeto; as relações entre substâncias, espécies e matrizes são bibliográficas.

## Alcance e controles

A biblioteca contém 40 compostos de literatura (32 voláteis e 8 comparadores de extratos), representados por 51 estruturas estereoquímicas. Quatro sistemas foram avaliados por redocking, com cinco sementes e critério definido antes dos resultados: primeira pose por Energy com RMSD pesado com simetria até 2 Å, no referencial fixo da proteína, em pelo menos 4/5 sementes. A referência de partida foi gerada independentemente da conformação cristalográfica, sem usar alinhamento do ligante no cálculo de RMSD.

- COX-2 humana, 5IKR, mefenamato/ID8: aprovado, 5/5; RMSD mediano 0,481 Å.
- Acetilcolinesterase humana, 4EY7 cadeia B, donepezila/E20 protonada +1: aprovado, 5/5; mediana 1,070 Å.
- Construção cristalográfica CB2 humano/T4-lisozima, 5ZTY, AM10257/9JU: aprovado, 5/5; mediana 0,537 Å. Estado ligado a antagonista; não demonstra agonismo de candidatos.
- Construção estável de 5-LOX humana, 6N2W cadeia B, NDGA/30Z e Fe2: reprovado, 0/5; mediana 3,549 Å. Poses próximas ao cristal foram amostradas, mas não ficaram em primeiro lugar por energia. Nenhum candidato foi triado nesse alvo sob este protocolo.

Nos três alvos aprovados, a biblioteca foi simulada nas sementes 42, 142 e 242, com até 20 poses por entrada. Repetições variam a busca, não a estrutura rígida da proteína. Escores de proteínas diferentes não são uma escala comum de afinidade. Os resultados não demonstram eficácia, atividade das misturas, seletividade, toxicidade ou segurança.

## Correção de identidade e biblioteca principal

A preparação inicial herdou 60 entradas do conjunto anterior ampliado. A conferência da identidade após ler a geometria 3D detectou uma divergência estereoquímica em `caryophyllene_oxide_conf03`. Seus arquivos e saídas foram preservados e colocados em quarentena.

Foi gerada uma geometria corrigida com identidade verificada e repetida nas três sementes nos alvos aprovados. Ela ficou praticamente igual a uma das duas conformações corretas já retidas, sem satisfazer o limiar de independência conformacional de 0,5 Å. Por isso, o teste corrigido ficou como diagnóstico separado, sem peso adicional na análise principal.

A biblioteca principal final contém **59 entradas conformacionais verificadas**, mantendo todos os 40 compostos e 51 estruturas. Todas as tabelas de síntese e contatos usam somente essa biblioteca. A auditoria completa identifica explicitamente as saídas em quarentena e o diagnóstico separado.

- `docking_inputs.csv` e `prepared/all_inputs.sdf`: lote histórico inicial de 60 entradas, incluindo a entrada que foi posteriormente colocada em quarentena.
- `primary_docking_inputs.csv` e `prepared/primary_inputs.sdf`: biblioteca principal curada de 59 entradas.
- `all_attempted_docking_inputs.csv`: 61 entradas tentadas, contando a entrada histórica problemática e o diagnóstico corrigido.
- `correction_docking_inputs.csv`: diagnóstico corrigido, fora do conjunto principal.
- `reports/input_identity_qc.csv` e `reports/input_curation.json`: identidades esperadas/observadas, motivos de exclusão e seleção da correção.

A existência desse problema na entrada herdada também limita a interpretação de sua linha no pacote anterior. As demais entradas foram conferidas contra as identidades declaradas. Não usar os scores da entrada em quarentena em médias, rankings ou conclusões sobre o óxido de cariofileno.

## Preparação e parâmetros

GNINA 1.3.3, executável oficial fixado por SHA256; CPU, quatro CPUs por processo, até dois processos simultâneos. `--cnn_scoring rescore --cnn fast --pose_sort_order Energy --exhaustiveness 16 --num_modes 20 --no_gpu`.

O modelo CNN fast foi selecionado para viabilizar a biblioteca ampliada. Seus valores não devem ser comparados numericamente aos do ensemble padrão da entrega anterior. A ordenação principal por energia foi definida antes dos resultados desta rodada. A comparação com a melhor CNNscore dentro das poses exportadas é secundária e restrita ao mesmo conjunto; não é outra busca independente.

Os receptores mantêm as coordenadas pesadas depositadas. Hidrogênios no modelo pH 7,4; sem reconstrução de alças ou relaxação pesada. Águas, glicanos e aditivos foram excluídos. A preparação da COX-2 é a mesma do piloto, com ambas as cadeias e porfirinas COH de cobalto. Na 5-LOX, ferro foi retido com aproximação empírica e a região que recobre o sítio permanece desordenada. CB2 mantém a construção mutante com fusão T4-lisozima, sem membrana ou lipídios.

As novas moléculas foram preparadas com ETKDGv3/MMFF; posições estereoquímicas indefinidas foram enumeradas. Estados do banco foram mantidos, exceto o carboxilato de ácido rotúndico. Não foi prevista a distribuição de microestados/pKa, uma limitação especialmente relevante para fenólicos. A estereoquímica do banco ou a enumeração não estabelecem o enantiômero nas amostras.

A auditoria verifica identidade química, coordenadas finitas, contatos não metálicos <2 Å, distância inferior a 70% da soma dos raios de van der Waals, energia positiva, saída da caixa >0,1 Å, alteração de ligação >0,15 Å e sobreposição intramolecular grave entre átomos separados por mais de três ligações. Comprimentos de ligação usam correspondência de átomos. Metais são registrados separadamente. Ausência de alerta não demonstra pose correta. Proximidades de resíduos até 4 Å não são automaticamente classificadas como ligações de hidrogênio ou mecanismos.

## Conteúdo

- `reports/`: PDF, síntese Markdown, sínteses por composto/estrutura, comparação de ordenação, contatos recorrentes e figuras PNG/SVG.
- `results/<alvo>/redocking/`: todos os controles, incluindo os reprovados.
- `results/<alvo>/docking/`: lotes históricos e diagnósticos corrigidos, com comandos, logs e hashes de conclusão.
- `results/<alvo>/primary_docking_poses.csv`: somente as poses da biblioteca principal.
- `results/<alvo>/pose_quality_audit.csv`: todas as poses, com elegibilidade e motivo de exclusão.
- `results/<alvo>/seed_stability_by_input.csv`: variação de energia e RMSD entre sementes para cada entrada principal.
- `structures.csv`: identidades químicas declaradas e metadados da preparação inicial; a contagem efetiva de conformações após a curadoria está em `primary_docking_inputs.csv`.
- `literature_species_links.csv`, `literature_sources.json`, `species_coverage.csv`: fontes, matrizes e cobertura botânica.
- `physicochemical_descriptors.csv`: descritores calculados; não são ensaios nem previsões validadas de toxicologia.
- `target_rationale.json`: justificativas, limitações de construções e contexto bibliográfico.
- `targets/<alvo>/view_redocking.pml`: inspeção no PyMOL sem alinhamento ou fit.
- `protocol.json`, `correction_protocol.json`, `method_checks.json`, `native_identity_checks.json`: parâmetros e verificações.
- `SHA256SUMS.txt`: inventário da entrega.

Eugenia catharinae está representada por comparadores de extratos; a composição completa do óleo continua pendente de acesso. Campomanesia littoralis permanece sem atribuição química nominal suficientemente confirmada. Não foi transferida automaticamente a composição genérica de C. xanthocarpa.

## Ambiente e reprodução no HOME

Extraia a pasta `Restinga_GNINA_Multialvo_V3` diretamente no seu HOME. O binário e o ambiente virtual não estão no pacote. A preparação baixa vários GB de binário/bibliotecas; não instala driver nem módulo de kernel e não exige sudo.

```bash
cd "$HOME/Restinga_GNINA_Multialvo_V3"
INSTALL_REPORT=1 bash 00_preparar_ambiente.sh
```

Se já possui uma cópia do mesmo binário, informe `GNINA_BIN` com caminho absoluto; o checksum ainda será conferido. O modo desta entrega é CPU. Para recriar somente as tabelas/figuras do material entregue:

```bash
bash run.sh audit --reuse-audited
bash run.sh report
```

Para uma **nova reprodução limpa**, gere uma pasta nova contendo somente as 59 entradas principais. A ferramenta recusa um destino já existente e não executa a simulação automaticamente:

```bash
.venv/bin/python prepare_clean_reproduction.py \
  --destination "$HOME/Restinga_GNINA_V3_Reproducao_Limpa"
```

Nessa pasta nova, prepare o ambiente ou use o Python do ambiente anterior e execute `run_v3.py validate`, depois `run_v3.py screen`, `audit_v3.py` e `synthesize_v3.py`. A validação reprovada continuará impedindo a triagem de candidatos naquele alvo. As entradas em quarentena e diagnósticas não são usadas na reprodução limpa. A nova pasta não contém resultados entregues, evitando sobrescritas.

Para **reproduzir historicamente os lotes iniciais e o diagnóstico da entrega**, a opção explícita de reprodução histórica é necessária, pois o pré-voo atual recusa a entrada em quarentena. Em uma cópia de trabalho, preserve a pasta de resultados com outro nome antes de recalcular:

```bash
bash run.sh validate --resume
bash run.sh screen --resume --historical-input-reproduction
.venv/bin/python run_input_correction.py --gnina "$GNINA_BIN" --resume
bash run.sh audit --reuse-audited
bash run.sh report
```

Defina `GNINA_BIN` antes desse bloco, ou substitua por um caminho absoluto. `--resume` reaproveita apenas arquivos concluídos, com conteúdo e parâmetros conferidos; caminhos locais podem mudar, mas hashes e parâmetros precisam coincidir. Saídas incompletas ou parâmetros diferentes são recusados. Os registros dos comandos originais são preservados.

Os scripts `prepare_targets.py` e `prepare_library.py` documentam a preparação original e recusam sobrescritas; `curate_input_identity.py` documenta a curadoria. Para calcular novamente a biblioteca principal, prefira a reprodução limpa. Não é necessário recriar o receptor para reproduzir a triagem com os arquivos preparados fornecidos.

Antes de priorizar moléculas para o TCC/dissertação, integrar GC-MS/índices de retenção/áreas por matriz e espécie, confirmar vouchers, resolver composição pendente e confrontar as hipóteses com ensaios específicos. O docking não permite concluir qual planta é mais potente ou mais segura.
