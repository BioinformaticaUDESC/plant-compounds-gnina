# Dados e fontes completos do estudo GNINA de sete espécies

Entrega de 30 de outubro de 2026. Os arquivos originais foram reunidos sem alterar coordenadas, buscas, escores ou conclusões dos controles. A biblioteca primária contém 40 compostos nomeados, 51 estruturas estereoquímicas e 59 entradas SDF. Foram recontadas 24.780 poses primárias de candidatos e 1.200 poses de controles nativos. Os arquivos diagnósticos e em quarentena permanecem identificados e não devem ser somados à análise primária.

## Arquivos solicitados e localização

| Conteúdo | Caminho a partir desta pasta |
|---|---|
| structures.csv | dados/Restinga_GNINA_Alvos_V4/structures.csv |
| primary_docking_inputs.csv | dados/Restinga_GNINA_Alvos_V4/primary_docking_inputs.csv |
| 59 entradas primárias SDF | Caminhos exatos na coluna input_path do manifesto; arquivos em dados/Restinga_GNINA_Alvos_V4/prepared/ligands |
| Configurações e protocolos | dados/Restinga_GNINA_Alvos_V4/protocol.json e protocol_round_b.json; dados/Restinga_GNINA_Multialvo_V3/protocol.json; targets/*/config.json |
| Saídas primárias V3 | dados/Restinga_GNINA_Multialvo_V3/results/COX2_5IKR, ACHE_4EY7 e CB2_5ZTY/docking; filtrar as entradas primárias pelo manifesto |
| Saídas primárias V4 | dados/Restinga_GNINA_Alvos_V4/results/SEH_3ANS, KEAP1_4XMB, PARP1_7KK4 e TYRS_1JIL/docking |
| Controles nativos | results/*/redocking em cada uma das duas pastas de dados |
| Sessões editáveis PyMOL | fontes_figuras/Figures_7_Plants_GNINA_V3/sessions; 12 controles nativos e 7 candidatos |
| Fontes Figuras 3, 4, S1 e S2 | render_pymol.py e compose_pymol_figures.py na pasta Figures_7_Plants_GNINA_V3 |
| Fontes Figuras 5 e S3–S7 | build_molecule_figures.py, figure_manifest.json e molecules/*.svg na mesma pasta; 51 SVG e 51 PNG |
| Figura 2 | fontes_figuras/Figura_2_reconstituida; ver a ressalva de procedência abaixo |
| Dados Suplementares S8 e S9 | informacoes_suplementares |
| Verificações anteriores | auditoria_verificada e arquivos originais verification_2026-10-27 das figuras |
| Respostas às perguntas | documentacao/Respostas_metodologicas.md |
| Racional por PDB | documentacao/Selecao_documentada_12_alvos.csv e os registros originais de seleção |

O lote SDF consolidado, saídas históricas adicionais e cópias de diagnóstico podem existir além das 59 entradas primárias. A unidade de inclusão é cada linha de primary_docking_inputs.csv e não a simples contagem de todos os arquivos *.sdf encontrados. O arquivo V3 comprimido também permanece em previous_delivery como original de procedência, mas foi expandido separadamente nesta entrega para acesso direto.

## Procedência da Figura 2

O script histórico de plotagem da Figura 2 não foi encontrado nos pacotes arquivados recuperados. rebuild_figure2.py é uma reconstrução nova, datada de 30/10/2026, e está identificada dessa forma no código e nos metadados. Ele lê as medianas de duas casas decimais efetivamente apresentadas na Tabela S1, calcula ranks por proteína com empates médios e confere os 224 valores contra a matriz anteriormente verificada. Os 224 valores coincidiram; não se afirma reprodução idêntica de pixels, fontes ou dimensões do gráfico antigo.

O PNG original incorporado ao manuscrito está preservado como Figure_2_Original_From_Manuscript.png, com identificação do arquivo de origem e hash. As novas saídas editáveis SVG e PDF e o PNG reconstruído estão em reconstructed. As demais fontes e as 19 sessões PyMOL são cópias inalteradas dos arquivos anteriormente entregues.

## Uso das fontes

Uma extração local pode ser mantida em /home/juan/Supplementary_Data_7_Plants_GNINA_Completo. Os scripts de figuras usam caminhos relativos à própria pasta. Para a reconstrução da Figura 2, executar o script rebuild_figure2.py com as dependências de requirements-figura2.txt. As versões originais do PyMOL e do RDKit constam de requirements-visualization.txt e visualization_environment.json na pasta de figuras. O ambiente da reconstrução não substitui o ambiente histórico de execução do GNINA.

Para conferir os ranks, o script aceita --input, --expected e --output-dir. Para inspecionar uma sessão, abrir o arquivo .pse correspondente no PyMOL; ela contém os objetos e a câmera. Para reproduzir as figuras moleculares, executar build_molecule_figures.py no ambiente registrado. render_pymol.py recompõe as vistas a partir das coordenadas; compose_pymol_figures.py monta as pranchas. Esses scripts podem sobrescrever suas saídas e manifestos; a reprodução deve ser feita numa cópia de trabalho, mantendo esta entrega original intacta.

Os comandos de docking históricos preservam caminhos absolutos originais como registros. Eles não são comandos portáveis prontos para execução na pasta atual. O executável GNINA não é redistribuído; sua identidade oficial e as dependências históricas estão documentadas no pacote. Esta entrega não executou novas buscas, dinâmica molecular ou ensaios biológicos.

## Verificações desta entrega

VERIFICACAO_COMPLETUDE.json registra as cópias originais preservadas, os 59 hashes de entradas conferidos, as contagens diretas de poses e sessões e a reconstrução dos ranks. COX2_verificacao_regra_da_caixa.json registra que a regra max(22 Å, extensão nativa + 12 Å) produz exatamente a caixa herdada de 22 × 22 × 22 Å, incluindo o centro. verificar_caixa_cox2.py reproduz somente essa conta geométrica com a biblioteca padrão de Python.

SHA256SUMS.txt cobre todos os arquivos deste pacote, exceto o próprio arquivo de checksums. Os checksums e registros históricos também foram preservados sem alterações dentro de suas pastas originais. Os documentos da IC em documentacao/comparacao_IC_Docking_OE são apenas exemplos para a comparação de ordenação e agregação de escores; não integram a triagem das sete plantas e não correspondem à auditoria integral daquele estudo.
