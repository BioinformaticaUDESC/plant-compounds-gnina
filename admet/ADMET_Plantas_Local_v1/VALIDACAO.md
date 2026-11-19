# Validação da entrega

Foi executada uma inferência completa no ambiente de preparação desta entrega, em Linux com Python 3.12 e CPU. Essa execução NÃO ocorreu na máquina de Juan. Seus arquivos estão em `validation_reference/`, com protocolo, ambiente, pesos identificados por SHA256, resultados e log próprios.

Foram verificadas 51 estruturas, 40 compostos e as identidades e hashes dos 59 SDF originais. O teste produziu 2.091 resultados primários (51 × 41 endpoints), 41 comparações adicionais para a neutralização do ácido rotúndico e os descritores RDKit. Todas as saídas de inferência eram finitas e os scores de classificação estavam em [0,1]. O resolvedor de dependências passou em `pip check`. A sintaxe dos scripts Python e Bash foi conferida.

Esse teste demonstra funcionamento e integridade das entradas e saídas; não mede acurácia experimental, calibração clínica ou domínio de aplicabilidade. As predições não são resultados laboratoriais.

O instalador Conda não foi executado nesta infraestrutura: a instalação de teste usou um ambiente Python isolado com as mesmas versões centrais e PyTorch CPU. O script Conda foi conferido quanto à sintaxe. A execução final do usuário poderá ter diferenças nas dependências transitivas; cada execução grava o ambiente resolvido completo.

Para repetir na máquina local, siga `README_PT.md`. Uma execução nova será criada em `results/`, separada dos dados de referência. O manuscrito não foi editado nesta entrega.
