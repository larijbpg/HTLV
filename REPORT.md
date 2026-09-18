# Relatório de Análise Filogenética e Variabilidade Genômica: Gene env (HTLV)

Projeto: `htlv-genomics-pipeline`
Autor(a): Larissa Jennifer Borba
Gene alvo: env (HTLV-1)

## Visão Geral
**Objetivo**
Este relatório tem como objetivo padronizar o conjunto de dados de sequências do gene *env*(Codifica as glicoproteínas da superfície do vírus que interagem com os receptores das células humanas para permitir a infecção.), realizar o alinhamento multiplo de sequencias (MSA) e garantir homologia posicional. Além disso, construir uma ávore filogenéticade de máxima verossimilhança para mapear clados geográficos e identificar sítios polimórficos e mutações de menor frequência (MAF)

## Metodologia e Pipeline
1. Aquisição de dados atraves do HTLV Database `data/raw/htlv_sequences.fasta` e a tabela de metadados correspondente em `data/raw/metadata.csv`

2. Correção do formato fasta `00_fix_fasta.py`: 
    Foi necessário criar esse script para fazer correções no arquivo bruto de sequencias, pois ao relizar o alinhamento via MAFFT, o arquivo final resultava em 0 bytes por erro de formatação. Com isso, foi realizado esse script para padronizar as sequencias de forma que obtivesse um arquivo com sequencias limpas.
        - padronização: apenas letras A,T,C,G / degraus do IUPAC em maiúsculas. Tambem, padroniza o cabeçalho das sequencias, de forma que cada linha comece com ">" , com é o formato FASTA.

3. Pre-processamento e validação das sequencias:
    3.1 Remoção de sequencias < 600 pb e > 1500 pb, visto que o gene env tem um tamanho entre600 e 1500 pares de base;
    3.2 Remove sequencias sem origem geográfica
    3.3 Remove sequencias repetidas
    3.4 Verifica qualidade das sequencias, de forma que é removido as sequencias que possuirem excesso de bases ambíguas (mais de 1%/5% de letras fora de A, T, C, G)

4. Alinhamento Multiplo de Sequencias `02_align_seq.py`: 
    MSA via MAFFT. Essa parte trata de organizar visualmente as sequencias do HTLV para que a mesma regiao do virus seja comparada exatamente no mesmo ponto em todas as amostras. O alinhamento multiplo insere gaps ("-") para que as sequencias fiquem emparelhadas.

5. Reconstrução e Visualização filogenética
    A reconstrução da árvore filogenética é realizada via **IQ-TREE 2** utilizando o método de Máxima Verossimilhança com suporte de nós por Ultrafast Bootstrap (UFBoot) e enraizamento pelo ponto médio (*midpoint rooting*). `03_run_phylo.py`
    Após testar a plotagem local(`04_plot_tree.py`) e via iTOL, optou-se exclusivamente pelas imagens geradas pelo **iTOL** para compor os resultados finais e figuras de publicação. Vale ressaltar que, para adicionar a origem geográfica das amostras, foi criado o script `05_make_itol_metadata.py`, ou seja, um arquivo baseado na padronização necessario do itol para que a arvore fosse "pintada". 

6. Análise de Variantes 
    Para essa análise foi criado o script `06_analyze_mutations.py` que primeiro verifica se o alinhamento foi correto, ou seja, se todas as sequencias possuem o mesmo tamanho (foi 4704 pb). Esse tamanho é justificado pela inserção de gaps ("-") após o alinhamento pelo MAFFT necessária para que regiões homologas fiquem emparelhadas entre sequências de tamanhos originalmente diferentes. 
    Em seguida, o script analisa todas as posições (colunas) do alinhamento entre as 2.339 sequencias aprovadas, comparando os alelos presentes em cada uma. Posições em que a unica diferença encontrada é um gap são descartadas da análise, já que representam sequências mais curtas ou parciais, e não uma variação biológica real. 
    Para as posições com mais d euma base distinta, o script conta a ocorrência de cada alelo gerando registros como:
         {"posição" : 305, "letras": {"A","C"}, "maf":0.0234}
    
A partir dessa contagem, é calculado o **MAF (Minor Allele Frequency)** — a frequência do alelo menos comum naquela posição. Esse valor é utilizado para diferenciar variações genéticas reais de possível ruído (erro de sequenciamento, amostras isoladas): posições em que o alelo minoritário aparece em menos de 1% das sequências são descartadas, permanecendo na análise final apenas as variações estatisticamente relevantes.

As posições aprovadas são exportadas para `results/tables/htlv_seq_mutations.csv`, e um gráfico de dispersão (Manhattan plot), posição x MAF, é gerado em `results/figures/dispersao_env.png`, permitindo visualizar as regiões do gene *env* com maior concentração de variabilidade genética.

## Resultados Obtidos

**Qualidade do Alinhamento e Sequências**
    Essa validação confirmou 0 sequências com tamanho incorreto entre as 2.339 aprovadas, todas com 4.704 pb (incluindo os gaps inseridos pelo MAFFT). Nenhuma sequência bruta ficou fora do filtro de tamanho definido para o gene "env", confirmando a consistência do pré-processamento realizado nas etapas anteriores.


**Reconstrução Filogenética**
    É possivel observar que a maioria das amostras tem origem geográfica não definida (cor roxo escuro).
![Árvore filogenética do gene env colorida por região geográfica](results/figures/htlv_tree_env.png)


**Análise de Mutações e Variabilidade**
    Das 4.704 posições do alinhamento, 2.962 (62,97%) apresentaram alguma variação entre as sequências antes da aplicação do filtro. Após a pad5ronização de maiúsculas e plicação de filtro do MAF (>=1%), esse numero foi reduzido para 403 posições (8,57%), consideradas mutações relevantes (uma faixa condizente com a diversidade genética esperada para uma amostra global do vírus).
    O Manhattan plot evidencia concentração de picos de MAF em determinada faixa de posições, ou regiões mais conservadas.
![Manhattan plot de MAF por posição no gene env](results/figures/dispersao_env.png)


