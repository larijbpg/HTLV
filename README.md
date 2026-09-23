# 🧬 HTLV Molecular Diversity & Phylogeography Pipeline
Este repositório contém um pipeline em Python para obtenção, pré-processamento e alinhamento de sequências do vírus T-linfotrópico humano (HTLV), seguido de análise de mutação e identificação de variantes, caracterização geográfica e classificação de subtipos por reconstrução filogenética.

O objetivo é automatizar o download, pré-processamento, alinhamento e análise comparativa de amostras de HTLV armazenadas em bancos de dados públicos de bioinformática.

# 📌 Sumário

* Visão Geral
* Estrutura do Projeto
* Pré-requisitos
* Instalação e Configuração
* Passo a passo da Execução
    1. Obtenção de Dados
    2. Pré-processamento e Validação
    3. Alinhamento de Sequências
    4. Filogenética 
    5. Análise de Substituições e Variabilidade
    6. Agrupamento por Região Geográfica
    7. Agrupamento por Subtipo
* Estrutura dos Dados de Entrada
* Resultados 

# 🛠 Visão Geral

O pipeline processa arquivos genômicos no formato FASTA acompanhados de metadados (como país de origem, tipo/subtipo e ano de isolamento) para:

1. Calcular as estatísticas descritivas das sequencias (tamanho, conteúdo GC)
2. Realizar alinhamento múltiplo de sequências (MSA)
3. Filogenética (IQ-TREE, iTOL)
4. Identificar posições polimórficas (mutações/SNPs) 
5. Cruzar mutações genéticas identificadas com a distribuição geográfica das sequências 
6. Analisar a distribuição de subtipos de HTLV-1 por continente

> **Nota sobre o escopo do projeto:** O pipeline foi inicialmente desenvolvido e executado utilizando o **genoma completo** do HTLV. Em uma fase posterior, o escopo foi refinado para focar exclusivamente no **gene env**, com um filtro de tamanho específico (600-1500 pb) aplicado desde a etapa de obtenção de dados. Por esse motivo, algumas seções deste documento (referentes a etapas executadas antes da mudança) mencionam números de sequências aprovadas da versão com genoma completo (ex: 257), enquanto os resultados mais recentes refletem a versão com o gene `env` (ex: 2341 sequências aprovadas de 5301).
Portanto: 
    - Total baixado: 5.302
    - Após limpeza: 5.302
    - Aprovadas - versão genoma completo: 257
    - Aprovadas - versão gene env (600-1500pb): 2.339

# 📁 Estrutura do Projeto

htlv-genomics-pipeline/
- data
    - raw (dados brutos)
    - processed (dados limpos e filtrados)
- htlv_env (ambiente virtual)
- results
    - aligments (sequencias alinhadas)
    - figures (figuras geradas a partir do projeto)
    - iqtree (arquivos gerados pelo IQ-TREE)
    - tables (tabelas geradas dos arquivos)
- scripts
    - 00_fix_fasta.py
    - 01_fetch_data.py
    - 02_align_seq.py
    - 03_run_phylo.py
    - 04_plot_tree.py
    - 05_make_itol_metadata.py
    - 06_analyze_mutations.py
    - 07_geo_mutations.py
    - 08_plot_subtypes.py
- README.md
- REPORT.md
- .gitignore

# ⚙️ Pré-requisitos
* Python 3+
* Ferramenta de alinhamento multiplo externa: MAFFT
* Para a execução dos scripts de alinhamento e filogenia sem a necessidade de instalação global no sistema, as ferramentas executáveis (como o MAFFT e IQ-TREE) foram centralizadas localmente na pasta `FerramentasBioinfo/` na raiz do ambiente de desenvolvimento. 

> **Nota:** Por conterem arquivos binários executáveis, os arquivos desta pasta estão ignorados pelo versionamento (`.gitignore`).

# 🚀 Instalação e Configuração

pip install biopython pandas matplotlib seaborn 

# 🔄 Passo a Passo da Execução

1. Obtenção de Dados

Baixe as sequências de HTLV do HTLV Database, NCBI Virus ou GenBank no formato FASTA.
Salve o arquivo de sequências em `data/raw/htlv_sequences.fasta` e a tabela de metadados correspondente em `data/raw/metadata.csv`

### 1.1 Correção de Formatação FASTA (Debugging de Arquivo)

**Problema Encontrado:** O MAFFT finalizava o alinhamento sem erros de terminal, mas gerava o arquivo final vazio (0 bytes). A causa raiz foi a presença de cabeçalhos (`>`) colados na mesma linha no arquivo baixado, quebrando o padrão FASTA.

**Investigação Sistemática Realizada:**
* **Teste de instalação:** Execução do MAFFT com um arquivo FASTA minimalista (2 sequências) para isolar falhas do sistema/executável.
* **Monitoramento de Recursos:** Acompanhamento de tempo de execução e uso da CPU no Gerenciador de Tarefas para identificar a interrupção precoce da leitura.

**Solução Automática (`scripts/00_fix_fasta.py`):**
* `conteudo.replace(">", "\n>")`: Garante quebra de linha antes de qualquer cabeçalho.
* `[linha for linha in ... if linha.strip() != ""]`: Filtragem de linhas em branco residuais via *list comprehension*.
* `"\n".join(linhas)`: Reestruturação completa da sintaxe do arquivo.

`python scripts/00_fix_fasta.py`

2. Pré-processamento e Validação

Execute a filtragem inicial para remover sequencias muito curtas, com bases ambiguas ou sem metadados geográficos/subtipo associados.

- Tamanho: Remover sequências menores que 600 pb e maiores que 1500 pb (tamanho referencia do gene env baseado em artigos científicos do NCBI)
- Qualidade de bases: Remover sequências com excesso de bases ambíguas (por exemplo, mais de 1% ou 5% de letras fora de A, C, G, T).
- Metadados Geográficos: Verificar se no cabeçalho/descrição do arquivo FASTA ou no registro do NCBI existe a indicação do país/região de origem (essencial para estudos de filogeografia e epidemiologia molecular).

`python scripts/01_fetch_data.py`

3. Alinhamento Multiplo de Sequências (MSA via MAFFT)

Realize o alinhamento multiplo das sequencias para alinhar regiões homologas (como genes gag, pol, env ou a região LTR). 
Ou seja, Essa parte trata de organizar visualmente as sequências do HTLV para que a mesma região do vírus seja comparada exatamente no mesmo ponto em todas as amostras.
Como o vírus sofre mutações, as sequências brutas vêm com tamanhos ligeiramente diferentes. O alinhamento múltiplo insere lacunas (gaps, marcados por traços -) nas sequências para que os genes correspondentes (regiões homólogas) fiquem emparelhados coluna por coluna.

**Execução:**
`python scripts/02_align_seq.py`

**Principais Aprendizados de Código nesta Etapa:**
* **Normalização de IDs:** Uso de `sequencia.id.split(".")[0].strip()` para isolar o Accession Number sem versões decimais, permitindo a busca cruzada com a tabela de metadados.
* **`Seq` vs. `str`:** Conversão de objetos Biopython para texto puro (`str(sequencia.seq).upper()`) garantindo compatibilidade com métodos nativos de string e padronização de bases em maiúsculas.
* **Resiliência do Pipeline:** Uso de blocos `try / except Exception as e` para captura de exceções sem derrubar o fluxo completo de processamento.

`python scripts/02_align_seq.py`

4. Filogenética e Visualização

A reconstrução da árvore filogenética é realizada via **IQ-TREE 2** utilizando o método de Máxima Verossimilhança com suporte de nós por Ultrafast Bootstrap (UFBoot) e enraizamento pelo ponto médio (*midpoint rooting*).

#### 4.1 Abordagens para Gerar/Visualizar a Árvore
Existem 3 métodos avaliados no pipeline:
1. **iTOL na web (Método Oficial):** Utilizado para renderização final e figuras de publicação.
2. **Script em Python (`scripts/04_plot_tree.py`):** Mantido no repositório para validação rápida local da árvore.
3. **FigTree:** Requer Java Runtime Environment (JRE) para interface gráfica (opcional).

> **Decisão de Design:** Após testar a plotagem local e via iTOL, optou-se exclusivamente pelas imagens geradas pelo **iTOL** para compor os resultados finais e figuras de publicação.

#### 4.2 Geração de Anotações para o iTOL

Para colorir a árvore no iTOL por região geográfica, o script faz a busca cruzada mantendo a equivalência de IDs do FASTA com os metadados do CSV.

**Execução:**
`python scripts/05_make_itol_metadata.py`

**Principais Aprendizados de Código nesta Etapa:**
* **Mapeamento de IDs (De-para):** Uso de dicionário (mapa_ids = dict(zip(ids_curtos, ids_longos))) para relacionar o ID limpo do CSV (ex: AB273635.1) ao ID longo mantido pelo IQ-TREE na árvore (ex: AB273635.1.522.undefined.-.9033).
* **Formatos de Anotação iTOL:** Estruturação de arquivo .txt do tipo DATASET_COLORSTRIP com mapeamento automático de paletas de cores (matplotlib.colors) por metadado.
* **Colormap:** Para evitar colisão de cores entre regiões geográficas, foi usado o colormap contínuo "turbo"

### 4.3 Correção de Leitura de Metadados e Compatibilização com a Árvore (Debugging de Dados)

**Problema Encontrado:** O script `05_make_itol_metadata.py` falhava ao ler `metadata.csv` com `pandas.errors.ParserError`, indicando incompatibilidade no número de campos por linha. Após corrigido, o arquivo de anotação gerado não coloriu nenhum galho da árvore no iTOL (avisos de "Couldn't find ID... in the tree" para centenas de IDs).

**Investigação Sistemática Realizada:**
* **Inspeção manual de linhas problemáticas:** Contagem de separadores por linha (`linha.count(",")`) para isolar registros malformados.
* **Verificação do separador real do arquivo:** Identificado que o CSV usa `;` como delimitador, não `,` (padrão do Pandas).
* **Conferência de nomes de colunas:** Uso de `df.columns` para confirmar a grafia exata (`Accession Number`, não `Acession Number`).
* **Comparação de IDs entre arquivos:** Verificado que `metadata.csv` contém todas as sequências baixadas originalmente (5301), enquanto a árvore filogenética (na versão inicial focada no genoma completo) continha apenas as sequências aprovadas pelo filtro de qualidade (257), causando incompatibilidade de IDs.

**Solução Automática:**
* `pd.read_csv(csv_input, sep=";", on_bad_lines="skip")`: define o separador correto e ignora linhas malformadas remanescentes.
* `SeqIO.parse(fasta_aprovado, "fasta")`: extrai os IDs das sequências aprovadas diretamente do FASTA já filtrado, usados como referência.
* `df[df[col_id].isin(ids_aprovados)]`: filtra o metadata, mantendo apenas as linhas correspondentes às sequências presentes na árvore.

5. Análise de Mutações e Frequência Alélica

Antes de comparar as sequências posição por posição, o script valida se o alinhamento gerado pelo MAFFT está correto: confirma se todas as sequências têm exatamente o mesmo tamanho (incluindo os gaps) e se nenhuma sequência bruta ficou fora do filtro de tamanho (600-1500 pb) definido na etapa de pré-processamento.

Com o alinhamento validado, o script percorre cada posição (coluna) do alinhamento e verifica se há mais de uma base entre as sequências naquele ponto. Posições onde a única diferença é um gap (`-`) são descartadas, já que representam sequências mais curtas ou parciais, não uma variação biológica real.

Para separar variações genuínas de possível ruído (erro de sequenciamento ou amostras isoladas), aplica-se um filtro de **MAF (Minor Allele Frequency)**: uma posição só é considerada mutação se o alelo menos frequente aparecer em pelo menos 1% das sequências analisadas. As posições aprovadas são exportadas para `results/tables/htlv_seq_mutations.csv`, com a posição, as bases encontradas e o valor de MAF calculado.

Por fim, o script gera um gráfico de dispersão (Manhattan plot) das posições mutadas, salvo em `results/figures/`, permitindo visualizar regiões do gene `env` com maior concentração de variabilidade genética.

**Execução:**
`python scripts/06_analyze_mutations.py`

**Principais Aprendizados de Código nesta Etapa:**
* **Padronização de bases:** Uso de `.upper()` para uniformizar maiúsculas/minúsculas nas sequências antes da contagem — sem essa correção, bases idênticas (ex: `C` e `c`) eram contabilizadas como alelos diferentes, inflando artificialmente a taxa de mutação detectada (de 62,97% para 8,57% após a correção).
* **Filtro de gap na comparação:** Uso de `set.discard("-")` para excluir gaps da contagem de variantes, evitando que sequências parciais (mais curtas que o alinhamento) fossem interpretadas como mutação.
* **Cálculo de MAF:** Contagem de ocorrência de cada base por posição (dicionário de frequências) para calcular a proporção do alelo minoritário, filtrando variações estatisticamente pouco relevantes (ex: 1-2 sequências divergentes em meio a milhares).
* **Manhattan plot:** Visualização de dispersão (posição x MAF) via Matplotlib para identificar visualmente regiões do gene com maior concentração de variabilidade.

6. Agrupamento por Região Geográfica

Para investigar se as mutações identificadas na etapa anterior estão associadas a alguma região geográfica específica, o script `07_geo_mutations.py` cruza as posições mutadas com os metadados de origem das sequências.

Como o alinhamento usa o id longo mantido pelo IQ-TREE (ex: `AB036346.1.31.undefined.-.1080`) e o metadata usa o Accession Number no formato curto (ex: `AB036346.1`), o id longo é padronizado (cortado nos dois primeiros segmentos) antes do cruzamento, seguindo a mesma lógica já usada no `05_make_itol_metadata.py`.

Para cada uma das posições mutadas, o script recupera a base de cada sequência naquela posição, junta com o continente de origem (via merge pelo id padronizado) e agrupa por Continent + base, contando quantas sequências de cada região apresentam cada variante. Esse processo é repetido para as 20 posições de maior MAF, e o resultado é visualizado como um heatmap (Continente x Posição), com a intensidade representando a contagem da base minoritária.

**Execução:**
`python scripts/07_geo_mutations.py`

**Principais Aprendizados de Código nesta Etapa:**
* **Reaproveitamento de cálculo fixo:** A relação id → região é calculada uma única vez (`df_id_regiao`), fora do loop das posições, evitando refazer o merge 403 vezes — já que essa relação não muda de posição para posição.
* **`zip()` para combinar listas paralelas:** Uso de `zip(seq_lista, lista_id)` para percorrer sequência e id já calculado ao mesmo tempo, sem recalcular o id dentro do loop.
* **`pivot_table` para montar a matriz do heatmap:** Transformação de uma tabela "longa" (posição, continente, base, contagem) em uma matriz (continente x posição) com `pivot_table(index=..., columns=..., values=..., fill_value=0)`.
* **Filtro de base minoritária:** Uso de `.idxmin()` por grupo de posição para isolar, em cada posição, apenas a contagem da base menos frequente — evitando que a base majoritária (presente na maioria das sequências, em todas as posições) mascarasse a variação real entre regiões.

7. Análise de Subtipos e Subgrupos por Continente

Complementando a análise geográfica, o script `08_plot_subtypes.py` investiga a distribuição dos subtipos (`Subtype`) de HTLV-1 entre os continentes, permitindo comparar os achados do pipeline com a literatura científica sobre a diversidade genética do vírus.

O script calcula a proporção de cada subtipo por continente (normalizando pelo total de sequências de cada região, para não enviesar o resultado pela quantidade de amostras) e gera um gráfico de barras empilhadas mostrando essa distribuição.

**Execução:**
`python scripts/08_plot_subtypes.py`

**Principais Aprendizados de Código nesta Etapa:**
* **Tabela cruzada com `groupby` + `unstack`:** Uso de `df.groupby(["Continent", "Subtype"]).size().unstack(fill_value=0)` para transformar contagens agrupadas em uma tabela no formato continente x subtipo, pronta para plotagem.
* **Normalização por linha:** Uso de `contagem.div(contagem.sum(axis=1), axis=0) * 100` para calcular a proporção percentual de cada subtipo dentro do total de cada continente, evitando que regiões com mais amostras (ex: América do Sul) dominassem visualmente o gráfico.
* **Gráfico de barras empilhadas com Seaborn/Matplotlib:** Uso de `plot(kind="bar", stacked=True)` com legenda posicionada fora da área de plotagem (`bbox_to_anchor`), evitando sobreposição com as barras.
