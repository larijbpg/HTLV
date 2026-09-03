# 🧬 HTLV Molecular Diversity & Phylogeography Pipeline

Este repositório contém um pipeline em Python para análise de diversidade genética, classificação de subtipos e caracterização geográfica de sequências do vírus T-linfotrópico humano (HTLV)

O objetivo é automatizar o dowload, pré-processamento, alinhamento e análise comparativa de amostras de HTLV armazenadas em bancos de dados públicos de bioinformática.

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
    6. Agrupamento por Região Geográfica e Subtipo
* Estrutura dos Dados de Entrada
* Resultados Esperados


# 🛠 Visão Geral

O pipeline processa arquivos genômicos no formato FASTA acompanhados de metadados (como país de origem, tipo/subtipo e ano de isolamento) para:

1. Calcular a estatística descritivas das sequencias (tamanho, conteúdo GC)
2. Realizar alinhamento multiplo de sequências (MSA)
3. Filogenética
4. Identificar posições polimórficas (mutações/SNPs) e calcular matrizes de identidade percentual.
5. Cruzar variações genéticas com a distribuição geográfica e subtipos de HTLV (HTLV-1, HTLV-2, etc.)

# 📁 Estrutura do Projeto

htlv-genomics-pipeline/
├── data/
│   ├── raw/
│   │   ├── htlv_sequences.fasta                    # Sequências brutas baixadas
│   │   ├── htlv_sequences_corrigido.fasta          # Sequências com quebras de linha ajustadas
│   │   └── metadata.csv                            # Tabela com metadados geográficos
│   └── processed/
│       ├── htlv_sequence_clean.fasta               # Sequências filtradas pós-validação
│       ├── htlv_sequence_aligned.fasta             # Sequências alinhadas pelo MAFFT
│       └── summary_sequences.csv                   # Resumo do processamento
├── results/
│   ├── alignments/                                 # Arquivos de alinhamento finais
│   ├── figures/                                    # Gráficos e árvores filogenéticas
│   └── tables/                                     # Tabelas estatísticas
├── scripts/
│   ├── 00_fix_fasta.py                             # Script de sanitização e correção do FASTA
│   ├── 01_fetch_data.py                            # Filtragem, validação e vínculo de metadados
│   └── 02_align_seq.py                             # Alinhamento múltiplo via MAFFT
├── README.md

# ⚙️ Pré-requisitos
* Python 3.8+
* Ferramenta de alinhamento multiplo externa (opcional, mas recomendada): MAFFT, Clustal Omega ou uso dos alinhadores internos do Biopython
* Para a execução dos scripts de alinhamento e filogenia sem a necessidade de instalação global no sistema, as ferramentas executáveis (como o MAFFT e IQ-TREE) foram centralizadas localmente na pasta `FerramentasBioinfo/` na raiz do ambiente de desenvolvimento. 

> **Nota:** Por conterem arquivos binários executáveis, os arquivos desta pasta estão ignorados pelo versionamento (`.gitignore`).

# 🚀 Instalação e Configuração

pip install biopython pandas matplotlib seaborn

# 🔄 Passo a Passo da Execução

1. Obtenção de Dados

Baixe as sequências de HTLV do HTLV Database, NCBI Virus ou GenBank no formato FASTA.
Salve o arquivo de sequências em `data/raw/htlv_sequences.fasta` e a tabela de metadados correspondente em `data/raw/metadata.csv`

**Alternativa via script:** Você pode usar a ferramenta Bio.Entrez para buscar diretamente do GenBank via terminal usando IDs de acesso.

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

- Tamanho: Remover sequências muito curtas (<8000bp), já que o genoma do HTLV-1 tem aprox. 9000bp.
- Qualidade de bases: Remover sequências com excesso de bases ambíguas (por exemplo, mais de 1% ou 5% de letras fora de A, C, G, T).
- Metadados Geográficos: Verificar se no cabeçalho/descrição do arquivo FASTA ou no registro do NCBI existe a indicação do país/região de origem (essencial para estudos de filogeografia e epidemiologia molecular).

`python scripts/01_fetch_data.py`

3. Alinhamento Multiplo de Sequências (MSA via MAFFT)

Realize o alinhamento multiplo das sequencias para alinhar regiões homologas (como genes gag, pol, env ou a região LTR). 
Ou seja, Essa parte trata de organizar visualmente as sequências do HTLV para que a mesma região do vírus seja comparada exatamente no mesmo ponto em todas as amostras.
Como o vírus sofre mutações, as sequências brutas vêm com tamanhos ligeiramente diferentes. O alinhamento múltiplo insere lacunas (gaps, marcados por traços -) nas sequências para que os genes correspondentes (regiões homólogas) fiquem emparelhados coluna por coluna.

**gag:** Codifica as proteínas estruturais do capsídeo e do nucleocapsídeo do vírus
**pol:** Codifica as enzimas essenciais para a replicação viral: a Transcriptase Reversa, Integrase e a Protease.
**env:** Codifica as glicoproteínas da superfície do vírus que interagem com os receptores das células humanas para permitir a infecção.
**região LTR:** São sequências não codificantes localizadas nas extremidades do genoma proviral que funcionam como promotoras e reguladoras da transcrição do vírus.

**Execução:**
`python scripts/02_align_seq.py`

**Principais Aprendizados de Código nesta Etapa:**
* **Normalização de IDs:** Uso de `sequencia.id.split(".")[0].strip()` para isolar o Accession Number sem versões decimais, permitindo a busca cruzada com a tabela de metadados.
* **`Seq` vs. `str`:** Conversão de objetos Biopython para texto puro (`str(sequencia.seq).upper()`) garantindo compatibilidade com métodos nativos de string e padronização de bases em maiúsculas.
* **Resiliência do Pipeline:** Uso de blocos `try / except Exception as e` para captura de exceções sem derrubar o fluxo completo de processamento.

`python scripts/02_align_seqs.py`

4. Filogenética e Visualização

A reconstrução da árvore filogenética é realizada via **IQ-TREE 3** utilizando o método de Máxima Verossimilhança com suporte de nós por Ultrafast Bootstrap (UFBoot) e enraizamento pelo ponto médio (*midpoint rooting*).

#### 4.1 Abordagens para Gerar/Visualizar a Árvore
Existem 3 métodos avaliados no pipeline:
1. **iTOL na web (Método Oficial):** Utilizado para renderização final e figuras de publicação.
2. **Script em Python (`scripts/04_plot_tree.py`):** Mantido no repositório para validação rápida local da árvore.
3. **FigTree:** Requer Java Runtime Environment (JRE) para interface gráfica (opcional).

> **Decisão de Design:** Após testar a plotagem local e via iTOL, optou-se exclusivamente pelas imagens geradas pelo **iTOL** para compor os resultados finais e figuras de publicação.

5. Análise de Substituições e Variabilidade
Calcule a matriz de identidade aos pares para identificar a divergencia nucleotidica entre as amostras ativas.

6. Agrupamento por Região Geográfica e Subtipo
Gere tabelas agregadas e gráficos de calor (heatmaps) que correlacionam o percetual de similaridade genética com o subtipo viral e o continente/pais de isolamento.
`python scripts/03_analyze.py`