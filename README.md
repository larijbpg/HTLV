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
    4. Análise de Substituições e Variabilidade
    5. Agrupamento por Região Geográfica e Subtipo
* Estrutura dos Dados de Entrada
* Resultados Esperados


# 🛠 Visão Geral

O pipeline processa arquivos genômicos no formato FASTA acompanhados de metadados (como país de origem, tipo/subtipo e ano de isolamento) para:

1. Calcular a estatística descritivas das sequencias (tamanho, conteúdo GC)
2. Realizar alinhamento multiplo de sequências (MSA)
3. Identificar posições polimórficas (mutações/SNPs) e calcular matrizes de identidade percentual.
4. Cruzar variações genéticas com a distribuição geográfica e subtipos de HTLV (HTLV-1, HTLV-2, etc.)

# 📁 Estrutura do Projeto

htlv-genomics-pipeline/
├── data/
│   ├── raw/              # Arquivos FASTA e metadados brutos baixados
│   └── processed/        # Sequências filtradas e limpas
├── results/
│   ├── alignments/       # Arquivos de alinhamento (.aln / .fasta)
│   ├── figures/          # Gráficos gerados (heatmaps, distribuições)
│   └── tables/           # Tabelas em CSV/TSV com os resultados
├── scripts/
│   ├── 01_fetch_data.py   # Download via Entrez/Biopython ou validação local
│   ├── 02_align_seqs.py   # Alinhamento local/global das sequências
│   └── 03_analyze.py      # Análise de similaridade, subtipos e geografia
├── requirements.txt      # Dependências do projeto
└── README.md             # Documentação do projeto

# ⚙️ Pré-requisitos
* Python 3.8+
* Ferramenta de alinhamento multiplo externa (opcional, mas recomendada): MAFFT, Clustal Omega ou uso dos alinhadores internos do Biopython

# 🚀 Instalação e Configuração

pip install biopython pandas matplotlib seaborn

# 🔄 Passo a Passo da Execução

1. Obtenção de Dados

Baixe as sequências de HTLV do HTLV Database, NCBI Virus ou GenBank no formato FASTA.
Salve o arquivo de sequências em `data/raw/htlv_sequences.fasta` e a tabela de metadados correspondente em `data/raw/metadata.csv`

**Alternativa via script:** Você pode usar a ferramenta Bio.Entrez para buscar diretamente do GenBank via terminal usando IDs de acesso.

2. Pré-processamento e Validação

Execute a filtragem inicial para remover sequencias muito curtas, com bases ambiguas ou sem metadados geográficos/subtipo associados.

- Tamanho: Remover sequências muito curtas (<8000bp), já que o genoma do HTLV-1 tem aprox. 9000bp.
- Qualidade de bases: Remover sequências com excesso de bases ambíguas (por exemplo, mais de 1% ou 5% de letras fora de A, C, G, T).
- Metadados Geográficos: Verificar se no cabeçalho/descrição do arquivo FASTA ou no registro do NCBI existe a indicação do país/região de origem (essencial para estudos de filogeografia e epidemiologia molecular).

`python scripts/01_fetch_data.py`

3. Alinhamento de Sequências

Realize o alinhamento multiplo das sequencias para alinhar regiões homologas (como genes gag, pol, env ou a região LTR). 
Ou seja, Essa parte trata de organizar visualmente as sequências do HTLV para que a mesma região do vírus seja comparada exatamente no mesmo ponto em todas as amostras.
Como o vírus sofre mutações, as sequências brutas vêm com tamanhos ligeiramente diferentes. O alinhamento múltiplo insere lacunas (gaps, marcados por traços -) nas sequências para que os genes correspondentes (regiões homólogas) fiquem emparelhados coluna por coluna.

**gag:** Codifica as proteínas estruturais do capsídeo e do nucleocapsídeo do vírus
**pol:** Codifica as enzimas essenciais para a replicação viral: a Transcriptase Reversa, Integrase e a Protease.
**env:** Codifica as glicoproteínas da superfície do vírus que interagem com os receptores das células humanas para permitir a infecção.
**região LTR:** São sequências não codificantes localizadas nas extremidades do genoma proviral que funcionam como promotoras e reguladoras da transcrição do vírus.

`python scripts/02_align_seqs.py`

4. Análise de Substituições e Variabilidade
Calcule a matriz de identidade aos pares para identificar a divergencia nucleotidica entre as amostras ativas.

5. Agrupamento por Região Geográfica e Subtipo
Gere tabelas agregadas e gráficos de calor (heatmaps) que correlacionam o percetual de similaridade genética com o subtipo viral e o continente/pais de isolamento.
`python scripts/03_analyze.py`