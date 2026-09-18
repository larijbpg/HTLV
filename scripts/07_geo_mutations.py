# 1. importar bibliotecas
# 2. caminho dos arquivos
# 3. cruzar dados de ids + letra naquela posição + região geográfica => id em comum nos arquivos
# 4. agrupar essa tabela por região + letra: revelando se alguma variação está concentrada geograficamente, ou espalhada por igual

from pathlib import Path
from Bio import SeqIO
import pandas as pd

# caminhos dos arquivos 
arq_mutacoes = Path("results/tables/htlv_seq_mutations.csv")
arq_id = Path("data/raw/metadata.csv")
arq_seq_al = Path("results/aligments/htlv_env_sequence_aligned.fasta")

# leitura dosarquivos
    # CSV com Pandas (pandas usa "," como separador mas o metadata está separado por ";")
    # FASTA com SeqIO.parse
df_mutacoes = pd.read_csv(arq_mutacoes) # pandas le o arquivo de mutações
df_metadata = pd.read_csv(arq_id, sep=";", on_bad_lines="skip")
seq_lista = []
for sequencia in SeqIO.parse(arq_seq_al, "fasta"):
    seq_lista.append({
        "id" : sequencia.id,
        "sequencia" : sequencia.seq.upper()
    })

# cruzar dados de ids + letra na posição i + região geográfica

# TESTE COM UMA POSIÇÃO APENAS
primeira_linha = df_mutacoes.iloc[0] # iloc[0] significa "pegue a linha pelo indice numerico 0"
posicao_teste = primeira_linha["posicao"]

lista_nova = []
for sequencia in seq_lista:
    # pegar o id longo (FASTA alinhado) → cortar/padronizar → comparar com o id curto que já existe pronto no metadata (sem precisar mudar nada lá)
    id_longo = sequencia["id"]
    partes = id_longo.split(".") 
    id_curto = ".".join(partes[:2])

    lista_nova.append({
        "id" : id_curto,
        "letra" : sequencia["sequencia"][posicao_teste]
    })
# agora eu quero juntar os dados da lista_nova (lista de dicionario) com os dados da tabela metadata (tabela). Para isso preciso que as duas sejam tabelas/df
df_lista_nova = pd.DataFrame(lista_nova) #pandas cria um dataframe da lista_nova

# agora posso dar o pd.merge() com as duas tabelas
# resultado = pd.merge(tabela_esquerda, tabela_direita, left_on="coluna_da_esquerda", right_on="coluna_da_direita")
    # tabela_esquerda e tabela_direita → as duas tabelas que você quer juntar
    # left_on → o nome da coluna de id na tabela da esquerda
    # right_on → o nome da coluna de id na tabela da direita (mesmo que o nome seja diferente!)
df_final = pd.merge(df_lista_nova, df_metadata, left_on="id", right_on="Accession Number")
print(df_final.head())

# agrupar
resultado = df_final.groupby(["Continent", "letra"]).size() # size(): vai contar quantas linhas tem em cada combinação
print(resultado)


