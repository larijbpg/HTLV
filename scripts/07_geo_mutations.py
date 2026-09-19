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

# # TESTE COM UMA POSIÇÃO APENAS
# primeira_linha = df_mutacoes.iloc[0] # iloc[0] significa "pegue a linha pelo indice numerico 0"
# posicao_teste = primeira_linha["posicao"]

# lista_nova = []
# for sequencia in seq_lista:
#     # pegar o id longo (FASTA alinhado) → cortar/padronizar → comparar com o id curto que já existe pronto no metadata (sem precisar mudar nada lá)
#     id_longo = sequencia["id"]
#     partes = id_longo.split(".") 
#     id_curto = ".".join(partes[:2])

#     lista_nova.append({
#         "id" : id_curto,
#         "letra" : sequencia["sequencia"][posicao_teste]
#     })

# # agora eu quero juntar os dados da lista_nova (lista de dicionario) com os dados da tabela metadata (tabela). Para isso preciso que as duas sejam tabelas/df
# df_lista_nova = pd.DataFrame(lista_nova) #pandas cria um dataframe da lista_nova

# # agora posso dar o pd.merge() com as duas tabelas
# # resultado = pd.merge(tabela_esquerda, tabela_direita, left_on="coluna_da_esquerda", right_on="coluna_da_direita")
#     # tabela_esquerda e tabela_direita → as duas tabelas que você quer juntar
#     # left_on → o nome da coluna de id na tabela da esquerda
#     # right_on → o nome da coluna de id na tabela da direita (mesmo que o nome seja diferente!)
# df_final = pd.merge(df_lista_nova, df_metadata, left_on="id", right_on="Accession Number")
# print(df_final.head())

# # agrupar
# resultado = df_final.groupby(["Continent", "letra"]).size() # size(): vai contar quantas linhas tem em cada combinação
# print(resultado)

# Sabendo que o id e a região de uma sequencia específica não muda conforme a posição atual da sequencia, vou fazer um df com o id + regiao fixo 
lista_id = []
for sequencia in seq_lista:
    id_longo = sequencia["id"]
    partes = id_longo.split(".")
    id_curto = ".".join(partes[:2])

    lista_id.append(id_curto)

# transformar essa lista em uma tabela (DataFrame)
tabela_ids = pd.DataFrame(lista_id, columns=["id_curto"]) # columns=["id_curto"]: escrevo o nome que quero para a coluna de ids 

# agora vou juntar elas com o .merge()
df_id_regiao = pd.merge(tabela_ids, df_metadata, left_on="id_curto", right_on="Accession Number")
print(df_id_regiao.head())

# Agora vou percorrer cada linha de df_mutacoes(403 posições)
resultados_todas_posicoes = []  # vou guardar aqui o resultado do groupby de cada posição

for indice, linha in df_mutacoes.iterrows(): 
    # iterrows(): percorre linha por linha de um df e devolve a cada volta o numero do indice da linha e a linha em si (parecido com o enumerate())
    posicao = linha["posicao"]

    letras_da_posicao = []  # vou guardar id_curto + letra dessa posição, pra cada sequência

    # zip(seq_lista, lista_id): percorre as duas listas ao mesmo tempo,
    # pegando a sequência e o id_curto correspondente a ela, na mesma ordem
    for sequencia, id_curto in zip(seq_lista, lista_id):
        letra = sequencia["sequencia"][posicao]
        letras_da_posicao.append({
            "id_curto": id_curto,
            "letra": letra
        })

    # transforma essa lista em tabela, pra poder juntar com a região
    df_letras = pd.DataFrame(letras_da_posicao)

    # junta com df_id_regiao (que já tem id_curto + Continent prontos)
    df_juntado = pd.merge(df_letras, df_id_regiao, on="id_curto")
    # quando as duas tabelas têm a coluna com o mesmo nome (nesse caso, as duas têm id_curto), você não precisa dizer left_on e right_on separadamente,
    # só usa on="id_curto", e o Pandas já sabe que é a mesma coluna nos dois lados.

    # agrupa por Continent + letra, contando quantas sequências em cada combinação
    contagem_regiao = df_juntado.groupby(["Continent", "letra"]).size()

    resultados_todas_posicoes.append({
        "posicao": posicao,
        "contagem_por_regiao": contagem_regiao
    })

print(f"Total de posições processadas: {len(resultados_todas_posicoes)}")
