import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

# 1. Usar o alinhamento que efetivamente gerou a árvore do IQ-TREE
fasta_aprovado = "results/aligments/htlv_env_sequence_aligned.fasta"

ids_longos = []
ids_curtos = []

for sequencia in SeqIO.parse(fasta_aprovado, "fasta"):
    full_id = sequencia.id
    # Extrai APENAS a letra e os números principais (ex: de "AB036355.1.40..." pega só "AB036355")
    short_id = full_id.split(".")[0].strip()
    
    ids_longos.append(full_id)
    ids_curtos.append(short_id)

# Cria a ponte: ID Curto (do CSV) -> ID Longo (da Árvore)
mapa_ids = dict(zip(ids_curtos, ids_longos))

# Configuração dos arquivos e colunas
csv_input = "data/raw/metadata.csv"
txt_output = "results/figures/itol_env_geo_strip.txt"

col_id = "Accession Number"
col_geo = "Geographic Origin"

# Carregar o CSV
df = pd.read_csv(csv_input, sep=";", on_bad_lines="skip")

# Garante que a coluna do CSV fique limpa igual ao short_id (pega só "AB036355")
df[col_id] = df[col_id].astype(str).str.strip().str.split(".").str[0]

# 3. Filtra usando os IDs limpos e recupera o ID longo que a árvore espera
df = df[df[col_id].isin(ids_curtos)].copy()
df["tree_id"] = df[col_id].map(mapa_ids)

# Remova linhas sem informação
df = df.dropna(subset=["tree_id", col_geo])

# mapear regiões para cores dinamicamente
regiao_unica = sorted(df[col_geo].unique())
# unique() filtra a coluna e extrai apenas os valores unicos(sem repetições). Se tiver 100 amostras do "Brasil", essa função reduz todas para apenas uma entrada "Brasil"
# sorted(...): organiza essa lista de regiões em ordem alfabética/ ordem ascendente
# unique_regions = salva esse resultado em uma variável para que o script possa criar a legenda e atribuir uma cor diferente para cada região listada
cmap = plt.get_cmap(
    "turbo"
)  # carrega um mapa de cores pré-definido do Matplotlib chamado "turbo", gera um numero ilimitado de cores distintas sem estipular um numero fixo

cores_por_regiao = (
    {}
)  # cria um dicionario vazio para armazenar o par de cada regiao geográfica com a sua respectiva cor em hexadecimal (ex: {"América do Sul": "#1f77b4"})

for i, regiao in enumerate(
    regiao_unica
):  # loop que percorre a lista de regiões únicas, o enumerate fornece ao mesmo tempo o índice numérico (i: 0, 1, 2...) e o nome da região (regio
    rgb = cmap(i / len(regiao_unica))[:3]
    # Em vez do resto da divisão (%), ele divide o índice pelo número total de regiões. Isso garante que cada país receba um tom 100% exclusivo.
    cores_por_regiao[regiao] = mcolors.to_hex(rgb)
    # Converte os valores RGB para formato de código hexadecimal (como #FF0000) e salva no dicionário associando à região correspondente.

# montar o arquivo no formato que o iTOL exige:
with open(txt_output, "w", encoding="utf-8") as f:
    f.write("DATASET_COLORSTRIP\n")
    f.write("SEPARATOR TAB\n")
    f.write(
        "DATASET_LABEL\tRegiao_Geografica\n"
    )  # \t insere um espaço de tabulação (Tab)
    f.write("COLOR\t#000000\n")
    f.write("STRIP_WIDTH\t25\n")

    # Legenda
    f.write("LEGEND_TITLE\tRegião Geográfica\n")
    f.write(
        f"LEGEND_SHAPES\t"
        + "\t".join(["1"] * len(regiao_unica))
        + "\n"
    )
    f.write(
        f"LEGEND_COLORS\t"
        + "\t".join([cores_por_regiao[r] for r in regiao_unica])
        + "\n"
    )
    f.write(
        f"LEGEND_LABELS\t"
        + "\t".join([str(r) for r in regiao_unica])
        + "\n"
    )

    # Dados das sequências
    f.write("DATA\n")  # dados das amostras
    for _, row in df.iterrows():  # _(underline) não vou usar o valor desse indice, ele descarta o numero da linha, row guarda as info da amostra atual
        # Passa por cada linha da sua planilha, pega a ID da sequência, descobre a região dela, busca a cor correspondente dessa região no dicionário e escreve uma linha no texto.
        seq_id = str(row["tree_id"]).strip()
        geo = row[col_geo]
        color = cores_por_regiao[geo]
        f.write(f"{seq_id}\t{color}\t{geo}\n")

print(f"Arquivo de anotação do iTOL gerado com sucesso em: {txt_output}")