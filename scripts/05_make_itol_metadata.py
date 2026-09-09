import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# configuração dos arquivos e colunas
csv_input = "data/raw/metadata.csv"
txt_output = "results/figures/itol_geo_strip.txt"

col_id = "Accession Number"
col_geo = "Geographic Origin"

# carregar o csv
df = pd.read_csv(csv_input, sep=";", on_bad_lines="skip")
# print(df.columns)
# remova linhas sem informação gráfica ou ID
df = df.dropna(subset=[col_id, col_geo])

# mapear regiões para cores dinamicamente
regiao_unica = sorted(df[col_geo].unique())
# unique() filtra a coluna e extrai apenas os valores unicos(sem repetições). Se tiver 100 amostras do "Brasil", essa função reduz todas para apenas uma entrada "Brasil"
# sorted(...): organiza essa lista de regiões em ordem alfabética/ ordem ascendente
# unique_regions = salva esse resultado em uma variável para que o script possa criar a legenda e atribuir uma cor diferente para cada região listada
cmap = plt.get_cmap("tab10") # carrega um mapa de cores pré-definido do Matplotlib chamado "tab10", que possui 10 cores bem distintas entre si.

cores_por_regiao = {} # cria um dicionario vazio para armazenar o par de cada regiao geográfica com a sua respectiva cor em hexadecimal (ex: {"América do Sul": "#1f77b4"})

for i, regiao in enumerate(regiao_unica): # loop que percorre a lista de regiões únicas, o enumerate fornece ao mesmo tempo o índice numérico (i: 0, 1, 2...) e o nome da região (region).
    rgb = cmap(i % 10) [:3] 
    # pega a cor correspondente ao indice i da paleta, extraindo apenas os 3 primeiros valores da tupla(red, green, blue). 
    # Os 10% faz o indice reinicar caso existam mais de 10 regioes
    cores_por_regiao[regiao] = mcolors.to_hex(rgb) 
    # Converte os valores RGB para formato de código hexadecimal (como #FF0000) e salva no dicionário associando à região correspondente.

# montar o arquivo no formato que o iTOL exige:
with open(txt_output, "w", encoding="utf-8") as f:
    f.write("DATASET_COLORSTRIP\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tRegiao_Geografica\n") # \t insere um espaço de tabulação (Tab)
    f.write("COLOR\t#000000\n")
    f.write("STRIP_WIDTH\t25\n")
    
    # Legenda
    f.write("LEGEND_TITLE\tRegião Geográfica\n")
    f.write(f"LEGEND_SHAPES\t" + "\t".join(["1"] * len(regiao_unica)) + "\n")
    f.write(f"LEGEND_COLORS\t" + "\t".join([cores_por_regiao[r] for r in regiao_unica]) + "\n")
    f.write(f"LEGEND_LABELS\t" + "\t".join([str(r) for r in regiao_unica]) + "\n")
    
    # Dados das sequências
    f.write("DATA\n") #dados das amostras 
    for _, row in df.iterrows(): # _(underline) não vou usar o valor desse indice, ele descarta o numero da linha, row guarda as info da amostra atual
    # Passa por cada linha da sua planilha, pega a ID da sequência, descobre a região dela, busca a cor correspondente dessa região no dicionário e escreve uma linha no texto.
        seq_id = str(row[col_id]).strip()
        geo = row[col_geo]
        color = cores_por_regiao[geo]
        f.write(f"{seq_id}\t{color}\t{geo}\n")

print(f"Arquivo de anotação do iTOL gerado com sucesso em: {txt_output}")
