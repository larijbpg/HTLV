# Este arquivo vai criar um gráfico para visualização dos subtipos e subgrupos

from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd 
import seaborn as sns

# caminhos dos arquivos
metadata_file = Path("data/raw/metadata.csv")
output_plot = Path("results/figures/subtypes_by_continent.png")

# leitura ds dados
df = pd.read_csv(metadata_file, sep=";", on_bad_lines="skip")

# limpeza dos dados => remove linhas onde Continent ou Subtype estejam ausentes ou com o "-"
df_clean = df[(df["Continent"] != "-") & 
              (df["Subtype"] != "-") &
              (df["Subtype"].notna())] 

# tabela cruzada (subtipos por continente)
contagem = df_clean.groupby(["Continent","Subtype"]).size().unstack(fill_value=0)
    # .size(): contagem de linhas
    # .unstack(fill_value=0): pega o resultado do group(que fica em formato de lista/tabela) e 
    # pivota a coluna Subtype para virar colunas do gráfico, preenchendo com 0 os subtipos que nao existirem em determinado continente
proporcao = contagem.div(contagem.sum(axis=1), axis=0) * 100 
    # .div(..., axis=0): divide a contagem de cada subtipo pelo total do continente (divisão por linha)
    # contagem.sum(axis=1): calcula o total de sequências por continente (soma na horizontal)
    # * 100 transforma para porcentagem

# plotagem para gráfico Matplotlib + Seaborn
sns.set_theme(style="whitegrid")
fig, ax = plt.subplots(figsize=(10,6))
# plt.subplots(): cria a estrutura da figura (fig) e os eixos do gráfico (ax) com dimensão 10x6 polegadas

# gráfico de barras empilhadas
# kind="bar": define o tipo como grafico de barras
# stacked=True: empilha as categorias(subtipos) em uma unica barra por continnte
# ax=ax: direciona o desenho para o eixo criado pelo matplotlib
# cmap="tab20": aplica a paleta de cores para diferenciar as categorias
# edgecolor="black": desenha o contorno preto ao redor de cada bloco das barras
proporcao.plot(kind="bar", stacked=True, ax=ax, cmap="tab20", edgecolor="black")

# ajuste de título e eixos
ax.set_title("Proporção de Subtipos de HTLV-1 por Continente")
ax.set_ylabel("Proporção (%)")
ax.set_xlabel("Continentes")
plt.xticks(rotation=45, ha="right")
# plt.xticks(): rotaciona os nomes dos continentes em 45 graus com alinhamento à direita

# legenda fora do gráfico para nao cobrir as barras
ax.legend(title="Subtipo", bbox_to_anchor=(1.05, 1), loc="upper left")
# bbox_to_anchor=(1.05, 1): posiciona a legenda à direita do gráfico (fora da área de plotagem)
# loc="upper left": alinha o canto superior esquerdo da caixa de legenda no ponto definido por bbox_to_anchor

# remove bordas topo/direita com seaborn
sns.despine()

# salvar grafico
plt.tight_layout()
# plt.tight_layout(): ajusta automaticamente as margens da figura para que rótulos e legendas não fiquem cortados

plt.savefig(output_plot, dpi=300)
print(f"Gráfico salvo em {output_plot}")