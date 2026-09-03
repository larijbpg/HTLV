import os
import matplotlib.pyplot as plt
from Bio import Phylo

# 1. Entrada de dados
input = "results/iqtree/iqtree.contree"
output = "results/figures/htlv_tree_python.png"

# garantir que a pasta de destino exista
os.makedirs("HTLV/results/figures", exist_ok=True) # exist_ok=True => evita que o Python lance um erro caso a pasta de destino já exista no diretório
print("Carregando a árvore filogenética...")
tree = Phylo.read(input, "newick") # "newick" => especifica o formato de arquivo da árvore filogenética que o Biopython deve ler

# Aplica enraizamento pelo ponto médio
tree.root_at_midpoint()

# configuração da figura
fig = plt.figure(figsize=(18,90), dpi=300) # figsize=(20, 55) => 20 de largura por 55 de altura, , dpi=300 => resolução da imagem em pontos por polegada (ideal)
axes = fig.add_subplot(1, 1, 1) 
# axes => área do gráfico onde a árvore será desenhada 
# add_subplot(1, 1, 1) => cria uma grade de 1 linha por 1 coluna e seleciona a primeira posição para plotar o gráfico

# quero exibir apenas bootstraps relevantes (>= 90)
# clade.confidence é o atributo do Biopython que armazena o valor de suporte estatístico (Bootstrap) de um nó(clade) interno da árvore filogenética
# clade.confidence = 100: Significa que aquele nó tem 100% de suporte no Bootstrap
for clade in tree.find_clades():
    if clade.confidence is not None and clade.confidence < 90:
        clade.confidence = None # oculta o numero se for menor que 90


# desenhar a arvore
Phylo.draw(
    tree,
    do_show=False, # evita que o Matplotlib abra uma janela interativa pop-up, permitindo salvar a imagem diretamente em arquivo
    axes=axes, # direciona a renderização da árvore para a área do gráfico definida anteriormente.
    show_confidence=True, # mostra os valores de bootstrap
    label_func=lambda leaf: str(leaf).split(".")[0] # simplifica os nomes das folhas no gráfico
)
# label_func_lambda: define uma função anônima inline para formatar o texto dos rótulos de cada nó da árvore.
# leaf: representa o objeto do nó folha (a sequência individual) que está sendo processado na função.
# str(leaf).split(".")[0]: converte o nome da folha em texto, corta na primeira ocorrência de ponto e pega apenas a parte anterior (ex: extrai AB273635 a partir de AB273635.1)

axes.set_title("Árvore Filogenética de Máxima Verossimilhança (HTLV - 257 Seq)", fontsize=14, fontweight = "bold")
axes.set_xlabel("Ditância Genética (substituições/sítio)", fontsize=6)

plt.tight_layout() # ajusta automaticamente o espaçamento dos elementos do gráfico para evitar sobreposição entre títulos e rótulos.
plt.savefig(output, dpi=300, bbox_inches="tight") #bbox_inches="tight" => recorta as margens brancas sobressalentes da imagem ao salvar, garantindo que nenhum rótulo seja cortado nas bordas.
print(f"Árvore gerada e salva com sucesso em {output}")


