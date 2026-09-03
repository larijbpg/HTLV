import os
import sys

# caminhos dos arquivos
entrada = "HTLV/data/processed/htlv_aligned.fasta"
saida = "HTLV/results/iqtree"
iqtree = r"C:\User\Larissa\Documents\Codigos\FerramentasBioinfo\iqtree.exe"

# quero garantir que ao rodar essas pasta vao existir e se não existir o os vai cria-las

os.makedirs("HTLV/results/iqtree", exist_ok=True)

if os.path.exists(entrada):
    print("Arquivo localizado!")
else:
    print(f"Erro. O arquivo {entrada} não foi encontrado.")
    sys.exit(1) 
    # interrompe o script aqui se não encontrar o arquivo de entrada
    # (1) pq é o código de saída que avisa o computador que paramos por causa de um erro
    #(0) indicaria que deu certo/sucesso

# valida o executável do IQ-TREE
if os.path.exists(iqtree):
    print("Executável do IQ-TREE localizado!")
else:
    print(f"Erro: O executável do IQ-TREE não foi encontrado em {iqtree}")
    sys.exit(1)

# Executar o IQ-TREE
# Executar pelo terminal, vou abrir meu ambiente virtual (bioinfo)
