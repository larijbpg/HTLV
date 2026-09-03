import os
import sys
import subprocess


# caminhos dos arquivos
entrada = "HTLV/data/processed/htlv_sequence_aligned.fasta"
saida = "HTLV/results/iqtree/iqtree"
iqtree = r"C:\Users\Larissa\Documents\Codigos\FerramentasBioinfo\iqtree3.exe"

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

# Executar o IQ-TREE => Executar pelo terminal, vou abrir meu ambiente virtual (bioinfo)
# 1. Montar a lista de argumentos para o terminal
cmd = [
    iqtree, # abre a lista
    "-s", entrada, # "s" => flag de entrada
    "-m", "MFP",   # "-m" => flag do modelo
    "-bb", "1000", # "-bb" => flag do bootstrap
    "-nt", "AUTO", # usa todos os 4 nucleos da CPU
    "-pre", saida
]

# 2. Executa o comando e aguarda a finalização
subprocess.run(cmd, check=True)
    # check=True é fundamental: se o IQ-TREE falhar por qualquer motivo durante a execução, o Python interrompe o script e avisa que ocorreu um erro.
