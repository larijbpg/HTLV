from pathlib import Path # portabilidade entre sistemas (Windows/LINUX)
from Bio import SeqIO
import pandas as pd

# Primeiro passo é definir os caminhos

arq_entrada = Path("results/aligments/htlv_env_sequence_aligned.fasta")
arq_saida = Path("results/tables/sequence_lengths.csv") # cria uma tabela
arq_limpo = Path("data/processed/htlv_env_sequence_clean.fasta")

# 1. Confirmar se o comprimento das sequencias estão corretos para a análise (se o filtro está correto e se elas tem o mesmo tamanho)

# Primeiro: vou conferir se tem sequencias fora do filtro e mostrar o id e o tamanho que elas tem (se existir)
tamanho_esperado = None
seq_erro = [] # crio uma lista para guardar as sequencia que estão com tamanho errado (se tiver alguma)

for sequencia in SeqIO.parse(arq_limpo, "fasta"): # para cada sequencia lida pelo SeqIO do arq_limpo:
    tamanho = len(sequencia.seq)
    if tamanho > 1500 or tamanho < 600:
        print(f"{sequencia.id}: {tamanho}pb")

# Segundo: vou ver se as sequencias tem o mesmo tamanho:

for sequencia in SeqIO.parse(arq_entrada, "fasta"):
    if tamanho_esperado == None:
        tamanho_esperado = len(sequencia.seq) # o tamanho da seq vai ser o tamanho da primeira sequencia lida
    else: 
        if len(sequencia.seq) != tamanho_esperado: # se não for None: se a sequencia for diferente do tamanho esperado:
            seq_erro.append({ # vou adicionar um dicionário nessa lista
                "id": sequencia.id, # vai pegar o id da sequencia
                "tamanho": len(sequencia.seq) # e vai pegar o tamanho dessa sequencia
            })
            print(f"Erro: A sequencia {sequencia.id} tem um tamanho:{len(sequencia.seq)}. O esperado era {tamanho_esperado}")  

# Terceiro: APENAS SE TIVER alguma sequencia dentro da lista seq_erro, vou criar um DataFrame com o pandas (para possíveis consultas futuras)
if seq_erro: 
    df = pd.DataFrame(seq_erro) 
    df.to_csv(arq_saida, index=False) # o df vai ser criado no mesmo caminho de "arq_saida"

print(f"Total de sequências com tamanho incorreto: {len(seq_erro)}")
print(f"O tamanho das sequências é: {tamanho_esperado}")

# Para ver se existe uma mutação eu vou comparar uma posição específica em todas as sequencias
# portanto vou criar uma lista de sequencias com seus ids (para consultas futuras)
seq_lista = []

for sequencia in SeqIO.parse(arq_entrada, "fasta"):
    seq_lista.append({
        "id": sequencia.id,
        "sequencia": sequencia.seq
    })
# Vou mostrar as primeiras 5 primeiras sequencias da lista de sequencias
print(seq_lista[0]["id"])

# Quarto: comparar cada posição do alinhamento entre todas as sequências, pra achar mutações

mutacoes = []  # vou guardar aqui as posições onde teve mais de uma letra (por exemplo: na coluna 0 todas as letras teriam que ser iguais)

# Vou percorrer cada posição do alinhamento de 0 até o tamanho esperado (4704), ou seja, percorrer colunas - pq sei que todas tem o mesmo tamanho
for i in range(0, tamanho_esperado): # para cada coluna numa escala de 0 a 4704
    letras_na_posicao = []  # letras de cada sequência só na posição i

    # para essa posição i(coluna), percorre todas as sequencias
    for sequencia in seq_lista:
        letra = sequencia["sequencia"][i] # a letra vai ser a posição i da sequencia (via chave do dicionario)
        letras_na_posicao.append(letra) # guarda essa letra na lista da posição atual

    letras_diferentes = set(letras_na_posicao)  # 'set' remove repetidos, só sobra o que é diferente
    letras_diferentes.discard("-")

    if len(letras_diferentes) > 1:  # se sobrou mais de uma letra diferente, é mutação
        mutacoes.append({
            "posicao": i,
            "letras": letras_diferentes
        })

print(f"Total de posições com mutação: {len(mutacoes)}")
porcentagem_mutacoes = (len(mutacoes) / tamanho_esperado * 100)
print(f"Porcentagem de mutações: {porcentagem_mutacoes:.2f}%")

# Avaliação da porcentagem de mutações
if porcentagem_mutacoes < 5:
    print(f"Porcentagem de mutações normal: {porcentagem_mutacoes:.2f}")
elif porcentagem_mutacoes >= 5 and porcentagem_mutacoes <= 30:
    print(f"Porcentagem de mutações considerável. Mas se for o caso de uma amostragem diversa, está razoável: {porcentagem_mutacoes:.2f}")
elif porcentagem_mutacoes > 30:
    print(f"Porcentagem de mutações alta: {porcentagem_mutacoes:.2f}. Revise o alinhamento ou o filtro de tamanho de sequências.")
