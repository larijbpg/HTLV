with open("HTLV/data/raw/htlv_sequences.fasta", "r", encoding="utf-8") as f:
    conteudo = f.read()

# troca qualquer símbolo >, inclusive o primeiríssimo do arquivo.
conteudo_corrigido = conteudo.replace(">", "\n>")
# Como ele insere um \n antes do primeiro >, a primeira linha do arquivo vira uma linha em branco. 
# Só que a sua linha seguinte (linhas = [linha for linha in ...]) remove imediatamente essa linha vazia gerada no começo.


# Remove linhas em branco extras que a quebra de linhas possa ter criado
linhas = [linha for linha in conteudo_corrigido.split("\n") if linha.strip() != ""]
# conteudo_corrigido.split("\n"): Pega todo o texto armazenado na variável conteudo_corrigido e o divide em uma lista de linhas, cortando o texto a cada quebra de linha (\n).
# Se a linha contiver apenas espaços (ou estiver zerada), ela vira "" e a condição avalia como falsa, descartando-a.
    # linha.strip() != "": É a condição de filtro. 
    # O método .strip() remove todos os espaços em branco, tabulações e quebras de linha invisíveis do início e do fim da linha. 
# Portanto, linhas vai receceber apenas as linhas do conteudo corrigido que nao estiverem vazias

conteudo_final = "\n".join(linhas) # junta todas as linhas da lista em um unico texto, colocando uma quebra de linha entre elas

with open("HTLV/data/raw/htlv_sequences_corrigido.fasta", "w", encoding="utf-8") as f:
    f.write(conteudo_final)

print("Arquivo corrigido salvo!")