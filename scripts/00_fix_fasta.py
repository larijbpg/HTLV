import re

# caminhos dos arquivos
# lê o arquivo bruto baixado do GenBank que contém as sequências do HTLV
with open("data/raw/htlv_sequences.fasta", "r", encoding="utf-8") as f:
    conteudo = f.read()

# troca qualquer símbolo >, inclusive o primeiríssimo do arquivo, garantindo a quebra de linha
conteudo_corrigido = conteudo.replace(">", "\n>")

# Remove linhas em branco extras que a quebra de linhas possa ter criado
# .split("\n"): divide o texto em uma lista a cada quebra de linha
# .strip() != "": filtra e descarta linhas vazias do resultado
linhas = [linha for linha in conteudo_corrigido.split("\n") if linha.strip() != ""]

linhas_limpas = []
for linha in linhas:
    # Se a linha for um cabeçalho FASTA (começa com >), mantém ela intacta
    if linha.startswith(">"):
        linhas_limpas.append(linha)
    else:
        # Se for linha de sequência, remove caracteres que não são nucleotídeos válidos (ex: números, Z, símbolos)
        # re.sub substitui tudo que NÃO for A, T, C, G, N ou degraus IUPAC por nada ("")
        seq_limpa = re.sub(r'[^ATCGNRYSWKMBDHVatcgnryswkmbdhv-]', '', linha)
        
        # garante que só adiciona a linha se sobrou algum nucleotídeo após a limpeza
        if seq_limpa:
            linhas_limpas.append(seq_limpa)

# junta todas as linhas da lista em um unico texto, colocando uma quebra de linha entre elas
conteudo_final = "\n".join(linhas_limpas)

# salva o resultado higienizado no arquivo corrigido pronto para o alinhamento
with open("data/raw/htlv_sequences_corrigido.fasta", "w", encoding="utf-8") as f:
    f.write(conteudo_final)

print("Arquivo corrigido e sanitizado com sucesso!")