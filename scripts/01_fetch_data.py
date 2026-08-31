from Bio import SeqIO
import pandas as pd
from Bio.SeqUtils import gc_fraction
import os

# Caminhos para os arquivos
fasta_file = "HTLV/data/raw/htlv_sequences.fasta"
output_clean_fasta = "HTLV/data/processed/htlv_sequence_clean.fasta"
output_report_csv = "HTLV/data/processed/pre-processing_report.csv"

os.makedirs("HTLV/data/processed", exist_ok=True)

sequences_info = []

# Leitura do arquivo FASTA usando Biopython
for sequencia in SeqIO.parse(fasta_file, "fasta"): # para cada sequencia na leitura e interpretação do arquivo fasta "fasta_file"
    sequences_info.append({ # adicione na lista sequence_info
        "ID": sequencia.id,
        "Tamanho_bp": len(sequencia.seq),
        "Descrição": sequencia.description,
        "Conteudo_GC": round(gc_fraction(sequencia.seq) * 100, 2)
    })

# Converte em DataFrame
df = pd.DataFrame(sequences_info)
print(f"Total de Sequencias lidas: {len(df)}")
print("\nPrimeiras sequencias:")
print(df.head())

# Salvar resumo dos dados
df.to_csv("HTLV/data/processed/summary_sequences.csv", index=False) 
# O parâmetro index=False instrui o Pandas a não salvar a coluna de índices numéricos no arquivo CSV final.
print("\nResumo Salvo em data/processed/summary_sequences.csv!")

# Limpeza de dados: sequencias curtas, bases ambíguas (N) e seq sem metadados geográficos
    # critérios de filtro:
tamanho_minimo = 8000 # genoma completo do HTLV 
seq_limpa = []
rejeitados = [] # serve para armazenar e rastrear o motivo exato pelo qual cada sequência foi descartada (ex: "Tamanho curto", "Muitas ambiguidades" ou "Sem país")

# ler arquivo de metadados para pegar a geografia (se existir)
try:
    df_meta = pd.read_csv("HTLV/data/raw/metadata.csv")
    ids_com_pais = set(df_meta[df_meta["country"].notna()]["Acession"]) 
    # set() função que cria um conjunto parecida com uma lista, mas nao aceita dados duplicados e busca extremamente rápida
except: # se ocorrer qualquer erro dentro do try, o python entra no except e cria a variavel ids_com_pais
    ids_com_pais = set()

# filtrar as sequencias
for sequencia in SeqIO.parse(fasta_file, "fasta"):
    seq_str = str(sequencia.seq).upper()
    tamanho = len(seq_str)
# validações
tem_tamanho = tamanho >= tamanho_minimo
tem_pouco_n = (seq_str.coun("N") / tamanho) <= 0.01 # max de 1% de Ns

# checa pais pelo ID no arquivo de metadados ou na descrição do FASTA
tem_pais = (sequencia.id in ids_com_pais) or ("country" in sequencia.description.lower())

if tem_tamanho and tem_pouco_n and tem_pais:
    seq_limpa.append(sequencia)

# salvar o arquivo FASTA limpo
SeqIO.write(seq_limpa, output_clean_fasta, "fasta")

print(f"\nLimpeza concluída! Sequências aprovadas: {len(seq_limpa)} de {len(df)}")
print(f"Arquivo limpo salvo em: {output_clean_fasta}")