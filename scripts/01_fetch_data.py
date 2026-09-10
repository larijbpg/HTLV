from Bio import SeqIO
import pandas as pd
from Bio.SeqUtils import gc_fraction
import os

# Caminhos para os arquivos
fasta_file = "data/raw/htlv_sequences_corrigido.fasta" # arquivo bruto corrigido
output_clean_fasta = "data/processed/htlv_env_sequence_clean.fasta" # arquivo limpo
output_report_csv = "data/processed/pre-processing_report_env.csv" # histórico e a justificativa do filtro de dados

os.makedirs("data/processed", exist_ok=True)

# Leitura do arquivo FASTA usando Biopython
sequences_info = []

for sequencia in SeqIO.parse(fasta_file, "fasta"): # para cada sequencia na leitura e interpretação do arquivo fasta "fasta_file"
    sequences_info.append({ # adicione na lista sequence_info
        "ID": sequencia.id,
        "Tamanho_bp": len(sequencia.seq),
        "Descrição": sequencia.description,
        "Conteudo_GC": round(gc_fraction(sequencia.seq) * 100, 2)
    })

# Converte em DataFrame - pandas
df = pd.DataFrame(sequences_info)
print(f"Total de Sequencias lidas: {len(df)}")
print("\nPrimeiras sequencias:")
print(df.head())

# Salvar resumo dos dados - pandas
df.to_csv("data/processed/summary_sequences.csv", index=False) 
# O parâmetro index=False instrui o Pandas a não salvar a coluna de índices numéricos no arquivo CSV final.
print("\nResumo Salvo em data/processed/summary_sequences.csv!")

# Limpeza das sequencias: sequencias curtas, bases ambíguas (N) e seq sem metadados geográficos
# critérios de filtro:
# Foco no gene ENV => possui entre 600 e 1500 pb
tamanho_minimo = 600  
tamanho_maximo = 1500
seq_limpa = []
rejeitados = [] # serve para armazenar e rastrear o motivo exato pelo qual cada sequência foi descartada (ex: "Tamanho curto", "Muitas ambiguidades" ou "Sem país")

# ler arquivo de metadados para pegar a geografia (se existir)
# primeiro preciso ver as informações que vieram desse arquivo
df_meta = pd.read_csv("data/raw/metadata.csv", sep=";")
print(df_meta.info())

try:
    df_meta = pd.read_csv("data/raw/metadata.csv", sep=";")

    # Extrai apenas a parte principal do ID (ex: "A36594" em vez de "A36594.1")
    ids_raw = df_meta[df_meta["Geographic Origin"].notna()]["Accession Number"].dropna().astype(str)
    ids_com_pais = set(ids_str.split(".")[0] for ids_str in ids_raw)
        # set() função que cria um conjunto parecida com uma lista, mas nao aceita dados duplicados e busca extremamente rápida
        # dropna(): remove os nulos (Nan ou None)
        # df_meta["Geographic Origin"].notna(): filtro booleano identificando quais linhas da colunas "Geographic Origin" não estão vazias
        # df_meta[...]: aplica esse filtro no DataFrame
        # ["Accession Number"]: seleciona apenas a coluna com o codigo de id da sequencia (Accession Number) dessas linhas filtradas
        # astype(str): Converte todos os IDs da coluna para o tipo texto (string).
   
# Se ocorrer qualquer erro durante a execução, o programa não trava abruptamente:
# ele captura a exceção, guarda os detalhes na variável 'e' e executa o bloco abaixo.
except Exception as e:
    print(f"Erro ao ler metadata.csv: {e}") 
    ids_com_pais = set()

# Inicializar contadores para o diagnóstico
total_sequencias = 0
qtd_maior_1500 = 0
qtd_menor_600 = 0
qtd_n_baixo = 0   # N <= 1%
qtd_n_alto = 0    # N > 1%
qtd_com_pais = 0
qtd_sem_pais = 0

# filtrar as sequencias
for sequencia in SeqIO.parse(fasta_file, "fasta"):
    total_sequencias +=1 # contador manual, é igual ====> total_seq = total_seq + 1
    seq_str = str(sequencia.seq).upper()
    # Quando a biblioteca Biopython lê um arquivo FASTA, ela não guarda a sequência como um texto comum. Ela cria um objeto especial chamado Seq.
    # Para usar métodos de string, é melhor transformar ela em string 
    tamanho = len(seq_str)
    id_limpo = sequencia.id.split(".")[0].strip()
    # sequencia.id acessa o id do cabeçalho FASTA lido pelo Biopython
    # aplit(".")[0]: corta o ID no ponto e pega apenas o primeiro pedaço(indice 0), descrtando o .1 que indica a versão da sequencia no GenBank

    # avaliação do tamanho para exclusão de sequências
    if tamanho > tamanho_maximo:
        qtd_maior_1500 +=1
    elif tamanho < tamanho_minimo: # não usei else para não colocar 800pb na lista de <600
        qtd_menor_600 +=1

    # avaliação do conteudo de "N"
    pct_n = (seq_str.count("N") / tamanho)
    if pct_n <= 0.01:
        qtd_n_baixo +=1
    else:
        qtd_n_alto +=1

    # avaliação da origem geografica
    tem_pais = (id_limpo in ids_com_pais) or ("geographic origin" in sequencia.description.lower())
    if tem_pais:
        qtd_com_pais +=1
    else:
        qtd_sem_pais +=1

    #Aprovação
    if (tamanho>=tamanho_minimo) and (tamanho<=tamanho_maximo) and (pct_n <= 0.01) and tem_pais:
        seq_limpa.append(sequencia)

# 4. Exibição dos Relatórios no Terminal
print("="*20,"RELATÓRIO DE DIAGNÓSTICO DO FASTA","="*20)
print(f"Total de sequências analisadas: {total_sequencias}")
print(f"Tamanho > 1500 bp: {qtd_maior_1500} | Tamanho < 600 bp: {qtd_menor_600}")
print(f"Ambiguidade N <= 1%: {qtd_n_baixo} | Ambiguidade N > 1%: {qtd_n_alto}")
print(f"Com informação de país: {qtd_com_pais} | Sem informação de país: {qtd_sem_pais}")
print("="*70)
   

# salvar o arquivo FASTA limpo
SeqIO.write(seq_limpa, output_clean_fasta, "fasta")

print(f"\nLimpeza concluída! Sequências aprovadas: {len(seq_limpa)} de {len(df_meta)}")
print(f"Arquivo limpo salvo em: {output_clean_fasta}")