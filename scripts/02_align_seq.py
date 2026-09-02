import subprocess
import os

# caminho das ferramentas centrais
mafft = r"C:\Users\Larissa\Documents\Codigos\FerramentasBioinfo\mafft.bat"

# caminhos dos arquivos de dados
entrada = "HTLV/data/processed/htlv_sequence_clean.fasta"
saida_alinhada = "HTLV/data/processed/htlv_sequence_aligned.fasta"

# Executar o Alinhamento Multiplo de Sequências (MSA)
print("Iniciando o Alinhamento das sequencias do HTLV com MAFFT...")


with open(saida_alinhada, "w", encoding="utf-8") as arquivo_saida:
    resultado = subprocess.run(
        [mafft, "--text", entrada],
        stdout=arquivo_saida,
        stderr=subprocess.PIPE,
        text=True
    )

print("ERRO:", resultado.stderr[:1000] if resultado.stderr else "(nenhum)")
print("Tamanho da saída:", os.path.getsize(saida_alinhada))

print(f"Alinhamento concluído com sucesso!")
print(f"Arquivo alinhado salvo em: {saida_alinhada}")