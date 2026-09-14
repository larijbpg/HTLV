import subprocess # chama ferramentas externas do sistema (MAFFT)
import os # lida com caminhos de arquivos e verifica o tamanho deles no disco

# caminho das ferramentas centrais
# no LINUX/WSL usamos direto o executavel do sistema, sem o .bat do Windows
mafft = "mafft"

# caminhos dos arquivos de dados do gene env
entrada = "data/processed/htlv_env_sequence_clean.fasta"
saida_alinhada = "results/aligments/htlv_env_sequence_aligned.fasta"

# Executar o Alinhamento Multiplo de Sequências (MSA)
print("Iniciando o Alinhamento das sequencias do HTLV com MAFFT...")

with open(saida_alinhada, "w", encoding="utf-8") as arquivo_saida:
    resultado = subprocess.run( # Python rode um programa externo do sistema operacional e salva as informações dessa execução dentro da variável "resultado"
        [mafft, "--auto", "--anysymbol", entrada], # "--auto": A instrução para o MAFFT escolher sozinho o algoritmo de alinhamento ideal com base no tamanho das sequências
        stdout = arquivo_saida, # stdout(saída padrão): Em vez de jogar as sequências alinhadas na tela do terminal, grave tudo direto dentro do arquivo de saída que abrimos
        stderr=subprocess.PIPE, # stderr(erro padrão): captura msg de aviso, progresso ou erro do MAFFT. PIPE: salva temporariamente caso quisermos checar depois se deu tudo certo
        text=True # garante que o texto gerado venha como string e não bytes
    )

print("ERRO:", resultado.stderr[:1000] if resultado.stderr else "(nenhum)")
# "resultado.stderr[:1000]": exibe os primeiros 1000 caracteres de erro/log se houver
# "if resultado.stderr": checa se existe alguma mensagem gravada dentro de stderr 
# Se o MAFFT tiver escrito qualquer coisa, essa primeira parte é exibida.
# Caso o stderr esteja completamente vazio, o Python simplesmente imprime a palavra "(nenhum)"
print("Tamanho da saída:", os.path.getsize(saida_alinhada))
# mostra no terminal o tamanho do arquivo gerado em bytes para validar que o alinhamento não salvou vazio
print(f"Alinhamento concluído com sucesso!")
print(f"Arquivo alinhado salvo em: {saida_alinhada}")