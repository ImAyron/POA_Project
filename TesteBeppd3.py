from Bio import SeqIO
import re
import pandas as pd

# Função para identificar e registrar os segmentos maiúsculos em uma sequência
def find_uppercase_segments(sequence):
    segments = []
    matches = re.finditer(r'[A-Z]+', sequence)  # Encontra todas as substrings maiúsculas

    for match in matches:
        segment = match.group()
        start = match.start()
        end = match.end() - 1  # Final da sequência maiúscula

        # Adiciona ao dicionário o segmento e suas posições
        segment_info = {
            'segment': segment,
            'start_position': start + 1,
            'end_position': end + 1
        }
        segments.append(segment_info)

    return segments

# Função para processar as sequências do arquivo FASTA e criar o DataFrame
def process_fasta_to_dataframe(input_file):
    data = []  # Lista para armazenar as informações das sequências

    # Lê o arquivo FASTA e processa cada sequência
    for record in SeqIO.parse(input_file, "fasta"):
        seq_str = str(record.seq)
        header = record.id  # Identificação do cabeçalho da sequência
        segments = find_uppercase_segments(seq_str)
        
        # Adiciona as informações de cada segmento ao conjunto de dados
        for segment_info in segments:
            data.append([header, segment_info['segment'], segment_info['start_position'], segment_info['end_position']])

    # Cria o DataFrame a partir dos dados coletados
    df = pd.DataFrame(data, columns=["Header", "Segment", "Start_Position", "End_Position"])
    return df


# Exemplo de uso
input_file = "/mnt/c/Users/ubira/OneDrive/Desktop/Atualização POA Project/Teste_A1/E2_CHIKV/E2_CHIKV_Bepipred3Results/Bcell_epitope_preds.fasta"  # Coloque o nome do seu arquivo aqui
df = process_fasta_to_dataframe(input_file)
print(df)  # Exibe o DataFrame com as informações
