import sys
from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

def read_fasta(fasta_file):
    """Lee un archivo FASTA y devuelve las secuencias en un diccionario."""
    sequences = {}
    with open(fasta_file, 'r') as f:
        sequence_id = None
        sequence_data = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if sequence_id is not None:
                    sequences[sequence_id] = ''.join(sequence_data).upper()
                sequence_id = line[1:]
                sequence_data = []
            else:
                sequence_data.append(line)
        if sequence_id is not None:
            sequences[sequence_id] = ''.join(sequence_data).upper()
    return sequences

def kmer_count(sequence, k):
    """Cuenta la frecuencia de k-mers en una secuencia."""
    counts = defaultdict(int)
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        counts[kmer] += 1
    return counts

def main(fasta_assembly, fasta_reference, k):
    sequences_assembly = read_fasta(fasta_assembly)
    sequences_reference = read_fasta(fasta_reference)
    
    total_counts_assembly = defaultdict(int)
    total_counts_reference = defaultdict(int)
    
    for seq_id, sequence in sequences_assembly.items():
        counts = kmer_count(sequence, k)
        for kmer, count in counts.items():
            total_counts_assembly[kmer] += count
    
    for seq_id, sequence in sequences_reference.items():
        counts = kmer_count(sequence, k)
        for kmer, count in counts.items():
            total_counts_reference[kmer] += count
    
    kmers = sorted(set(total_counts_assembly.keys()).union(set(total_counts_reference.keys())))
    assembly_counts = [total_counts_assembly[kmer] for kmer in kmers]
    reference_counts = [total_counts_reference[kmer] for kmer in kmers]
    
    # Calcular correlación de Pearson
    correlation, _ = pearsonr(assembly_counts, reference_counts)
    
    # Crear un DataFrame de pandas para los k-mers y sus frecuencias
    df = pd.DataFrame({
        'K-mer': kmers,
        'Assembly': assembly_counts,
        'Reference': reference_counts
    })
    
    # Visualizar los resultados con matplotlib
    plt.figure(figsize=(10, 8))
    bar_width = 0.4
    index = np.arange(len(kmers))
    
    plt.barh(index, df["Reference"], bar_width, color='red', label='Reference')
    plt.barh(index + bar_width, df["Assembly"], bar_width, color='blue', label='Assembly')
    
    plt.xlabel('Frequency')
    plt.ylabel('K-mer')
    plt.xticks(fontsize=8)
    plt.yticks(index + bar_width / 2, df["K-mer"], fontsize=8)
    plt.title('K-mer Frequency')
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Añadir MSE y correlación de Pearson a la leyenda
    plt.legend(title=f'Pearson Correlation: {correlation:.4f}\nTotal k-mers: {len(kmers)}\nTotal unique k-mers: {len(df[df["Assembly"] != df["Reference"]])}')
    
    plt.savefig("kmerFreq_comparison.png")
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python kmer_count.py <assembly_fasta_file> <reference_fasta_file> <k>")
        sys.exit(1)
    
    fasta_assembly = sys.argv[1]
    fasta_reference = sys.argv[2]
    k = int(sys.argv[3])
    main(fasta_assembly, fasta_reference, k)
