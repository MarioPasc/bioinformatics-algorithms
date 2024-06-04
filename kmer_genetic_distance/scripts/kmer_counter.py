import numpy as np
import pandas as pd
from typing import List
from Bio import SeqIO
import os
from tqdm import tqdm
from collections import Counter

class BioKMerAnalysis:
    def __init__(self, input_files: List[str], output_dir: str, k:int) -> None:
        self.input_files = input_files
        self.output_dir = output_dir
        self.k = k

    def read_fna(self, file_path: str) -> List[str]:
        sequences = []
        for record in SeqIO.parse(file_path, "fasta"):
            sequences.append(str(record.seq))
        return sequences

    def calculate_raw_frequencies(self, sequences: List[str]) -> pd.Series:
        kmers = [''.join(seq[i:i+self.k]) for seq in sequences for i in range(len(seq)-self.k+1)]
        frequency = Counter(kmers)
        return pd.Series(frequency)

    def normalize_frequencies(self, raw_frequencies: pd.Series, lengths: List[int]) -> pd.DataFrame:
        """
        Normaliza las frecuencias de k-mers utilizando la norma N = L - k + 1.
        
        Args:
            raw_frequencies (pd.Series): Serie de pandas con las frecuencias de k-mers.
            lengths (List[int]): Lista de longitudes de las secuencias de ADN.
        
        Returns:
            pd.DataFrame: DataFrame de pandas con las frecuencias normalizadas de k-mers.
        """
        norm_frequencies = raw_frequencies.copy()
        total_possible_kmers = sum([l - self.k + 1 for l in lengths])
        norm_frequencies /= total_possible_kmers
        return pd.DataFrame({'kmer': norm_frequencies.index, 'frequency': norm_frequencies.values})

    def save_to_csv(self, data: pd.DataFrame, filename: str) -> None:
        data.to_csv(filename, index=False)

    def run(self) -> None:
        for file_path, name in tqdm(zip(self.input_files, ["Escherichia_Coli", "Klebsiella_Pneumoniae"]), desc="Processing files", total=len(self.input_files)):
            sequences = self.read_fna(file_path)
            lengths = [len(seq) for seq in sequences]
            raw_frequencies = self.calculate_raw_frequencies(sequences)
            norm_frequencies = self.normalize_frequencies(raw_frequencies, lengths)
            output_file = f"{self.output_dir}/{name}_k{self.k}.csv"
            self.save_to_csv(norm_frequencies, output_file)

def main() -> int:
    base_path_genomic_data = "/home/mariopasc/C++/reference_genomes"
    input_files = [os.path.join(base_path_genomic_data, "Escherichia_Coli", "GCF_000005845.2_ASM584v2_genomic.fna"),
                   os.path.join(base_path_genomic_data, "Klebsiella_Pneumoniae", "GCF_000240185.2_ASM24018v2_genomic.fna")]
    output_dir = os.path.join(os.getcwd(), "kmer_genetic_distance", "data")

    analysis = BioKMerAnalysis(input_files, output_dir, k=2)
    analysis.run()

    return 0

if __name__ == "__main__":
    main()
