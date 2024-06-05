from Bio import AlignIO
from collections import Counter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from scipy.spatial.distance import euclidean, cityblock
from scipy.stats import pearsonr
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
import os
from typing import Dict


class KmerAnalysis:
    def __init__(self, aligned_file: str) -> None:
        self.alignment = AlignIO.read(aligned_file, "fasta")
        self.genome1 = str(self.alignment[0].seq)
        self.genome2 = str(self.alignment[1].seq)

    def kmer_frequencies(self, sequence: str, k: int) -> Counter:
        kmer_counts = Counter([sequence[i:i+k] for i in range(len(sequence) - k + 1)])
        return kmer_counts

    def normalize_frequencies(self, kmer_freqs: Counter, sequence_length: int, k: int) -> Dict[str, float]:
        normalization_factor = sequence_length - k + 1
        normalized_freqs = {kmer: freq / normalization_factor for kmer, freq in kmer_freqs.items()}
        return normalized_freqs

    def calculate_kmer_frequencies(self, k: int) -> pd.DataFrame:
        freqs1 = self.kmer_frequencies(self.genome1, k)
        freqs2 = self.kmer_frequencies(self.genome2, k)
        
        norm_freqs1 = self.normalize_frequencies(freqs1, len(self.genome1), k)
        norm_freqs2 = self.normalize_frequencies(freqs2, len(self.genome2), k)
        
        all_kmers = set(freqs1.keys()).union(set(freqs2.keys()))
        
        data = []
        for kmer in tqdm(all_kmers, desc="Calculating k-mer frequencies"):
            data.append({
                "Kmer": kmer,
                "Frequency": freqs1.get(kmer, 0),
                "Normalized Frequency": norm_freqs1.get(kmer, 0),
                "Organism": "Helicobacter pylori"
            })
            data.append({
                "Kmer": kmer,
                "Frequency": freqs2.get(kmer, 0),
                "Normalized Frequency": norm_freqs2.get(kmer, 0),
                "Organism": "Neisseria gonorrhoeae"
            })
        
        df = pd.DataFrame(data)
        return df

    def save_frequencies_to_csv(self, df: pd.DataFrame, filename: str) -> None:
        df.to_csv(filename, index=False)

    def plot_heatmap(self, df: pd.DataFrame) -> None:
        pivot_df = df.pivot(index='Kmer', columns='Organism', values='Normalized Frequency')
        plt.figure(figsize=(10, 8))
        sns.heatmap(pivot_df, cmap='viridis')
        plt.title('Heatmap of k-mer Frequencies')
        plt.show()

    def plot_dendrogram(self, df: pd.DataFrame) -> None:
        pivot_df = df.pivot(index='Kmer', columns='Organism', values='Normalized Frequency').fillna(0)
        linked = linkage(pivot_df.T, 'ward')
        
        plt.figure(figsize=(10, 8))
        dendrogram(linked, labels=pivot_df.columns, leaf_rotation=90)
        plt.title('Dendrogram of k-mer Clustering')
        plt.show()

    def plot_pca(self, df: pd.DataFrame) -> None:
        pivot_df = df.pivot(index='Kmer', columns='Organism', values='Normalized Frequency').fillna(0)
        
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(pivot_df.T)
        
        plt.figure(figsize=(10, 8))
        plt.scatter(pca_result[:, 0], pca_result[:, 1])
        for i, organism in enumerate(pivot_df.columns):
            plt.text(pca_result[i, 0], pca_result[i, 1], organism)
        plt.title('PCA of k-mer Frequencies')
        plt.xlabel('PCA Component 1')
        plt.ylabel('PCA Component 2')
        plt.show()

    def plot_mds(self, df: pd.DataFrame) -> None:
        pivot_df = df.pivot(index='Kmer', columns='Organism', values='Normalized Frequency').fillna(0)
        
        mds = MDS(n_components=2, dissimilarity='precomputed')
        dist_matrix = np.zeros((pivot_df.shape[1], pivot_df.shape[1]))
        for i in range(pivot_df.shape[1]):
            for j in range(pivot_df.shape[1]):
                dist_matrix[i, j] = euclidean(pivot_df.iloc[:, i], pivot_df.iloc[:, j])
        mds_result = mds.fit_transform(dist_matrix)
        
        plt.figure(figsize=(10, 8))
        plt.scatter(mds_result[:, 0], mds_result[:, 1])
        for i, organism in enumerate(pivot_df.columns):
            plt.text(mds_result[i, 0], mds_result[i, 1], organism)
        plt.title('MDS of k-mer Frequencies')
        plt.xlabel('MDS Dimension 1')
        plt.ylabel('MDS Dimension 2')
        plt.show()

    def calculate_distances(self, df: pd.DataFrame) -> pd.DataFrame:
        pivot_df = df.pivot(index='Kmer', columns='Organism', values='Normalized Frequency').fillna(0)
        manhattan_distance = cityblock(pivot_df.iloc[:, 0], pivot_df.iloc[:, 1])
        euclidean_distance = euclidean(pivot_df.iloc[:, 0], pivot_df.iloc[:, 1])
        correlation, _ = pearsonr(pivot_df.iloc[:, 0], pivot_df.iloc[:, 1])
        
        return pd.DataFrame({
            'Distance': ['Manhattan', 'Euclidean', 'Correlation'],
            'Value': [manhattan_distance, euclidean_distance, correlation]
        })

    def save_distances_to_csv(self, df: pd.DataFrame, filename: str) -> None:
        df.to_csv(filename, index=False)


def main() -> int:
    # Uso de la clase
    aligned_file = "/home/mariopasc/C++/genomes_used/Alignment/aligned_genomes.fasta"
    analysis = KmerAnalysis(aligned_file)

    # Calcular frecuencias de k-mers para diferentes valores de k
    k_values = [2]
    for k in k_values:
        df = analysis.calculate_kmer_frequencies(k)
        analysis.save_frequencies_to_csv(df, f"./kmer_genetic_distance/data/kmer_frequencies_k{k}.csv")
        
        # Generar visualizaciones
        analysis.plot_heatmap(df)
        analysis.plot_dendrogram(df)
        analysis.plot_pca(df)
        analysis.plot_mds(df)

    # Calcular y guardar distancias para k=2
    df_k2 = analysis.calculate_kmer_frequencies(2)
    distances_df = analysis.calculate_distances(df_k2)
    analysis.save_distances_to_csv(distances_df, "./kmer_genetic_distance/data/distances_k2.csv")

    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
