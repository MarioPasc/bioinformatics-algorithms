import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List, Tuple
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.manifold import MDS
import os

class GeneticDistanceVisualizer:
    def __init__(self, kmer_data_path: str, distance_data_path: str, output_dir: str, k: int, measure: str):
        self.kmer_data = pd.read_csv(kmer_data_path)
        self.distance_data = pd.read_csv(distance_data_path)
        self.organisms = ["E_coli", "S_flexneri", "S_enterica", "B_subtilis"]
        self.output_dir = os.path.join(output_dir, f'k{k}')
        self.k = k
        self.measure = measure

        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)

    def plot_heatmap(self):
        distance_matrix = self._create_distance_matrix(self.measure)
        plt.figure(figsize=(8, 6))
        sns.heatmap(distance_matrix, annot=True, fmt=".2f", cmap="viridis", xticklabels=self.organisms, yticklabels=self.organisms)
        plt.title(f"Heatmap of Genetic Distances (k={self.k}, Measure={self.measure})")
        heatmap_path = os.path.join(self.output_dir, f'heatmap_{self.measure}.png')
        plt.savefig(heatmap_path)
        plt.close()

    def plot_dendrogram(self):
        distance_matrix = self._create_distance_matrix(self.measure)
        linked = linkage(distance_matrix, 'single')
        plt.figure(figsize=(8, 6))
        dendrogram(linked, labels=self.organisms, orientation='top', distance_sort='ascending')
        plt.title(f"Dendrogram of Genetic Distances (k={self.k}, Measure={self.measure})")
        dendrogram_path = os.path.join(self.output_dir, f'dendrogram_{self.measure}.png')
        plt.savefig(dendrogram_path)
        plt.close()

    def plot_mds(self):
        distance_matrix = self._create_distance_matrix(self.measure)
        mds = MDS(n_components=2, dissimilarity="precomputed")
        mds_results = mds.fit_transform(distance_matrix)
        plt.figure(figsize=(8, 6))
        plt.scatter(mds_results[:, 0], mds_results[:, 1], label=self.measure)

        for i, organism in enumerate(self.organisms):
            plt.text(mds_results[i, 0], mds_results[i, 1], organism)

        plt.title(f"MDS Plot of Genetic Distances (k={self.k}, Measure={self.measure})")
        plt.xlabel("MDS1")
        plt.ylabel("MDS2")
        mds_path = os.path.join(self.output_dir, f'mds_{self.measure}.png')
        plt.savefig(mds_path)
        plt.close()

    def plot_combined_mds(self, k: int, measures: List[str]):
        plt.figure(figsize=(10, 8))
        colors = ['r', 'b']
        
        for i, measure in enumerate(measures):
            distance_matrix = self._create_distance_matrix(measure)
            mds = MDS(n_components=2, dissimilarity="precomputed")
            mds_results = mds.fit_transform(distance_matrix)
            plt.scatter(mds_results[:, 0], mds_results[:, 1], c=colors[i], label=f'{measure}', alpha=0.5)

            for j, organism in enumerate(self.organisms):
                plt.text(mds_results[j, 0], mds_results[j, 1], organism, color=colors[i])
        
        plt.title(f"Combined MDS Plot of Genetic Distances (k={k})")
        plt.xlabel("MDS1")
        plt.ylabel("MDS2")
        plt.legend()
        combined_mds_path = os.path.join(self.output_dir, f'combined_mds_k{k}.png')
        plt.savefig(combined_mds_path)
        plt.close()

    def plot_kmer_frequencies(self):
        plt.figure(figsize=(10, 6))
        for organism in self.organisms:
            plt.plot(self.kmer_data['kmer'], self.kmer_data[organism], label=organism)
        plt.xlabel('K-mer')
        plt.ylabel('Frequency')
        plt.title(f'K-mer Frequencies for Each Species (k={self.k})')
        plt.xticks(rotation=90, fontsize=8)  # Rotate x-axis labels 90 degrees and set font size to 8
        if (self.k==4):
            plt.xticks(rotation=90, fontsize=2)
        plt.legend()
        kmer_freq_path = os.path.join(self.output_dir, 'kmer_frequencies.png')
        plt.savefig(kmer_freq_path)
        plt.close()

    def plot_clustered_heatmap(self):
        plt.figure(figsize=(10, 8))
        kmer_matrix = self.kmer_data.drop(columns=['k']).set_index('kmer')
        kmer_matrix = kmer_matrix.T
        sns.clustermap(kmer_matrix, cmap="viridis", col_cluster=True, row_cluster=True)
        plt.title(f'Clustered Heatmap of K-mer Frequencies (k={self.k})')
        clustered_heatmap_path = os.path.join(self.output_dir, 'clustered_heatmap.png')
        plt.savefig(clustered_heatmap_path)
        plt.close()

    def _create_distance_matrix(self, measure: str) -> pd.DataFrame:
        filtered_data = self.distance_data[self.distance_data['Measure'] == measure]
        matrix = pd.DataFrame(index=self.organisms, columns=self.organisms)
        for _, row in filtered_data.iterrows():
            matrix.at[row['Org_1'], row['Org_2']] = row['Value']
            matrix.at[row['Org_2'], row['Org_1']] = row['Value']
        return matrix.fillna(0).astype(float)

    def plot_distance_comparison(self, all_distances: List[Tuple[str, List[pd.DataFrame]]], ks: List[int]):
        fig, axes = plt.subplots(2, 1, figsize=(12, 12))
        
        # Subplot 0: Heatmap of average distances for each k and each measure
        measures = [distance[0] for distance in all_distances]
        avg_distances = []
        for measure, distances in all_distances:
            avg_distances.append([dist.mean().mean() for dist in distances])
        heatmap_data = pd.DataFrame(avg_distances, columns=[f'k={k}' for k in ks], index=measures)
        sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="viridis", ax=axes[0])
        axes[0].set_title('Average Genetic Distances for Different k Values and Measures')

        # Subplot 1: Line plot of distances between species for different k
        for measure, distances in all_distances:
            for i, org1 in enumerate(self.organisms):
                for j, org2 in enumerate(self.organisms):
                    if i < j:
                        dist_values = [dist.at[org1, org2] for dist in distances]
                        axes[1].plot(ks, dist_values, label=f'{org1} - {org2} ({measure})')
        axes[1].set_xlabel('k')
        axes[1].set_ylabel('Distance')
        axes[1].set_title('Comparison of Genetic Distances for Different k Values')
        axes[1].legend()

        comparison_path = os.path.join(self.output_dir, '..', 'distance_comparison.png')
        plt.tight_layout()
        plt.savefig(comparison_path)
        plt.close()

def main():
    base_path = './kmer_genetic_distance/data'
    output_dir = os.path.join(base_path, 'figures')
    ks = [2, 3, 4]
    measures = ["Euclidean Distance", "Manhattan Distance", "Pearson Correlation"]
    all_distances = []

    for measure in measures:
        distances_for_measure = []
        for k in ks:
            kmer_data_path = os.path.join(base_path, f'norm_kmer_freq_k{k}.csv')
            distance_data_path = os.path.join(base_path, f'distances_k{k}.csv')
            visualizer = GeneticDistanceVisualizer(kmer_data_path, distance_data_path, output_dir, k, measure)
            visualizer.plot_heatmap()
            visualizer.plot_dendrogram()
            visualizer.plot_mds()
            visualizer.plot_kmer_frequencies()
            visualizer.plot_clustered_heatmap()
            distances_for_measure.append(visualizer._create_distance_matrix(measure))
        all_distances.append((measure, distances_for_measure))

    # Create comparison plot
    comparison_visualizer = GeneticDistanceVisualizer(kmer_data_path, distance_data_path, output_dir, k=2, measure='Euclidean Distance')
    comparison_visualizer.plot_distance_comparison(all_distances, ks)

    # Plot combined MDS for Euclidean and Manhattan distances
    for k in ks:
        combined_visualizer = GeneticDistanceVisualizer(kmer_data_path, distance_data_path, output_dir, k, measure='Euclidean Distance')
        combined_visualizer.plot_combined_mds(k, ['Euclidean Distance', 'Manhattan Distance'])

if __name__ == "__main__":
    main()
