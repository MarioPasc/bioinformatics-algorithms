from Bio import SeqIO
from collections import Counter
from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
import csv
from scipy.spatial.distance import euclidean, cityblock
from scipy.stats import pearsonr
import os

class KmerAnalyzer:
    def __init__(self, k_values: List[int], file_paths: Dict[str, str], output_dir: str) -> None:
        self.k_values = k_values
        self.file_paths = file_paths
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

    def read_fasta(self, file: str) -> List[str]:
        sequences = []
        for record in SeqIO.parse(file, "fasta"):
            sequences.append(str(record.seq))
        return sequences

    def count_kmers(self, sequence: str, k: int) -> Counter:
        return Counter([sequence[i:i+k] for i in range(len(sequence) - k + 1)])

    def normalize_frequencies(self, kmer_counts: Counter, sequence_length: int, k: int) -> Dict[str, float]:
        norm_factor = sequence_length - k + 1
        return {k: v / norm_factor for k, v in kmer_counts.items() if v / norm_factor >= 10e-5}

    def plot_kmer_frequencies(self, kmer_frequencies: Dict[str, float], title: str, filename: str) -> None:
        k_mers = list(kmer_frequencies.keys())
        frequencies = list(kmer_frequencies.values())
        plt.figure(figsize=(10, 5))
        plt.bar(k_mers, frequencies, color='blue')
        plt.xlabel('k-mers')
        plt.ylabel('Frequency')
        plt.title(title)
        plt.xticks(rotation=90)
        plt.savefig(filename)
        plt.close()

    def calculate_distances(self, freq1: Dict[str, float], freq2: Dict[str, float]) -> Tuple[float, float, float]:
        kmers = set(freq1.keys()).union(set(freq2.keys()))
        vec1 = [freq1.get(k, 0) for k in kmers]
        vec2 = [freq2.get(k, 0) for k in kmers]
        euclid_dist = euclidean(vec1, vec2)
        manhattan_dist = cityblock(vec1, vec2)
        correlation, _ = pearsonr(vec1, vec2)
        return euclid_dist, manhattan_dist, correlation

    def save_kmer_frequencies(self, results: Dict[int, Dict[str, Dict[str, float]]]) -> None:
        for k, k_results in results.items():
            all_kmers = set()
            for freqs in k_results.values():
                all_kmers.update(freqs.keys())

            with open(os.path.join(self.output_dir, f"norm_kmer_freq_k{k}.csv"), 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                header = ["k", "kmer"] + list(k_results.keys())
                writer.writerow(header)
                for kmer in all_kmers:
                    row = [k, kmer] + [k_results[org].get(kmer, 0) for org in k_results]
                    writer.writerow(row)

    def save_distances(self, distances: Dict[int, Dict[str, Dict[str, float]]]) -> None:
        for k, k_distances in distances.items():
            with open(os.path.join(self.output_dir, f"distances_k{k}.csv"), 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(["Org_1", "Org_2", "Measure", "Value"])
                for pair, measures in k_distances.items():
                    org1, org2 = pair.split("_vs_")
                    for measure, value in measures.items():
                        writer.writerow([org1, org2, measure, value])

    def run_analysis(self) -> None:
        results = {k: {} for k in self.k_values}
        distances = {k: {} for k in self.k_values}
        for name, file in self.file_paths.items():
            sequences = self.read_fasta(file)
            sequence = sequences[0]
            sequence_length = len(sequence)
            for k in self.k_values:
                kmer_counts = self.count_kmers(sequence, k)
                normalized_kmers = self.normalize_frequencies(kmer_counts, sequence_length, k)
                results[k][name] = normalized_kmers
                plot_filename = os.path.join(self.output_dir, f"{name}_kmers_k{k}.png")
                self.plot_kmer_frequencies(normalized_kmers, f"k-mer Frequencies for {name} (k={k})", plot_filename)

        # Save k-mer frequencies to separate CSV files for each k
        self.save_kmer_frequencies(results)

        # Calculate and save distances between the organisms for each k
        for k in self.k_values:
            keys = list(results[k].keys())
            for i in range(len(keys)):
                for j in range(i+1, len(keys)):
                    org1, org2 = keys[i], keys[j]
                    euclid_dist, manhattan_dist, correlation = self.calculate_distances(results[k][org1], results[k][org2])
                    distances[k][f"{org1}_vs_{org2}"] = {
                        "Euclidean Distance": euclid_dist,
                        "Manhattan Distance": manhattan_dist,
                        "Pearson Correlation": correlation
                    }
            self.save_distances(distances)

def main() -> None:
    k_values = [2, 3, 4]  # List of k values
    base_path = "/home/mariopasc/C++/genomes_used"
    file_paths = {
        "E_coli": os.path.join(base_path, "Escherichia_Coli", "GCF_000005845.2_ASM584v2_genomic.fna"),
        "S_flexneri": os.path.join(base_path, "Shigella_flexneri", "GCF_000006925.2_ASM692v2_genomic.fna"),
        "S_enterica": os.path.join(base_path, "Salmonella_Enterica", "GCF_000006945.2_ASM694v2_genomic.fna"),
        "B_subtilis": os.path.join(base_path, "Bacillus_subtilis", "GCF_000009045.1_ASM904v1_genomic.fna")
    }
    output_dir = "./kmer_genetic_distance/data"

    analyzer = KmerAnalyzer(k_values, file_paths, output_dir)
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
