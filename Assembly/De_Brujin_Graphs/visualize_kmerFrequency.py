import sys
import pandas as pd
import matplotlib.pyplot as plt

def main():
    input_data = sys.stdin.read()
    data = input_data.strip().split("\n")
    kmer_data = [line.split(", ") for line in data]
    kmer_info = [(item[0].split(": ")[1], int(item[1].split(": ")[1])) for item in kmer_data]

    df = pd.DataFrame(kmer_info, columns=["K-mer", "Frequency"])
    df = df.sort_values(by="K-mer", ascending=True)

    plt.figure(figsize=(10, 8))
    plt.barh(df["K-mer"], df["Frequency"], color='skyblue')
    plt.xlabel('Frequency')
    plt.ylabel('K-mer')
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=2)
    plt.title('K-mer Frequency')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.savefig("images/kmerFreq.png")
    plt.show()

    """
    input_data = sys.stdin.read()
    data = input_data.strip().split("\n")
    kmer_data = [line.split(", ") for line in data if ", " in line and len(line.split(", ")) == 2]
    kmer_info = [(item[0].split(": ")[1], int(item[1].split(": ")[1])) for item in kmer_data]

    # Creamos un DataFrame con los datos.
    df = pd.DataFrame(kmer_info, columns=["K-mer", "Frequency"])

    # Obtenemos los 50 k-mers más frecuentes.
    top_kmers = df.nlargest(20, 'Frequency')
    top_kmers = top_kmers.sort_values(by="Frequency", ascending=True)
    
    # Obtenemos los 50 k-mers menos frecuentes.
    bottom_kmers = df.nsmallest(20, 'Frequency')
    bottom_kmers = bottom_kmers.sort_values(by="Frequency", ascending=False)

    # Configuramos el tamaño de la figura para acomodar dos subplots.
    plt.figure(figsize=(12, 20))

    # Primer subplot para los 50 k-mers más frecuentes.
    plt.subplot(2, 1, 1)
    plt.barh(top_kmers["K-mer"], top_kmers["Frequency"], color='skyblue')
    plt.xlabel('Frequency')
    plt.ylabel('K-mer')
    plt.title('Top 50 Most Frequent K-mers')
    plt.tight_layout()

    # Segundo subplot para los 50 k-mers menos frecuentes.
    plt.subplot(2, 1, 2)
    plt.barh(bottom_kmers["K-mer"], bottom_kmers["Frequency"], color='salmon')
    plt.xlabel('Frequency')
    plt.ylabel('K-mer')
    plt.title('Top 50 Least Frequent K-mers')
    plt.tight_layout()

    # Guardamos la figura en el directorio images.
    plt.savefig("images/kmerFreq.png")
    """

if __name__ == "__main__":
    main()


