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
    plt.title('K-mer Frequency')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.savefig("images/kmerFreq.png")
    plt.show()

if __name__ == "__main__":
    main()


