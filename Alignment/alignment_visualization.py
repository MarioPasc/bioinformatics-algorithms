import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.colors as mcolors

def read_matrix():
    matrix = []
    line = sys.stdin.readline()
    while line.strip():
        row = list(map(int, line.split()))
        matrix.append(row)
        line = sys.stdin.readline()
    return np.array(matrix)

def read_trace_matrix():
    matrix = []
    line = sys.stdin.readline()
    while line.strip():
        row = list(line.strip().split())
        matrix.append(row)
        line = sys.stdin.readline()
    return np.array(matrix)


def main():
    # Leer "Score Matrix:"
    sys.stdin.readline()  # Leer y descartar la línea "Score Matrix:"
    score_matrix = read_matrix()
    
    # Leer y descartar líneas vacías y "Trace Matrix:"
    sys.stdin.readline()  # Leer y descartar línea en blanco
    trace_matrix = read_trace_matrix()
    
    # Leer y descartar líneas hasta "Alignment:"
    while "Alignment:" not in sys.stdin.readline():
        continue
    
    # Leer alineamientos
    sequence_a = sys.stdin.readline().split(': ')[1].strip()
    sequence_b = sys.stdin.readline().split(': ')[1].strip()
    
    # Visualización de la matriz de puntuación
    plt.figure(figsize=(10, 8))
    sns.heatmap(score_matrix, annot=True, fmt="d", cmap="viridis")
    plt.title("Score Matrix")
    plt.tight_layout()
    plt.savefig('/home/mariopasc/C++/bioinformatics-algorithms/images/score_matrix.png')
    plt.close()
    
    # Define the colors for each trace direction
    trace_directions = ['U', 'D', 'L']
    trace_colors = ['red', 'green', 'blue']  # Different color for each direction

    # Create a color map with a color for each direction
    cmap = mcolors.ListedColormap(trace_colors)

    # Create an array of the same shape as trace_matrix with numeric values
    # for each trace direction.
    trace_matrix_numeric = np.zeros_like(trace_matrix, dtype=int)
    for index, direction in enumerate(trace_directions):
        trace_matrix_numeric[trace_matrix == direction] = index

    # Plot the heatmap using the numeric matrix and the custom colormap
    plt.figure(figsize=(10, 8))
    sns.heatmap(trace_matrix_numeric, annot=trace_matrix, fmt='', cmap=cmap, cbar=False)
    plt.title("Trace Matrix")
    plt.tight_layout() 
    plt.savefig('/home/mariopasc/C++/bioinformatics-algorithms/images/trace_matrix.png')
    plt.close()
    
    print("Alignment:")
    print("Sequence A:", sequence_a)
    print("Sequence B:", sequence_b)

if __name__ == "__main__":
    main()
