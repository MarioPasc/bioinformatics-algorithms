import sys
from Bio import Phylo
import matplotlib.pyplot as plt

def parse_newick_from_output(file_path):
    """
    Parse the Newick string from a specified file that contains it in a specific section.
    """
    with open(file_path, 'r') as file:
        in_newick_section = False
        newick_string = ""

        for line in file:
            if line.strip() == "Árbol filogenético en formato Newick:":
                in_newick_section = True
            elif in_newick_section:
                if line.strip().endswith(';'):
                    newick_string = line.strip()
                    break

        return newick_string

def plot_tree(newick_string):
    """
    Plot a phylogenetic tree from a Newick string and save it to a specified file path as a PNG.
    """
    from io import StringIO

    handle = StringIO(newick_string)
    tree = Phylo.read(handle, "newick")

    # Configuramos la figura antes de dibujar el árbol
    fig = plt.figure(figsize=(10, 5))  # Ajusta el tamaño según sea necesario
    axes = fig.add_subplot(1, 1, 1)

    # Dibujamos el árbol en los ejes configurados
    Phylo.draw(tree, axes=axes)

    # Guardamos la figura completa, asegurando que todo se haya dibujado
    plt.savefig("/home/mariopasc/C++/bioinformatics-algorithms/images/phylogenetic_tree.png")
    plt.close()  # Cerramos la figura para liberar memoria

if __name__ == "__main__":
    file_path = sys.argv[1]
    newick_string = parse_newick_from_output(file_path)
    if newick_string:
        print("Newick string found:", newick_string)
        plot_tree(newick_string)
    else:
        print("No Newick string found in the output.")
