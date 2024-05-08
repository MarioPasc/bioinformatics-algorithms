#include "neighbour_joining.h"
#include <iostream>
#include <vector>
#include <string>

int main() {
    // Lista de secuencias a alinear
    std::vector<std::string> sequences = {
        "ATGCGA",
        "ATGGACGA",
        "ATCCA",
        "ATGGA",
        "CTGCAA"
    };

    // Creación del objeto NeighbourJoining
    NeighbourJoining nj(sequences);

    // Cálculo del alineamiento múltiple
    std::string result = nj.multiple_alignment();
    
    // Impresión del árbol guía y de la matriz de distancias para debugging
    std::cout << "Guide Tree Root: " << nj.get_root()->sequence << std::endl;
    std::cout << "Multiple Alignment Result:" << std::endl;
    std::cout << result << std::endl;

    // Impresión de la matriz de distancias
    std::cout << "Distance Matrix:" << std::endl;
    auto distance_matrix = nj.get_distance_matrix();
    for (const auto& row : distance_matrix) {
        for (int dist : row) {
            std::cout << dist << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
