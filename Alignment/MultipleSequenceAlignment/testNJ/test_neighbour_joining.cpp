#include <iostream>
#include <vector>
#include <string>
#include "neighbour_joining.h"

int main() {
    // Secuencias de prueba
    std::vector<std::string> sequences = {
        "ATGCGA",
        "ATGGACGA",
        "ATCCA",
        "ATGGA",
        "CTGCAA"
    };

    // Crear objeto NeighbourJoining
    NeighbourJoining nj(sequences);

    // Calcular la matriz de distancias
    nj.calculate_distance_matrix();

    // Mostrar la matriz de distancias
    std::cout << "Matriz de distancias:" << std::endl;
    nj.print_distance_matrix();

    // Encontrar el par de nodos con la distancia más pequeña
    std::cout << "Par de nodos con la distancia más pequeña:" << std::endl;
    nj.find_smallest_distance_node();

    return 0;
}
