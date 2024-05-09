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

    // Mostrar la matriz de distancias original
    std::cout << "Matriz de distancias original:" << std::endl;
    nj.print_distance_matrix();
    
    // Fusionar los nodos con la distancia más pequeña
    nj.join_smallest_distance_nodes();

    // Mostrar la siguiente matriz de distancias
    std::cout << "Siguiente matriz de distancias:" << std::endl;
    nj.print_distance_matrix();
    
    /////////////////////////////////////////////////////////////////

    // Fusionar los nodos con la distancia más pequeña
    nj.join_smallest_distance_nodes();

    // Mostrar la siguiente matriz de distancias
    std::cout << "Siguiente matriz de distancias:" << std::endl;
    nj.print_distance_matrix();

    /////////////////////////////////////////////////////////////////

    // Fusionar los nodos con la distancia más pequeña
    nj.join_smallest_distance_nodes();

    // Mostrar la siguiente matriz de distancias
    std::cout << "Siguiente matriz de distancias:" << std::endl;
    nj.print_distance_matrix();
    return 0;
}