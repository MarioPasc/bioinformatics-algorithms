#include <iostream>
#include <vector>
#include <string>
#include "neighbour_joining.h"
#include <unordered_map>

int test_node_fusion() {
    // Secuencias de prueba en formato unordered_map
    std::unordered_map<std::string, std::string> sequences = {
        {"SeqA", "ATGCGA"},
        {"SeqB", "ATGGACGA"},
        {"SeqC", "ATCCA"},
        {"SeqD", "ATGGA"},
        {"SeqE", "CTGCAA"}
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

int main() {
    // Secuencias de prueba en formato unordered_map
    std::unordered_map<std::string, std::string> sequences = {
        {"SeqA", "ATGCGA"},
        {"SeqB", "ATGGACGA"},
        {"SeqC", "ATCCA"},
        {"SeqD", "ATGGA"}
    };

    // Crear objeto NeighbourJoining
    NeighbourJoining nj(sequences);

    // Mostrar la matriz de distancias original
    std::cout << "Matriz de distancias original:" << std::endl;
    nj.print_distance_matrix();

    nj.build_tree();
    return 0;
}