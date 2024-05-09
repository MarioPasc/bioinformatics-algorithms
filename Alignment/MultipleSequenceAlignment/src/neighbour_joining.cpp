#include "neighbour_joining.h"
#include <iostream>
#include <limits>

NeighbourJoining::NeighbourJoining(const std::vector<std::string>& sequences) : sequences(sequences) {
    int num_sequences = sequences.size();
    distance_matrix = std::make_unique<std::vector<std::vector<int>>>(num_sequences, std::vector<int>(num_sequences, 0));
    nodes.resize(num_sequences);
    for (int i = 0; i < num_sequences; ++i) {
        nodes[i].id = i;
        nodes[i].left_child = -1;
        nodes[i].right_child = -1;
    }
}
void NeighbourJoining::calculate_distance_matrix() {
    int num_sequences = sequences.size();
    distance_matrix = std::make_unique<std::vector<std::vector<int>>>(num_sequences, std::vector<int>(num_sequences, 0));
    for (int i = 0; i < num_sequences; ++i) {
        for (int j = i + 1; j < num_sequences; ++j) {
            NeedlemanWunsch nw(sequences[i], sequences[j], 3, -1, -2);
            nw.align();
            int alignment_score = nw.get_alignment_score();
            (*distance_matrix)[i][j] = alignment_score;
            (*distance_matrix)[j][i] = alignment_score;
        }
    }
}

void NeighbourJoining::print_distance_matrix() const {
    int num_sequences = distance_matrix->size();
    for (int i = 0; i < num_sequences; ++i) {
        for (int j = 0; j < num_sequences; ++j) {
            std::cout << (*distance_matrix)[i][j] << " ";
        }
        std::cout << std::endl;
    }
    NeighbourJoining::find_smallest_distance_node();
}

void NeighbourJoining::find_smallest_distance_node() const {
    int num_sequences = distance_matrix->size();
    int min_distance = std::numeric_limits<int>::max();
    int min_i = -1;
    int min_j = -1;
    for (int i = 0; i < num_sequences; ++i) {
        for (int j = i + 1; j < num_sequences; ++j) {
            if ((*distance_matrix)[i][j] < min_distance) {
                min_distance = (*distance_matrix)[i][j];
                min_i = i;
                min_j = j;
            }
        }
    }
    std::cout << "Smallest distance: " << min_distance << " between nodes " << min_i << " and " << min_j << std::endl;
}

void NeighbourJoining::join_smallest_distance_nodes() {
    int num_sequences = sequences.size();
    int min_distance = std::numeric_limits<int>::max();
    int min_i = -1;
    int min_j = -1;

    // Encontrar el par de nodos con la distancia más pequeña
    for (int i = 0; i < num_sequences; ++i) {
        for (int j = i + 1; j < num_sequences; ++j) {
            if ((*distance_matrix)[i][j] < min_distance) {
                min_distance = (*distance_matrix)[i][j];
                min_i = i;
                min_j = j;
            }
        }
    }

    // Crear un nuevo nodo que represente la fusión de los nodos con la distancia más pequeña
    Node new_node;
    new_node.id = nodes.size();
    new_node.left_child = min_i;
    new_node.right_child = min_j;
    nodes.push_back(new_node);

    // Crear una nueva matriz de distancias con la dimensión correspondiente a la cantidad de nodos activos
    int new_num_sequences = num_sequences - 1;
    auto new_distance_matrix = std::make_unique<std::vector<std::vector<int>>>(new_num_sequences, std::vector<int>(new_num_sequences, 0));

    // Rellenar la nueva matriz de distancias con los nuevos valores
    int new_index = 0;
    for (int i = 0; i < num_sequences; ++i) {
        if (i != min_i && i != min_j) {
            int new_jndex = 0;
            for (int j = 0; j < num_sequences; ++j) {
                if (j != min_i && j != min_j) {
                    (*new_distance_matrix)[new_index][new_jndex] = (*distance_matrix)[i][j];
                    new_jndex++;
                }
            }
            int distance = ((*distance_matrix)[min_i][i] + (*distance_matrix)[min_j][i] - (*distance_matrix)[min_i][min_j]) / 2;
            (*new_distance_matrix)[new_index][new_num_sequences - 1] = distance;
            (*new_distance_matrix)[new_num_sequences - 1][new_index] = distance;
            new_index++;
        }
    }

    // Asignar la nueva matriz de distancias utilizando std::move
    distance_matrix = std::move(new_distance_matrix);
}