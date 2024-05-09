#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>  // Add this line
#include "../include/neighbour_joining.h"
#include <limits>
#include <algorithm>

NeighbourJoining::NeighbourJoining(const std::unordered_map<std::string, std::string>& sequence_map) {
    int num_sequences = sequence_map.size();
    distance_matrix = std::make_unique<std::vector<std::vector<int>>>(num_sequences, std::vector<int>(num_sequences, 0));
    nodes.reserve(num_sequences);

    for (const auto& entry : sequence_map) {
        const std::string& seq_id = entry.first;
        const std::string& sequence = entry.second;

        Node* leaf_node = new Node();
        leaf_node->id = seq_id;
        leaf_node->left_child = nullptr;
        leaf_node->right_child = nullptr;
        leaf_node->sequence = sequence;
        nodes.push_back(leaf_node);

        sequences.push_back(sequence);
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
    std::cout << "Current distance matrix:" << std::endl;
    for (int i = 0; i < num_sequences; ++i) {
        for (int j = 0; j < num_sequences; ++j) {
            std::cout << (*distance_matrix)[i][j] << " ";
        }
        std::cout << std::endl;
    }
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
    if (min_i != -1 && min_j != -1) {
        std::cout << "Smallest distance: " << min_distance << " between nodes " << nodes[min_i]->id << " and " << nodes[min_j]->id << std::endl;
    }
}

void NeighbourJoining::join_smallest_distance_nodes() {
    int num_nodes = distance_matrix->size();
    if (num_nodes <= 1) return;  // No hay suficientes nodos para proceder

    std::pair<int, int> smallest_distance_pair = find_smallest_distance_pair(num_nodes);
    int min_i = smallest_distance_pair.first;
    int min_j = smallest_distance_pair.second;

    if (min_i == -1 || min_j == -1) {
        // No se encontraron pares válidos para fusionar
        return;
    }

    // Crear un nuevo nodo y actualizar la matriz de distancias
    Node* new_node = create_new_node(min_i, min_j);
    nodes.push_back(new_node);  // Añadir el nuevo nodo a la lista
    // Contar solo los nodos activos
    int new_num_nodes = std::count_if(nodes.begin(), nodes.end(), [](const Node* n) { return n->active; });
    // Crear una nueva matriz de distancias
    auto new_distance_matrix = create_new_distance_matrix(num_nodes, new_num_nodes, min_i, min_j);
    distance_matrix = std::move(new_distance_matrix);

    // Imprimir el estado de los nodos después de la fusión
    for (int i = 0; i < nodes.size(); ++i) {
        if (nodes[i]->active) {
            std::cout << "Node Active: " << nodes[i]->id << std::endl;
        } else {
            std::cout << "Node Inactive: " << nodes[i]->id << std::endl;
        }
    }
}


std::pair<int, int> NeighbourJoining::find_smallest_distance_pair(int num_nodes) {
    int min_i = -1;
    int min_j = -1;
    int min_distance = 1000000;

    // El problema es que si añadimos nodos, num_nodes aumenta, por lo que realmente tendríamos que comprobar
    // con toda la lista de nodos.
    // El problema es que al usar nodes.size() estamos accediendo a zonas fuera de la matriz, ya que existe la 
    // posición 4 en el vector de nodos, pero la matriz es 3x3

    for (int i = 0; i < nodes.size(); ++i) {
        if (!nodes[i]->active) continue;  // Saltar nodos inactivos

        for (int j = i + 1; j < nodes.size(); ++j) {
            if (!nodes[j]->active) continue;  // Saltar nodos inactivos

            int distance = (*distance_matrix)[i][j];
            if (distance < min_distance) {
                min_distance = distance;
                min_i = i;
                min_j = j;
            }
        }
    }
    if (min_i == -1 || min_j == -1) {
        // No se encontraron pares válidos para fusionar
        return std::make_pair(-1, -1);
    }

    std::cout << "Smallest Distance Pair: (" << min_i << ", " << min_j << ") with distance " << min_distance << std::endl;
    return std::make_pair(min_i, min_j);
}


NeighbourJoining::Node* NeighbourJoining::create_new_node(int min_i, int min_j) {
    Node* new_node = new Node();
    new_node->id = nodes[min_i]->id + "-" + nodes[min_j]->id;
    new_node->left_child = nodes[min_i];
    new_node->right_child = nodes[min_j];
    new_node->sequence = "";
    new_node->active = true;

    // Desactivar los nodos antiguos
    nodes[min_i]->active = false;
    nodes[min_j]->active = false;

    return new_node;
}

// Asegurando el cálculo correcto de las nuevas distancias
std::unique_ptr<std::vector<std::vector<int>>> NeighbourJoining::create_new_distance_matrix(int num_nodes, int new_num_nodes, int min_i, int min_j) {
    auto new_distance_matrix = std::make_unique<std::vector<std::vector<int>>>(new_num_nodes, std::vector<int>(new_num_nodes, 0));
    std::vector<int> active_indices;

    // Recolectar todos los índices activos, excepto los dos que se van a fusionar
    for (int i = 0; i < num_nodes; ++i) {
        if (nodes[i]->active && i != min_i && i != min_j) {
            active_indices.push_back(i);
        }
    }

    // Añadir el índice del nuevo nodo al final para representarlo como último elemento en la matriz
    active_indices.push_back(num_nodes);  // Nuevo nodo es el último en la lista de nodos

    // Rellenar la nueva matriz de distancias
    for (int i = 0; i < new_num_nodes; ++i) {
        for (int j = 0; j < new_num_nodes; ++j) {
            if (i == new_num_nodes - 1 || j == new_num_nodes - 1) {
                // Calculando las distancias para el nuevo nodo fusionado
                if (i < new_num_nodes - 1 || j < new_num_nodes - 1) {
                    int other_index = active_indices[(i == new_num_nodes - 1) ? j : i];
                    (*new_distance_matrix)[i][j] = ((*distance_matrix)[other_index][min_i] + (*distance_matrix)[other_index][min_j] - (*distance_matrix)[min_i][min_j]) / 2;
                    (*new_distance_matrix)[j][i] = (*new_distance_matrix)[i][j];  // Asegurar simetría
                }
            } else {
                // Copiar las distancias existentes entre otros nodos activos
                (*new_distance_matrix)[i][j] = (*distance_matrix)[active_indices[i]][active_indices[j]];
            }
        }
    }
    return new_distance_matrix;
}

// Verificar la condición de parada en el método build_tree
void NeighbourJoining::build_tree() {
    calculate_distance_matrix();
    std::cout << "Matriz de distancias inicial en build_tree:" << std::endl;
    print_distance_matrix();

    while (distance_matrix->size() > 1) {
        join_smallest_distance_nodes();
        std::cout << "Matriz de distancias después de una fusión:" << std::endl;
        for (int i = 0; i < nodes.size(); ++i) {
            if (nodes[i]->active) {
                std::cout << "Node Active: " << nodes[i]->id << std::endl;
            } else {
                std::cout << "Node Inactive: " << nodes[i]->id << std::endl;
            }
        }
        print_distance_matrix();
    }

    std::cout << "Construcción del árbol completada." << std::endl;
    if (!nodes.empty()) {
        Node* root = nodes.back();  // Asumiendo que el último nodo es la raíz
        std::string alignment_order = get_alignment_order(root);
        std::cout << "Orden de alineamiento: " << alignment_order << std::endl;
    }
}

std::string NeighbourJoining::get_alignment_order(Node* node) {
    if (node == nullptr) {
        return "";
    }
    std::string left_order = get_alignment_order(node->left_child);
    std::string right_order = get_alignment_order(node->right_child);

    if (node->left_child == nullptr && node->right_child == nullptr) {
        return node->id;
    }

    return left_order + (left_order.empty() ? "" : "->") + right_order;
}


