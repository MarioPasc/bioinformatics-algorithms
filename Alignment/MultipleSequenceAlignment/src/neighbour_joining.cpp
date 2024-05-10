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
        //std::cout << "Smallest distance: " << min_distance << " between nodes " << nodes[min_i]->id << " and " << nodes[min_j]->id << std::endl;
    }
}

void NeighbourJoining::join_smallest_distance_nodes() {
    int num_nodes = distance_matrix->size();
    if (num_nodes <= 1) return;  // No hay suficientes nodos para proceder
    update_active_nodes();

    std::pair<int, int> smallest_distance_pair = find_smallest_distance_pair();
    int min_i = smallest_distance_pair.first;
    int min_j = smallest_distance_pair.second;

    if (min_i == -1 || min_j == -1) {
        // No se encontraron pares válidos para fusionar
        return;
    }

    // Crear un nuevo nodo y actualizar la matriz de distancias
    create_new_node(min_i, min_j);


    // Contar solo los nodos activos
    int new_num_nodes = std::count_if(nodes.begin(), nodes.end(), [](const Node* n) { return n->active; });
    // Crear una nueva matriz de distancias
    auto new_distance_matrix = create_new_distance_matrix(min_i, min_j);
    distance_matrix = std::move(new_distance_matrix);
    /*
    // Imprimir el estado de los nodos después de la fusión
    for (int i = 0; i < nodes.size(); ++i) {
        if (nodes[i]->active) {
            std::cout << "Node Active: " << nodes[i]->id << std::endl;
        } else {
            std::cout << "Node Inactive: " << nodes[i]->id << std::endl;
        }
    }
    */
}

void NeighbourJoining::update_active_nodes() {
    active_nodes.clear();  // Limpiar la lista anterior de nodos activos
    for (Node* node : nodes) {
        if (node->active) {
            active_nodes.push_back(node);  // Añadir solo nodos activos a la lista
        }
    }
    /*
    for (Node* node: active_nodes) {
        std::cout << "Active Node: " << node->id << std::endl;
    }
    */
}

std::pair<int, int> NeighbourJoining::find_smallest_distance_pair() {
    int min_i = -1;
    int min_j = -1;
    int min_distance = std::numeric_limits<int>::max();
    
    for (int i = 0; i < active_nodes.size(); ++i) {
        for (int j = i + 1; j < active_nodes.size(); ++j) {
            if (active_nodes[i] && active_nodes[j]) { // Verificar que los punteros son válidos
                int distance = (*distance_matrix)[i][j];
                if (distance < min_distance) {
                    min_distance = distance;
                    min_i = i;
                    min_j = j;
                }
            }
        }
    }

    if (min_i == -1 || min_j == -1) {
        std::cout << "No valid pairs found to merge." << std::endl;
        return std::make_pair(-1, -1);
    }

    //std::cout << "Smallest Distance Pair: (" << active_nodes[min_i]->id << ", " << active_nodes[min_j]->id << ") with distance " << min_distance << std::endl;
    //std::cout << "Smallest Distance Pair: (" << min_i << ", " << min_j << ") with distance " << min_distance << std::endl;
    return std::make_pair(min_i, min_j);
}



void NeighbourJoining::create_new_node(int min_i, int min_j) {
    // Hay que usar la lista de nodos activos para la creación del nuevo. 
    Node* new_node = new Node();
    new_node->id = active_nodes[min_i]->id + "-" + active_nodes[min_j]->id;
    new_node->left_child = active_nodes[min_i];
    new_node->right_child = active_nodes[min_j];
    new_node->sequence = "";
    new_node->active = true;
    new_node->depth = std::max(new_node->left_child->depth, new_node->right_child->depth) + 1;

    // Desactivar los nodos antiguos
    new_node->left_child->active = false;
    new_node->right_child->active = false;
    nodes.push_back(new_node);  // Añadir el nuevo nodo a la lista
    update_active_nodes(); // Ahora que hemos creado el nuevo nodo, actualizamos la lista.
}

std::unique_ptr<std::vector<std::vector<int>>> NeighbourJoining::create_new_distance_matrix(int min_i, int min_j) {
    int new_num_nodes = active_nodes.size();
    auto new_distance_matrix = std::make_unique<std::vector<std::vector<int>>>(new_num_nodes, std::vector<int>(new_num_nodes, 0));

    // Rellenar la nueva matriz de distancias
    for (int i = 0; i < new_num_nodes - 1; ++i) {  // -1 porque el último nodo será el nuevo nodo fusionado
        for (int j = 0; j < new_num_nodes - 1; ++j) {
            if (i == min_i || i == min_j || j == min_i || j == min_j) {
                // No calcular distancias que involucren directamente a los nodos fusionados aquí
                continue;
            }
            // Copiar las distancias existentes entre otros nodos activos
            (*new_distance_matrix)[i][j] = (*distance_matrix)[i][j];
        }
    }

    // Calcular las distancias del nuevo nodo fusionado
    int new_node_index = new_num_nodes - 1;
    for (int i = 0; i < new_num_nodes - 1; ++i) {
        // Utilizando la fórmula de Neighbour Joining para calcular la distancia al nuevo nodo
        (*new_distance_matrix)[i][new_node_index] = ((*distance_matrix)[i][min_i] + (*distance_matrix)[i][min_j] - (*distance_matrix)[min_i][min_j]) / 2;
        (*new_distance_matrix)[new_node_index][i] = (*new_distance_matrix)[i][new_node_index];  // Mantener la simetría
    }

    return new_distance_matrix;
}

/*
Esta es la función principal, que va creando el árbol
*/
void NeighbourJoining::build_tree() {
    calculate_distance_matrix();
    while (distance_matrix->size() > 1) {
        join_smallest_distance_nodes();
    }

    std::cout << "Construcción del árbol completada." << std::endl;
    if (!nodes.empty()) {
        Node* root = nodes.back();  // Asumiendo que el último nodo es la raíz
        std::cout << "Árbol filogenético:" << std::endl;
        print_tree(root, "", false);
        std::cout << "Árbol filogenético en formato Newick:" << std::endl;
        std::string newick_format = generate_newick_format(root);
        std::cout << newick_format << ";" << std::endl;  // Añade el punto y coma final necesario en formato Newick
        std::cout << "===================" << std::endl;
        std::string alignment = align_sequences();
        std::cout << "Alineamiento final:" << std::endl;
        std::cout << alignment << std::endl;
    }
}


/*
Conseguir el orden del alineamiento
*/
std::vector<NeighbourJoining::Node*> NeighbourJoining::get_alignment_order(Node* node) {
    if (node == nullptr) {
        return {};
    }

    std::vector<Node*> left_order = get_alignment_order(node->left_child);
    std::vector<Node*> right_order = get_alignment_order(node->right_child);

    if (node->left_child == nullptr && node->right_child == nullptr) {
        return {node};
    }

    left_order.insert(left_order.end(), right_order.begin(), right_order.end());
    return left_order;
}

/*
Imprimir el árbol
*/
void NeighbourJoining::print_tree(Node* node, std::string prefix, bool is_left) {
    if (node == nullptr) {
        return;
    }

    std::cout << prefix;
    std::cout << (is_left? "├──" : "└──");
    std::cout << node->id << std::endl;

    print_tree(node->left_child, prefix + (is_left? "│   " : "    "), true);
    print_tree(node->right_child, prefix + (is_left? "│   " : "    "), false);
}

/*
Realizar el alineamiento final, en proceso de realización...
*/
std::string NeighbourJoining::align_sequences() {
    std::vector<Node*> order = get_alignment_order(nodes.back()); // Asumiendo que el último nodo es la raíz
    std::string alignment = order[0]->sequence; // Inicia el alineamiento con la secuencia del primer nodo

    for (size_t i = 1; i < order.size(); ++i) {
        std::string next_sequence = order[i]->sequence;
        NeedlemanWunsch nw(alignment, next_sequence, 3, -1, -2); // Ajusta los parámetros según sea necesario
        nw.align();
        auto aligned_sequences = nw.get_alignment();
        alignment = aligned_sequences.first; // Usa la secuencia alineada de 'alignment' como la nueva secuencia base
        std::cout << alignment << std::endl;
    }

    return alignment;
}

/*
Función para generar la representación Newick del árbol
*/
std::string NeighbourJoining::generate_newick_format(Node* node) {
    if (node == nullptr) {
        return "";
    }

    std::string result;
    if (node->left_child || node->right_child) {
        result += "(";
        if (node->left_child) {
            result += generate_newick_format(node->left_child);
        }
        if (node->right_child) {
            result += ",";
            result += generate_newick_format(node->right_child);
        }
        result += ")";
    }
    result += node->id;
    // Usa la profundidad dividida por 10 como una estimación de la distancia
    result += ":" + std::to_string(node->depth / 10.0);
    return result;
}




