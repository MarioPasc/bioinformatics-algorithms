#include "neighbour_joining.h"
#include <limits>
#include <algorithm>
#include <tuple>

NeighbourJoining::NeighbourJoining(const std::vector<std::string>& sequences) : sequences(sequences), root(nullptr) {
    calculate_distance_matrix();
    build_guide_tree();
}

void NeighbourJoining::calculate_distance_matrix() {
    int n = sequences.size();
    distance_matrix.resize(n, std::vector<int>(n, 0));

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int distance = calculate_distance(sequences[i], sequences[j]);
            distance_matrix[i][j] = distance;
            distance_matrix[j][i] = distance;
        }
    }
}

int NeighbourJoining::calculate_distance(const std::string& seq1, const std::string& seq2) {
    NeedlemanWunsch nw(seq1, seq2, 1, -1, -1);
    nw.align();
    return nw.get_alignment_score();
}

void NeighbourJoining::build_guide_tree() {
    std::vector<Node*> nodes;
    for (int i = 0; i < sequences.size(); ++i) {
        nodes.push_back(new Node{sequences[i], i, nullptr, nullptr});
    }

    while (nodes.size() > 1) {
        Node* node_i;
        Node* node_j;
        std::tie(node_i, node_j) = find_nearest_neighbours(nodes);
        join_neighbours(node_i, node_j, nodes);
    }

    root = nodes[0];
}

std::pair<NeighbourJoining::Node*, NeighbourJoining::Node*> NeighbourJoining::find_nearest_neighbours(std::vector<Node*>& nodes) {
    int n = nodes.size();
    int min_distance = std::numeric_limits<int>::max();
    Node* node_i = nullptr;
    Node* node_j = nullptr;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int distance = distance_matrix[nodes[i]->id][nodes[j]->id];
            if (distance < min_distance) {
                min_distance = distance;
                node_i = nodes[i];
                node_j = nodes[j];
            }
        }
    }

    return {node_i, node_j};
}

void NeighbourJoining::join_neighbours(Node* node_i, Node* node_j, std::vector<Node*>& nodes) {
    // Crear un nuevo nodo que combine dos nodos más cercanos
    Node* new_node = new Node{"", static_cast<int>(nodes.size()), node_i, node_j};
    nodes.push_back(new_node);

    // Actualizar la matriz de distancias para incluir el nuevo nodo
    update_distance_matrix(node_i->id, node_j->id);

    // Eliminar los nodos antiguos de la lista
    auto remove_node = [&](Node* node) {
        auto it = std::remove(nodes.begin(), nodes.end(), node);
        nodes.erase(it, nodes.end());
    };
    remove_node(node_i);
    remove_node(node_j);
}

void NeighbourJoining::update_distance_matrix(int i, int j) {
    int n = distance_matrix.size();
    std::vector<std::vector<int>> new_distance_matrix;

    // Crear una nueva matriz de distancias sin las filas y columnas i, j
    for (int k = 0; k < n; ++k) {
        if (k != i && k != j) {
            std::vector<int> new_row;
            for (int l = 0; l < n; ++l) {
                if (l != i && l != j) {
                    new_row.push_back(distance_matrix[k][l]);
                }
            }
            new_distance_matrix.push_back(new_row);
        }
    }

    // Agregar la nueva fila y columna para el nuevo nodo combinado
    std::vector<int> new_distances;
    for (int k = 0; k < new_distance_matrix.size(); ++k) {
        int distance = (distance_matrix[i][k] + distance_matrix[j][k] - distance_matrix[i][j]) / 2;
        new_distances.push_back(distance);
        new_distance_matrix[k].push_back(distance);  // Añadir nueva columna al final de cada fila
    }
    new_distances.push_back(0);  // Distancia del nuevo nodo a sí mismo
    new_distance_matrix.push_back(new_distances); // Añadir la nueva fila

    // Reemplazar la matriz de distancias vieja con la nueva
    distance_matrix = new_distance_matrix;
}



std::string NeighbourJoining::multiple_alignment() {
    return align(root);
}

std::string NeighbourJoining::align(Node* node) {
    if (node->left == nullptr && node->right == nullptr) {
        return node->sequence;
    }

    std::string left_alignment = align(node->left);
    std::string right_alignment = align(node->right);

    NeedlemanWunsch nw(left_alignment, right_alignment, 1, -1, -1);
    nw.align();
    auto [aligned_left, aligned_right] = nw.get_alignment();

    return aligned_left + "\n" + aligned_right;
}

NeighbourJoining::~NeighbourJoining() {
    delete_tree(root);
}

void NeighbourJoining::delete_tree(Node* node) {
    if (node == nullptr) {
        return;
    }
    delete_tree(node->left);
    delete_tree(node->right);
    delete node;
}

NeighbourJoining::Node* NeighbourJoining::get_root() const {
    return root;
}

std::vector<std::vector<int>> NeighbourJoining::get_distance_matrix() const {
    return distance_matrix;
}