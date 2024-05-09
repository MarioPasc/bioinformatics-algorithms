#include "neighbour_joining.h"
#include <iostream>
#include <limits>

NeighbourJoining::NeighbourJoining(const std::vector<std::string>& sequences) : sequences(sequences) {
    int num_sequences = sequences.size();
    distance_matrix = std::make_unique<std::vector<std::vector<int>>>(num_sequences, std::vector<int>(num_sequences, 0));
    nodes.resize(num_sequences);
    for (int i = 0; i < num_sequences; ++i) {
        nodes[i].id = i;
        nodes[i].left_child = nullptr;
        nodes[i].right_child = nullptr;
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
    int num_nodes = distance_matrix->size();
    std::pair<int, int> smallest_distance_pair = find_smallest_distance_pair(num_nodes);
    int min_i = smallest_distance_pair.first;
    int min_j = smallest_distance_pair.second;

    Node new_node = create_new_node(min_i, min_j);
    nodes.push_back(new_node);

    int new_num_nodes = num_nodes - 1;
    auto new_distance_matrix = create_new_distance_matrix(num_nodes, new_num_nodes, min_i, min_j);

    distance_matrix = std::move(new_distance_matrix);
}

std::pair<int, int> NeighbourJoining::find_smallest_distance_pair(int num_nodes) {
    int min_distance = std::numeric_limits<int>::max();
    int min_i = -1;
    int min_j = -1;

    for (int i = 0; i < num_nodes; ++i) {
        for (int j = i + 1; j < num_nodes; ++j) {
            if ((*distance_matrix)[i][j] < min_distance) {
                min_distance = (*distance_matrix)[i][j];
                min_i = i;
                min_j = j;
            }
        }
    }

    return std::make_pair(min_i, min_j);
}

NeighbourJoining::Node NeighbourJoining::create_new_node(int min_i, int min_j) {
    Node new_node;
    new_node.id = nodes.size();
    new_node.left_child = &nodes[min_i];
    new_node.right_child = &nodes[min_j];
    return new_node;
}

std::unique_ptr<std::vector<std::vector<int>>> NeighbourJoining::create_new_distance_matrix(int num_nodes, int new_num_nodes, int min_i, int min_j) {
    auto new_distance_matrix = std::make_unique<std::vector<std::vector<int>>>(new_num_nodes, std::vector<int>(new_num_nodes, 0));

    int new_index = 0;
    for (int i = 0; i < num_nodes; ++i) {
        if (i != min_i && i != min_j) {
            int new_jndex = 0;
            for (int j = 0; j < num_nodes; ++j) {
                if (j != min_i && j != min_j) {
                    (*new_distance_matrix)[new_index][new_jndex] = (*distance_matrix)[i][j];
                    new_jndex++;
                }
            }
            int distance = ((*distance_matrix)[i][min_i] + (*distance_matrix)[i][min_j] - (*distance_matrix)[min_i][min_j]) / 2;
            (*new_distance_matrix)[new_index][new_num_nodes - 1] = distance;
            (*new_distance_matrix)[new_num_nodes - 1][new_index] = distance;
            new_index++;
        }
    }

    return new_distance_matrix;
}