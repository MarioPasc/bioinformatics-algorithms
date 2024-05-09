#include "neighbour_joining.h"
#include <iostream>
#include <limits>

NeighbourJoining::NeighbourJoining(const std::vector<std::string>& sequences) : sequences(sequences) {
    int num_sequences = sequences.size();
    distance_matrix.resize(num_sequences, std::vector<int>(num_sequences, 0));
    nodes.resize(num_sequences);
    for (int i = 0; i < num_sequences; ++i) {
        nodes[i].id = i;
        nodes[i].left_child = -1;
        nodes[i].right_child = -1;
    }
}

void NeighbourJoining::calculate_distance_matrix() {
    int num_sequences = sequences.size();
    for (int i = 0; i < num_sequences; ++i) {
        for (int j = i + 1; j < num_sequences; ++j) {
            NeedlemanWunsch nw(sequences[i], sequences[j], 3, -1, -2);
            nw.align();
            int alignment_score = nw.get_alignment_score();
            distance_matrix[i][j] = alignment_score;
            distance_matrix[j][i] = alignment_score;
        }
    }
}

void NeighbourJoining::print_distance_matrix() const {
    int num_sequences = sequences.size();
    for (int i = 0; i < num_sequences; ++i) {
        for (int j = 0; j < num_sequences; ++j) {
            std::cout << distance_matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void NeighbourJoining::find_smallest_distance_node() const {
    int num_sequences = sequences.size();
    int min_distance = std::numeric_limits<int>::max();
    int min_i = -1;
    int min_j = -1;
    for (int i = 0; i < num_sequences; ++i) {
        for (int j = i + 1; j < num_sequences; ++j) {
            if (distance_matrix[i][j] < min_distance) {
                min_distance = distance_matrix[i][j];
                min_i = i;
                min_j = j;
            }
        }
    }
    std::cout << "Smallest distance: " << min_distance << " between nodes " << min_i << " and " << min_j << std::endl;
}
