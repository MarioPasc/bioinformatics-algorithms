#ifndef NEIGHBOUR_JOINING_H
#define NEIGHBOUR_JOINING_H

#include <vector>
#include <string>
#include "needleman_wunsch.h"
#include <memory>
class NeighbourJoining {
public:
    NeighbourJoining(const std::vector<std::string>& sequences);

    void calculate_distance_matrix();
    void print_distance_matrix() const;
    void find_smallest_distance_node() const;
    void join_smallest_distance_nodes();

private:

    struct Node {
        int id;
        Node* left_child;
        Node* right_child;
        std::string sequence;
    };

    std::vector<std::string> sequences;
    std::unique_ptr<std::vector<std::vector<int>>> distance_matrix;
    std::vector<Node> nodes;
    std::pair<int, int> find_smallest_distance_pair(int num_nodes);
    Node create_new_node(int min_i, int min_j);
    std::unique_ptr<std::vector<std::vector<int>>> create_new_distance_matrix(int num_nodes, int new_num_nodes, int min_i, int min_j);

};
#endif // NEIGHBOUR_JOINING_H
