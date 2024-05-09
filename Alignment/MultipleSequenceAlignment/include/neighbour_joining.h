#ifndef NEIGHBOUR_JOINING_H
#define NEIGHBOUR_JOINING_H

#include <vector>
#include <string>
#include "needleman_wunsch.h"
#include <memory>
#include <unordered_map>

class NeighbourJoining {
public:
    struct Node {
        std::string id;
        Node* left_child = nullptr;
        Node* right_child = nullptr;
        std::string sequence;
        bool active = true;
    };

    NeighbourJoining(const std::unordered_map<std::string, std::string>& sequence_map);

    void calculate_distance_matrix();
    void print_distance_matrix() const;
    void find_smallest_distance_node() const;
    void join_smallest_distance_nodes();
    void build_tree();
    std::string get_alignment_order(NeighbourJoining::Node* node);
private:
    std::vector<std::string> sequences;
    std::unique_ptr<std::vector<std::vector<int>>> distance_matrix;
    std::vector<Node*> nodes;  // Solo un vector de punteros

    std::pair<int, int> find_smallest_distance_pair(int num_nodes);
    NeighbourJoining::Node* create_new_node(int min_i, int min_j);
    std::unique_ptr<std::vector<std::vector<int>>> create_new_distance_matrix(int num_nodes, int new_num_nodes, int min_i, int min_j);
};

#endif 
// NEIGHBOUR_JOINING_H

