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
        int depth = 0;
        bool active = true;
    };

    NeighbourJoining(const std::unordered_map<std::string, std::string>& sequence_map);

    void calculate_distance_matrix();
    void print_distance_matrix() const;
    void find_smallest_distance_node() const;
    void join_smallest_distance_nodes();
    void build_tree();
    std::vector<Node*> get_alignment_order(Node* node);
private:
    std::vector<std::string> sequences;
    std::unique_ptr<std::vector<std::vector<int>>> distance_matrix;
    std::vector<Node*> nodes;  
    std::vector<Node*> active_nodes;
    void update_active_nodes();
    std::pair<int, int> find_smallest_distance_pair();
    void create_new_node(int min_i, int min_j);
    std::unique_ptr<std::vector<std::vector<int>>> create_new_distance_matrix(int min_i, int min_j);
    void print_tree(Node* node, std::string prefix = "", bool is_left = false);
    std::string align_sequences();
    std::string generate_newick_format(Node* node);
};

#endif 
// NEIGHBOUR_JOINING_H

