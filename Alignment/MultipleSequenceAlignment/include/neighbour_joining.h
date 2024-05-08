#ifndef NEIGHBOUR_JOINING_H
#define NEIGHBOUR_JOINING_H

#include <vector>
#include <string>
#include "needleman_wunsch.h"
#include <tuple>

class NeighbourJoining {
private:
struct Node;

public:
    NeighbourJoining(const std::vector<std::string>& sequences);

    void calculate_distance_matrix();
    void build_guide_tree();
    std::string multiple_alignment();
    Node* get_root() const;
    std::vector<std::vector<int>> get_distance_matrix() const;
    ~NeighbourJoining();


private:
    struct Node {
        std::string sequence;
        int id;
        Node* left;
        Node* right;
    };

    std::vector<std::string> sequences;
    std::vector<std::vector<int>> distance_matrix;
    Node* root;

    void delete_tree(Node* node);
    std::pair<Node*, Node*> find_nearest_neighbours(std::vector<Node*>& nodes);
    void join_neighbours(NeighbourJoining::Node *node_i, NeighbourJoining::Node *node_j, std::vector<NeighbourJoining::Node*>& nodes);
    int calculate_distance(const std::string& seq1, const std::string& seq2);
    void update_distance_matrix(int i, int j);
    std::string align(Node* node);
};

#endif // NEIGHBOUR_JOINING_H