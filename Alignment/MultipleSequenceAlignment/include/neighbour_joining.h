#ifndef NEIGHBOUR_JOINING_H
#define NEIGHBOUR_JOINING_H

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include "needleman_wunsch.h"

class NeighbourJoining {
public:
    struct Node {
        std::string id;  // Identifier for the node, typically used for leaf nodes
        Node* left_child = nullptr;  // Pointer to the left child in the tree
        Node* right_child = nullptr;  // Pointer to the right child in the tree
        std::string sequence;  // The DNA or protein sequence represented by this node
        int depth = 0;  // Depth of the node in the tree, useful for visual representation
        bool active = true;  // Flag to indicate if the node is active in the current context
    };

    NeighbourJoining(const std::unordered_map<std::string, std::string>& sequence_map);
    void calculate_distance_matrix();  // Computes the pairwise distance matrix using Needleman-Wunsch
    void print_distance_matrix() const;  // Outputs the current distance matrix to the console
    void join_smallest_distance_nodes();  // Merges the two nodes with the smallest distance
    void build_tree();  // Constructs the phylogenetic tree by iteratively merging nodes
    std::vector<Node*> get_alignment_order(Node* node);  // Returns the order of nodes for alignment

private:
    std::vector<std::string> sequences;  // Stores the original sequences
    std::unique_ptr<std::vector<std::vector<int>>> distance_matrix;  // Matrix of distances between sequences
    std::vector<Node*> nodes;  // Vector of nodes corresponding to the sequences
    std::vector<Node*> active_nodes;  // Nodes that are currently active and part of the ongoing tree construction
    void update_active_nodes();  // Updates the list of active nodes after a merging event
    std::pair<int, int> find_smallest_distance_pair();  // Finds the indices of the closest pair of nodes
    void create_new_node(int min_i, int min_j);  // Creates a new node by merging two nodes
    std::unique_ptr<std::vector<std::vector<int>>> create_new_distance_matrix(int min_i, int min_j);  // Creates a new distance matrix post node merging
    void print_tree(Node* node, std::string prefix = "", bool is_left = false);  // Helper function to print the tree for debugging
    std::string align_sequences();  // Function to perform the final sequence alignment
    std::string generate_newick_format(Node* node);  // Generates a Newick format string for the tree
};

#endif // NEIGHBOUR_JOINING_H
