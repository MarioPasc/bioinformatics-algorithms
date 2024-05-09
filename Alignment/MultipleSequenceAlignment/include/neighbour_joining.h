#ifndef NEIGHBOUR_JOINING_H
#define NEIGHBOUR_JOINING_H

#include <vector>
#include <string>
#include "needleman_wunsch.h"
class NeighbourJoining {
public:
    struct Node {
        int id;
        int left_child;
        int right_child;
    };

    NeighbourJoining(const std::vector<std::string>& sequences);

    void calculate_distance_matrix();
    void print_distance_matrix() const;
    void find_smallest_distance_node() const;
    void join_smallest_distance_nodes();

private:
    std::vector<std::string> sequences;
    std::vector<std::vector<int>> distance_matrix;
    std::vector<Node> nodes;
};

#endif // NEIGHBOUR_JOINING_H
