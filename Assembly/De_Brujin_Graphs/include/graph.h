#ifndef GRAPH_H
#define GRAPH_H

#include <unordered_map>
#include <vector>
#include <string>

struct Node {
    std::string kmer;
    std::vector<Node*> neighbors;
};

std::unordered_map<std::string, Node*> buildGraph(const std::vector<std::string>& reads, int k);
std::string findEulerianPath(std::unordered_map<std::string, Node*>& graph);

#endif
