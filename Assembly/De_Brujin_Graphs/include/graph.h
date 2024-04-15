#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <unordered_map>

struct Node {
    std::string kmer;
    std::vector<Node*> neighbors;
    bool passed = false;
};

std::unordered_map<std::string, Node*> buildGraph(const std::vector<std::string>& reads, int k);
void printGraph(const std::unordered_map<std::string, Node*>& graph);
std::unordered_map<std::string, int> getKmerFrequency(const std::string& sequence, int k);
void findEulerianPath(Node* startNode, Node* exitNode, std::vector<std::string>& path);
Node* findStartNode(std::unordered_map<std::string, Node*>& graph)
#endif
