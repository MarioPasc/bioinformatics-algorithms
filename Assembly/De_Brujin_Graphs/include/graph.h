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
std::unordered_map<std::string, int> getKmerFrequency(const std::string& sequence, int k);
Node* findStartingNode(const std::unordered_map<std::string, Node*>& graph);
void printGraph(const std::unordered_map<std::string, Node*>& graph);
void Fleury(Node* start, Node* originalStart = nullptr); // Nueva declaración

#endif // GRAPH_H