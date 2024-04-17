// Graph.h

#ifndef GRAPH_H
#define GRAPH_H

#include <unordered_map>
#include <vector>
#include <string>

struct Node; 

// Estructura para representar las aristas del grafo
struct Edge {
    Node* from;   
    Node* to;     
    bool passed;

    Edge(Node* f, Node* t);
};

// Estructura para representar los nodos del grafo
struct Node {
    std::string kmer;         
    std::vector<Edge*> edges;

    Node(std::string k);
    ~Node();
};

void addEdge(Node* fromNode, Node* toNode, std::unordered_map<std::string, Node*>& graph);

std::unordered_map<std::string, Node*> buildGraph(const std::vector<std::string>& reads, int k);

void printGraph(const std::unordered_map<std::string, Node*>& graph);

std::vector<Node*> fleuryAlgorithm(std::unordered_map<std::string, Node*>& graph);

void getKmerFrequency(const std::vector<std::string>& reads, std::unordered_map<std::string, Node*>& graph);

#endif // GRAPH_H
