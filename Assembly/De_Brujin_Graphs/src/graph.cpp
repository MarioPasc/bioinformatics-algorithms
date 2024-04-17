#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include "graph.h"

using namespace std;

struct Node;

// Structures used to create the graph
Edge::Edge(Node* f, Node* t) : from(f), to(t), passed(false) {}

// Implementación de Node
Node::Node(std::string k) : kmer(std::move(k)) {}

Node::~Node() {
    for (Edge* edge : edges) {
        delete edge;
    }
}

// Graph generation from the given reads

void addEdge(Node* fromNode, Node* toNode, std::unordered_map<std::string, Node*>& graph) {
    // The 'find_if' function takes three parameters: the start iterator, the end iterator, and a predicate function.
    // The predicate function here is a lambda that captures 'toNode' by reference and checks if any of the edges
    // in the 'fromNode' edges list points to 'toNode'.
    auto it = std::find_if(fromNode->edges.begin(), fromNode->edges.end(),
        [&toNode](Edge* edge) { return edge->to == toNode; });

    if (it == fromNode->edges.end()) { // Si no se encuentra, añadir nueva arista
        Edge* edge = new Edge(fromNode, toNode);
        fromNode->edges.push_back(edge);
    }
}

std::unordered_map<std::string, Node*> buildGraph(const std::vector<std::string>& reads, int k) {
    std::unordered_map<std::string, Node*> graph;
    int k1 = k-1;
    for (const std::string& read : reads) {
        for (size_t i = 0; i <= read.length() - k1; ++i) {
            std::string k1mer = read.substr(i, k1);
            if (graph.find(k1mer) == graph.end()) {
                graph[k1mer] = new Node(k1mer);
            }

            if (i > 0) {
                std::string prev_k1mer = read.substr(i - 1, k1);
                if (graph.find(prev_k1mer) == graph.end()) {
                    graph[prev_k1mer] = new Node(prev_k1mer);
                }
                addEdge(graph[prev_k1mer], graph[k1mer], graph);
            }
        }
    }

    return graph;
}

void printGraph(const std::unordered_map<std::string, Node*>& graph) {
    for (const auto& pair : graph) {
        std::cout << "Node " << pair.first << " has edges to: ";
        for (const Edge* edge : pair.second->edges) {
            std::cout << edge->to->kmer << " ";
        }
        std::cout << std::endl;
    }
}

std::vector<Node*> fleuryAlgorithm(std::unordered_map<std::string, Node*>& graph) {
    std::vector<Node*> path;
    // Get the starting node
    Node* current = graph["AG"];

    while (true) {
        bool foundEdge = false;
        Edge* nextEdge = nullptr;

        // Look for an unexplored edge
        for (Edge* edge : current->edges) {
            // We wont check if the edge forms a bridge since we are not getting
            // really complicated reads as inputs 
            if (!edge->passed) { 
                nextEdge = edge;
                foundEdge = true;
                // Exit the for loop
                break;
            }
        }

        // If we havent found any edge, we've finished
        if (!foundEdge) {break;}

        // Mark this edge as passed by the algorithm
        nextEdge->passed = true;
        // Add the current node
        path.push_back(current);
        // Next current node will be the one that the edge points to
        current = nextEdge->to;

    }
    // Add the last node
    if (!path.empty() && path.back() != current) {
        path.push_back(current);
    }
    return path;

}



