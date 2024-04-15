#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>

using namespace std;

struct Node {
    string kmer;
    vector<Node*> neighbors;
    bool passed = false;
};

unordered_map<string, Node*> buildGraph(const vector<string>& reads, int k) {
    unordered_map<string, Node*> graph;

    for (const string& read : reads) {
        for (int i = 0; i <= read.length() - k + 1; i++) {
            string kmer = read.substr(i, k - 1);
            
            // Check if the current k-mer is already in the graph
            if (graph.find(kmer) == graph.end()) {
                // If not, create a new node with the current k-mer and add it to the graph
                graph[kmer] = new Node{kmer, {}};
            }
            
            // Check if there is a previous k-mer
            if (i > 0) {
                string prev_kmer = read.substr(i - 1, k - 1);
                // Add an edge from the previous k-mer to the current k-mer
                // this is a pointer from the previous k-mer to the actual
                graph[prev_kmer]->neighbors.push_back(graph[kmer]);
            }
        }
    }
    // Return the constructed graph
    return graph;
}

std::unordered_map<string, int> getKmerFrequency(const std::string& sequence, int k) {
    std::unordered_map<string, int> kmerFrequency;

    for (int i = 0; i <= sequence.length() - k; ++i) {
        string kmer = sequence.substr(i ,k-1);
        
        kmerFrequency[kmer]++;
    }

    return kmerFrequency;
}


void printGraph(const std::unordered_map<std::string, Node*>& graph) {
    for (const auto& pair : graph) {
        const Node* node = pair.second;
        std::cout << "Nodo: " << node->kmer << std::endl;
        std::cout << "Vecinos: ";
        for (const Node* neighbor : node->neighbors) {
            std::cout << neighbor->kmer << " ";
        }
        std::cout << std::endl << std::endl;
    }
}

Node* findStartNode(unordered_map<string, Node*>& graph) {
    string path;
    unordered_map<Node*, int> inDegree, outDegree;

    // Calculate the in-degree and out-degree of each node
    for (auto& pair : graph) {
        Node* node = pair.second;
        outDegree[node] = node->neighbors.size();
        for (Node* neighbor : node->neighbors) {
            inDegree[neighbor]++;
        }
    }

    // Find the starting node (a node with out-degree greater than in-degree)
    Node* startNode = nullptr;
    for (auto& pair : graph) {
        Node* node = pair.second;
        if (outDegree[node] > inDegree[node]) {
            startNode = node;
            break;
        }
    }
    return startNode;
}

void findEulerianPath(Node* startNode, Node* exitNode, std::vector<std::string>& path) {
    // Check the start node as passed
    startNode->passed = true;
    // Base case: Actual node has no neighbours o
    // Caso base: si el nodo actual no tiene vecinos o el único vecino es el nodo de salida
    if (startNode->neighbors.empty() || (startNode->neighbors.size() == 1 && startNode->neighbors[0] == exitNode)) {
        // Agregar el k-mer del nodo actual al camino de Euler
        path.push_back(startNode->kmer);
        return;
    }

    // Recorrer los vecinos del nodo actual
    for (Node* neighbor : startNode->neighbors) {
        if (!neighbor->passed) {
            findEulerianPath(neighbor, exitNode, path);
        }
    }

    // Agregar el k-mer del nodo actual al camino de Euler
    path.push_back(startNode->kmer);
}
