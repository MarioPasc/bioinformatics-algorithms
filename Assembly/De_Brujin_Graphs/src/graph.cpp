#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>

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

Node* findStartingNode(const unordered_map<string, Node*>& graph) {
    Node* startNode = nullptr;
    int maxDegree = 0;  // Track the maximum degree

    for (const auto& pair : graph) {
        int degree = pair.second->neighbors.size();
        if (degree > maxDegree) {
            maxDegree = degree;
            startNode = pair.second;
        }
    }
    return startNode;
}

int countEdges(const std::unordered_map<std::string, Node*>& graph) {
    int edgeCount = 0;

    for (const auto& pair : graph) {
        const Node* node = pair.second;
        for (const Node* neighbor : node->neighbors) {
            if (!neighbor->passed) {
                edgeCount++;
            }
        }
    }

    return edgeCount;
}

void Fleury(Node* start, std::unordered_map<std::string, Node*>& graph, Node* originalStart = nullptr) {
    if (originalStart == nullptr) {
        originalStart = start; // Establece el nodo de inicio original en la primera llamada
    }

    cout << start->kmer << endl; // Imprime el k-mer del nodo actual

    // Marca el nodo actual como visitado
    start->passed = true;

    // Verifica si se han recorrido todas las aristas
    int remainingEdges = countEdges(graph);

    if (remainingEdges == 0) {
        if (start == originalStart) {
            cout << "Se ha encontrado un ciclo de Euler." << endl;
        } else {
            cout << "Se ha encontrado un camino de Euler." << endl;
        }
        return;
    } 

    while (!start->neighbors.empty()) {
        Node* neighbor = start->neighbors.back();
        start->neighbors.pop_back();

        if (!neighbor->passed) {
            // Llama recursivamente a Fleury con el vecino como nuevo nodo de inicio
            Fleury(neighbor, graph, originalStart);
        }
    }
}
