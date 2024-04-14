#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>

using namespace std;

struct Node {
    string kmer;
    vector<Node*> neighbors;
};


unordered_map<string, Node*> buildGraph(const vector<string>& reads, int k) {
    unordered_map<string, Node*> graph;

    for (const string& read : reads) {
        for (int i = 0; i <= read.length() - k; i++) {
            string kmer = read.substr(i, k);
            if (graph.find(kmer) == graph.end()) {
                graph[kmer] = new Node{kmer, {}};
            }
            if (i > 0) {
                string prev_kmer = read.substr(i - 1, k);
                graph[prev_kmer]->neighbors.push_back(graph[kmer]);
            }
        }
    }
    return graph;
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

string findEulerianPath(unordered_map<string, Node*>& graph) {
    string path;
    unordered_map<Node*, int> inDegree, outDegree;

    // Calcular los grados de entrada y salida de cada nodo
    for (auto& pair : graph) {
        Node* node = pair.second;
        outDegree[node] = node->neighbors.size();
        for (Node* neighbor : node->neighbors) {
            inDegree[neighbor]++;
        }
    }

    // Encontrar el nodo de inicio (con grado de salida mayor que grado de entrada)
    Node* startNode = nullptr;
    for (auto& pair : graph) {
        Node* node = pair.second;
        if (outDegree[node] > inDegree[node]) {
            startNode = node;
            break;
        }
    }

    // Realizar el recorrido del camino de Euler
    Node* currentNode = startNode;
    while (outDegree[currentNode] > 0) {
        path += currentNode->kmer[0]; // equivalente a (*currentNode).kmer[0]
        Node* nextNode = currentNode->neighbors.back();
        currentNode->neighbors.pop_back();
        outDegree[currentNode]--;
        currentNode = nextNode;
    }
    path += currentNode->kmer;

    return path;
}