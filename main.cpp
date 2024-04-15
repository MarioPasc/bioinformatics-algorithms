#include <iostream>
#include <vector>
#include <string>
#include "graph.h" // Asegúrate de que la ruta sea correcta

int main() {
    //std::vector<std::string> reads = {"ATG", "TGC", "GCA", "CAT", "ATA"};
    std::vector<std::string> reads = {"ATGCTAGCAC"};
    int k = 3; // Tamaño del k-mer

    std::unordered_map<std::string, Node*> graph = buildGraph(reads, k);

    Node* startNode = findStartingNode(graph);
    if (startNode == nullptr) {
        std::cout << "No se encontró un nodo de inicio adecuado." << std::endl;
        return 1; // Termina el programa si no se encuentra un nodo de inicio
    } else {
        std::cout << "Nodo de inicio: " << startNode->kmer << std::endl;
    }
    std::cout << "Inicio del algoritmo de Fleury" << std::endl;
    Fleury(startNode);
    return 0;
}

