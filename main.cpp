#include <iostream>
#include <vector>
#include <string>
#include "graph.h" 
#include <fstream>
#include <vector>
#include <string>
#include <kseq++/seqio.hpp>

using namespace klibpp;

// Function to read FASTQ files and return a vector of sequences
std::vector<std::string> readFastqSequences(const std::string& filename) {
    std::vector<std::string> sequences;
    KSeq record;

    // Corrección: Convertir std::string a const char* usando c_str()
    SeqStreamIn iss(filename.c_str());  // Usa c_str() para convertir std::string a const char*

    // Leer cada registro del archivo FASTQ
    while (iss >> record) {
        sequences.push_back(record.seq);  // Almacenar solo la secuencia
    }

    return sequences;
}


void runTest(const std::vector<std::string>& reads, int k, const std::string& testName) {
    std::cout << "Test: " << testName << std::endl;
    
    // Construir el grafo
    std::unordered_map<std::string, Node*> graph = buildGraph(reads, k);

    // Imprimir el grafo
    std::cout << "Graph structure:" << std::endl;
    printGraph(graph);

    // Imprimir la frecuencia de cada k-mero
    getKmerFrequency(reads, graph);


    // Encontrar el circuito Euleriano si existe
    std::vector<Node*> circuit = fleuryAlgorithm(graph);
    
    std::cout << "Eulerian Circuit:";
    if (circuit.empty()) {
        std::cout << " No Eulerian Circuit found.";
    } else {
        for (Node* node : circuit) {
            std::cout << " " << node->kmer << " ->";
        }
    }
    std::cout << " END" << std::endl << std::endl;

    // Limpiar la memoria
    for (auto& pair : graph) {
        delete pair.second;
    }
}

void runTests() {
    std::vector<std::string> readsPerfectEulerian = {"AGT", "GTA", "TAG", "AGT"};
    std::vector<std::string> readsNonEulerian = {"AGT", "GTC", "TCA", "CAT"};
    std::vector<std::string> readsEulerianWithDeadEnds = {"AGT", "GTA", "TAG", "AGC", "GCT"};
    int k = 3;
    runTest(readsPerfectEulerian, k, "Perfect Eulerian Cycle");
    runTest(readsNonEulerian, k, "No Eulerian Cycle");
    runTest(readsEulerianWithDeadEnds, k, "Eulerian Cycle with Extras");
}

void calculateKmerFrequencyFastq(const std::string& fastqPath) {
    std::vector<std::string> sequences = readFastqSequences(fastqPath);

    std::unordered_map<std::string, Node*> graph = buildGraph(sequences, 3);

    getKmerFrequency(sequences, graph);
}

int main() {

    calculateKmerFrequencyFastq("/home/mariopasc/Downloads/DataSet_Trabajo_Casa/DataSet_Trabajo_Casa/ERR103404_2.fastq.gz");

    return 0;
}

