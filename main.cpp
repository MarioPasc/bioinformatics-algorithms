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

    // Verificar si el archivo se abrió correctamente
    if (!iss) {
        std::cerr << "Error al abrir el archivo: " << filename << std::endl;
        return sequences;
    }

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

    std::unordered_map<std::string, Node*> graph = buildGraph(sequences, 5);

    getKmerFrequency(sequences, graph);
}

int main_assembly(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Uso: " << argv[0] << " <kmerfreq|testFleury> <ruta_archivo|modo_prueba>" << std::endl;
        return 1;
    }

    std::string mode = argv[1];
    std::string input = argv[2];

    if (mode == "kmerfreq") {
        // Llama directamente a la función con la ruta proporcionada
        calculateKmerFrequencyFastq(input);

    } else if (mode == "testFleury") {
        std::vector<std::string> reads;
        int k = 3;

        if (input == "readsPerfectEulerian") {
            reads = {"AGT", "GTA", "TAG", "AGT"};
            runTest(reads, k, "Perfect Eulerian Cycle");
        } else if (input == "readsNonEulerian") {
            reads = {"AGT", "GTC", "TCA", "CAT"};
            runTest(reads, k, "No Eulerian Cycle");
        } else if (input == "readsEulerianWithDeadEnds") {
            reads = {"AGT", "GTA", "TAC", "ACT", "CTG", "TGA", "ACG", "GAG", "CTT"};
            runTest(reads, k, "Eulerian Cycle with Extras");
        } else if (input == "laboratory"){
            reads = {"ATGCTAGCAC"};
            runTest(reads, k, "Assembly Lab Reads");
        } else {
            std::cout << "Modo de prueba inválido." << std::endl;
            return 1;
        }
    } else {
        std::cout << "Modo inválido. Use 'kmerfreq' o 'testFleury'." << std::endl;
        return 1;
    }
}

int main() {
    

    return 0;
}

