#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include "assembly.h"
#include "graph.h"
using namespace std;

int main() {
    std::vector<std::string> reads = {"ATGCTAGCAC"};
    int k = 3;

    assembly(reads, k); 

    for (const string& read: reads) {
        std::unordered_map<std::string, int> kmerFrequency = getKmerFrequency(read, k);
        // Print the k-mers and their frequencies
        for (const auto& pair : kmerFrequency) {
            std::cout << "K-mer: " << pair.first << ", Frequency: " << pair.second << std::endl;
        }
    }
    return 0;
}

