#include <iostream>
#include "assembly.h"
#include "graph.h"

using namespace std;

void assembly(const vector<string>& reads, int k) {
    unordered_map<string, Node*> graph = buildGraph(reads, k);
    printGraph(graph);

    
}