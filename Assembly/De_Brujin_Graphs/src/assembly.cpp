#include <iostream>
#include "assembly.h"
#include "graph.h"

using namespace std;

void assembly(const vector<string>& reads, int k) {
    unordered_map<string, Node*> graph = buildGraph(reads, k);
    printGraph(graph);

    string eulerianPath = findEulerianPath(graph);
    cout << "Camino de Euler: " << eulerianPath << endl;
    cout << "Alineación: " << eulerianPath.substr(0, k) << eulerianPath.substr(k) << endl;
}