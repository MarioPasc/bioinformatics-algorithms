#include <iostream>
#include <vector>
#include <string>
#include "assembly.h"

int main() {
    std::vector<std::string> reads = {"ATGCTAGCAC"};
    int k = 4;

    assembly(reads, k);

    return 0;
}
