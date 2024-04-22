#include "needleman_wunsch.h"
#include <iostream>

int main() {
    // Secuencias de prueba.
    std::string seq1 = "TGGCATTCCGA";
    std::string seq2 = "GCCAATGAC";

    // Crear una instancia de NeedlemanWunsch.
    // Asumiendo que las penalizaciones de mismatch y gap son las mismas.
    NeedlemanWunsch nw(seq2, seq1, 3, -1, -2);

    // Ejecutar el algoritmo de alineamiento.
    nw.align();

    // Imprimir las matrices de puntuación y traza.
    std::cout << "Score Matrix:" << std::endl;
    nw.print_score_matrix();
    std::cout << std::endl << "Trace Matrix:" << std::endl;
    nw.print_trace_matrix();

    // Obtener y mostrar el alineamiento resultante.
    auto alignment = nw.get_alignment();
    std::cout << std::endl << "Alignment:" << std::endl;
    std::cout << "Sequence A: " << alignment.first << std::endl;
    std::cout << "Sequence B: " << alignment.second << std::endl;

    return 0;
}
