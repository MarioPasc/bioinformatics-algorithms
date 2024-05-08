// needleman_wunsch.h

#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include <vector>
#include <string>

class NeedlemanWunsch {
public:
    // Constructor que toma las secuencias y los parámetros de puntuación.
    NeedlemanWunsch(const std::string& seq_a, const std::string& seq_b, 
                    int match_score, int mismatch_penalty, int gap_penalty);

    // Ejecutar el algoritmo de Needleman-Wunsch.
    void align();

    // Obtener el alineamiento óptimo tras ejecutar align().
    std::pair<std::string, std::string> get_alignment() const;
    void print_score_matrix() const;
    void print_trace_matrix() const;
    int get_alignment_score() const;

private:
    // Métodos para el algoritmo.
    void initialize_matrices();
    void calculate_scores_and_traces();
    void traceback_alignment();

    // Secuencias a alinear.
    std::string sequence_a;
    std::string sequence_b;

    // Matrices de puntuación y traza.
    std::vector<std::vector<int>> score_matrix;
    std::vector<std::vector<char>> trace_matrix;

    // Puntuaciones y penalizaciones.
    int match;
    int mismatch;
    int gap;

    // Alineamientos resultantes.
    std::string aligned_a;
    std::string aligned_b;
};

#endif // NEEDLEMAN_WUNSCH_H
