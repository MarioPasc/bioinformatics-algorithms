// needleman_wunsch.cpp
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include "smith_waterman.h"
using namespace std;

SmithWaterman::SmithWaterman(const std::string& seq_a, const std::string& seq_b, 
                                 int match_score, int mismatch_penalty, int gap_penalty): 
    sequence_a(seq_a),
    sequence_b(seq_b),
    match(match_score),
    mismatch(mismatch_penalty),
    gap(gap_penalty),
    score_matrix(seq_a.length() + 1, std::vector<int>(seq_b.length() + 1)),
    trace_matrix(seq_a.length() + 1, std::vector<char>(seq_b.length() + 1)) {
    // Inicializar las matrices con sus valores iniciales.
    initialize_matrices();
}

void SmithWaterman::initialize_matrices() {
    // Inicializar la primera fila.
    for (size_t j = 0; j < score_matrix[0].size(); ++j) {
        // score_matrix[0][j] = j * gap;
        score_matrix[0][j] = 0;
        trace_matrix[0][j] = 'L';  // Indica que viene de la izquierda (gap en secuencia A).
    }

    // Inicializar la primera columna.
    for (size_t i = 0; i < score_matrix.size(); ++i) {
        //score_matrix[i][0] = i * gap;
        score_matrix[i][0] = 0;
        trace_matrix[i][0] = 'U';  // Indica que viene de arriba (gap en secuencia B).
    }
}

void SmithWaterman::align() {
    calculate_scores_and_traces();

    traceback_alignment();
}

void SmithWaterman::calculate_scores_and_traces() {
    // Define la puntuación de similaridad basada en una matriz de puntuación.
    auto s = [this](char a, char b) -> int {
        static const std::unordered_map<char, std::unordered_map<char, int>> score_map{
            {'A', {{'A', match}, {'C', mismatch}, {'G', -mismatch}, {'T', mismatch}}},
            {'C', {{'A', mismatch}, {'C', match}, {'G', mismatch}, {'T', -mismatch}}},
            {'G', {{'A', -mismatch}, {'C', mismatch}, {'G', match}, {'T', mismatch}}},
            {'T', {{'A', mismatch}, {'C', -mismatch}, {'G', mismatch}, {'T', match}}},
        };
        return score_map.at(a).at(b);
    };

    for (size_t i = 1; i < score_matrix.size(); ++i) {
        for (size_t j = 1; j < score_matrix[i].size(); ++j) {
            int match_score = score_matrix[i - 1][j - 1] + s(sequence_a[i - 1], sequence_b[j - 1]);
            int delete_score = score_matrix[i - 1][j] + gap;
            int insert_score = score_matrix[i][j - 1] + gap;

            int max_score = std::max({0, match_score, delete_score, insert_score});  

            score_matrix[i][j] = max_score;

            // Actualizar la matriz de trazas.
            if (max_score == match_score) {
                trace_matrix[i][j] = 'D';  // Diagonal
            } else if (max_score == delete_score) {
                trace_matrix[i][j] = 'U';  // Arriba
            } else {
                trace_matrix[i][j] = 'L';  // Izquierda
            }
        }
    }
}

void SmithWaterman::traceback_alignment() {
    std::string alignA;
    std::string alignB;

    // Encontrar el valor máximo en la matriz para comenzar el trazado
    int max_value = 0;
    size_t max_i = 0, max_j = 0;
    for (size_t i = 0; i < score_matrix.size(); ++i) {
        for (size_t j = 0; j < score_matrix[i].size(); ++j) {
            if (score_matrix[i][j] > max_value) {
                max_value = score_matrix[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    // Iniciar el trazado desde la posición del valor máximo
    size_t i = max_i, j = max_j;

    // Continuar el trazado hasta llegar a un valor cero
    while (i > 0 && j > 0 && score_matrix[i][j] > 0) {
        if (trace_matrix[i][j] == 'D') {  // Diagonal
            alignA = sequence_a[i - 1] + alignA;
            alignB = sequence_b[j - 1] + alignB;
            --i;
            --j;
        } else if (trace_matrix[i][j] == 'U') {  // Arriba
            alignA = sequence_a[i - 1] + alignA;
            alignB = '-' + alignB;  // Indica un gap en B.
            --i;
        } else {  // Izquierda
            alignA = '-' + alignA;  // Indica un gap en A.
            alignB = sequence_b[j - 1] + alignB;
            --j;
        }
    }

    // Asignar los alineamientos a las variables miembro de la clase.
    aligned_a = alignA;
    aligned_b = alignB;
}


std::pair<std::string, std::string> SmithWaterman::get_alignment() const {
    return {aligned_a, aligned_b};
}


void SmithWaterman::print_score_matrix() const {
    for (const auto& row : score_matrix) {
        for (const auto& cell : row) {
            std::cout << std::setw(4) << cell;
        }
        std::cout << std::endl;
    }
}

void SmithWaterman::print_trace_matrix() const {
    for (const auto& row : trace_matrix) {
        for (const auto& cell : row) {
            std::cout << std::setw(4) << cell;
        }
        std::cout << std::endl;
    }
}

