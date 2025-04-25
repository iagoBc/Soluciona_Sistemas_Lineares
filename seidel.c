#include "gauss.h"

int gauss_seidel(real_t **A, real_t *b, real_t *x, int_t n) {
    real_t x_antigo[n];

    // Chute inicial
    for (int_t i = 0; i < n; i++) {
        x[i] = 0.0;
    }

    int iteracao = 0;
    while (iteracao < MAX_IT) {
        // Salva o vetor x anterior
        for (int_t i = 0; i < n; i++) {
            x_antigo[i] = x[i];
        }

        // Iteração de Gauss-Seidel
        for (int_t i = 0; i < n; i++) {
            real_t soma = 0.0;
            for (int_t k = 0; k < n; k++) {
                if (k != i) soma += A[i][k] * x[k];
            }
            x[i] = (b[i] - soma) / A[i][i];
        }

        // Cálculo da norma máxima do erro
        real_t erro_max = 0.0;
        for (int_t i = 0; i < n; i++) {
            real_t erro_i = fabs(x[i] - x_antigo[i]);
            if (erro_i > erro_max) erro_max = erro_i;
        }

        if (erro_max < TOL) break;

        iteracao++;
    }

    return iteracao;
}