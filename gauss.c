#include "gauss.h"

int_t encontraMAX(real_t **A, int_t i, int_t n){
    real_t maior = A[i][i];
    int_t k = i;

    for(int j = 1; j < n; j++){
        if(A[j][i] > maior){
            maior = A[j][i];
            k = j;
        } 
    }

    return k; // Retorna a linha do maior
}

void trocaLinha(real_t **A, real_t *b, int_t i, int_t n, int_t iPivo){
    real_t pivo = A[0][iPivo];
    real_t *aux = (real_t*) malloc(sizeof(real_t) * n);
    int_t B = b[i];

    for(int j = 0; j < n; j++) aux[j] = A[i][j]; // Aux recebe a linha i
    for(int j = 0; j < n; j++) A[i][j] = A[iPivo][j]; // Coloca a linha do pivo no lugar da linha i
    for(int j = 0; j < n; j++) A[iPivo][j] = aux[j]; // Coloca a linha i no lugar da linha do pivo

    b[i] = b[iPivo];
    b[iPivo] = B;

}


void eliminacaoGauss(real_t **A, real_t *b, int_t n){
    for(int_t i = 0; i < n; i++){
        int_t iPivo = encontraMAX(A, i, n);
        if(i != iPivo) trocaLinha(A, b, i, n, iPivo);
        for(int_t k = i + 1; k < n; k++){ //Começa na linha seguinte de i
            real_t m = A[k][i] / A[i][i]; //Calcula o valor necessario para zerar embaixo
            A[k][i] = 0.0;
            for(int_t j = i + 1; j < n; j++) A[k][j] -= A[i][j] * m; // Calcula o valor dos coeficientes da linha de baixo, coluna a coluna
            b[k] -= b[i] * m; // Calcula o valor da costante da linha debaixo
        }
    }
}