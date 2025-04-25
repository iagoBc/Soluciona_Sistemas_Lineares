#include "gauss.h"
#include "seidel.h"

int main() {
    fesetround(FE_DOWNWARD); // Define arredondamento para baixo
    SL sist1;
    SL sist2;
    real_t *r1, *r2;

    scanf("%ld", &sist1.ordem); // Ordem do sistema
    sist2.ordem = sist1.ordem;

    // Aloca os ponteiros para linhas
    sist1.matriz = (real_t**) malloc(sizeof(real_t*) * sist1.ordem);
    sist2.matriz = (real_t**) malloc(sizeof(real_t*) * sist2.ordem);

    // Aloca os vetores
    sist1.B = (real_t*) malloc(sizeof(real_t) * sist1.ordem);
    sist1.x = (real_t*) malloc(sizeof(real_t) * sist1.ordem);
    sist2.B = (real_t*) malloc(sizeof(real_t) * sist2.ordem);
    sist2.x = (real_t*) malloc(sizeof(real_t) * sist2.ordem);
    r1 = (real_t*) malloc(sizeof(real_t) * sist1.ordem);
    r2 = (real_t*) malloc(sizeof(real_t) * sist2.ordem);

    // Verifica alocação
    if (sist1.matriz == NULL || sist1.B == NULL || sist1.x == NULL ||
        sist2.matriz == NULL || sist2.B == NULL || sist2.x == NULL) {
        printf("Erro ao alocar memoria\n");
        return 1;
    }

    // Aloca cada linha da matriz
    for (int i = 0; i < sist1.ordem; i++) {
        sist1.matriz[i] = (real_t*) malloc(sizeof(real_t) * sist1.ordem);
        sist2.matriz[i] = (real_t*) malloc(sizeof(real_t) * sist2.ordem);
    }

    // Lê os coeficientes e termos independentes
    for (int i = 0; i < sist1.ordem; i++) {
        for (int k = 0; k < (sist1.ordem + 1); k++) {
            if (k == sist1.ordem) {
                scanf("%lf", &sist1.B[i]);
                sist2.B[i] = sist1.B[i];
                r1[i] = sist1.B[i];
                r2[i] = sist1.B[i];
            } 
            
            else {
                scanf("%lf", &sist1.matriz[i][k]);
                sist2.matriz[i][k] = sist1.matriz[i][k];
            }
        }
    }

    printf("\n");

    int iteracoes = gauss_seidel(sist2.matriz, sist2.B, sist2.x, sist2.ordem);
    eliminacaoGauss(sist1.matriz, sist1.B, sist1.ordem);
    retrossubs(sist1.matriz, sist1.B, sist1.x, sist1.ordem);

    for (int i = 0; i < sist1.ordem; i++) {
        real_t soma1 = 0.0;
        real_t soma2 = 0.0;
        for (int k = 0; k < sist1.ordem; k++) {
            soma1 += sist1.matriz[i][k] * sist1.x[k];
            soma2 += sist2.matriz[i][k] * sist2.x[k];
        }
        r1[i] = sist1.B[i] - soma1;
        r2[i] = sist2.B[i] - soma2;
    }
     

    printf("EG:\n");
    for (int i = 0; i < sist1.ordem; i++) {
        printf("%.12e ", sist1.x[i]);
    }

    printf("\n");
    for (int i = 0; i < sist1.ordem; i++) {
        printf("%.12e ", r1[i]);
    }

    printf("\n\nGS [ %d iterações ]:\n", iteracoes);
    for (int i = 0; i < sist2.ordem; i++) {
        printf("%.12e ", sist2.x[i]);
    }

    printf("\n");
    
    for (int i = 0; i < sist2.ordem; i++) {
        printf("%.12e ", r2[i]);
    }

    printf("\n");

    // Liberação de memória
    for (int i = 0; i < sist1.ordem; i++) {
        free(sist1.matriz[i]);
        free(sist2.matriz[i]);
    }
    free(sist1.matriz);
    free(sist2.matriz);
    free(sist1.B);
    free(sist1.x);
    free(sist2.B);
    free(sist2.x);

    return 0;
}
