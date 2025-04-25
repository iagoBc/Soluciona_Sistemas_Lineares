#ifndef __GAUSS_H__
#define __GAUSS_H__
#include "utils.h"

int_t encontraMAX(real_t **A, int_t i, int_t n);
void trocaLinha(real_t **A, real_t *b, int_t i, int_t n, int_t iPivo);
void eliminacaoGauss(real_t **A, real_t *b, int_t n);

#endif