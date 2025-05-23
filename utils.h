#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <stdint.h>
#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fenv.h>  // Biblioteca para controle do arredondamento

#pragma STDC FENV_ACCESS ON 

// Valor absoluto de um número. Alternativa ao uso da função 'fabs()'
#define ABS(num)  ((num) < 0.0 ? -(num) : (num))

#define MAX_IT 50
#define TOL 1e-4 
// real_t: tipo usado para representar valores em ponto flutuante
typedef double real_t;

// int_t: tipo usado para representar valores em inteiros
typedef int64_t int_t;

// string_t: tipo usado para representar ponteiros para char/strings
typedef char * string_t;

// rtime_t: tipo usado para representar valores de tempo em ponto flutuante
typedef double rtime_t;

// Número máximo de dígitos em um número
#define numDigits(n)  6  // ( (int) log10(n) + 1 )

// Macro para verificar se valor 'n' é potência de 2 ou não
#define isPot2(n) (n && !(n & (n - 1)))

typedef struct {
  real_t **matriz;
  real_t *B, *x;
  int_t ordem;
} SL;

// Funções
rtime_t timestamp(void);
string_t markerName(string_t baseName, int n);
void retrossubs(real_t **A, real_t *b, real_t *x, int_t n);

#endif // __UTILS_H__

