/*
* titan_avx2.h - Interfaccia Vettoriale AVX2 per TITAN-KEM
*/
#ifndef TITAN_AVX2_H
#define TITAN_AVX2_H

#include <stdint.h>

#define TITAN_N 256
typedef int16_t poly_avx[TITAN_N];

#ifdef __cplusplus
extern "C" {
#endif

void titan_pointwise_avx2(poly_avx c, const poly_avx a, const poly_avx b);
void titan_poly_add_avx2(poly_avx c, const poly_avx a, const poly_avx b);
void titan_poly_sub_avx2(poly_avx c, const poly_avx a, const poly_avx b);
void titan_ntt_avx2(poly_avx r);
void titan_ntt_inverse_avx2(poly_avx r);
void titan_cbd2_avx2(const uint8_t *noise, poly_avx r);
void titan_shake128x4_matrixA(uint8_t out[4][672], const uint8_t seed[32]);

#ifdef __cplusplus
}
#endif

#endif
