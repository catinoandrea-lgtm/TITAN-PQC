#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "titan_avx2.h"

#define Q 3329

static double now_us(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e6 + (double)ts.tv_nsec * 1e-3;
}

static int mod_q(int x) {
    x %= Q;
    if (x < 0) x += Q;
    return x;
}

static int roundtrip_test(void) {
    poly_avx a, orig;
    for (int i = 0; i < 256; ++i) {
        a[i] = (int16_t)((i * 17 + 123) % Q);
        orig[i] = a[i];
    }
    titan_ntt_avx2(a);
    titan_ntt_inverse_avx2(a);
    for (int i = 0; i < 256; ++i) {
        if (mod_q(a[i]) != mod_q(orig[i])) return 0;
    }
    return 1;
}

int main(void) {
    poly_avx a, b, c;
    uint8_t seed[32] = {0x01, 0x02, 0x03, 0x04};
    uint8_t outA[4][672];
    const int ITERS = 200000;

    for (int i = 0; i < 256; ++i) {
        a[i] = (int16_t)(i % Q);
        b[i] = (int16_t)((i * 3) % Q);
        c[i] = 0;
    }

    printf("==========================================================\n");
    printf(" TITAN-KEM-512 | OFFICIAL AVX2 BENCHMARK v2\n");
    printf("==========================================================\n");
    printf("Forward-Inverse Roundtrip : %s\n", roundtrip_test() ? "PASS" : "FAIL");

    double t0 = now_us();
    for(int i = 0; i < ITERS; i++) titan_ntt_avx2(a);
    printf("NTT Forward              : %.3f us\n", (now_us() - t0) / ITERS);

    t0 = now_us();
    for(int i = 0; i < ITERS; i++) titan_ntt_inverse_avx2(a);
    printf("NTT Inverse              : %.3f us\n", (now_us() - t0) / ITERS);

    t0 = now_us();
    for(int i = 0; i < ITERS; i++) titan_pointwise_avx2(c, a, b);
    printf("Pointwise Mul            : %.3f us\n", (now_us() - t0) / ITERS);

    t0 = now_us();
    for(int i = 0; i < ITERS / 10; i++) titan_shake128x4_matrixA(outA, seed);
    printf("Keccakx4 Matrix A        : %.3f us\n", (now_us() - t0) / (ITERS / 10));

    printf("==========================================================\n");
    return 0;
}
