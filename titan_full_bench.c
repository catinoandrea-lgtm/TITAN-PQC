#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// Includiamo l'header del tuo KEM originale
#include "kem.h"

static double now_us(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e6 + (double)ts.tv_nsec * 1e-3;
}

int main(void) {
    KemPublicKey pk; KemSecretKey sk; KemCiphertext ct; KemSharedSecret ss1, ss2;
    uint8_t seed[32] = {0x42};
    const int N = 1000;
    int16_t poly[256] = {1};

    printf("==========================================================\n");
    printf("  TITAN-KEM-512 | OFFICIAL ARM BENCHMARK (Cxxdroid / M1)\n");
    printf("==========================================================\n");

    // --- TEST CORRETTEZZA KEM ---
    kem_keygen(&pk, &sk, seed);
    kem_encaps(&ct, &ss1, &pk, seed);
    kem_decaps(&ss2, &ct, &sk);

    if(memcmp(ss1.data, ss2.data, 32) == 0) {
        printf("[0] Correctness Test: [PASS] (Shared secrets match)\n\n");
    } else {
        printf("[0] Correctness Test: [FAIL] (Mismatch!)\n\n");
        return 1;
    }

    // --- KEM PROTOCOL BENCHMARK ---
    printf("[A] SUSTAINED BENCHMARK (N=%d)\n", N);
    
    double t0 = now_us();
    for(int i=0; i<N; i++) kem_keygen(&pk, &sk, seed);
    printf("    Keygen (Avg) : %7.2f us\n", (now_us()-t0)/N);

    t0 = now_us();
    for(int i=0; i<N; i++) kem_encaps(&ct, &ss1, &pk, seed);
    printf("    Encaps (Avg) : %7.2f us\n", (now_us()-t0)/N);

    t0 = now_us();
    for(int i=0; i<N; i++) kem_decaps(&ss2, &ct, &sk);
    printf("    Decaps (Avg) : %7.2f us\n", (now_us()-t0)/N);

    printf("\n[B] ATOMIC OPERATIONS (Math Engine)\n");
    
    // --- FORWARD NTT BENCHMARK ---
    t0 = now_us();
    for(int i=0; i<N*10; i++) ntt_forward(poly);
    printf("    NTT Forward  : %7.2f us\n", (now_us()-t0)/(N*10));

    // --- INVERSE NTT BENCHMARK ---
    t0 = now_us();
    for(int i=0; i<N*10; i++) ntt_inverse(poly);
    printf("    NTT Inverse  : %7.2f us\n", (now_us()-t0)/(N*10));

    printf("==========================================================\n");
    return 0;
}
