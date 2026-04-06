#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <sys/utsname.h>
#include "kem.h"

static double now_us(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e6 + (double)ts.tv_nsec * 1e-3;
}

int main(void) {
    KemPublicKey pk; KemSecretKey sk; KemCiphertext ct;
    KemSharedSecret ss1, ss2;
    uint8_t seed[64] = {0x42};
    uint8_t rnd[32]  = {0x55};
    const int N = 1000;
    int16_t poly[256] = {1};
    struct utsname un; uname(&un);

    printf("==========================================================\n");
    printf("  TITAN ML-KEM-512 v2.3 | ARM NEON Benchmark\n");
    printf("  Platform: %s\n", un.machine);
    printf("==========================================================\n");

    kem_keygen(&pk, &sk, seed);
    kem_encaps(&ct, &ss1, &pk, rnd);
    kem_decaps(&ss2, &ct, &sk);

    if(memcmp(ss1.data, ss2.data, 32) == 0)
        printf("[OK] Correctness: PASS (shared secrets match)\n\n");
    else {
        printf("[!!] Correctness: FAIL\n\n");
        return 1;
    }

    /* Warmup */
    for(int i=0; i<100; i++) kem_keygen(&pk, &sk, seed);

    printf("[A] KEM PROTOCOL (N=%d)\n", N);
    double t0 = now_us();
    for(int i=0; i<N; i++) kem_keygen(&pk, &sk, seed);
    double kg = (now_us()-t0)/N;
    printf("    Keygen : %7.2f us\n", kg);

    t0 = now_us();
    for(int i=0; i<N; i++) kem_encaps(&ct, &ss1, &pk, rnd);
    double en = (now_us()-t0)/N;
    printf("    Encaps : %7.2f us\n", en);

    t0 = now_us();
    for(int i=0; i<N; i++) kem_decaps(&ss2, &ct, &sk);
    double de = (now_us()-t0)/N;
    printf("    Decaps : %7.2f us\n", de);

    printf("    --------------------------------\n");
    printf("    Total  : %7.2f us (full handshake)\n\n", kg+en+de);

    printf("[B] NTT ENGINE (N=%d)\n", N*10);
    t0 = now_us();
    for(int i=0; i<N*10; i++) ntt_forward(poly);
    printf("    NTT Fwd: %7.2f us\n", (now_us()-t0)/(N*10));

    t0 = now_us();
    for(int i=0; i<N*10; i++) ntt_inverse(poly);
    printf("    NTT Inv: %7.2f us\n", (now_us()-t0)/(N*10));

    printf("==========================================================\n");
    return 0;
}
