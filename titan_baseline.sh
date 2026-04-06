#!/bin/bash
# titan_baseline.sh — Compare TITAN vs pqcrystals/kyber reference
# Runs on x86_64 Linux (Ubuntu, Colab, GitHub CI)
#
# Usage: bash titan_baseline.sh
#
set -e

echo "=========================================================="
echo " TITAN vs pqcrystals/kyber — Baseline Comparison"
echo "=========================================================="
echo ""

# ── 1. Build TITAN (scalar) ──────────────────────────────────
echo "[1/4] Building TITAN (scalar C99)..."
gcc -O3 -std=gnu11 -Wall titan_kat_test.c kem32.c keccakf1600.c -o titan_kat -lm
echo "      TITAN KAT built OK"

# Build TITAN scalar benchmark
cat > titan_scalar_bench.c << 'EOFBENCH'
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "kem32.h"

static double now_us(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e6 + (double)ts.tv_nsec * 1e-3;
}

int main(void) {
    Kem32PublicKey pk; Kem32SecretKey sk;
    Kem32Ciphertext ct; Kem32SharedSecret ss1, ss2;
    uint8_t seed[64] = {0x42};
    uint8_t rnd[32] = {0x55};
    const int N = 10000;

    /* Warmup */
    for(int i=0; i<100; i++) kem32_keygen(&pk, &sk, seed);

    double t0 = now_us();
    for(int i=0; i<N; i++) kem32_keygen(&pk, &sk, seed);
    double kg = (now_us()-t0)/N;

    t0 = now_us();
    for(int i=0; i<N; i++) kem32_encaps(&ct, &ss1, &pk, rnd);
    double en = (now_us()-t0)/N;

    t0 = now_us();
    for(int i=0; i<N; i++) kem32_decaps(&ss2, &ct, &sk);
    double de = (now_us()-t0)/N;

    printf("TITAN_SCALAR keygen=%.2f encaps=%.2f decaps=%.2f total=%.2f\n",
           kg, en, de, kg+en+de);
    return 0;
}
EOFBENCH

gcc -O3 -std=gnu11 titan_scalar_bench.c kem32.c keccakf1600.c -o titan_scalar_bench -lm
echo "      TITAN scalar bench built OK"

# ── 2. Download & build pqcrystals reference ─────────────────
echo ""
echo "[2/4] Downloading pqcrystals/kyber reference..."
if [ ! -d "kyber-ref" ]; then
    git clone --depth 1 https://github.com/pq-crystals/kyber.git kyber-ref 2>/dev/null
fi
echo "      Downloaded OK"

echo "      Building reference (ref implementation)..."
cd kyber-ref/ref

# Build reference benchmark
cat > bench_compare.c << 'EOFREF'
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "api.h"
#include "kem.h"
#include "randombytes.h"

static double now_us(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e6 + (double)ts.tv_nsec * 1e-3;
}

int main(void) {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t ct[CRYPTO_CIPHERTEXTBYTES];
    uint8_t ss1[CRYPTO_BYTES], ss2[CRYPTO_BYTES];
    const int N = 10000;

    /* Warmup */
    for(int i=0; i<100; i++) crypto_kem_keypair(pk, sk);

    double t0 = now_us();
    for(int i=0; i<N; i++) crypto_kem_keypair(pk, sk);
    double kg = (now_us()-t0)/N;

    t0 = now_us();
    for(int i=0; i<N; i++) crypto_kem_enc(ct, ss1, pk);
    double en = (now_us()-t0)/N;

    t0 = now_us();
    for(int i=0; i<N; i++) crypto_kem_dec(ss2, ct, sk);
    double de = (now_us()-t0)/N;

    printf("REF_C keygen=%.2f encaps=%.2f decaps=%.2f total=%.2f\n",
           kg, en, de, kg+en+de);
    return 0;
}
EOFREF

# Compile reference (kyber512)
gcc -O3 -std=gnu11 -DKYBER_K=2 \
    bench_compare.c kem.c indcpa.c polyvec.c poly.c ntt.c cbd.c reduce.c \
    verify.c symmetric-shake.c fips202.c randombytes.c \
    -o ../../ref_bench -lm 2>/dev/null

cd ../..
echo "      Reference bench built OK"

# ── 3. Run both ──────────────────────────────────────────────
echo ""
echo "[3/4] Running benchmarks (N=10000)..."
echo ""

echo "  --- TITAN KAT ---"
./titan_kat
echo ""

echo "  --- Performance ---"
TITAN_RESULT=$(./titan_scalar_bench)
REF_RESULT=$(./ref_bench)

echo "  $TITAN_RESULT"
echo "  $REF_RESULT"

# ── 4. Parse and compare ─────────────────────────────────────
echo ""
echo "[4/4] Comparison"
echo "=========================================================="

T_KG=$(echo $TITAN_RESULT | grep -oP 'keygen=\K[0-9.]+')
T_EN=$(echo $TITAN_RESULT | grep -oP 'encaps=\K[0-9.]+')
T_DE=$(echo $TITAN_RESULT | grep -oP 'decaps=\K[0-9.]+')
T_TOT=$(echo $TITAN_RESULT | grep -oP 'total=\K[0-9.]+')

R_KG=$(echo $REF_RESULT | grep -oP 'keygen=\K[0-9.]+')
R_EN=$(echo $REF_RESULT | grep -oP 'encaps=\K[0-9.]+')
R_DE=$(echo $REF_RESULT | grep -oP 'decaps=\K[0-9.]+')
R_TOT=$(echo $REF_RESULT | grep -oP 'total=\K[0-9.]+')

printf "  %-20s %12s %12s %8s\n" "Operation" "TITAN (us)" "Ref C (us)" "Ratio"
printf "  %-20s %12s %12s %8s\n" "--------------------" "----------" "----------" "------"
printf "  %-20s %12s %12s %8.2fx\n" "Keygen" "$T_KG" "$R_KG" "$(echo "$R_KG / $T_KG" | bc -l)"
printf "  %-20s %12s %12s %8.2fx\n" "Encaps" "$T_EN" "$R_EN" "$(echo "$R_EN / $T_EN" | bc -l)"
printf "  %-20s %12s %12s %8.2fx\n" "Decaps" "$T_DE" "$R_DE" "$(echo "$R_DE / $T_DE" | bc -l)"
printf "  %-20s %12s %12s %8.2fx\n" "TOTAL" "$T_TOT" "$R_TOT" "$(echo "$R_TOT / $T_TOT" | bc -l)"
echo "=========================================================="
echo "  Ratio > 1.0 = TITAN faster, < 1.0 = Reference faster"
echo "=========================================================="
