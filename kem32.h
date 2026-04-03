
#pragma once
#include <stdint.h>
#include <stddef.h>
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

#define KEM32_K         2
#define KEM32_N         256
#define KEM32_Q         3329
#define KEM32_PK_BYTES  800
#define KEM32_SK_BYTES  1632
#define KEM32_CT_BYTES  768
#define KEM32_SS_BYTES  32
#define KEM32_SYM_BYTES 32

// Public key
typedef struct {
    uint8_t data[KEM32_PK_BYTES];
} Kem32PublicKey;

// Secret key
typedef struct {
    uint8_t data[KEM32_SK_BYTES];
} Kem32SecretKey;

// Ciphertext
typedef struct {
    uint8_t data[KEM32_CT_BYTES];
} Kem32Ciphertext;

// Shared secret
typedef struct {
    uint8_t data[KEM32_SS_BYTES];
} Kem32SharedSecret;

// Function prototypes
void kem32_keygen (Kem32PublicKey *pk, Kem32SecretKey *sk, const uint8_t seed[KEM32_SYM_BYTES]);
void kem32_encaps (Kem32Ciphertext *ct, Kem32SharedSecret *ss,
                   const Kem32PublicKey *pk, const uint8_t rnd[KEM32_SYM_BYTES]);
int  kem32_decaps (Kem32SharedSecret *ss,
                   const Kem32Ciphertext *ct, const Kem32SecretKey *sk);

// NTT interface matching your test
void ntt32_forward      (int16_t r[KEM32_N]);
void ntt32_inverse      (int16_t r[KEM32_N]);
void ntt32_pointwise_mul(int16_t *c, const int16_t *a, const int16_t *b);

#ifdef __cplusplus
}
#endif