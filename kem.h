#pragma once
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define KEM_K         2
#define KEM_N         256
#define KEM_Q         3329
#define KEM_PK_BYTES  800
#define KEM_SK_BYTES  1632
#define KEM_CT_BYTES  768
#define KEM_SS_BYTES  32
#define KEM_SYM_BYTES 32
#define KEM_KEYGEN_BYTES 64   /* FIPS 203: d[32] || z[32] */

typedef struct { uint8_t data[KEM_PK_BYTES];  } KemPublicKey;
typedef struct { uint8_t data[KEM_SK_BYTES];  } KemSecretKey;
typedef struct { uint8_t data[KEM_CT_BYTES];  } KemCiphertext;
typedef struct { uint8_t data[KEM_SS_BYTES];  } KemSharedSecret;

void kem_keygen (KemPublicKey *pk, KemSecretKey *sk,
                 const uint8_t seed[KEM_KEYGEN_BYTES]);
void kem_encaps (KemCiphertext *ct, KemSharedSecret *ss,
                 const KemPublicKey *pk,
                 const uint8_t randomness[KEM_SYM_BYTES]);
int  kem_decaps (KemSharedSecret *ss, const KemCiphertext *ct,
                 const KemSecretKey *sk);

void ntt_forward       (int16_t poly[KEM_N]);
void ntt_inverse       (int16_t poly[KEM_N]);
void ntt_pointwise_mul (int16_t *c, const int16_t *a, const int16_t *b);

#ifdef __cplusplus
}
#endif
