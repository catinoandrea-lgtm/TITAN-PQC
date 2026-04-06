/*
 * kem_bridge.c – Mappa kem_* -> kem32_* (scalare C99)
 * NON compilare insieme a kem.c
 * Compile: clang -std=c99 -O2 -lm
 *   test_titan32.c kem_bridge.c kem32.c keccakf1600.c sha256.c -o test_scalar
 */
#include "kem.h"
#include "kem32.h"

void kem_keygen(KemPublicKey *pk, KemSecretKey *sk,
                const uint8_t seed[64]) {
    kem32_keygen((Kem32PublicKey *)pk,(Kem32SecretKey *)sk,seed);
}
void kem_encaps(KemCiphertext *ct, KemSharedSecret *ss,
                const KemPublicKey *pk, const uint8_t rnd[32]) {
    kem32_encaps((Kem32Ciphertext *)ct,(Kem32SharedSecret *)ss,
                 (const Kem32PublicKey *)pk,rnd);
}
int kem_decaps(KemSharedSecret *ss, const KemCiphertext *ct,
               const KemSecretKey *sk) {
    return kem32_decaps((Kem32SharedSecret *)ss,
                        (const Kem32Ciphertext *)ct,
                        (const Kem32SecretKey *)sk);
}
/* alias: alcune versioni del test chiamano ntt_forward */
void ntt_forward(int16_t r[256]) { ntt32_forward(r); }
