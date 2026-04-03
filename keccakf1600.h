#pragma once
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint64_t st[25];
    uint32_t rate;
    uint32_t pos;
    uint8_t  pad;
    uint8_t  squeezing;
} KeccakState;

void keccakf1600      (uint64_t st[25]);
void shake128         (uint8_t *out, size_t olen, const uint8_t *in, size_t ilen);
void shake256         (uint8_t *out, size_t olen, const uint8_t *in, size_t ilen);
void sha3_256         (uint8_t out[32], const uint8_t *in, size_t len);
void sha3_512         (uint8_t out[64], const uint8_t *in, size_t len);
void shake128_init    (KeccakState *s);
void shake256_init    (KeccakState *s);
void shake128_absorb  (KeccakState *s, const uint8_t *in, size_t len);
void shake256_absorb  (KeccakState *s, const uint8_t *in, size_t len);
void shake128_finalize(KeccakState *s);
void shake256_finalize(KeccakState *s);
void shake128_squeeze (KeccakState *s, uint8_t *out, size_t len);
void shake256_squeeze (KeccakState *s, uint8_t *out, size_t len);

#ifdef __cplusplus
}
#endif
