#pragma once
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct { uint32_t state[8]; uint64_t count; uint8_t buf[64]; } SHA256_CTX;

void sha256_init  (SHA256_CTX *ctx);
void sha256_update(SHA256_CTX *ctx, const uint8_t *in, size_t len);
void sha256_final (SHA256_CTX *ctx, uint8_t out[32]);
void sha256       (uint8_t out[32], const uint8_t *in, size_t len);

#ifdef __cplusplus
}
#endif
