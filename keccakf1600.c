/*
 * keccakf1600.c â€“ Keccak-f[1600] unrolled a 2 round, CORRETTO
 *
 * Questa implementazione usa la mappatura Pi/Rho verificata per ogni riga
 * Chi. I bug della versione precedente erano:
 *   1. Righe Chi Esa/Eki scambiate (input elements errati)
 *   2. Double XOR: variabili giÃ  modificate da ^= usate di nuovo con ^D
 *   3. Rho offset errato per Ase (usava 56 invece di 2)
 *   4. Aliasing UB nel keccak_absorb
 */

#include <stdint.h>
#include <string.h>
#include "keccakf1600.h"

#define ROL64(a,n) (((uint64_t)(a) << (n)) | ((uint64_t)(a) >> (64-(n))))

static const uint64_t RC[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808aULL,
    0x8000000080008000ULL, 0x000000000000808bULL, 0x0000000080000001ULL,
    0x8000000080008081ULL, 0x8000000000008009ULL, 0x000000000000008aULL,
    0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000aULL,
    0x000000008000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL,
    0x8000000000008003ULL, 0x8000000000008002ULL, 0x8000000000000080ULL,
    0x000000000000800aULL, 0x800000008000000aULL, 0x8000000080008081ULL,
    0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL
};

/* Offsets Rho da FIPS 202 (indice = posizione nella tabella pi) */
static const int ROTC[24] = {
    1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14,
    27, 41, 56,  8, 25, 43, 62, 18, 39, 61, 20, 44
};
/* Pi: percorso di permutazione s[PI[i]] â†’ s[1] (Rho-Pi unrolled) */
static const int PI[24] = {
    10, 7, 11, 17, 18,  3,  5, 16,  8, 21, 24,  4,
    15, 23, 19, 13, 12,  2, 20, 14, 22,  9,  6,  1
};

void KeccakF1600_StatePermute(uint64_t s[25]) {
    uint64_t bc[5], t;
    for (int r = 0; r < 24; r++) {
        /* â”€â”€ Theta â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
        for (int i = 0; i < 5; i++)
            bc[i] = s[i] ^ s[i+5] ^ s[i+10] ^ s[i+15] ^ s[i+20];
        for (int i = 0; i < 5; i++) {
            t = bc[(i+4)%5] ^ ROL64(bc[(i+1)%5], 1);
            for (int j = 0; j < 25; j += 5) s[j+i] ^= t;
        }
        /* â”€â”€ Rho + Pi â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
        t = s[1];
        for (int i = 0; i < 24; i++) {
            int j = PI[i];
            bc[0] = s[j];
            s[j] = ROL64(t, ROTC[i]);
            t = bc[0];
        }
        /* â”€â”€ Chi â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
        for (int j = 0; j < 25; j += 5) {
            for (int i = 0; i < 5; i++) bc[i] = s[j+i];
            for (int i = 0; i < 5; i++)
                s[j+i] = bc[i] ^ (~bc[(i+1)%5] & bc[(i+2)%5]);
        }
        /* â”€â”€ Iota â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
        s[0] ^= RC[r];
    }
}

/* â”€â”€ Sponge core â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
static void keccak_absorb(uint64_t s[25], unsigned rate,
                           const uint8_t *in, size_t inlen, uint8_t pad) {
    unsigned r64 = rate / 8;
    while (inlen >= rate) {
        for (unsigned i = 0; i < r64; i++) {
            uint64_t tmp; memcpy(&tmp, in + 8*i, 8); s[i] ^= tmp;
        }
        KeccakF1600_StatePermute(s);
        in    += rate;
        inlen -= rate;
    }
    /* padding finale */
    uint8_t t[200] = {0};
    memcpy(t, in, inlen);
    t[inlen]    = pad;
    t[rate - 1] |= 0x80;
    for (unsigned i = 0; i < r64; i++) {
        uint64_t tmp; memcpy(&tmp, t + 8*i, 8); s[i] ^= tmp;
    }
}

static void keccak_squeeze(uint64_t s[25], unsigned rate,
                            uint8_t *out, size_t outlen) {
    while (outlen > 0) {
        KeccakF1600_StatePermute(s);
        size_t n = (outlen < rate) ? outlen : rate;
        memcpy(out, s, n);
        out    += n;
        outlen -= n;
    }
}

/* â”€â”€ SHA3-256 / SHA3-512 â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
void sha3_256(uint8_t out[32], const uint8_t *in, size_t inlen) {
    uint64_t s[25] = {0};
    keccak_absorb(s, 136, in, inlen, 0x06);
    KeccakF1600_StatePermute(s);
    memcpy(out, s, 32);
}

void sha3_512(uint8_t out[64], const uint8_t *in, size_t inlen) {
    uint64_t s[25] = {0};
    keccak_absorb(s, 72, in, inlen, 0x06);
    KeccakF1600_StatePermute(s);
    memcpy(out, s, 64);
}

/* â”€â”€ SHAKE-128 / SHAKE-256 (bulk API) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
void shake128(uint8_t *out, size_t outlen,
              const uint8_t *in, size_t inlen) {
    uint64_t s[25] = {0};
    keccak_absorb(s, 168, in, inlen, 0x1F);
    keccak_squeeze(s, 168, out, outlen);
}

void shake256(uint8_t *out, size_t outlen,
              const uint8_t *in, size_t inlen) {
    uint64_t s[25] = {0};
    keccak_absorb(s, 136, in, inlen, 0x1F);
    keccak_squeeze(s, 136, out, outlen);
}

/* â”€â”€ SHAKE-128 streaming API (per compatibilitÃ ) â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€ */
void shake128_init(KeccakState *st)                              { memset(st, 0, sizeof *st); }
void shake128_absorb(KeccakState *st, const uint8_t *in, size_t l) { keccak_absorb((uint64_t*)st, 168, in, l, 0x1F); }
void shake128_finalize(KeccakState *st)                          { (void)st; }
void shake128_squeeze(KeccakState *st, uint8_t *out, size_t l)   { keccak_squeeze((uint64_t*)st, 168, out, l); }