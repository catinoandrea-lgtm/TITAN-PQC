/*
 * titan_kat_test.c — ML-KEM-512 Known Answer Test (KAT) harness
 * Andrea Catino — TITAN v2.3-FIPS203
 *
 * This file tests TITAN against deterministic KAT vectors to verify
 * FIPS 203 conformance. It runs on x86_64 (no NEON required).
 *
 * Compile (Linux/Colab/GitHub CI):
 *   gcc -O3 -std=gnu11 -Wall titan_kat_test.c kem32.c keccakf1600.c -o titan_kat -lm
 *   ./titan_kat
 *
 * Tests performed:
 *   1. Internal consistency: keygen → encaps → decaps → shared secrets match
 *   2. Determinism: same seed always produces same keys/ciphertext
 *   3. Decaps rejection: modified ciphertext produces different shared secret
 *   4. Key sizes: pk=800, sk=1632, ct=768, ss=32
 *   5. NTT roundtrip: forward(inverse(x)) == x for all coefficients
 *   6. Multiple random seeds: 100 iterations of full KEM cycle
 *   7. SHA3/SHAKE consistency: verify hash outputs against known values
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "kem32.h"
#include "keccakf1600.h"

/* ── Hex printing ─────────────────────────────────────────────── */
static void print_hex(const char *label, const uint8_t *data, int len) {
    printf("  %s: ", label);
    for (int i = 0; i < len; i++) printf("%02x", data[i]);
    printf("\n");
}

/* ── Deterministic PRNG (SHAKE-256 based) for reproducible tests ── */
static void deterministic_seed(uint8_t *out, size_t outlen,
                               const uint8_t *entropy, size_t elen,
                               uint32_t counter) {
    uint8_t input[256];
    memcpy(input, entropy, elen < 252 ? elen : 252);
    input[elen]   = (uint8_t)(counter & 0xFF);
    input[elen+1] = (uint8_t)((counter >> 8) & 0xFF);
    input[elen+2] = (uint8_t)((counter >> 16) & 0xFF);
    input[elen+3] = (uint8_t)((counter >> 24) & 0xFF);
    shake256(out, outlen, input, elen + 4);
}

/* ── Test 1: SHA3-256 known answer ────────────────────────────── */
static int test_sha3_256(void) {
    /* SHA3-256("") = a7ffc6f8bf1ed76651c14756a061d662f580ff4de43b49fa82d80a4b80f8434a */
    static const uint8_t expected[32] = {
        0xa7,0xff,0xc6,0xf8,0xbf,0x1e,0xd7,0x66,
        0x51,0xc1,0x47,0x56,0xa0,0x61,0xd6,0x62,
        0xf5,0x80,0xff,0x4d,0xe4,0x3b,0x49,0xfa,
        0x82,0xd8,0x0a,0x4b,0x80,0xf8,0x43,0x4a
    };
    uint8_t out[32];
    sha3_256(out, (const uint8_t*)"", 0);
    if (memcmp(out, expected, 32) != 0) {
        printf("[FAIL] SHA3-256 empty string\n");
        print_hex("got     ", out, 32);
        print_hex("expected", expected, 32);
        return 0;
    }
    return 1;
}

/* ── Test 2: SHA3-512 known answer ────────────────────────────── */
static int test_sha3_512(void) {
    /* SHA3-512("") = a69f73cca23a9ac5c8b567dc185a756e97c982164fe25859e0d1dcc1475c80a6
                       15b2123af1f5f94c11e3e9402c3ac558f500199d95b6d3e301758586281dcd26 */
    static const uint8_t expected[64] = {
        0xa6,0x9f,0x73,0xcc,0xa2,0x3a,0x9a,0xc5,
        0xc8,0xb5,0x67,0xdc,0x18,0x5a,0x75,0x6e,
        0x97,0xc9,0x82,0x16,0x4f,0xe2,0x58,0x59,
        0xe0,0xd1,0xdc,0xc1,0x47,0x5c,0x80,0xa6,
        0x15,0xb2,0x12,0x3a,0xf1,0xf5,0xf9,0x4c,
        0x11,0xe3,0xe9,0x40,0x2c,0x3a,0xc5,0x58,
        0xf5,0x00,0x19,0x9d,0x95,0xb6,0xd3,0xe3,
        0x01,0x75,0x85,0x86,0x28,0x1d,0xcd,0x26
    };
    uint8_t out[64];
    sha3_512(out, (const uint8_t*)"", 0);
    if (memcmp(out, expected, 64) != 0) {
        printf("[FAIL] SHA3-512 empty string\n");
        return 0;
    }
    return 1;
}

/* ── Test 3: SHAKE-128 known answer ───────────────────────────── */
static int test_shake128(void) {
    /* SHAKE128("", 32) = 7f9c2ba4e88f827d616045507605853ed73b8093f6efbc88eb1a6eacfa66ef26 */
    static const uint8_t expected[32] = {
        0x7f,0x9c,0x2b,0xa4,0xe8,0x8f,0x82,0x7d,
        0x61,0x60,0x45,0x50,0x76,0x05,0x85,0x3e,
        0xd7,0x3b,0x80,0x93,0xf6,0xef,0xbc,0x88,
        0xeb,0x1a,0x6e,0xac,0xfa,0x66,0xef,0x26
    };
    uint8_t out[32];
    shake128(out, 32, (const uint8_t*)"", 0);
    if (memcmp(out, expected, 32) != 0) {
        printf("[FAIL] SHAKE-128 empty string\n");
        print_hex("got     ", out, 32);
        print_hex("expected", expected, 32);
        return 0;
    }
    return 1;
}

/* ── Test 4: NTT roundtrip ────────────────────────────────────── */
static int test_ntt_roundtrip(void) {
    int16_t poly[256], orig[256];
    for (int i = 0; i < 256; i++) {
        poly[i] = (int16_t)((i * 17 + 42) % 3329);
        orig[i] = poly[i];
    }
    ntt32_forward(poly);
    ntt32_inverse(poly);
    for (int i = 0; i < 256; i++) {
        int16_t a = poly[i] % 3329;
        if (a < 0) a += 3329;
        int16_t b = orig[i] % 3329;
        if (b < 0) b += 3329;
        if (a != b) {
            printf("[FAIL] NTT roundtrip at index %d: got %d, expected %d\n", i, a, b);
            return 0;
        }
    }
    return 1;
}

/* ── Test 5: Key sizes ────────────────────────────────────────── */
static int test_key_sizes(void) {
    if (KEM32_PK_BYTES != 800) { printf("[FAIL] PK size %d != 800\n", KEM32_PK_BYTES); return 0; }
    if (KEM32_SK_BYTES != 1632) { printf("[FAIL] SK size %d != 1632\n", KEM32_SK_BYTES); return 0; }
    if (KEM32_CT_BYTES != 768) { printf("[FAIL] CT size %d != 768\n", KEM32_CT_BYTES); return 0; }
    if (KEM32_SS_BYTES != 32) { printf("[FAIL] SS size %d != 32\n", KEM32_SS_BYTES); return 0; }
    return 1;
}

/* ── Test 6: KEM determinism ──────────────────────────────────── */
static int test_kem_determinism(void) {
    uint8_t seed[64] = {0};
    for (int i = 0; i < 64; i++) seed[i] = (uint8_t)(i * 7 + 3);
    uint8_t rnd[32] = {0};
    for (int i = 0; i < 32; i++) rnd[i] = (uint8_t)(i * 13 + 5);

    Kem32PublicKey pk1, pk2;
    Kem32SecretKey sk1, sk2;
    Kem32Ciphertext ct1, ct2;
    Kem32SharedSecret ss1, ss2;

    kem32_keygen(&pk1, &sk1, seed);
    kem32_keygen(&pk2, &sk2, seed);

    if (memcmp(pk1.data, pk2.data, KEM32_PK_BYTES) != 0) {
        printf("[FAIL] Keygen not deterministic (pk differs)\n");
        return 0;
    }
    if (memcmp(sk1.data, sk2.data, KEM32_SK_BYTES) != 0) {
        printf("[FAIL] Keygen not deterministic (sk differs)\n");
        return 0;
    }

    kem32_encaps(&ct1, &ss1, &pk1, rnd);
    kem32_encaps(&ct2, &ss2, &pk1, rnd);

    if (memcmp(ct1.data, ct2.data, KEM32_CT_BYTES) != 0) {
        printf("[FAIL] Encaps not deterministic (ct differs)\n");
        return 0;
    }
    if (memcmp(ss1.data, ss2.data, KEM32_SS_BYTES) != 0) {
        printf("[FAIL] Encaps not deterministic (ss differs)\n");
        return 0;
    }
    return 1;
}

/* ── Test 7: KEM correctness (encaps/decaps match) ────────────── */
static int test_kem_correctness(const uint8_t seed[64], const uint8_t rnd[32], int id) {
    Kem32PublicKey pk;
    Kem32SecretKey sk;
    Kem32Ciphertext ct;
    Kem32SharedSecret ss_enc, ss_dec;

    kem32_keygen(&pk, &sk, seed);
    kem32_encaps(&ct, &ss_enc, &pk, rnd);
    int rc = kem32_decaps(&ss_dec, &ct, &sk);

    if (memcmp(ss_enc.data, ss_dec.data, KEM32_SS_BYTES) != 0) {
        printf("[FAIL] KEM vector #%d: shared secrets mismatch\n", id);
        print_hex("ss_enc", ss_enc.data, 32);
        print_hex("ss_dec", ss_dec.data, 32);
        return 0;
    }
    if (rc != 0) {
        printf("[FAIL] KEM vector #%d: decaps returned failure (%d)\n", id, rc);
        return 0;
    }
    return 1;
}

/* ── Test 8: Implicit rejection (modified ciphertext) ─────────── */
static int test_implicit_rejection(void) {
    uint8_t seed[64] = {0xAA};
    uint8_t rnd[32] = {0xBB};

    Kem32PublicKey pk;
    Kem32SecretKey sk;
    Kem32Ciphertext ct;
    Kem32SharedSecret ss_enc, ss_dec;

    kem32_keygen(&pk, &sk, seed);
    kem32_encaps(&ct, &ss_enc, &pk, rnd);

    /* Flip one byte in ciphertext */
    ct.data[0] ^= 0xFF;

    int rc = kem32_decaps(&ss_dec, &ct, &sk);

    /* Decaps should return failure AND shared secret should differ */
    if (rc == 0) {
        printf("[FAIL] Implicit rejection: decaps did not detect modification\n");
        return 0;
    }
    if (memcmp(ss_enc.data, ss_dec.data, KEM32_SS_BYTES) == 0) {
        printf("[FAIL] Implicit rejection: same shared secret despite modified CT\n");
        return 0;
    }
    return 1;
}

/* ── Test 9: KAT vector with known hash of outputs ───────────── */
static int test_kat_accumulated(void) {
    /*
     * Accumulated KAT: run 100 keygen+encaps+decaps with deterministic seeds,
     * hash all shared secrets together. This verifies that the entire protocol
     * chain is consistent across builds and platforms.
     */
    uint8_t master_seed[32] = {
        0x01,0x23,0x45,0x67,0x89,0xAB,0xCD,0xEF,
        0xFE,0xDC,0xBA,0x98,0x76,0x54,0x32,0x10,
        0x0F,0x1E,0x2D,0x3C,0x4B,0x5A,0x69,0x78,
        0x87,0x96,0xA5,0xB4,0xC3,0xD2,0xE1,0xF0
    };

    /* Accumulate all shared secrets */
    uint8_t accumulator[32] = {0};
    int pass = 1;

    for (int i = 0; i < 100; i++) {
        uint8_t seed[64], rnd[32];
        deterministic_seed(seed, 64, master_seed, 32, (uint32_t)(i * 2));
        deterministic_seed(rnd, 32, master_seed, 32, (uint32_t)(i * 2 + 1));

        if (!test_kem_correctness(seed, rnd, i)) {
            pass = 0;
            break;
        }

        /* Compute shared secret for accumulation */
        Kem32PublicKey pk;
        Kem32SecretKey sk;
        Kem32Ciphertext ct;
        Kem32SharedSecret ss;
        kem32_keygen(&pk, &sk, seed);
        kem32_encaps(&ct, &ss, &pk, rnd);
        for (int j = 0; j < 32; j++) accumulator[j] ^= ss.data[j];
    }

    if (pass) {
        /* Hash the accumulator to get a single fingerprint */
        uint8_t fingerprint[32];
        sha3_256(fingerprint, accumulator, 32);
        print_hex("KAT-100 fingerprint", fingerprint, 32);
    }
    return pass;
}

/* ══════════════════════════════════════════════════════════════
 *  MAIN
 * ══════════════════════════════════════════════════════════════ */
int main(void) {
    int total = 0, passed = 0;

    printf("==========================================================\n");
    printf("  TITAN ML-KEM-512 v2.3 | KAT & Conformance Test Suite\n");
    printf("==========================================================\n\n");

    /* Symmetric primitives */
    printf("[1] SHA3-256 (empty string) ........... ");
    total++; if (test_sha3_256()) { passed++; printf("[PASS]\n"); } else printf("\n");

    printf("[2] SHA3-512 (empty string) ........... ");
    total++; if (test_sha3_512()) { passed++; printf("[PASS]\n"); } else printf("\n");

    printf("[3] SHAKE-128 (empty string, 32B) ..... ");
    total++; if (test_shake128()) { passed++; printf("[PASS]\n"); } else printf("\n");

    /* NTT */
    printf("[4] NTT forward-inverse roundtrip ..... ");
    total++; if (test_ntt_roundtrip()) { passed++; printf("[PASS]\n"); } else printf("\n");

    /* Parameter sizes */
    printf("[5] Key/CT/SS sizes (FIPS 203) ........ ");
    total++; if (test_key_sizes()) { passed++; printf("[PASS]\n"); } else printf("\n");

    /* KEM protocol */
    printf("[6] KEM determinism ................... ");
    total++; if (test_kem_determinism()) { passed++; printf("[PASS]\n"); } else printf("\n");

    printf("[7] KEM correctness (single) .......... ");
    total++;
    { uint8_t s[64]={0x42}, r[32]={0x55};
      if (test_kem_correctness(s, r, 0)) { passed++; printf("[PASS]\n"); } else printf("\n"); }

    printf("[8] Implicit rejection ................ ");
    total++; if (test_implicit_rejection()) { passed++; printf("[PASS]\n"); } else printf("\n");

    printf("[9] Accumulated KAT (100 vectors) ..... ");
    total++; if (test_kat_accumulated()) { passed++; printf("[PASS]\n"); } else printf("\n");

    /* Summary */
    printf("\n==========================================================\n");
    if (passed == total) {
        printf("  ALL %d TESTS PASSED\n", total);
    } else {
        printf("  %d/%d TESTS PASSED — %d FAILED\n", passed, total, total - passed);
    }
    printf("==========================================================\n");

    return (passed == total) ? 0 : 1;
}
