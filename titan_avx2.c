#include "titan_avx2.h"
#include <immintrin.h>
#include <string.h>

#define Q 3329
#define QINV ((int16_t)(-3327))
#define BARR 20159
#define INV_N_MONT ((int16_t)512)

extern const int16_t ZETAS[128];

static inline int16_t fqmul_scalar(int16_t a, int16_t b) {
    int32_t p = (int32_t)a * b;
    int16_t t = (int16_t)((int16_t)p * QINV);
    return (int16_t)((p - (int32_t)t * Q) >> 16);
}

static inline int16_t barrett_scalar(int16_t a) {
    int16_t t = (int16_t)(((int32_t)BARR * a + (1 << 25)) >> 26);
    return (int16_t)(a - t * (int16_t)Q);
}

static inline __m256i fqmul_avx2(__m256i a, __m256i b) {
    __m256i lo = _mm256_mullo_epi16(a, b);
    __m256i hi = _mm256_mulhi_epi16(a, b);
    __m256i t = _mm256_mullo_epi16(lo, _mm256_set1_epi16(QINV));
    __m256i sub = _mm256_mulhi_epi16(t, _mm256_set1_epi16(Q));
    return _mm256_sub_epi16(hi, sub);
}

void titan_pointwise_avx2(poly_avx c, const poly_avx a, const poly_avx b) {
    for (int i = 0; i < TITAN_N; i += 16) {
        __m256i va = _mm256_loadu_si256((const __m256i*)(a + i));
        __m256i vb = _mm256_loadu_si256((const __m256i*)(b + i));
        _mm256_storeu_si256((__m256i*)(c + i), fqmul_avx2(va, vb));
    }
}

void titan_poly_add_avx2(poly_avx c, const poly_avx a, const poly_avx b) {
    for (int i = 0; i < TITAN_N; i += 16) {
        __m256i va = _mm256_loadu_si256((const __m256i*)(a + i));
        __m256i vb = _mm256_loadu_si256((const __m256i*)(b + i));
        _mm256_storeu_si256((__m256i*)(c + i), _mm256_add_epi16(va, vb));
    }
}

void titan_poly_sub_avx2(poly_avx c, const poly_avx a, const poly_avx b) {
    for (int i = 0; i < TITAN_N; i += 16) {
        __m256i va = _mm256_loadu_si256((const __m256i*)(a + i));
        __m256i vb = _mm256_loadu_si256((const __m256i*)(b + i));
        _mm256_storeu_si256((__m256i*)(c + i), _mm256_sub_epi16(va, vb));
    }
}

void titan_cbd2_avx2(const uint8_t *noise, poly_avx r) {
    for (int i = 0; i < 256; i += 16) {
        __m128i n8 = _mm_loadl_epi64((const __m128i*)(noise + i/2));
        __m128i n16 = _mm_cvtepu8_epi16(n8);
        __m128i a = _mm_and_si128(n16, _mm_set1_epi16(0x55));
        __m128i b = _mm_and_si128(_mm_srli_epi16(n16, 1), _mm_set1_epi16(0x55));
        __m128i sum = _mm_add_epi16(a, b);
        __m128i c0 = _mm_sub_epi16(_mm_and_si128(sum, _mm_set1_epi16(3)), _mm_and_si128(_mm_srli_epi16(sum, 2), _mm_set1_epi16(3)));
        __m128i c1 = _mm_sub_epi16(_mm_and_si128(_mm_srli_epi16(sum, 4), _mm_set1_epi16(3)), _mm_and_si128(_mm_srli_epi16(sum, 6), _mm_set1_epi16(3)));
        _mm_storeu_si128((__m128i*)(r + i), _mm_unpacklo_epi16(c0, c1));
        _mm_storeu_si128((__m128i*)(r + i + 8), _mm_unpackhi_epi16(c0, c1));
    }
}

void titan_ntt_avx2(poly_avx r) {
    int k = 1;
    for (int len = 128; len >= 16; len >>= 1) {
        for (int start = 0; start < 256; start += 2 * len) {
            __m256i zv = _mm256_set1_epi16(ZETAS[k++]);
            for (int j = start; j < start + len; j += 16) {
                __m256i pj = _mm256_loadu_si256((__m256i*)&r[j]);
                __m256i pjl = _mm256_loadu_si256((__m256i*)&r[j + len]);
                __m256i t = fqmul_avx2(zv, pjl);
                _mm256_storeu_si256((__m256i*)&r[j + len], _mm256_sub_epi16(pj, t));
                _mm256_storeu_si256((__m256i*)&r[j], _mm256_add_epi16(pj, t));
            }
        }
    }
    for (int len = 8; len >= 2; len >>= 1) {
        for (int start = 0; start < 256; start += 2 * len) {
            int16_t z = ZETAS[k++];
            for (int j = start; j < start + len; j++) {
                int16_t t = fqmul_scalar(z, r[j + len]);
                r[j + len] = barrett_scalar((int16_t)(r[j] - t));
                r[j] = barrett_scalar((int16_t)(r[j] + t));
            }
        }
    }
}

void titan_ntt_inverse_avx2(poly_avx r) {
    int k = 127;
    for (int len = 2; len <= 8; len <<= 1) {
        for (int start = 0; start < 256; start += 2 * len) {
            int16_t z = ZETAS[k--];
            for (int j = start; j < start + len; j++) {
                int16_t t = r[j];
                r[j] = barrett_scalar((int16_t)(t + r[j + len]));
                r[j + len] = fqmul_scalar(z, (int16_t)(r[j + len] - t));
            }
        }
    }
    for (int len = 16; len <= 128; len <<= 1) {
        for (int start = 0; start < 256; start += 2 * len) {
            int16_t z = ZETAS[k--];
            __m256i vz = _mm256_set1_epi16(z);
            for (int j = start; j < start + len; j += 16) {
                __m256i r0 = _mm256_loadu_si256((__m256i*)&r[j]);
                __m256i r1 = _mm256_loadu_si256((__m256i*)&r[j + len]);
                __m256i sum = _mm256_add_epi16(r0, r1);
                __m256i diff = _mm256_sub_epi16(r1, r0);
                __m256i prod = fqmul_avx2(vz, diff);
                _mm256_storeu_si256((__m256i*)&r[j], sum);
                _mm256_storeu_si256((__m256i*)&r[j + len], prod);
            }
        }
    }
    __m256i vinv = _mm256_set1_epi16(INV_N_MONT);
    for (int j = 0; j < 256; j += 16) {
        __m256i v = _mm256_loadu_si256((__m256i*)&r[j]);
        _mm256_storeu_si256((__m256i*)&r[j], fqmul_avx2(v, vinv));
    }
    for (int j = 0; j < 256; ++j) r[j] = barrett_scalar(r[j]);
}

#define ROL(a,n) _mm256_xor_si256(_mm256_slli_epi64((a), (n)), _mm256_srli_epi64((a), 64-(n)))
static const uint64_t RC[24]={0x1,0x8082,0x800000000000808a,0x8000000080008000,0x808b,0x80000001,0x8000000080008081,0x8000000000008009,0x8a,0x88,0x80008009,0x8000000a,0x8000808b,0x800000000000008b,0x8000000000008089,0x8000000000008003,0x8000000000008002,0x8000000000000080,0x800a,0x800000008000000a,0x8000000080008081,0x8000000000008080,0x80000001,0x8000000080008008};
static const int PI[24]={10,7,11,17,18,3,5,16,8,21,24,4,15,23,19,13,12,2,20,14,22,9,6,1};
static const int RT[24]={1,3,6,10,15,21,28,36,45,55,2,14,27,41,56,8,25,43,62,18,39,61,20,44};

void titan_shake128x4_matrixA(uint8_t out[4][672], const uint8_t seed[32]) {
    __m256i s[25] = {0};
    for(int i=0; i<4; i++) {
        uint64_t sd; memcpy(&sd, seed + i*8, 8);
        s[i] = _mm256_set1_epi64x(sd);
    }
    s[4] = _mm256_set_epi64x((0x1FLL<<16)|(1<<8)|1, (0x1FLL<<16)|(0<<8)|1, (0x1FLL<<16)|(1<<8)|0, (0x1FLL<<16)|(0<<8)|0);
    s[20] = _mm256_set1_epi64x(0x80ULL << 56);
    for(int r=0; r<24; r++){
        __m256i bc[5], t, b0;
        for(int i=0; i<5; i++) bc[i] = _mm256_xor_si256(s[i], _mm256_xor_si256(s[i+5], _mm256_xor_si256(s[i+10], _mm256_xor_si256(s[i+15], s[i+20]))));
        for(int i=0; i<5; i++){
            t = _mm256_xor_si256(bc[(i+4)%5], ROL(bc[(i+1)%5], 1));
            for(int j=0; j<25; j+=5) s[j+i] = _mm256_xor_si256(s[j+i], t);
        }
        t = s[1];
        for(int i=0; i<24; i++){ int j = PI[i]; b0 = s[j]; s[j] = ROL(t, RT[i]); t = b0; }
        for(int j=0; j<25; j+=5){
            for(int i=0; i<5; i++) bc[i] = s[j+i];
            for(int i=0; i<5; i++) s[j+i] = _mm256_xor_si256(bc[i], _mm256_andnot_si256(bc[(i+1)%5], bc[(i+2)%5]));
        }
        s[0] = _mm256_xor_si256(s[0], _mm256_set1_epi64x(RC[r]));
    }
    for(int i=0; i<21; i++) {
        uint64_t l[4]; _mm256_storeu_si256((__m256i*)l, s[i]);
        for(int k=0; k<4; k++) memcpy(out[k] + i*8, &l[k], 8);
    }
}
