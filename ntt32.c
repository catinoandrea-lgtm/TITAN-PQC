/*
 * ntt32.c â€“ NTT NEON ibrido per TITAN-KEM-512
 * Tratto dalla relazione tecnica (giÃ  corretto nella logica NTT)
 * Compile: clang -std=c99 -O3 -march=armv8-a+simd
 */
#include <stdint.h>
#include <arm_neon.h>

#define Q          3329
#define QINV       ((int16_t)(-3327))
#define BARR       20159
#define INV_N_MONT ((int16_t)512)

static inline int16_t barrett_scalar(int16_t a) {
    int16_t t = (int16_t)(((int32_t)BARR * a + (1 << 25)) >> 26);
    return (int16_t)(a - t * (int16_t)Q);
}

static const int16_t ZETAS[128] = {
    -1044, -758, -359,-1517, 1493, 1422,  287,  202, -171,  622, 1577,  182,  962,-1202,-1474, 1468,
      573,-1325,  264,  383, -829, 1458,-1602, -130, -681, 1017,  732,  608,-1542,  411, -205,-1571,
     1223,  652, -552, 1015,-1293, 1491, -282,-1544,  516,   -8, -320, -666,-1618,-1162,  126, 1469,
     -853,  -90, -271,  830,  107,-1421, -247, -951, -398,  961,-1508, -725,  448,-1065,  677,-1275,
    -1103,  430,  555,  843,-1251,  871, 1550,  105,  422,  587,  177, -235, -291, -460, 1574, 1653,
     -246,  778, 1159, -147, -777, 1483, -602, 1119,-1590,  644, -872,  349,  418,  329, -156,  -75,
      817, 1097,  603,  610, 1322,-1285,-1465,  384,-1215, -136, 1218,-1335, -874,  220,-1187,-1659,
    -1185,-1530,-1278,  794,-1510, -854, -870,  478, -108, -308,  996,  991,  958,-1460, 1522, 1628
};

static inline int16x8_t fqmul_neon(int16x8_t a, int16x8_t b) {
    const int16x8_t v_q    = vdupq_n_s16(Q);
    const int16x8_t v_qinv = vdupq_n_s16(QINV);
    int32x4_t pl = vmull_s16(vget_low_s16(a), vget_low_s16(b));
    int32x4_t ph = vmull_high_s16(a, b);
    int16x8_t p  = vcombine_s16(vmovn_s32(pl), vmovn_s32(ph));
    int16x8_t t  = vmulq_s16(p, v_qinv);
    pl = vmlsl_s16(pl, vget_low_s16(t), vget_low_s16(v_q));
    ph = vmlsl_high_s16(ph, t, v_q);
    return vcombine_s16(vshrn_n_s32(pl, 16), vshrn_n_s32(ph, 16));
}

static inline int16x8_t barrett_neon(int16x8_t a) {
    const int16x8_t v_barr = vdupq_n_s16(BARR);
    const int16x8_t v_q    = vdupq_n_s16(Q);
    int32x4_t al = vshrq_n_s32(vmull_s16(vget_low_s16(a), vget_low_s16(v_barr)), 26);
    int32x4_t ah = vshrq_n_s32(vmull_high_s16(a, v_barr), 26);
    return vsubq_s16(a, vmulq_s16(vcombine_s16(vmovn_s32(al), vmovn_s32(ah)), v_q));
}

void ntt32_forward(int16_t r[256]) {
    int k = 1, len, s, j;
    /* Strati NEON: len = 128, 64, 32, 16, 8 */
    for (len = 128; len >= 8; len >>= 1) {
        for (s = 0; s < 256; s += 2 * len) {
            int16x8_t vz = vdupq_n_s16(ZETAS[k++]);
            for (j = s; j < s + len; j += 8) {
                int16x8_t r0 = vld1q_s16(&r[j]), r1 = vld1q_s16(&r[j+len]);
                int16x8_t t  = fqmul_neon(vz, r1);
                vst1q_s16(&r[j],     vaddq_s16(r0, t));
                vst1q_s16(&r[j+len], vsubq_s16(r0, t));
            }
        }
        /* Barrett dopo ogni layer NEON per contenere i coefficienti */
        for (int i = 0; i < 256; i += 8)
            vst1q_s16(&r[i], barrett_neon(vld1q_s16(&r[i])));
    }
    /* Strati scalari: len = 4, 2 — con Barrett */
    for (; len >= 2; len >>= 1) {
        for (s = 0; s < 256; s += 2 * len) {
            int16_t z = ZETAS[k++];
            for (j = s; j < s + len; j++) {
                int32_t p  = (int32_t)z * r[j+len];
                int16_t t  = (int16_t)((p - (int32_t)((int16_t)((int16_t)p * QINV)) * Q) >> 16);
                r[j+len] = barrett_scalar((int16_t)(r[j] - t));
                r[j]     = barrett_scalar((int16_t)(r[j] + t));
            }
        }
    }
}

void ntt32_inverse(int16_t r[256]) {
    int k = 127, len, s, j;
    /* Strati scalari: len = 2, 4 — con Barrett */
    for (len = 2; len <= 4; len <<= 1) {
        for (s = 0; s < 256; s += 2 * len) {
            int16_t z = ZETAS[k--];
            for (j = s; j < s + len; j++) {
                int16_t t  = r[j];
                r[j]       = barrett_scalar((int16_t)(t + r[j+len]));
                int32_t p  = (int32_t)z * (r[j+len] - t);
                r[j+len]   = (int16_t)((p - (int32_t)((int16_t)((int16_t)p * QINV)) * Q) >> 16);
            }
        }
    }
    /* Strati NEON: len = 8..128 — Barrett sulla somma */
    for (len = 8; len <= 128; len <<= 1) {
        for (s = 0; s < 256; s += 2 * len) {
            int16x8_t vz = vdupq_n_s16(ZETAS[k--]);
            for (j = s; j < s + len; j += 8) {
                int16x8_t r0 = vld1q_s16(&r[j]), r1 = vld1q_s16(&r[j+len]);
                vst1q_s16(&r[j],     barrett_neon(vaddq_s16(r0, r1)));
                vst1q_s16(&r[j+len], fqmul_neon(vz, vsubq_s16(r1, r0)));
            }
        }
    }
    /* Normalizzazione finale: * INV_N_MONT */
    int16x8_t vinv = vdupq_n_s16(INV_N_MONT);
    for (j = 0; j < 256; j += 8)
        vst1q_s16(&r[j], fqmul_neon(vld1q_s16(&r[j]), vinv));
}

void ntt32_pointwise_mul(int16_t *c, const int16_t *a, const int16_t *b) {
    for (int i = 0; i < 256; i += 8)
        vst1q_s16(&c[i], fqmul_neon(vld1q_s16(&a[i]), vld1q_s16(&b[i])));
}