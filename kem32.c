/*
 * kem32.c  –  TITAN-KEM-512 scalar (C99, no NEON)
 * Andrea Catino – Independent Researcher, Italy
 *
 * v2.3  – Streaming gen_matrix (incremental SHAKE-128 squeeze)
 *
 * Compile: clang -std=c99 -O2 -lm
 *   titan_full_bench.c kem_bridge.c kem32.c keccakf1600.c -o titan_scalar
 */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "kem32.h"
#include "keccakf1600.h"

/* ── costanti ────────────────────────────────────────────────────── */
#define Q           3329
#define QINV        ((int16_t)(-3327))
#define BARR        20159
#define R2          ((int16_t)1353)
#define INV_N_MONT  ((int16_t)512)
#define SHAKE128_RATE 168
#define XOF_BUFLEN   (4 * SHAKE128_RATE)   /* 672 byte */

/* ── aritmetica scalare ──────────────────────────────────────────── */
static inline int16_t fqmul(int16_t a, int16_t b) {
    int32_t p = (int32_t)a * b;
    int16_t t = (int16_t)((int16_t)p * QINV);
    return (int16_t)((p - (int32_t)t * Q) >> 16);
}
static inline int16_t barrett(int16_t a) {
    int16_t t = (int16_t)(((int32_t)BARR * a + (1 << 25)) >> 26);
    return (int16_t)(a - t * (int16_t)Q);
}

/* ── zetas (Montgomery-form, identici a kem.c) ───────────────────── */
static const int16_t ZETAS[128] = {
    -1044, -758, -359,-1517, 1493, 1422,  287,  202,
     -171,  622, 1577,  182,  962,-1202,-1474, 1468,
      573,-1325,  264,  383, -829, 1458,-1602, -130,
     -681, 1017,  732,  608,-1542,  411, -205,-1571,
     1223,  652, -552, 1015,-1293, 1491, -282,-1544,
      516,   -8, -320, -666,-1618,-1162,  126, 1469,
     -853,  -90, -271,  830,  107,-1421, -247, -951,
     -398,  961,-1508, -725,  448,-1065,  677,-1275,
    -1103,  430,  555,  843,-1251,  871, 1550,  105,
      422,  587,  177, -235, -291, -460, 1574, 1653,
     -246,  778, 1159, -147, -777, 1483, -602, 1119,
    -1590,  644, -872,  349,  418,  329, -156,  -75,
      817, 1097,  603,  610, 1322,-1285,-1465,  384,
    -1215, -136, 1218,-1335, -874,  220,-1187,-1659,
    -1185,-1530,-1278,  794,-1510, -854, -870,  478,
     -108, -308,  996,  991,  958,-1460, 1522, 1628
};

/* ── NTT (scalar, esposta come ntt32_*) ──────────────────────────── */
void ntt32_forward(int16_t r[256]) {
    int k=1, len, s, j;
    for (len=128; len>=2; len>>=1)
        for (s=0; s<256; s+=2*len) {
            int16_t z = ZETAS[k++];
            for (j=s; j<s+len; j++) {
                int16_t t = fqmul(z, r[j+len]);
                r[j+len] = barrett((int16_t)(r[j]-t));
                r[j]     = barrett((int16_t)(r[j]+t));
            }
        }
}
void ntt32_inverse(int16_t r[256]) {
    int k=127, len, s, j;
    for (len=2; len<=128; len<<=1)
        for (s=0; s<256; s+=2*len) {
            int16_t z = ZETAS[k--];
            for (j=s; j<s+len; j++) {
                int16_t t = r[j];
                r[j]     = barrett((int16_t)(t+r[j+len]));
                r[j+len] = fqmul(z,(int16_t)(r[j+len]-t));
            }
        }
    for (j=0; j<256; j++) r[j] = fqmul(r[j], INV_N_MONT);
}
void ntt32_pointwise_mul(int16_t *c, const int16_t *a, const int16_t *b) {
    int i;
    for (i=0; i<256; i++) c[i] = fqmul(a[i], b[i]);
}

/* ── packing / compressione ──────────────────────────────────────── */
static void poly_pack(uint8_t *r, const int16_t p[256]) {
    int i;
    for (i=0; i<256; i+=2) {
        int16_t a=p[i], b=p[i+1];
        while(a<0)a+=Q; while(a>=Q)a-=Q;
        while(b<0)b+=Q; while(b>=Q)b-=Q;
        r[3*(i/2)]  =(uint8_t)a;
        r[3*(i/2)+1]=(uint8_t)((a>>8)|(b<<4));
        r[3*(i/2)+2]=(uint8_t)(b>>4);
    }
}
static void poly_unpack(int16_t p[256], const uint8_t *r) {
    int i;
    for (i=0; i<256; i+=2) {
        p[i]  =(int16_t)(r[3*(i/2)]|((int16_t)(r[3*(i/2)+1]&0x0F)<<8));
        p[i+1]=(int16_t)((r[3*(i/2)+1]>>4)|((int16_t)r[3*(i/2)+2]<<4));
    }
}
static void poly_compress10(uint8_t *r, const int16_t p[256]) {
    int i, j;
    for (i=0; i<256; i+=4) {
        uint16_t t[4];
        for (j=0; j<4; j++) {
            int16_t x=p[i+j];
            while(x<0)x+=Q; while(x>=Q)x-=Q;
            t[j]=(uint16_t)(((uint32_t)x*1024+Q/2)/Q & 0x3FF);
        }
        r[5*(i/4)]  =(uint8_t)t[0];
        r[5*(i/4)+1]=(uint8_t)((t[0]>>8)|(t[1]<<2));
        r[5*(i/4)+2]=(uint8_t)((t[1]>>6)|(t[2]<<4));
        r[5*(i/4)+3]=(uint8_t)((t[2]>>4)|(t[3]<<6));
        r[5*(i/4)+4]=(uint8_t)(t[3]>>2);
    }
}
static void poly_decompress10(int16_t p[256], const uint8_t *r) {
    int i, j;
    for (i=0; i<256; i+=4) {
        uint16_t t[4];
        t[0]=(uint16_t)(r[5*(i/4)]  |((uint16_t)(r[5*(i/4)+1]&0x03)<<8));
        t[1]=(uint16_t)((r[5*(i/4)+1]>>2)|((uint16_t)(r[5*(i/4)+2]&0x0F)<<6));
        t[2]=(uint16_t)((r[5*(i/4)+2]>>4)|((uint16_t)(r[5*(i/4)+3]&0x3F)<<4));
        t[3]=(uint16_t)((r[5*(i/4)+3]>>6)|((uint16_t)r[5*(i/4)+4]<<2));
        for (j=0; j<4; j++) p[i+j]=(int16_t)(((uint32_t)t[j]*Q+512)>>10);
    }
}
static void poly_compress4(uint8_t *r, const int16_t p[256]) {
    int i;
    for (i=0; i<256; i+=2) {
        int16_t x0=p[i], x1=p[i+1];
        while(x0<0)x0+=Q; while(x0>=Q)x0-=Q;
        while(x1<0)x1+=Q; while(x1>=Q)x1-=Q;
        r[i/2]=(uint8_t)((((uint32_t)x0*16+Q/2)/Q&0x0F)|
                         ((((uint32_t)x1*16+Q/2)/Q&0x0F)<<4));
    }
}
static void poly_decompress4(int16_t p[256], const uint8_t *r) {
    int i;
    for (i=0; i<256; i+=2) {
        p[i]  =(int16_t)(((r[i/2]&0x0F)*Q+8)>>4);
        p[i+1]=(int16_t)(((r[i/2]>>4)*Q+8)>>4);
    }
}

/* ── Montgomery helpers ──────────────────────────────────────────── */
static void poly_tomont(int16_t r[256]) {
    int i;
    for (i=0; i<256; i++) r[i]=fqmul(r[i], R2);
}
static void poly_frommont(int16_t r[256]) {
    int i;
    for (i=0; i<256; i++) r[i]=fqmul(r[i], (int16_t)1);
}

/* ── CBD eta=2 ───────────────────────────────────────────────────── */
static void cbd2(int16_t r[256], const uint8_t *noise) {
    int l, m;
    for (l=0; l<256; l+=8) {
        uint32_t bv; memcpy(&bv, noise+l/2, 4);
        uint32_t d=(bv&0x55555555u)+((bv>>1)&0x55555555u);
        for (m=0; m<8; m++)
            r[l+m]=(int16_t)(((d>>(4*m))&3u)-((d>>(4*m+2))&3u));
    }
}

/* ── gen_matrix: SHAKE-128 streaming (rejection sampling robusto) ── */
static void gen_matrix32(int16_t A[2][2][256],
                         const uint8_t rho[32], int transposed) {
    int i, j, ctr, off;
    for (i=0; i<2; i++) for (j=0; j<2; j++) {
        uint8_t seed[34];
        memcpy(seed, rho, 32);
        seed[32] = transposed ? (uint8_t)j : (uint8_t)i;
        seed[33] = transposed ? (uint8_t)i : (uint8_t)j;
        KeccakState xof;
        shake128_init(&xof);
        shake128_absorb(&xof, seed, 34);
        shake128_finalize(&xof);
        ctr=0;
        while (ctr<256) {
            uint8_t xbuf[SHAKE128_RATE];
            shake128_squeeze(&xof, xbuf, SHAKE128_RATE);
            off=0;
            while (ctr<256 && off<=SHAKE128_RATE-3) {
                int16_t d1=(int16_t)(xbuf[off]|((int16_t)(xbuf[off+1]&0x0F)<<8));
                int16_t d2=(int16_t)((xbuf[off+1]>>4)|((int16_t)xbuf[off+2]<<4));
                off+=3;
                if (d1<Q) A[i][j][ctr++]=d1;
                if (ctr<256 && d2<Q) A[i][j][ctr++]=d2;
            }
        }
        poly_tomont(A[i][j]);
    }
}

/* ── pke_enc (interno) ───────────────────────────────────────────── */
static void pke_enc(Kem32Ciphertext *ct, const Kem32PublicKey *pk,
                    const uint8_t m[32], const uint8_t r_seed[32]) {
    int16_t that[2][256],rvhat[2][256],u[2][256],v[256];
    int16_t A[2][2][256],e1[2][256],e2[256];
    uint8_t rho[32];
    int i, j, l;
    memcpy(rho, pk->data+768, 32);
    for (i=0; i<2; i++) {
        poly_unpack(that[i], pk->data+i*384);
        poly_tomont(that[i]); ntt32_forward(that[i]);
    }
    gen_matrix32(A, rho, 1);   /* A^T */
    /* r (nonce 0,1) */
    for (i=0; i<2; i++) {
        uint8_t ext[33],ns[128];
        memcpy(ext,r_seed,32); ext[32]=(uint8_t)i;
        shake256(ns,128,ext,33); cbd2(rvhat[i],ns);
        poly_tomont(rvhat[i]); ntt32_forward(rvhat[i]);
    }
    /* e1 (nonce 2,3) */
    for (i=0; i<2; i++) {
        uint8_t ext[33],ns[128];
        memcpy(ext,r_seed,32); ext[32]=(uint8_t)(2+i);
        shake256(ns,128,ext,33); cbd2(e1[i],ns); poly_tomont(e1[i]);
    }
    /* e2 (nonce 4) */
    { uint8_t ext[33],ns[128];
      memcpy(ext,r_seed,32); ext[32]=4;
      shake256(ns,128,ext,33); cbd2(e2,ns); poly_tomont(e2); }
    /* u = INTT(A^T*r) + e1 */
    for (i=0; i<2; i++) {
        int16_t tmp[256];
        memset(u[i],0,512);
        for (j=0; j<2; j++) {
            ntt32_pointwise_mul(tmp,A[i][j],rvhat[j]);
            for (l=0; l<256; l++) u[i][l]=barrett((int16_t)(u[i][l]+tmp[l]));
        }
        ntt32_inverse(u[i]);
        for (l=0; l<256; l++) u[i][l]=barrett((int16_t)(u[i][l]+e1[i][l]));
        poly_frommont(u[i]);
    }
    /* v = INTT(t^T*r) + e2 + m */
    memset(v,0,512);
    for (j=0; j<2; j++) {
        int16_t tmp[256];
        ntt32_pointwise_mul(tmp,that[j],rvhat[j]);
        for (l=0; l<256; l++) v[l]=barrett((int16_t)(v[l]+tmp[l]));
    }
    ntt32_inverse(v);
    for (l=0; l<256; l++) {
        int16_t mbit=(int16_t)(((m[l/8]>>(l%8))&1)*1665);
        v[l]=barrett((int16_t)(v[l]+e2[l]+fqmul(mbit,R2)));
    }
    poly_frommont(v);
    { uint8_t *p=ct->data;
      for (i=0; i<2; i++) { poly_compress10(p,u[i]); p+=320; }
      poly_compress4(p,v); }
}

/* ── KEYGEN  (FIPS 203 §7.1: seed = d[32] || z[32]) ────────────── */
void kem32_keygen(Kem32PublicKey *pk, Kem32SecretKey *sk,
                  const uint8_t seed[64]) {
    const uint8_t *d=seed, *z=seed+32;
    uint8_t rho[32],sigma[32],buf[64];
    int16_t A[2][2][256],s[2][256],e[2][256],shat[2][256],t[2][256];
    int i, j, l;
    sha3_512(buf,d,32);
    memcpy(rho,buf,32); memcpy(sigma,buf+32,32);
    gen_matrix32(A, rho, 0);   /* A */
    /* s (nonce 0,1) */
    for (i=0; i<2; i++) {
        uint8_t ext[33],ns[128];
        memcpy(ext,sigma,32); ext[32]=(uint8_t)i;
        shake256(ns,128,ext,33); cbd2(s[i],ns);
        poly_tomont(s[i]); memcpy(shat[i],s[i],512); ntt32_forward(shat[i]);
    }
    /* e (nonce 2,3) */
    for (i=0; i<2; i++) {
        uint8_t ext[33],ns[128];
        memcpy(ext,sigma,32); ext[32]=(uint8_t)(2+i);
        shake256(ns,128,ext,33); cbd2(e[i],ns); poly_tomont(e[i]);
    }
    /* t = INTT(A*s) + e */
    for (i=0; i<2; i++) {
        int16_t tmp[256];
        memset(t[i],0,512);
        for (j=0; j<2; j++) {
            ntt32_pointwise_mul(tmp,A[i][j],shat[j]);
            for (l=0; l<256; l++) t[i][l]=barrett((int16_t)(t[i][l]+tmp[l]));
        }
        ntt32_inverse(t[i]);
        for (l=0; l<256; l++) t[i][l]=barrett((int16_t)(t[i][l]+e[i][l]));
        poly_frommont(t[i]);
    }
    { uint8_t *p=pk->data;
      for (i=0; i<2; i++) { poly_pack(p,t[i]); p+=384; }
      memcpy(p,rho,32); }
    { uint8_t *p=sk->data;
      for (i=0; i<2; i++) { poly_pack(p,shat[i]); p+=384; }
      memcpy(p,pk->data,800); p+=800;
      sha3_256(p,pk->data,800); p+=32;
      memcpy(p,z,32); }  /* z copiato direttamente (FIPS 203) */
}

/* ── ENCAPS ──────────────────────────────────────────────────────── */
void kem32_encaps(Kem32Ciphertext *ct, Kem32SharedSecret *ss,
                  const Kem32PublicKey *pk, const uint8_t rnd[32]) {
    uint8_t m[32],hpk[32],cin[64],kr[64];
    sha3_256(m,rnd,32); sha3_256(hpk,pk->data,800);
    memcpy(cin,m,32); memcpy(cin+32,hpk,32); sha3_512(kr,cin,64);
    pke_enc(ct,pk,m,kr+32);
    { uint8_t hct[32],ssin[64];
      sha3_256(hct,ct->data,768);
      memcpy(ssin,kr,32); memcpy(ssin+32,hct,32);
      sha3_256(ss->data,ssin,64); }
}

/* ── DECAPS ──────────────────────────────────────────────────────── */
int kem32_decaps(Kem32SharedSecret *ss, const Kem32Ciphertext *ct,
                 const Kem32SecretKey *sk) {
    const uint8_t *skp=sk->data;
    int16_t shat[2][256],u[2][256],v[256],w[256];
    uint8_t m[32];
    int i, l;
    for (i=0; i<2; i++) { poly_unpack(shat[i],skp); skp+=384; }
    const uint8_t *pkb=skp;
    const uint8_t *hpk=skp+800;
    const uint8_t *z  =hpk+32;
    { const uint8_t *p=ct->data;
      for (i=0; i<2; i++) { poly_decompress10(u[i],p); poly_tomont(u[i]); p+=320; }
      poly_decompress4(v,p); poly_tomont(v); }
    memset(w,0,512);
    for (i=0; i<2; i++) {
        int16_t uhat[256],tmp[256];
        memcpy(uhat,u[i],512); ntt32_forward(uhat);
        ntt32_pointwise_mul(tmp,shat[i],uhat);
        for (l=0; l<256; l++) w[l]=barrett((int16_t)(w[l]+tmp[l]));
    }
    ntt32_inverse(w);
    for (l=0; l<256; l++) w[l]=barrett((int16_t)(v[l]-w[l]));
    poly_frommont(w);
    memset(m,0,32);
    for (l=0; l<256; l++) {
        int16_t x=w[l]; while(x<0)x+=Q; while(x>=Q)x-=Q;
        int16_t d=(int16_t)(x-1665); if(d<0)d=(int16_t)-d;
        if(d<832) m[l/8]|=(uint8_t)(1<<(l%8));
    }
    { Kem32Ciphertext ct2; Kem32PublicKey pk_t;
      uint8_t kr[64],cin[64];
      memcpy(pk_t.data,pkb,800);
      memcpy(cin,m,32); memcpy(cin+32,hpk,32); sha3_512(kr,cin,64);
      pke_enc(&ct2,&pk_t,m,kr+32);
      uint64_t diff=0;
      for (i=0; i<768; i++) diff|=(uint64_t)(ct->data[i]^ct2.data[i]);
      uint8_t fail=(uint8_t)((-diff)>>63);
      uint8_t mask=(uint8_t)(-(int8_t)fail);
      for (i=0; i<32; i++) kr[i]^=(uint8_t)(mask&(kr[i]^z[i]));
      { uint8_t hct[32],ssin[64];
        sha3_256(hct,ct->data,768);
        memcpy(ssin,kr,32); memcpy(ssin+32,hct,32);
        sha3_256(ss->data,ssin,64); }
      return -(int)fail; }
}
