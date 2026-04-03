/*
 * kem.c – TITAN-KEM-512  (ARM NEON, ML-KEM-512 compatible)
 * Andrea Catino – Independent Researcher, Italy
 *
 * v2.2 – FIX CRITICO PERFORMANCE:
 *   - Matrice A: shake128_squeeze(3 byte/volta) → shake128(buf,672,...)
 *     Riduzione: 315 permutazioni Keccak/poly → 4 permutazioni/poly
 *     Speedup atteso keygen/encaps/decaps: ~15-25x
 *   - e[] aggiunto in keygen (nonce 2,3) – correzione crittografica
 *   - poly_* static (no conflitti linker)
 *   - barrett() nell'accumulo NTT
 *
 * Compile: clang -std=c99 -O3 -march=armv8-a+simd -lm
 *   test_titan32.c kem.c ntt32.c keccakf1600.c sha256.c -o test_titan
 */
#include "kem.h"
#include "keccakf1600.h"
#include "sha256.h"
#include <string.h>
#include <stdint.h>
#include <arm_neon.h>

#define Q           3329
#define QINV        ((int16_t)(-3327))
#define BARR        20159
#define R2          ((int16_t)1353)
#define INV_N_MONT  ((int16_t)512)

/* 4 blocchi SHAKE-128 = 672 byte → ~364 coefficienti validi attesi > 256 */
#define SHAKE128_RATE   168
#define XOF_BUFLEN      (4 * SHAKE128_RATE)     /* 672 bytes */

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

/* ── NEON arithmetic ─────────────────────────────────────────────── */
static inline int16x8_t fqmul_neon(int16x8_t a, int16x8_t b) {
    const int16x8_t vq=vdupq_n_s16(Q),vqinv=vdupq_n_s16(QINV);
    int32x4_t pl=vmull_s16(vget_low_s16(a),vget_low_s16(b));
    int32x4_t ph=vmull_high_s16(a,b);
    int16x8_t p=vcombine_s16(vmovn_s32(pl),vmovn_s32(ph));
    int16x8_t t=vmulq_s16(p,vqinv);
    pl=vmlsl_s16(pl,vget_low_s16(t),vget_low_s16(vq));
    ph=vmlsl_high_s16(ph,t,vq);
    return vcombine_s16(vshrn_n_s32(pl,16),vshrn_n_s32(ph,16));
}
static inline int16x8_t barrett_neon(int16x8_t a) {
    const int16x8_t vb=vdupq_n_s16(BARR),vq=vdupq_n_s16(Q);
    int32x4_t al=vshrq_n_s32(vmull_s16(vget_low_s16(a),vget_low_s16(vb)),26);
    int32x4_t ah=vshrq_n_s32(vmull_high_s16(a,vb),26);
    return vsubq_s16(a,vmulq_s16(vcombine_s16(vmovn_s32(al),vmovn_s32(ah)),vq));
}
static inline int16_t fqmul(int16_t a,int16_t b){
    int32_t p=(int32_t)a*b; int16_t t=(int16_t)((int16_t)p*QINV);
    return (int16_t)((p-(int32_t)t*Q)>>16);
}
static inline int16_t barrett(int16_t a){
    int16_t t=(int16_t)(((int32_t)BARR*a+(1<<25))>>26);
    return (int16_t)(a-t*(int16_t)Q);
}

/* ── NTT forward/inverse (ibrido NEON+scalar) ───────────────────── */
void ntt_forward(int16_t r[256]) {
    int k=1,len,s,j;
    for(len=128;len>=8;len>>=1){
        for(s=0;s<256;s+=2*len){
            int16x8_t vz=vdupq_n_s16(ZETAS[k++]);
            for(j=s;j<s+len;j+=8){
                int16x8_t r0=vld1q_s16(&r[j]),r1=vld1q_s16(&r[j+len]);
                int16x8_t t=fqmul_neon(vz,r1);
                vst1q_s16(&r[j],    vaddq_s16(r0,t));
                vst1q_s16(&r[j+len],vsubq_s16(r0,t));
            }
        }
        for(int i=0;i<256;i+=8)
            vst1q_s16(&r[i],barrett_neon(vld1q_s16(&r[i])));
    }
    for(;len>=2;len>>=1)
        for(s=0;s<256;s+=2*len){
            int16_t z=ZETAS[k++];
            for(j=s;j<s+len;j++){
                int16_t t=fqmul(z,r[j+len]);
                r[j+len]=barrett((int16_t)(r[j]-t));
                r[j]    =barrett((int16_t)(r[j]+t));
            }
        }
}
void ntt_inverse(int16_t r[256]) {
    int k=127,len,s,j;
    for(len=2;len<=4;len<<=1)
        for(s=0;s<256;s+=2*len){
            int16_t z=ZETAS[k--];
            for(j=s;j<s+len;j++){
                int16_t t=r[j];
                r[j]    =barrett((int16_t)(t+r[j+len]));
                r[j+len]=fqmul(z,(int16_t)(r[j+len]-t));
            }
        }
    for(len=8;len<=128;len<<=1)
        for(s=0;s<256;s+=2*len){
            int16x8_t vz=vdupq_n_s16(ZETAS[k--]);
            for(j=s;j<s+len;j+=8){
                int16x8_t r0=vld1q_s16(&r[j]),r1=vld1q_s16(&r[j+len]);
                vst1q_s16(&r[j],    barrett_neon(vaddq_s16(r0,r1)));
                vst1q_s16(&r[j+len],fqmul_neon(vz,vsubq_s16(r1,r0)));
            }
        }
    int16x8_t vi=vdupq_n_s16(INV_N_MONT);
    for(j=0;j<256;j+=8) vst1q_s16(&r[j],fqmul_neon(vld1q_s16(&r[j]),vi));
}
void ntt_pointwise_mul(int16_t *c,const int16_t *a,const int16_t *b){
    for(int i=0;i<256;i+=8)
        vst1q_s16(&c[i],fqmul_neon(vld1q_s16(&a[i]),vld1q_s16(&b[i])));
}

/* ── Montgomery helpers ──────────────────────────────────────────── */
static void poly_tomont(int16_t r[256]){
    int16x8_t vr2=vdupq_n_s16(R2);
    for(int i=0;i<256;i+=8) vst1q_s16(&r[i],fqmul_neon(vld1q_s16(&r[i]),vr2));
}
static void poly_frommont(int16_t r[256]){
    int16x8_t v1=vdupq_n_s16(1);
    for(int i=0;i<256;i+=8) vst1q_s16(&r[i],fqmul_neon(vld1q_s16(&r[i]),v1));
}

/* ── CBD eta=2 ───────────────────────────────────────────────────── */
static void cbd2(int16_t r[256],const uint8_t *noise){
    for(int l=0;l<256;l+=8){
        uint32_t b; memcpy(&b,noise+(l/2),4);
        uint32_t d=(b&0x55555555u)+((b>>1)&0x55555555u);
        for(int m=0;m<8;m++)
            r[l+m]=(int16_t)(((d>>(4*m))&3u)-((d>>(4*m+2))&3u));
    }
}

/* ── Packing / Compression ──────────────────────────────────────── */
static void poly_pack(uint8_t *r,const int16_t p[256]){
    for(int i=0;i<256;i+=2){
        int16_t a=p[i],b=p[i+1];
        while(a<0)a+=Q; while(a>=Q)a-=Q;
        while(b<0)b+=Q; while(b>=Q)b-=Q;
        r[3*(i/2)]  =(uint8_t)a;
        r[3*(i/2)+1]=(uint8_t)((a>>8)|(b<<4));
        r[3*(i/2)+2]=(uint8_t)(b>>4);
    }
}
static void poly_unpack(int16_t p[256],const uint8_t *r){
    for(int i=0;i<256;i+=2){
        p[i]  =(int16_t)(r[3*(i/2)]|((int16_t)(r[3*(i/2)+1]&0x0F)<<8));
        p[i+1]=(int16_t)((r[3*(i/2)+1]>>4)|((int16_t)r[3*(i/2)+2]<<4));
    }
}
static void poly_compress10(uint8_t *r,const int16_t p[256]){
    for(int i=0;i<256;i+=4){
        uint16_t t[4];
        for(int j=0;j<4;j++){
            int16_t x=p[i+j];
            while(x<0)x+=Q; while(x>=Q)x-=Q;
            t[j]=(uint16_t)(((uint32_t)x*1024+Q/2)/Q&0x3FF);
        }
        r[5*(i/4)]  =(uint8_t)t[0];
        r[5*(i/4)+1]=(uint8_t)((t[0]>>8)|(t[1]<<2));
        r[5*(i/4)+2]=(uint8_t)((t[1]>>6)|(t[2]<<4));
        r[5*(i/4)+3]=(uint8_t)((t[2]>>4)|(t[3]<<6));
        r[5*(i/4)+4]=(uint8_t)(t[3]>>2);
    }
}
static void poly_decompress10(int16_t p[256],const uint8_t *r){
    for(int i=0;i<256;i+=4){
        uint16_t t[4];
        t[0]=(uint16_t)(r[5*(i/4)]|((uint16_t)(r[5*(i/4)+1]&0x03)<<8));
        t[1]=(uint16_t)((r[5*(i/4)+1]>>2)|((uint16_t)(r[5*(i/4)+2]&0x0F)<<6));
        t[2]=(uint16_t)((r[5*(i/4)+2]>>4)|((uint16_t)(r[5*(i/4)+3]&0x3F)<<4));
        t[3]=(uint16_t)((r[5*(i/4)+3]>>6)|((uint16_t)r[5*(i/4)+4]<<2));
        for(int j=0;j<4;j++) p[i+j]=(int16_t)(((uint32_t)t[j]*Q+512)>>10);
    }
}
static void poly_compress4(uint8_t *r,const int16_t p[256]){
    for(int i=0;i<256;i+=2){
        int16_t x0=p[i],x1=p[i+1];
        while(x0<0)x0+=Q; while(x0>=Q)x0-=Q;
        while(x1<0)x1+=Q; while(x1>=Q)x1-=Q;
        r[i/2]=(uint8_t)((((uint32_t)x0*16+Q/2)/Q&0x0F)|
                         ((((uint32_t)x1*16+Q/2)/Q&0x0F)<<4));
    }
}
static void poly_decompress4(int16_t p[256],const uint8_t *r){
    for(int i=0;i<256;i+=2){
        p[i]  =(int16_t)(((r[i/2]&0x0F)*Q+8)>>4);
        p[i+1]=(int16_t)(((r[i/2]>>4)*Q+8)>>4);
    }
}

/* ── gen_matrix: bulk SHAKE-128, 672 byte per poly ──────────────── */
/*   FIX CRITICO: una sola chiamata shake128 = 4 permutazioni/poly   */
/*   vs. ~315 permutazioni/poly con streaming 3-byte-at-a-time       */
static void gen_matrix(int16_t A[2][2][256],
                       const uint8_t rho[32], int transposed){
    for(int i=0;i<2;i++) for(int j=0;j<2;j++){
        uint8_t seed[34];
        memcpy(seed, rho, 32);
        seed[32] = transposed ? (uint8_t)j : (uint8_t)i;
        seed[33] = transposed ? (uint8_t)i : (uint8_t)j;
        uint8_t buf[XOF_BUFLEN];
        shake128(buf, XOF_BUFLEN, seed, 34);  /* 4 perm totali */
        int ctr=0, off=0;
        while(ctr<256 && off<=XOF_BUFLEN-3){
            int16_t d1=(int16_t)(buf[off]|((int16_t)(buf[off+1]&0x0F)<<8));
            int16_t d2=(int16_t)((buf[off+1]>>4)|((int16_t)buf[off+2]<<4));
            off+=3;
            if(d1<Q) A[i][j][ctr++]=d1;
            if(ctr<256 && d2<Q) A[i][j][ctr++]=d2;
        }
        poly_tomont(A[i][j]);
    }
}

/* ── PKE encrypt ────────────────────────────────────────────────── */
static void pke_enc(KemCiphertext *ct,const KemPublicKey *pk,
                    const uint8_t msg[32],const uint8_t r_seed[32]){
    int16_t that[2][256],rv[2][256],u[2][256],v[256];
    int16_t A[2][2][256],e1[2][256],e2[256];
    uint8_t rho[32]; memcpy(rho,pk->data+768,32);
    for(int i=0;i<2;i++){poly_unpack(that[i],pk->data+i*384);
                         poly_tomont(that[i]); ntt_forward(that[i]);}
    /* A^T: transposed=1 */
    gen_matrix(A, rho, 1);
    /* r (nonce 0,1) */
    for(int i=0;i<2;i++){
        uint8_t ext[33],ns[128]; memcpy(ext,r_seed,32); ext[32]=(uint8_t)i;
        shake256(ns,128,ext,33); cbd2(rv[i],ns); poly_tomont(rv[i]); ntt_forward(rv[i]);
    }
    /* e1 (nonce 2,3) */
    for(int i=0;i<2;i++){
        uint8_t ext[33],ns[128]; memcpy(ext,r_seed,32); ext[32]=(uint8_t)(2+i);
        shake256(ns,128,ext,33); cbd2(e1[i],ns); poly_tomont(e1[i]);
    }
    /* e2 (nonce 4) */
    { uint8_t ext[33],ns[128]; memcpy(ext,r_seed,32); ext[32]=4;
      shake256(ns,128,ext,33); cbd2(e2,ns); poly_tomont(e2); }
    /* u = INTT(A^T * r) + e1 */
    for(int i=0;i<2;i++){
        memset(u[i],0,512);
        for(int j=0;j<2;j++){int16_t tmp[256]; ntt_pointwise_mul(tmp,A[i][j],rv[j]);
            for(int l=0;l<256;l++) u[i][l]=barrett((int16_t)(u[i][l]+tmp[l]));}
        ntt_inverse(u[i]);
        for(int l=0;l<256;l++) u[i][l]=barrett((int16_t)(u[i][l]+e1[i][l]));
        poly_frommont(u[i]);
    }
    /* v = INTT(t^T * r) + e2 + m */
    memset(v,0,512);
    for(int j=0;j<2;j++){int16_t tmp[256]; ntt_pointwise_mul(tmp,that[j],rv[j]);
        for(int l=0;l<256;l++) v[l]=barrett((int16_t)(v[l]+tmp[l]));}
    ntt_inverse(v);
    for(int l=0;l<256;l++){
        int16_t mb=(int16_t)(((msg[l/8]>>(l%8))&1)*1665);
        v[l]=barrett((int16_t)(v[l]+e2[l]+fqmul(mb,R2)));
    }
    poly_frommont(v);
    uint8_t *ctp=ct->data;
    for(int i=0;i<2;i++){poly_compress10(ctp,u[i]); ctp+=320;}
    poly_compress4(ctp,v);
}

/* ── KEYGEN ─────────────────────────────────────────────────────── */
void kem_keygen(KemPublicKey *pk,KemSecretKey *sk,const uint8_t seed[32]){
    uint8_t rho[32],sigma[32],buf[64];
    int16_t A[2][2][256],s[2][256],e[2][256],shat[2][256],t[2][256];
    sha3_512(buf,seed,32); memcpy(rho,buf,32); memcpy(sigma,buf+32,32);
    gen_matrix(A, rho, 0);          /* A (non-transposed) */
    /* s (nonce 0,1) */
    for(int i=0;i<2;i++){
        uint8_t ext[33],ns[128]; memcpy(ext,sigma,32); ext[32]=(uint8_t)i;
        shake256(ns,128,ext,33); cbd2(s[i],ns);
        poly_tomont(s[i]); memcpy(shat[i],s[i],512); ntt_forward(shat[i]);
    }
    /* e (nonce 2,3) – CORREZIONE CRITTOGRAFICA */
    for(int i=0;i<2;i++){
        uint8_t ext[33],ns[128]; memcpy(ext,sigma,32); ext[32]=(uint8_t)(2+i);
        shake256(ns,128,ext,33); cbd2(e[i],ns); poly_tomont(e[i]);
    }
    /* t = INTT(A*s) + e */
    for(int i=0;i<2;i++){
        memset(t[i],0,512);
        for(int j=0;j<2;j++){int16_t tmp[256]; ntt_pointwise_mul(tmp,A[i][j],shat[j]);
            for(int l=0;l<256;l++) t[i][l]=barrett((int16_t)(t[i][l]+tmp[l]));}
        ntt_inverse(t[i]);
        for(int l=0;l<256;l++) t[i][l]=barrett((int16_t)(t[i][l]+e[i][l]));
        poly_frommont(t[i]);
    }
    uint8_t *pkp=pk->data;
    for(int i=0;i<2;i++){poly_pack(pkp,t[i]); pkp+=384;} memcpy(pkp,rho,32);
    uint8_t *skp=sk->data;
    for(int i=0;i<2;i++){poly_pack(skp,shat[i]); skp+=384;}
    memcpy(skp,pk->data,KEM_PK_BYTES); skp+=KEM_PK_BYTES;
    sha3_256(skp,pk->data,KEM_PK_BYTES); skp+=32; sha256(skp,seed,32);
}

/* ── ENCAPS ─────────────────────────────────────────────────────── */
void kem_encaps(KemCiphertext *ct,KemSharedSecret *ss,
                const KemPublicKey *pk,const uint8_t rnd[32]){
    uint8_t m[32],hpk[32],cin[64],kr[64],hct[32],ssin[64];
    sha3_256(m,rnd,32); sha3_256(hpk,pk->data,KEM_PK_BYTES);
    memcpy(cin,m,32); memcpy(cin+32,hpk,32); sha3_512(kr,cin,64);
    pke_enc(ct,pk,m,kr+32);
    sha3_256(hct,ct->data,KEM_CT_BYTES);
    memcpy(ssin,kr,32); memcpy(ssin+32,hct,32); sha3_256(ss->data,ssin,64);
}

/* ── DECAPS ─────────────────────────────────────────────────────── */
int kem_decaps(KemSharedSecret *ss,const KemCiphertext *ct,
               const KemSecretKey *sk){
    const uint8_t *skp=sk->data;
    int16_t shat[2][256],u[2][256],v[256],w[256]; uint8_t m[32];
    for(int i=0;i<2;i++){poly_unpack(shat[i],skp); skp+=384;}
    const uint8_t *pkb=skp, *hpk=skp+KEM_PK_BYTES, *z=hpk+32;
    const uint8_t *ctp=ct->data;
    for(int i=0;i<2;i++){poly_decompress10(u[i],ctp); poly_tomont(u[i]); ctp+=320;}
    poly_decompress4(v,ctp); poly_tomont(v);
    memset(w,0,512);
    for(int i=0;i<2;i++){
        int16_t uh[256],tmp[256]; memcpy(uh,u[i],512);
        ntt_forward(uh); ntt_pointwise_mul(tmp,shat[i],uh);
        for(int l=0;l<256;l++) w[l]=barrett((int16_t)(w[l]+tmp[l]));
    }
    ntt_inverse(w);
    for(int l=0;l<256;l++) w[l]=barrett((int16_t)(v[l]-w[l]));
    poly_frommont(w);
    memset(m,0,32);
    for(int l=0;l<256;l++){
        int16_t x=w[l]; while(x<0)x+=Q; while(x>=Q)x-=Q;
        int16_t d=(int16_t)(x-1665); if(d<0)d=(int16_t)-d;
        if(d<832) m[l/8]|=(uint8_t)(1<<(l%8));
    }
    KemCiphertext ct2; KemPublicKey pk_t; uint8_t kr[64],cin[64];
    memcpy(pk_t.data,pkb,KEM_PK_BYTES);
    memcpy(cin,m,32); memcpy(cin+32,hpk,32); sha3_512(kr,cin,64);
    pke_enc(&ct2,&pk_t,m,kr+32);
    uint64_t diff=0;
    for(int i=0;i<KEM_CT_BYTES;i++) diff|=(uint64_t)(ct->data[i]^ct2.data[i]);
    uint8_t fail=(uint8_t)((-diff)>>63), mask=(uint8_t)(-(int8_t)fail);
    for(int i=0;i<32;i++) kr[i]^=(uint8_t)(mask&(kr[i]^z[i]));
    uint8_t hct[32],ssin[64];
    sha3_256(hct,ct->data,KEM_CT_BYTES);
    memcpy(ssin,kr,32); memcpy(ssin+32,hct,32); sha3_256(ss->data,ssin,64);
    return -(int)fail;
}
