/*
 * titan_ios.c — TITAN ML-KEM-512 (single file for iOS / Apple Silicon)
 * Andrea Catino — v2.3-FIPS203
 *
 * === a-Shell (iPhone/iPad) ===
 *   clang titan_ios.c -O3 -std=c99 -o titan_ios && ./titan_ios
 *
 * === Xcode (macOS terminal, Apple Silicon native) ===
 *   clang titan_ios.c -O3 -std=c99 -mcpu=apple-m1 -o titan_m1 && ./titan_m1
 *
 * === Xcode iOS project ===
 *   Add this file to your project, call titan_main() from Swift/ObjC.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <arm_neon.h>

/* ══════════════════════════════════════════════════════════════════
 *  KECCAK-f[1600]
 * ══════════════════════════════════════════════════════════════════ */
#define ROL64(a,n) (((uint64_t)(a)<<(n))|((uint64_t)(a)>>(64-(n))))

static const uint64_t RC[24] = {
    0x0000000000000001ULL,0x0000000000008082ULL,0x800000000000808aULL,
    0x8000000080008000ULL,0x000000000000808bULL,0x0000000080000001ULL,
    0x8000000080008081ULL,0x8000000000008009ULL,0x000000000000008aULL,
    0x0000000000000088ULL,0x0000000080008009ULL,0x000000008000000aULL,
    0x000000008000808bULL,0x800000000000008bULL,0x8000000000008089ULL,
    0x8000000000008003ULL,0x8000000000008002ULL,0x8000000000000080ULL,
    0x000000000000800aULL,0x800000008000000aULL,0x8000000080008081ULL,
    0x8000000000008080ULL,0x0000000080000001ULL,0x8000000080008008ULL
};
static const int ROTC[24]={1,3,6,10,15,21,28,36,45,55,2,14,27,41,56,8,25,43,62,18,39,61,20,44};
static const int PI[24]={10,7,11,17,18,3,5,16,8,21,24,4,15,23,19,13,12,2,20,14,22,9,6,1};

static void KeccakF1600(uint64_t s[25]) {
    uint64_t bc[5],t;
    for(int r=0;r<24;r++){
        for(int i=0;i<5;i++) bc[i]=s[i]^s[i+5]^s[i+10]^s[i+15]^s[i+20];
        for(int i=0;i<5;i++){t=bc[(i+4)%5]^ROL64(bc[(i+1)%5],1);for(int j=0;j<25;j+=5)s[j+i]^=t;}
        t=s[1]; for(int i=0;i<24;i++){int j=PI[i];bc[0]=s[j];s[j]=ROL64(t,ROTC[i]);t=bc[0];}
        for(int j=0;j<25;j+=5){for(int i=0;i<5;i++)bc[i]=s[j+i];for(int i=0;i<5;i++)s[j+i]=bc[i]^(~bc[(i+1)%5]&bc[(i+2)%5]);}
        s[0]^=RC[r];
    }
}
static void keccak_absorb(uint64_t s[25],unsigned rate,const uint8_t*in,size_t inlen,uint8_t pad){
    unsigned r64=rate/8;
    while(inlen>=rate){for(unsigned i=0;i<r64;i++){uint64_t tmp;memcpy(&tmp,in+8*i,8);s[i]^=tmp;}KeccakF1600(s);in+=rate;inlen-=rate;}
    uint8_t t[200]={0};memcpy(t,in,inlen);t[inlen]=pad;t[rate-1]|=0x80;
    for(unsigned i=0;i<r64;i++){uint64_t tmp;memcpy(&tmp,t+8*i,8);s[i]^=tmp;}
}
static void keccak_squeeze(uint64_t s[25],unsigned rate,uint8_t*out,size_t outlen){
    while(outlen>0){KeccakF1600(s);size_t n=(outlen<rate)?outlen:rate;memcpy(out,s,n);out+=n;outlen-=n;}
}
static void sha3_256(uint8_t out[32],const uint8_t*in,size_t inlen){uint64_t s[25]={0};keccak_absorb(s,136,in,inlen,0x06);KeccakF1600(s);memcpy(out,s,32);}
static void sha3_512(uint8_t out[64],const uint8_t*in,size_t inlen){uint64_t s[25]={0};keccak_absorb(s,72,in,inlen,0x06);KeccakF1600(s);memcpy(out,s,64);}
static void shake256(uint8_t*out,size_t outlen,const uint8_t*in,size_t inlen){uint64_t s[25]={0};keccak_absorb(s,136,in,inlen,0x1F);keccak_squeeze(s,136,out,outlen);}

/* Streaming SHAKE-128 */
typedef struct { uint64_t s[25]; } KeccakState;
static void shake128_init(KeccakState*st){memset(st,0,sizeof*st);}
static void shake128_absorb(KeccakState*st,const uint8_t*in,size_t l){keccak_absorb((uint64_t*)st,168,in,l,0x1F);}
static void shake128_finalize(KeccakState*st){(void)st;}
static void shake128_squeeze(KeccakState*st,uint8_t*out,size_t l){keccak_squeeze((uint64_t*)st,168,out,l);}

/* ══════════════════════════════════════════════════════════════════
 *  ML-KEM-512 (TITAN engine)
 * ══════════════════════════════════════════════════════════════════ */
#define Q           3329
#define QINV        ((int16_t)(-3327))
#define BARR        20159
#define R2          ((int16_t)1353)
#define INV_N_MONT  ((int16_t)512)
#define SHAKE128_RATE 168
#define KEM_PK_BYTES  800
#define KEM_SK_BYTES  1632
#define KEM_CT_BYTES  768
#define KEM_SS_BYTES  32

typedef struct { uint8_t data[KEM_PK_BYTES];  } KemPublicKey;
typedef struct { uint8_t data[KEM_SK_BYTES];  } KemSecretKey;
typedef struct { uint8_t data[KEM_CT_BYTES];  } KemCiphertext;
typedef struct { uint8_t data[KEM_SS_BYTES];  } KemSharedSecret;

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

/* ── NTT forward/inverse ─────────────────────────────────────────── */
static void ntt_forward(int16_t r[256]) {
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
        for(int i=0;i<256;i+=8) vst1q_s16(&r[i],barrett_neon(vld1q_s16(&r[i])));
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
static void ntt_inverse(int16_t r[256]) {
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
static void ntt_pointwise_mul(int16_t *c,const int16_t *a,const int16_t *b){
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
            int16_t x=p[i+j]; while(x<0)x+=Q; while(x>=Q)x-=Q;
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

/* ── gen_matrix: streaming SHAKE-128 ─────────────────────────────── */
static void gen_matrix(int16_t A[2][2][256],
                       const uint8_t rho[32], int transposed){
    for(int i=0;i<2;i++) for(int j=0;j<2;j++){
        uint8_t seed[34];
        memcpy(seed, rho, 32);
        seed[32] = transposed ? (uint8_t)j : (uint8_t)i;
        seed[33] = transposed ? (uint8_t)i : (uint8_t)j;
        KeccakState xof;
        shake128_init(&xof);
        shake128_absorb(&xof, seed, 34);
        shake128_finalize(&xof);
        int ctr=0;
        while(ctr<256){
            uint8_t buf[SHAKE128_RATE];
            shake128_squeeze(&xof, buf, SHAKE128_RATE);
            int off=0;
            while(ctr<256 && off<=SHAKE128_RATE-3){
                int16_t d1=(int16_t)(buf[off]|((int16_t)(buf[off+1]&0x0F)<<8));
                int16_t d2=(int16_t)((buf[off+1]>>4)|((int16_t)buf[off+2]<<4));
                off+=3;
                if(d1<Q) A[i][j][ctr++]=d1;
                if(ctr<256 && d2<Q) A[i][j][ctr++]=d2;
            }
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
    gen_matrix(A, rho, 1);
    for(int i=0;i<2;i++){
        uint8_t ext[33],ns[128]; memcpy(ext,r_seed,32); ext[32]=(uint8_t)i;
        shake256(ns,128,ext,33); cbd2(rv[i],ns); poly_tomont(rv[i]); ntt_forward(rv[i]);
    }
    for(int i=0;i<2;i++){
        uint8_t ext[33],ns[128]; memcpy(ext,r_seed,32); ext[32]=(uint8_t)(2+i);
        shake256(ns,128,ext,33); cbd2(e1[i],ns); poly_tomont(e1[i]);
    }
    { uint8_t ext[33],ns[128]; memcpy(ext,r_seed,32); ext[32]=4;
      shake256(ns,128,ext,33); cbd2(e2,ns); poly_tomont(e2); }
    for(int i=0;i<2;i++){
        memset(u[i],0,512);
        for(int j=0;j<2;j++){int16_t tmp[256]; ntt_pointwise_mul(tmp,A[i][j],rv[j]);
            for(int l=0;l<256;l++) u[i][l]=barrett((int16_t)(u[i][l]+tmp[l]));}
        ntt_inverse(u[i]);
        for(int l=0;l<256;l++) u[i][l]=barrett((int16_t)(u[i][l]+e1[i][l]));
        poly_frommont(u[i]);
    }
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

/* ── KEYGEN (FIPS 203: seed = d[32] || z[32]) ────────────────────── */
static void kem_keygen(KemPublicKey *pk,KemSecretKey *sk,const uint8_t seed[64]){
    const uint8_t *d=seed, *z=seed+32;
    uint8_t rho[32],sigma[32],buf[64];
    int16_t A[2][2][256],s[2][256],e[2][256],shat[2][256],t[2][256];
    sha3_512(buf,d,32); memcpy(rho,buf,32); memcpy(sigma,buf+32,32);
    gen_matrix(A, rho, 0);
    for(int i=0;i<2;i++){
        uint8_t ext[33],ns[128]; memcpy(ext,sigma,32); ext[32]=(uint8_t)i;
        shake256(ns,128,ext,33); cbd2(s[i],ns);
        poly_tomont(s[i]); memcpy(shat[i],s[i],512); ntt_forward(shat[i]);
    }
    for(int i=0;i<2;i++){
        uint8_t ext[33],ns[128]; memcpy(ext,sigma,32); ext[32]=(uint8_t)(2+i);
        shake256(ns,128,ext,33); cbd2(e[i],ns); poly_tomont(e[i]);
    }
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
    sha3_256(skp,pk->data,KEM_PK_BYTES); skp+=32;
    memcpy(skp,z,32);
}

/* ── ENCAPS ─────────────────────────────────────────────────────── */
static void kem_encaps(KemCiphertext *ct,KemSharedSecret *ss,
                const KemPublicKey *pk,const uint8_t rnd[32]){
    uint8_t m[32],hpk[32],cin[64],kr[64],hct[32],ssin[64];
    sha3_256(m,rnd,32); sha3_256(hpk,pk->data,KEM_PK_BYTES);
    memcpy(cin,m,32); memcpy(cin+32,hpk,32); sha3_512(kr,cin,64);
    pke_enc(ct,pk,m,kr+32);
    sha3_256(hct,ct->data,KEM_CT_BYTES);
    memcpy(ssin,kr,32); memcpy(ssin+32,hct,32); sha3_256(ss->data,ssin,64);
}

/* ── DECAPS ─────────────────────────────────────────────────────── */
static int kem_decaps(KemSharedSecret *ss,const KemCiphertext *ct,
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

/* ══════════════════════════════════════════════════════════════════
 *  BENCHMARK — with platform detection
 * ══════════════════════════════════════════════════════════════════ */
#include <sys/utsname.h>

static double now_us(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e6 + (double)ts.tv_nsec * 1e-3;
}

static const char* detect_platform(void) {
    struct utsname un;
    if (uname(&un) == 0) return un.machine;
    return "unknown";
}

/* Entry point — also callable as titan_main() from Swift/ObjC */
int titan_main(void);
int titan_main(void) {
    KemPublicKey pk; KemSecretKey sk; KemCiphertext ct;
    KemSharedSecret ss1, ss2;
    uint8_t seed[64] = {0x42};
    uint8_t rnd[32]  = {0x55};
    const int N = 1000;
    int16_t poly[256] = {1};

    printf("==========================================================\n");
    printf("  TITAN ML-KEM-512 v2.3 | ARM NEON Native Benchmark\n");
    printf("  Platform: %s\n", detect_platform());
    printf("==========================================================\n");

    kem_keygen(&pk, &sk, seed);
    kem_encaps(&ct, &ss1, &pk, rnd);
    kem_decaps(&ss2, &ct, &sk);

    if(memcmp(ss1.data, ss2.data, 32) == 0)
        printf("[OK] Correctness: PASS (shared secrets match)\n\n");
    else {
        printf("[!!] Correctness: FAIL\n\n");
        return 1;
    }

    /* Warmup */
    for(int i=0; i<100; i++) kem_keygen(&pk, &sk, seed);

    printf("[A] KEM PROTOCOL (N=%d)\n", N);
    double t0 = now_us();
    for(int i=0; i<N; i++) kem_keygen(&pk, &sk, seed);
    double kg = (now_us()-t0)/N;
    printf("    Keygen : %7.2f us\n", kg);

    t0 = now_us();
    for(int i=0; i<N; i++) kem_encaps(&ct, &ss1, &pk, rnd);
    double en = (now_us()-t0)/N;
    printf("    Encaps : %7.2f us\n", en);

    t0 = now_us();
    for(int i=0; i<N; i++) kem_decaps(&ss2, &ct, &sk);
    double de = (now_us()-t0)/N;
    printf("    Decaps : %7.2f us\n", de);

    printf("    --------------------------------\n");
    printf("    Total  : %7.2f us (full handshake)\n\n", kg+en+de);

    printf("[B] NTT ENGINE (N=%d)\n", N*10);
    t0 = now_us();
    for(int i=0; i<N*10; i++) ntt_forward(poly);
    printf("    NTT Fwd: %7.2f us\n", (now_us()-t0)/(N*10));

    t0 = now_us();
    for(int i=0; i<N*10; i++) ntt_inverse(poly);
    printf("    NTT Inv: %7.2f us\n", (now_us()-t0)/(N*10));

    printf("==========================================================\n");
    return 0;
}

int main(void) { return titan_main(); }
