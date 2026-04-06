# TITAN: Hardware-Accelerated ML-KEM-512 for ARM64 NEON and AVX2

[![Ubuntu AVX2 Build](https://github.com/catinoandrea-lgtm/TITAN-PQC/actions/workflows/ubuntu_avx2_test.yml/badge.svg)](https://github.com/catinoandrea-lgtm/TITAN-PQC/actions)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/catinoandrea-lgtm/TITAN-PQC/blob/master/titan_pqc.ipynb)

## Overview

**TITAN** is an optimized implementation of **ML-KEM-512** (FIPS 203) targeting ARM64 NEON and x86_64 AVX2 architectures. The primary contribution is a high-performance NTT/INTT engine using 16-bit Montgomery multiplication with hybrid NEON/scalar butterfly scheduling, designed for deployment on mobile and IoT devices.

The NTT engine is portable and can be reused for other lattice-based schemes (ML-DSA, FN-DSA, future NIST standards) without modification.

### Engineering Contributions
- **Hybrid NEON/scalar NTT**: layers with stride >= 8 use ARM NEON 8-wide butterflies; final layers (stride 4, 2) fall back to scalar with Barrett reduction at every butterfly, avoiding register pressure on narrow strides.
- **16-bit Montgomery core**: `fqmul` computes `a * b * R^{-1} mod Q` using `vmull` + `vmlsl` + `vshrn` (6 NEON instructions / 8 elements). Constants derived via Hensel lifting: `Q_INV = -3327 (mod 2^16)`, `R^2 mod Q = 1353`.
- **Streaming matrix generation**: SHAKE-128 incremental squeeze (168 B/block) with guaranteed termination for rejection sampling.
- **Constant-time decapsulation**: bitwise OR comparison + arithmetic masking for implicit rejection (no branches on secret data).
- **Cross-platform**: single codebase compiles natively on Apple Silicon, Android ARM64, and x86_64 via architecture-detecting Makefile.

### ML-KEM-512 Parameters (FIPS 203)
| Parameter | Value |
| :--- | :--- |
| Module rank `k` | 2 |
| Polynomial degree `n` | 256 |
| Modulus `q` | 3329 |
| Secret distribution `eta` | CBD(2) |
| `d_u` / `d_v` | 10 / 4 bits |
| Public key | 800 bytes |
| Secret key | 1632 bytes |
| Ciphertext | 768 bytes |
| Shared secret | 32 bytes |

## Performance Benchmarks

### 1. Full KEM Protocol | Apple M1 (ARM64 NEON)
*Clang -O3, N=1000 iterations, clock_gettime(CLOCK_MONOTONIC).*

| Operation | Time (us) |
| :--- | :---: |
| **Keygen** | **9.51** |
| **Encaps** | **12.36** |
| **Decaps** | **12.57** |
| NTT Forward | 0.28 |
| NTT Inverse | 0.24 |

### 2. Full KEM Protocol | Android ARM64 (Cortex-A55/A78)
*CxxDroid on mobile SoC.*

| Operation | Time (us) |
| :--- | :---: |
| **Keygen** | 67.86 |
| **Encaps** | 90.13 |
| **Decaps** | 99.64 |

### 3. AVX2 Math Engine | x86_64 (Google Colab)
*GCC -O3 -mavx2 -mfma, N=200000 iterations.*

| Operation | Time (us) |
| :--- | :---: |
| NTT Roundtrip | **PASS** |
| NTT Forward | 3.795 |
| NTT Inverse | 1.813 |
| Pointwise Mul | 0.015 |
| Keccakx4 | 1.609 |

## Repository Structure
| File | Description |
| :--- | :--- |
| `kem.c` / `kem.h` | Full ML-KEM-512 (keygen, encaps, decaps) — ARM NEON |
| `kem32.c` / `kem32.h` | Pure scalar C99 fallback — identical protocol |
| `kem_bridge.c` | API adapter: `kem_*` -> `kem32_*` for scalar builds |
| `ntt32.c` | Hybrid NEON/scalar NTT engine |
| `titan_avx2.c/.h` | AVX2 vector math engine |
| `keccakf1600.c/.h` | Keccak-f[1600] + SHA3/SHAKE (bulk and streaming) |
| `sha256.c/.h` | SHA-256 |
| `titan_constants.c` | Shared zeta table (AVX2 builds) |
| `titan_full_bench.c` | KEM protocol benchmark |
| `titan_bench.c` | Math engine benchmark (AVX2) |
| `Makefile` | Universal build (auto-detects arch) |
| `titan_pqc.ipynb` | Google Colab notebook |

## Build and Run

```bash
git clone https://github.com/catinoandrea-lgtm/TITAN-PQC.git
cd TITAN-PQC
make clean && make

./titan_m1_bench        # Apple Silicon
./titan_android_bench   # Android / Linux ARM64
./titan_avx2_bench      # x86_64
```

## Implementation Notes

### Montgomery Domain Tracking
All polynomial arithmetic operates in Montgomery domain. Explicit `poly_tomont` (multiply by R^2) and `poly_frommont` (multiply by 1, extracting from Montgomery form) gates are placed at encode/decode boundaries. The NTT zetas are pre-multiplied by R.

### Barrett Reduction
`v = 20159` (approximation of `2^26 / Q`). Applied at every butterfly level in both forward and inverse NTT to guarantee coefficient bounds of `|a| < Q` throughout the computation.

### FIPS 203 Conformance
- KeyGen accepts 64 bytes of randomness: `d[32] || z[32]` per FIPS 203 section 7.1
- Symmetric primitives: SHA3-256, SHA3-512, SHAKE-128, SHAKE-256 (Keccak-f[1600])
- IND-CCA2 transform with Fujisaki-Okamoto implicit rejection
- SHA-256 is included but not used in the KEM protocol

### Known Limitations
- Not yet validated against NIST KAT vectors (planned for next release)
- Zeta table is duplicated across source files (not centralized via shared header)
- AVX2 Keccakx4 demonstrates single-block vectorized squeeze only

## License
MIT License — Copyright (c) 2026 Andrea Catino

### Citation
```
A. Catino, "TITAN: Hardware-Accelerated NTT Engine for Post-Quantum
Cryptography on ARM64 NEON and AVX2," 2026.
https://github.com/catinoandrea-lgtm/TITAN-PQC
```
