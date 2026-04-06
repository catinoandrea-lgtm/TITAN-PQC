# TITAN: Hardware-Accelerated ML-KEM-512 for ARM64 NEON and AVX2

[![KAT Tests](https://github.com/catinoandrea-lgtm/TITAN-PQC/actions/workflows/kat_test.yml/badge.svg)](https://github.com/catinoandrea-lgtm/TITAN-PQC/actions/workflows/kat_test.yml)
[![AVX2 Build](https://github.com/catinoandrea-lgtm/TITAN-PQC/actions/workflows/ubuntu_avx2_test.yml/badge.svg)](https://github.com/catinoandrea-lgtm/TITAN-PQC/actions/workflows/ubuntu_avx2_test.yml)
[![Baseline](https://github.com/catinoandrea-lgtm/TITAN-PQC/actions/workflows/baseline_test.yml/badge.svg)](https://github.com/catinoandrea-lgtm/TITAN-PQC/actions/workflows/baseline_test.yml)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/catinoandrea-lgtm/TITAN-PQC/blob/master/titan_pqc.ipynb)

## Overview

**TITAN** is an optimized implementation of **ML-KEM-512** (FIPS 203) featuring a high-performance NTT engine targeting ARM64 NEON and x86_64 AVX2. The primary contribution is a hybrid NEON/scalar NTT with 16-bit Montgomery multiplication, designed for mobile and IoT deployment.

### Engineering Contributions
- **Hybrid NEON/scalar NTT**: 8-wide NEON butterflies for stride >= 8, scalar with Barrett reduction for stride 4 and 2
- **Montgomery core**: `vmull` + `vmlsl` + `vshrn` (6 NEON instructions / 8 elements), constants via Hensel lifting
- **Streaming matrix generation**: incremental SHAKE-128 squeeze with guaranteed termination
- **Constant-time decapsulation**: bitwise OR comparison + arithmetic masking (no secret-dependent branches)
- **Cross-platform**: single codebase for Apple Silicon, Android ARM64, x86_64, with auto-detecting Makefile

### ML-KEM-512 Parameters (FIPS 203)
| Parameter | Value |
| :--- | :--- |
| Module rank k / Degree n / Modulus q | 2 / 256 / 3329 |
| Secret distribution / d_u / d_v | CBD(2) / 10 bits / 4 bits |
| Public key / Secret key / Ciphertext / Shared secret | 800 / 1632 / 768 / 32 bytes |

## Performance Benchmarks

### Full KEM Protocol (ARM64 NEON)

| Platform | Keygen (us) | Encaps (us) | Decaps (us) | Handshake (us) |
| :--- | :---: | :---: | :---: | :---: |
| **Apple M1** (GitHub CI) | 8.13 | 11.55 | 11.92 | **31.60** |
| **Samsung S22 Ultra** (Cortex-X2) | 26.16 | 34.04 | 28.60 | **88.79** |
| **Android mid-range** (A55/A78) | 63.87 | 84.88 | 94.68 | **243.43** |

### NTT Engine (all platforms)

| Platform | NTT Fwd (us) | NTT Inv (us) | Pointwise (us) |
| :--- | :---: | :---: | :---: |
| Apple M1 (NEON) | 0.29 | 0.25 | -- |
| Samsung S22 Ultra (NEON) | 0.87 | 0.51 | -- |
| Android mid-range (NEON) | 4.15 | 1.88 | -- |
| AMD EPYC 7742 (AVX2) | 2.40 | 1.32 | 0.020 |
| Google Colab (AVX2) | 6.08 | 1.92 | 0.016 |

### KAT Conformance

Cross-platform deterministic fingerprint (100 accumulated vectors):
```
KAT-100: 807c873a3ad4a624b16bfcd44018d8cf0949da547561fe99daf5039cdf1342bd
```
Verified identical on: x86_64 Ubuntu, Android ARM64 (Samsung S22 Ultra, mid-range device).

## Repository Structure

| File | Description |
| :--- | :--- |
| `kem.c` / `kem.h` | ML-KEM-512 with ARM NEON |
| `kem32.c` / `kem32.h` | Pure scalar C99 fallback |
| `kem_bridge.c` | API adapter `kem_*` -> `kem32_*` |
| `ntt32.c` | Hybrid NEON/scalar NTT engine |
| `titan_avx2.c` / `.h` | AVX2 vector math engine |
| `keccakf1600.c` / `.h` | Keccak-f[1600] + SHA3/SHAKE |
| `sha256.c` / `.h` | SHA-256 (utility, not used in KEM) |
| `titan_constants.c` | Shared zeta table (AVX2) |
| `titan_full_bench.c` | KEM protocol benchmark (multi-file ARM) |
| `titan_bench.c` | Math engine benchmark (AVX2) |
| `titan_kat_test.c` | KAT conformance test suite (9 tests) |
| `titan_baseline.sh` | Baseline comparison vs pqcrystals/kyber ref |
| `titan_cxxdroid.c` | Single-file build for Android (CxxDroid) |
| `titan_ios.c` | Single-file build for iOS (a-Shell) / macOS |
| `Makefile` | Universal build (auto-detects architecture) |
| `titan_pqc.ipynb` | Google Colab notebook |

## Build and Run

### Benchmarks
```bash
git clone https://github.com/catinoandrea-lgtm/TITAN-PQC.git
cd TITAN-PQC
make clean && make       # auto-detects: ARM64/Apple/x86

./titan_m1_bench         # Apple Silicon
./titan_android_bench    # Android / Linux ARM64
./titan_avx2_bench       # x86_64
```

### KAT Tests (any platform)
```bash
make kat && ./titan_kat
```

### Baseline Comparison (x86_64, requires git + bc)
```bash
bash titan_baseline.sh
```

### Mobile builds
```bash
# Android (CxxDroid) — open titan_cxxdroid.c, set language to C, run
# iOS (a-Shell) —
clang titan_ios.c -O3 -std=c99 -o titan && ./titan
```

## KAT Test Suite

The `titan_kat_test.c` file runs 9 conformance tests:

| # | Test | Description |
| :--- | :--- | :--- |
| 1 | SHA3-256 | Empty string against NIST known answer |
| 2 | SHA3-512 | Empty string against NIST known answer |
| 3 | SHAKE-128 | Empty string (32B output) against NIST known answer |
| 4 | NTT roundtrip | forward(inverse(x)) == x for all coefficients |
| 5 | Key sizes | pk=800, sk=1632, ct=768, ss=32 per FIPS 203 |
| 6 | Determinism | Same seed produces identical keys and ciphertexts |
| 7 | Correctness | encaps/decaps shared secrets match |
| 8 | Implicit rejection | Modified ciphertext produces different shared secret |
| 9 | Accumulated KAT | 100 deterministic vectors with cross-platform fingerprint |

## Implementation Notes

### Montgomery Arithmetic
`Q_INV = -3327 (mod 2^16)`, `R^2 mod Q = 1353`, `N^{-1} * R mod Q = 512`. Barrett reduction `v = 20159` applied at every butterfly level.

### FIPS 203 Conformance
- KeyGen: 64 bytes input `d[32] || z[32]`, z stored directly in secret key
- All symmetric ops: SHA3-256/SHA3-512/SHAKE-128/SHAKE-256 (Keccak-f[1600])
- IND-CCA2 with Fujisaki-Okamoto implicit rejection

### Known Limitations
- Not yet validated against NIST ACVP test vectors (accumulated self-test only)
- Zeta table duplicated across source files
- AVX2 Keccakx4 is single-block squeeze only

## License
MIT License — Copyright (c) 2026 Andrea Catino

### Citation
```
A. Catino, "TITAN: Hardware-Accelerated NTT Engine for Post-Quantum
Cryptography on ARM64 NEON and AVX2," 2026.
https://github.com/catinoandrea-lgtm/TITAN-PQC
```
