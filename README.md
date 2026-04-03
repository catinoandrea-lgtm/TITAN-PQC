# TITAN-PQC# TITAN-KEM-512: Multi-Architecture Accelerated Implementation

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/catinoandrea-lgtm/TITAN-PQC/blob/master/titan_pqc.ipynb)

This repository contains the official artifact for the hardware-accelerated implementation of the **TITAN-KEM-512** protocol. The code is highly optimized for **ARM64 (NEON)** for the full Key Encapsulation Mechanism, and includes a high-performance **x86_64 (AVX2/FMA)** vector math engine.

## 📂 Repository Structure
- **Core KEM (ARM/Scalar):** `kem.c`, `kem.h`, `ntt32.c`
- **AVX2 Vector Engine:** `titan_avx2.c`, `titan_avx2.h` (Includes 256-bit NTT, INTT, Pointwise, Keccakx4)
- **Symmetric Primitives:** `keccakf1600.c`, `sha256.c`
- **Automation:** Universal `Makefile` and Google Colab notebook (`titan_pqc.ipynb`).

---

## 🚀 Official Performance Benchmarks

### 1. Full KEM Protocol | Apple Silicon (M1) - ARM64 NEON
*Tested natively on Apple M1 processor (macOS). Compiled with Clang.*
| Operation | Time (µs) |
| :--- | :---: |
| **Keygen (Avg)** | **9.51** |
| **Encaps (Avg)** | **12.36** |
| **Decaps (Avg)** | **12.57** |
| NTT Forward (NEON) | 0.28 |
| NTT Inverse (NEON) | 0.24 |

### 2. Full KEM Protocol | Android Mobile (ARMv8 A55/A78)
*Tested on mobile environment representing low-power IoT/Edge devices.*
| Operation | Time (µs) |
| :--- | :---: |
| **Keygen (Avg)** | 67.86 |
| **Encaps (Avg)** | 90.13 |
| **Decaps (Avg)** | 99.64 |

### 3. Math Engine | x86_64 Server (Intel Xeon / AMD) - AVX2
*Tested on Google Colab Ubuntu instances with AVX2 & FMA enabled.*
| Operation | Time (µs) | Note |
| :--- | :---: | :--- |
| **Forward-Inverse Roundtrip** | **PASS** | 100% Mathematical correctness |
| Keccakx4 (Matrix A) | 1.609 | 4 parallel SHAKE128 streams |
| NTT Forward (AVX2) | 3.795 | 256-bit registers |
| NTT Inverse (AVX2) | 1.813 | 256-bit registers |
| Pointwise Mul | 0.015 | ~15 nanoseconds |

---

## 🛠 How to Build & Run
The provided `Makefile` automatically detects your architecture.

```bash
# 1. Clone the repository
git clone [https://github.com/catinoandrea-lgtm/TITAN-PQC.git](https://github.com/catinoandrea-lgtm/TITAN-PQC.git)
cd TITAN-PQC

# 2. Compile
make clean && make

# 3. Run the benchmarks
# On Apple M1:
./titan_m1_bench

# On x86_64 (Linux/Colab):
./titan_avx2_bench

# On Android/Linux ARM:
./titan_android_bench


## 📜 License
This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.
Copyright (c) 2026 Andrea Catino

### Academic Use
If you use this code in your research, please cite our work. This artifact is released to foster research in Post-Quantum Cryptography and hardware acceleration.
