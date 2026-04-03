CC = gcc
CFLAGS = -O3 -std=gnu11 -Wall -Wextra

# Rilevamento Architettura e OS
UNAME_M := $(shell uname -m)
UNAME_S := $(shell uname -s)

# --- CONFIGURAZIONE x86_64 (Server AMD/Intel e Google Colab) ---
ifeq ($(UNAME_M),x86_64)
    CFLAGS += -mavx2 -mfma -march=native
    SRCS = titan_bench.c titan_avx2.c titan_constants.c
    TARGET = titan_avx2_bench

# --- CONFIGURAZIONE APPLE SILICON (M1/M2/M3) ---
else ifeq ($(UNAME_S),Darwin)
    CC = clang
    CFLAGS += -mcpu=apple-m1
    SRCS = titan_full_bench.c kem.c ntt32.c keccakf1600.c sha256.c
    TARGET = titan_m1_bench

# --- CONFIGURAZIONE ANDROID / LINUX ARM64 ---
else ifeq ($(UNAME_M),aarch64)
    CC = clang
    CFLAGS += -march=armv8-a+simd
    SRCS = titan_full_bench.c kem.c ntt32.c keccakf1600.c sha256.c
    TARGET = titan_android_bench
endif

all: $(TARGET)

$(TARGET):
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET) -lm
	@echo "=> Build completata! Esegui: ./$(TARGET)"

clean:
	rm -f titan_avx2_bench titan_m1_bench titan_android_bench
