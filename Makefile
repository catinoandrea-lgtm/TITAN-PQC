CC = clang
CFLAGS = -O3 -std=c99 -Wall -Wextra

#  M1 vs Android
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

ifeq ($(UNAME_S),Darwin)
    CFLAGS += -mcpu=apple-m1
    TARGET = titan_m1_bench
else ifeq ($(UNAME_M),aarch64)
    CFLAGS += -march=armv8-a+simd
    TARGET = titan_android_bench
endif


SRCS = titan_full_bench.c kem.c ntt32.c keccakf1600.c sha256.c
OBJS = $(SRCS:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET) -lm
	@echo "=> Build completata! Esegui: ./$(TARGET)"

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o titan_m1_bench titan_android_bench
