# FilterCoeff - host build (GCC/Clang)
# Targets: all, lib, test, ref, bench, examples, clean

CC      ?= gcc
CFLAGS  ?= -O2 -std=c99 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion
LDLIBS  += -lm

INCLUDE := -Iinclude

SRC := \
	src/fce_math.c \
	src/fce_fir.c \
	src/fce_iir.c \
	src/fce_quant.c \
	src/fce_validate.c \
	src/fce_sim.c \
	src/fce_export.c \
	src/fce_export_stdio.c \
	src/fce_generate.c

OBJ := $(SRC:.c=.o)
LIB := libfiltercoeff.a

TEST_SRC := $(wildcard tests/test_*.c)
TEST_OBJ := $(TEST_SRC:.c=.o)

EXAMPLE_SRC := $(wildcard examples/example_*.c)
EXAMPLE_BIN := $(EXAMPLE_SRC:.c=)

.PHONY: all lib test ref bench examples clean

all: lib

lib: $(LIB)

$(LIB): $(OBJ)
	$(AR) rcs $@ $^

%.o: %.c include/filtercoeff.h include/filtercoeff_config.h src/fce_internal.h
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

test: $(LIB) $(TEST_OBJ)
	$(CC) $(CFLAGS) $(INCLUDE) $(TEST_OBJ) $(LIB) $(LDLIBS) -o tests/run_tests
	./tests/run_tests

ref: $(LIB)
	$(MAKE) -C tools/reference run

bench: $(LIB)
	$(CC) $(CFLAGS) $(INCLUDE) bench/bench_design.c $(LIB) $(LDLIBS) -o bench/bench_design
	./bench/bench_design

examples: $(LIB) $(EXAMPLE_BIN)

$(EXAMPLE_BIN): %: %.c $(LIB)
	$(CC) $(CFLAGS) $(INCLUDE) $< $(LIB) $(LDLIBS) -o $@

clean:
	rm -f $(OBJ) $(LIB) tests/run_tests tests/*.o
	rm -f $(EXAMPLE_BIN) tools/reference/fce_dump bench/bench_design
