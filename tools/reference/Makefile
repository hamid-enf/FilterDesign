# Reference validation vs SciPy
#   make            -> build the dump tool + run the comparison
#   make dump       -> build only the C dump tool
# Requires: python3 with numpy + scipy (development only, NOT a firmware dep)

CC      ?= gcc
CFLAGS  ?= -O2 -std=c99 -Wall
INCLUDE := -I../../include
LIB     := ../../libfiltercoeff.a

# Python with numpy+scipy: allow override with PY=...; auto-detect venvs.
PY ?= $(shell \
    python3 -c "import scipy" 2>/dev/null && echo python3 || \
    (command -v .venv/bin/python3 >/dev/null 2>&1 && echo .venv/bin/python3) || \
    (command -v /home/user/.venv-ref/bin/python3 >/dev/null 2>&1 && echo /home/user/.venv-ref/bin/python3) || \
    echo python3)

.PHONY: all dump run clean

all: dump
	$(PY) compare_reference.py ./fce_dump

dump: fce_dump

fce_dump: fce_dump.c $(LIB)
	$(CC) $(CFLAGS) $(INCLUDE) fce_dump.c $(LIB) -lm -o fce_dump

run: all

clean:
	rm -f fce_dump
