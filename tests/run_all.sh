#!/bin/sh
# Build and run every test binary.
set -e
cd "$(dirname "$0")"
CC="${CC:-gcc}"
CFLAGS="${CFLAGS:--O2 -std=c99 -Wall -Wextra -Wpedantic}"
for t in test_core test_math test_fir test_iir test_quant test_validate \
         test_export test_edge test_sim; do
    echo "=== $t ==="
    $CC $CFLAGS -I../include -I../src -I. $t.c ../libfiltercoeff.a -lm -o /tmp/fce_$t
    /tmp/fce_$t
done
echo "ALL TESTS PASSED"
