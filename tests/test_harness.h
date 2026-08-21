/*
 * test_harness.h - minimal test framework (no external dependencies).
 */
#ifndef TEST_HARNESS_H
#define TEST_HARNESS_H

#include <stdio.h>
#include <math.h>
#include <string.h>

extern int g_tests_run;
extern int g_tests_failed;

#define TEST_ASSERT(cond) \
    do { \
        g_tests_run++; \
        if (!(cond)) { \
            g_tests_failed++; \
            printf("  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        } \
    } while (0)

#define TEST_NEAR(a, b, tol) \
    do { \
        double _a = (a), _b = (b), _t = (tol); \
        g_tests_run++; \
        if (!(fabs(_a - _b) <= _t)) { \
            g_tests_failed++; \
            printf("  FAIL %s:%d: |%.17g - %.17g| = %.3g > %.3g\n", \
                   __FILE__, __LINE__, _a, _b, fabs(_a - _b), _t); \
        } \
    } while (0)

#define TEST_OK(st) \
    do { \
        fce_status_t _st = (st); \
        g_tests_run++; \
        if (_st != FCE_OK) { \
            g_tests_failed++; \
            printf("  FAIL %s:%d: status=%d (%s)\n", __FILE__, __LINE__, \
                   (int)_st, fce_status_str(_st)); \
        } \
    } while (0)

#endif /* TEST_HARNESS_H */
