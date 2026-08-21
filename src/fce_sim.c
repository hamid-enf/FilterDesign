/*
 * fce_sim.c - light-weight signal simulation.
 *
 * This module exists ONLY to sanity-check generated coefficients
 * (impulse response matches coefficients, sine amplitude matches the
 * frequency response, ...). It is intentionally NOT a production filter
 * runtime: no state-of-the-art optimizations, no DMA, no streaming.
 */
#include "fce_internal.h"

#if FCE_ENABLE_SIMULATION

void fce_sim_fir(const double* h, uint16_t n,
                 const double* x, double* y, uint32_t nx)
{
    uint32_t i, j;
    for (i = 0; i < nx; i++)
    {
        double acc = 0.0;
        for (j = 0; j < n && j <= i; j++)
            acc += h[j] * x[i - j];
        y[i] = acc;
    }
}

void fce_sim_sos(const double* sos, uint16_t n_sections,
                 const double* x, double* y, uint32_t nx,
                 double* state /* 2 per section: w1, w2 (transposed DF-II) */)
{
    uint32_t i, s;
    for (i = 0; i < nx; i++)
    {
        double v = x[i];
        for (s = 0; s < n_sections; s++)
        {
            const double* c = sos + 5u * s;
            double w1 = state[2u * s];
            double w2 = state[2u * s + 1u];
            double wn = v - c[3] * w1 - c[4] * w2; /* a0 = 1 */
            double out = c[0] * wn + c[1] * w1 + c[2] * w2;
            state[2u * s] = wn;
            state[2u * s + 1u] = w1;
            v = out;
        }
        y[i] = v;
    }
}

uint32_t fce_sim_signal(fce_sim_signal_t kind, double a, double b,
                        double fs, double* out, uint32_t n)
{
    uint32_t i;
    if (out == NULL || n == 0u)
        return 0u;

    switch (kind)
    {
    case FCE_SIM_SINE:
        for (i = 0; i < n; i++)
            out[i] = b * sin(2.0 * FCE_PI * a * (double)i / fs);
        return n;

    case FCE_SIM_MULTITONE:
        for (i = 0; i < n; i++)
        {
            double t = (double)i / fs;
            out[i] = b * (sin(2.0 * FCE_PI * a * t)
                        + sin(2.0 * FCE_PI * b * t)
                        + sin(2.0 * FCE_PI * (2.0 * a + b) * t));
        }
        return n;

    case FCE_SIM_IMPULSE:
        for (i = 0; i < n; i++)
            out[i] = (i == (uint32_t)a) ? b : 0.0;
        return n;

    case FCE_SIM_STEP:
        for (i = 0; i < n; i++)
            out[i] = (i >= (uint32_t)a) ? b : 0.0;
        return n;

    case FCE_SIM_WHITE_NOISE:
    {
        uint32_t seed = (uint32_t)a;
        for (i = 0; i < n; i++)
        {
            /* xorshift32 */
            seed ^= seed << 13;
            seed ^= seed >> 17;
            seed ^= seed << 5;
            out[i] = b * (2.0 * (double)seed / 4294967295.0 - 1.0);
        }
        return n;
    }

    case FCE_SIM_CHIRP:
    {
        double f0 = a, f1 = b;
        double k = (f1 - f0) / ((double)n / fs);
        for (i = 0; i < n; i++)
        {
            double t = (double)i / fs;
            out[i] = sin(2.0 * FCE_PI * (f0 * t + 0.5 * k * t * t));
        }
        return n;
    }

    case FCE_SIM_DC_PLUS_NOISE:
    {
        uint32_t seed = 12345u;
        for (i = 0; i < n; i++)
        {
            seed ^= seed << 13;
            seed ^= seed >> 17;
            seed ^= seed << 5;
            out[i] = a + b * (2.0 * (double)seed / 4294967295.0 - 1.0);
        }
        return n;
    }

    default:
        return 0u;
    }
}

#endif /* FCE_ENABLE_SIMULATION */
