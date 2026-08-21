#include "fce_internal.h"
#include <stdio.h>

/* extern prototypes from fce_iir.c */
fce_status_t proto_ellip(uint32_t n, double rp, double rs,
                         fce_cplx_t* z, uint32_t* nz,
                         fce_cplx_t* p, uint32_t* np, double* k);
fce_status_t proto_bessel(uint32_t n, fce_cplx_t* p, double* k);

int main(void)
{
    /* elliptic K */
    printf("ellipk(0)    = %.17g (expect 1.5707963267948966)\n", fce_ellipk(0.0));
    printf("ellipk(0.5)  = %.17g (expect 1.854074677301372)\n", fce_ellipk(0.5));
    printf("ellipk(0.99) = %.17g\n", fce_ellipk(0.99));

    /* ellipj */
    {
        double sn, cn, dn;
        fce_ellipj(0.5, 0.5, &sn, &cn, &dn);
        printf("ellipj(0.5,0.5) sn=%.17g cn=%.17g dn=%.17g\n", sn, cn, dn);
        fce_ellipj(1.0, 0.25, &sn, &cn, &dn);
        printf("ellipj(1,0.25)  sn=%.17g cn=%.17g dn=%.17g\n", sn, cn, dn);
    }

    /* arc_jac_sn */
    {
        fce_cplx_t w = fce_arc_jac_sn(fce_cx(0.0, 2.0), 0.5);
        printf("asn(2j,0.5) = %.17g%+.17gj\n", w.re, w.im);
        w = fce_arc_jac_sn(fce_cx(0.5, 0.0), 0.25);
        printf("asn(0.5,0.25) = %.17g%+.17gj\n", w.re, w.im);
    }

    /* ellipdeg */
    printf("ellipdeg(5, 0.0458) = %.17g\n", fce_ellipdeg(5, 0.0458));

    /* proto ellip order 5, rp=0.5, rs=60 */
    {
        fce_cplx_t z[16], p[16];
        uint32_t nz, np;
        double k;
        uint32_t i;
        fce_status_t st = proto_ellip(5, 0.5, 60.0, z, &nz, p, &np, &k);
        printf("proto_ellip(5,0.5,60) st=%d nz=%u np=%u k=%.17g\n",
               (int)st, nz, np, k);
        for (i = 0; i < np; i++)
            printf("  p[%u] = %.17g %+.17gj  |p|=%.17g\n", i, p[i].re, p[i].im,
                   fce_cx_abs(p[i]));
        for (i = 0; i < nz; i++)
            printf("  z[%u] = %.17g %+.17gj\n", i, z[i].re, z[i].im);
    }

    /* proto bessel order 4 */
    {
        fce_cplx_t p[16];
        double k;
        uint32_t i;
        fce_status_t st = proto_bessel(4, p, &k);
        printf("proto_bessel(4) st=%d k=%.17g\n", (int)st, k);
        for (i = 0; i < 4; i++)
            printf("  p[%u] = %.17g %+.17gj  |p|=%.17g\n", i, p[i].re, p[i].im,
                   fce_cx_abs(p[i]));
    }
    return 0;
}
