/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka and the R Development Core Team
 *  Copyright (C) 2000-2001 The R Development Core Team
 *  based on AS243 (C) 1989 Royal Statistical Society
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 */

/*  Algorithm AS 243  Lenth,R.V. (1989). Appl. Statist., Vol.38, 185-189.
 *  ----------------
 *  Cumulative probability at t of the non-central t-distribution
 *  with df degrees of freedom (may be fractional) and non-centrality
 *  parameter delta.
 *
 *  NOTE
 *
 *    Requires the following auxiliary routines:
 *
 *	lgammafn(x)	- log gamma function
 *	pbeta(x, a, b)	- incomplete beta function
 *	pnorm(x)	- normal distribution function
 *
 *  CONSTANTS
 *
 *    M_SQRT_2dPI  = 1/ {gamma(1.5) * sqrt(2)} = sqrt(2 / pi)
 *    M_LN_SQRT_PI = ln(sqrt(pi)) = ln(pi)/2
 */

#include "nmath.h"
#include "dpq.h"

/*----------- DEBUGGING -------------
 *
 *	make CFLAGS='-DDEBUG_pnt -g -I/usr/local/include -I../include'

 * -- Feb.3, 1999; M.Maechler:
    - For 't > delta > 20' (or so)	the result is completely WRONG!
 */

double pnt(double t, double df, double delta, int lower_tail, int log_p)
{
    double a, albeta, b, del, errbd, geven, godd, lambda, p, q, rxb, s, tnc, tt, x, xeven, xodd;
    int it, negdel;

    /* note - itrmax and errmax may be changed to suit one's needs. */

    const int itrmax = 1000;
    const static double errmax = 1.e-12;

    if (df <= 0.)
        ML_ERR_return_NAN;

    if (!R_FINITE(t))
        return (t < 0) ? R_DT_0 : R_DT_1;
    if (t >= 0.)
    {
        negdel = FALSE;
        tt = t;
        del = delta;
    }
    else
    {
        negdel = TRUE;
        tt = -t;
        del = -delta;
    }

    if (df > 4e5 || del * del > 2 * M_LN2 * (-(DBL_MIN_EXP)))
    {
        /*-- 2nd part: if del > 37.62, then p=0 below
          FIXME: test should depend on `df', `tt' AND `del' ! */
        /* Approx. from	 Abramowitz & Stegun 26.7.10 (p.949) */
        s = 1. / (4. * df);

        return pnorm(tt * (1. - s), del, sqrt(1. + tt * tt * 2. * s), lower_tail != negdel, log_p);
    }

    /* initialize twin series */
    /* Guenther, J. (1978). Statist. Computn. Simuln. vol.6, 199. */

    x = t * t;
    x = x / (x + df); /* in [0,1) */
#ifdef DEBUG_pnt
    REprintf("pnt(t=%7g, df=%7g, delta=%7g) ==> x= %10g:", t, df, delta, x);
#endif
    if (x > 0.)
    { /* <==>  t != 0 */
        lambda = del * del;
        p = .5 * exp(-.5 * lambda);
#ifdef DEBUG_pnt
        REprintf("\t p=%10g\n", p);
#endif
        if (p == 0.)
        { /* underflow! */

            /*========== really use an other algorithm for this case !!! */
            ML_ERROR(ME_UNDERFLOW);
            ML_ERROR(ME_RANGE); /* |delta| too large */
            return R_DT_0;
        }
#ifdef DEBUG_pnt
        REprintf("it  1e5*(godd,  geven)       p	 q	    s	 "
                 /* 1.3 1..4..7.9 1..4..7.9  1..4..7.9 1..4..7.9 1..4..7.9_ */
                 "	  pnt(*)      errbd\n");
        /* 1..4..7..0..3..6 1..4..7.9*/
#endif
        q = M_SQRT_2dPI * p * del;
        s = .5 - p;
        a = .5;
        b = .5 * df;
        rxb = pow(1. - x, b);
        albeta = M_LN_SQRT_PI + lgammafn(b) - lgammafn(.5 + b);
        xodd = pbeta(x, a, b, /*lower*/ TRUE, /*log_p*/ FALSE);
        godd = 2. * rxb * exp(a * log(x) - albeta);
        xeven = 1. - rxb;
        geven = b * x * rxb;
        tnc = p * xodd + q * xeven;

        /* repeat until convergence or iteration limit */
        for (it = 1; it <= itrmax; it++)
        {
            a += 1.;
            xodd -= godd;
            xeven -= geven;
            godd *= x * (a + b - 1.) / a;
            geven *= x * (a + b - .5) / (a + .5);
            p *= lambda / (2 * it);
            q *= lambda / (2 * it + 1);
            tnc += p * xodd + q * xeven;
            s -= p;
            if (s <= 0.)
            { /* happens e.g. for (t,df,delta)=(40,10,38.5), after 799 it.*/
                ML_ERROR(ME_PRECISION);
#ifdef DEBUG_pnt
                REprintf("s = %#14.7g < 0 !!! ---> non-convergence!!\n", s);
#endif
                goto finis;
            }
            errbd = 2. * s * (xodd - godd);
#ifdef DEBUG_pnt
            REprintf("%3d %#9.4g %#9.4g	 %#9.4g %#9.4g %#9.4g %#14.10g %#9.4g\n", it, 1e5 * godd, 1e5 * geven, p, q, s,
                     tnc, errbd);
#endif
            if (errbd < errmax)
                goto finis; /*convergence*/
        }
        /* non-convergence:*/
        ML_ERROR(ME_PRECISION);
    }
    else
    { /* x = t = 0 */
        tnc = 0.;
    }
finis:
    tnc += pnorm(-del, 0., 1., /*lower*/ TRUE, /*log_p*/ FALSE);

    lower_tail = lower_tail != negdel; /* xor */
    return R_DT_val(tnc);
}
