/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998--1999  The R Core Team
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
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *

 * Reference:
 * Cran, G. W., K. J. Martin and G. E. Thomas (1977).
 *	Remark AS R19 and Algorithm AS 109,
 *	Applied Statistics, 26(1), 111-114.
 * Remark AS R83 (v.39, 309-310) and the correction (v.40(1) p.236)
 *	have been incorporated in this version.
 */

#include "Mathlib.h"

#define zero 0.0

/* set the exponent of accu to -2r-2 for r digits of accuracy */
#ifdef OLD
#define acu 1.0e-32
#define lower 0.0001
#define upper 0.9999

#else /*---- NEW ---- -- still fails for p = 1e11, q=.5*/

#define fpu 3e-308
/* acu_min:  Minimal value for accuracy 'acu' which will depend on (a,p);
         acu_min >= fpu ! */
#define acu_min 1e-300
#define lower fpu
#define upper 1 - 2.22e-16

#endif

#define const1 2.30753
#define const2 0.27061
#define const3 0.99229
#define const4 0.04481

static volatile double xtrunc;

double qbeta(double alpha, double p, double q)
{
    int swap_tail, i_pb, i_inn;
    double a, adj, logbeta, g, h, pp, prev, qq, r, s, t, tx, w, y, yprev;
    double acu;
    volatile double xinbta;

    /* define accuracy and initialize */

    xinbta = alpha;

    /* test for admissibility of parameters */

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(q) || ISNAN(alpha))
        return p + q + alpha;
#endif
    if (p < zero || q < zero || alpha < zero || alpha > 1)
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }
    if (alpha == zero || alpha == 1)
        return alpha;

    logbeta = lbeta(p, q);

    /* change tail if necessary;  afterwards   0 < a <= 1/2	 */
    if (alpha <= 0.5)
    {
        a = alpha;
        pp = p;
        qq = q;
        swap_tail = 0;
    }
    else
    { /* change tail, swap  p <-> q :*/
        a = 1 - alpha;
        pp = q;
        qq = p;
        swap_tail = 1;
    }

    /* calculate the initial approximation */

    r = sqrt(-log(a * a));
    y = r - (const1 + const2 * r) / (1 + (const3 + const4 * r) * r);
    if (pp > 1 && qq > 1)
    {
        r = (y * y - 3) / 6;
        s = 1 / (pp + pp - 1);
        t = 1 / (qq + qq - 1);
        h = 2 / (s + t);
        w = y * sqrt(h + r) / h - (t - s) * (r + 5 / 6 - 2 / (3 * h));
        xinbta = pp / (pp + qq * exp(w + w));
    }
    else
    {
        r = qq + qq;
        t = 1 / (9 * qq);
        t = r * pow(1 - t + y * sqrt(t), 3);
        if (t <= zero)
            xinbta = 1 - exp((log((1 - a) * qq) + logbeta) / qq);
        else
        {
            t = (4 * pp + r - 2) / t;
            if (t <= 1)
                xinbta = exp((log(a * pp) + logbeta) / pp);
            else
                xinbta = 1 - 2 / (t + 1);
        }
    }

    /* solve for x by a modified newton-raphson method, */
    /* using the function pbeta_raw */

    r = 1 - pp;
    t = 1 - qq;
    yprev = zero;
    adj = 1;
    if (xinbta < lower)
        xinbta = lower;
    else if (xinbta > upper)
        xinbta = upper;

    /* Desired accuracy should depend on  (a,p)
     * This is from Remark .. on AS 109, adapted.
     * However, it's not clear if this is "optimal" for IEEE double prec.

     * acu = fmax2(acu_min, pow(10., -25. - 5./(pp * pp) - 1./(a * a)));

     * NEW: 'acu' accuracy NOT for squared adjustment, but simple;
     * ---- i.e.,  "new acu" = sqrt(old acu)

     */
    acu = fmax2(acu_min, pow(10., -13 - 2.5 / (pp * pp) - 0.5 / (a * a)));
    tx = prev = zero; /* keep -Wall happy */

    for (i_pb = 0; i_pb < 1000; i_pb++)
    {
        y = pbeta_raw(xinbta, pp, qq);
        /* y = pbeta_raw2(xinbta, pp, qq, logbeta) -- to SAVE CPU; */
#ifdef IEEE_754
        if (!FINITE(y))
#else
        if (errno)
#endif
        {
            ML_ERROR(ME_DOMAIN);
            return ML_NAN;
        }
        y = (y - a) * exp(logbeta + r * log(xinbta) + t * log(1 - xinbta));
        if (y * yprev <= zero)
            prev = fmax2(fabs(adj), fpu);
        g = 1;
        for (i_inn = 0; i_inn < 1000; i_inn++)
        {
            adj = g * y;
            if (fabs(adj) < prev)
            {
                tx = xinbta - adj; /* trial new x */
                if (tx >= zero && tx <= 1)
                {
                    if (prev <= acu)
                        goto L_converged;
                    if (fabs(y) <= acu)
                        goto L_converged;
                    if (tx != zero && tx != 1)
                        break;
                }
            }
            g /= 3;
        }
        xtrunc = tx; /* this prevents trouble with excess FPU */
                     /* precision on some machines. */
        if (xtrunc == xinbta)
            goto L_converged;
        xinbta = tx;
        yprev = y;
    }
    /*-- NOT converged: Iteration count --*/
    ML_ERROR(ME_PRECISION);

L_converged:
    if (swap_tail)
        xinbta = 1 - xinbta;
    return xinbta;
}
