/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000 The R Development Core Team
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
 *
 *  DESCRIPTION
 *
 *	The "Student" t distribution quantile function.
 *
 *  NOTES
 *
 *	This is a C translation of the Fortran routine given in:
 *	Algorithm 396: Student's t-quantiles by
 *	G.W. Hill CACM 13(10), 619-620, October 1970
 */

#include "Mathlib.h"
#include "dpq.h"

const static double eps = 1.e-12;

double qt(double p, double ndf, int lower_tail, int log_p)
{
    double a, b, c, d, prob, P, q, x, y;
    int neg;

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(ndf))
        return p + ndf;
#endif
    if (p == R_DT_0)
        return ML_NEGINF;
    if (p == R_DT_1)
        return ML_POSINF;
    R_Q_P01_check(p);
    if (ndf < 1)
        ML_ERR_return_NAN;

    /* FIXME: This test should depend on  ndf  AND p  !!
     * -----  and in fact should be replaced by
     * something like Abramowitz & Stegun 26.7.5 (p.949)
     */
    if (ndf > 1e20)
        return qnorm(p, 0., 1., lower_tail, log_p);

    if (log_p)
    { /* FIXME: *can* do better in the case where P = 2*p;
       *	  using	 qnorm(p, ... , log_p) .. */
        p = exp(p);
    }

    if ((lower_tail && p > 0.5) || (!lower_tail && p < 0.5))
    {
        neg = 0;
        P = 2 * R_D_Cval(p);
    }
    else
    {
        neg = 1;
        P = 2 * R_D_Lval(p);
    }

    if (fabs(ndf - 2) < eps)
    { /* df ~= 2 */
        q = sqrt(2 / (P * (2 - P)) - 2);
    }
    else if (ndf < 1 + eps)
    { /* df ~= 1 */
        prob = P * M_PI_2;
        q = cos(prob) / sin(prob);
    }
    else
    { /*-- usual case;  including, e.g.,  df = 1.1 */
        a = 1 / (ndf - 0.5);
        b = 48 / (a * a);
        c = ((20700 * a / b - 98) * a - 16) * a + 96.36;
        d = ((94.5 / (b + c) - 3) / b + 1) * sqrt(a * M_PI_2) * ndf;
        y = pow(d * P, 2 / ndf);

        if (y > 0.05 + a)
        {
            /* Asymptotic inverse expansion about normal */
            x = qnorm(0.5 * P, 0.0, 1.0, /*lower_tail*/ LTRUE, /*log_p*/ LFALSE);
            y = x * x;
            if (ndf < 5)
                c += 0.3 * (ndf - 4.5) * (x + 0.6);
            c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c;
            y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1) * x;
            y = a * y * y;
            if (y > 0.002) /* FIXME: This cutoff is machine-precision dependent*/
                y = exp(y) - 1;
            else
            { /* Taylor of	 e^y -1 : */
                y = (0.5 * y + 1) * y;
            }
        }
        else
        {
            y = ((1 / (((ndf + 6) / (ndf * y) - 0.089 * d - 0.822) * (ndf + 2) * 3) + 0.5 / (ndf + 4)) * y - 1) *
                    (ndf + 1) / (ndf + 2) +
                1 / y;
        }
        q = sqrt(ndf * y);
    }
    if (neg)
        q = -q;
    return q;
}
