/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
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
 */

#include "Mathlib.h"

double pnchisq(double x, double n, double lambda)
{
    double df, df1, val, lambda2, c, t, term;
    double igamma();
    double acc = 1.0e-12;

    n = floor(n + 0.5);
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n) || ISNAN(lambda))
        return x + n + lambda;
    if (!finite(n) || !finite(lambda))
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }
#endif
    if (n <= 0 || lambda < 0)
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }
    if (x <= 0)
        return 0;
    df = n;
    df1 = 0.5 * df;
    x = 0.5 * x;
    val = pgamma(x, df1, 1.0);
    lambda2 = 0.5 * lambda;
    c = 1;
    t = 0;
    do
    {
        t = t + 1;
        c = c * lambda2 / t;
        df1 = df1 + 1.0;
        term = c * pgamma(x, df1, 1.0);
        val = val + term;
    } while (term >= acc);
    return val * exp(-lambda2);
}
