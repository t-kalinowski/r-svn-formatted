/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 1999 R Core Team
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
 *  SYNOPSIS
 *
 *    #include "Mathlib.h"
 *    double phyper(double x, double NR, double NB, double n);
 *
 *  DESCRIPTION
 *
 *    The distribution function of the hypergeometric distribution.
 */

#include "Mathlib.h"

double phyper(double x, double NR, double NB, double n)
{
    /* Sample of  n balls from  NR red  and  NB black ones;  x are red */

    /* basically the same code is used also in  ./qhyper.c -- keep in sync! */
    double N, xstart, xend, xr, xb, sum, term;
    int small_N;
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(NR) || ISNAN(NB) || ISNAN(n))
        return x + NR + NB + n;
    if (!finite(x) || !finite(NR) || !finite(NB) || !finite(n))
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }
#endif

    x = floor(x);
    NR = floor(NR + 0.5);
    NB = floor(NB + 0.5);
    N = NR + NB;
    n = floor(n + 0.5);
    if (NR < 0 || NB < 0 || n < 0 || n > N)
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }
    small_N = (N < 1000); /* won't have underflow in product below */
    xstart = fmax2(0, n - NB);
    xend = fmin2(n, NR);
    if (x < xstart)
        return 0.0;
    if (x >= xend)
        return 1.0;
    xr = xstart;
    xb = n - xr;
    /* if N is small,  term := product.ratio( bin.coef );
       otherwise work with its logarithm to protect against underflow */
    term = lfastchoose(NR, xr) + lfastchoose(NB, xb) - lfastchoose(N, n);
    if (small_N)
        term = exp(term);
    NR -= xr;
    NB -= xb;
    sum = 0.0;
    while (xr <= x)
    {
        sum += (small_N ? term : exp(term));
        xr++;
        NB++;
        if (small_N)
            term *= (NR / xr) * (xb / NB);
        else
            term += log((NR / xr) * (xb / NB));
        xb--;
        NR--;
    }
    return sum;
}
