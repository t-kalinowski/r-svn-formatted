/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
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
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double choose(double n, double k);
 *    double lchoose(double n, double k);
 *    (and private)
 *    double fastchoose(double n, double k);
 *    double lfastchoose(double n, double k);
 *
 *  DESCRIPTION
 *
 *    Binomial coefficients.
 */
#include "nmath.h"

double lfastchoose(double n, double k)
{
    return -log(n + 1.) - lbeta(n - k + 1., k + 1.);
    /* the same (but less stable):
     * == lgammafn(n + 1.0) - lgammafn(k + 1.0) - lgammafn(n - k + 1.0); */
}

double fastchoose(double n, double k)
{
    return exp(lfastchoose(n, k));
}

double lchoose(double n, double k)
{
    n = floor(n + 0.5);
    k = floor(k + 0.5);
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(n) || ISNAN(k))
        return n + k;
#endif
    if (n < 0)
        ML_ERR_return_NAN;
    if (k < 0 || n < k)
        return ML_NEGINF;

    return lfastchoose(n, k);
}

double choose(double n, double k)
{
    n = floor(n + 0.5);
    k = floor(k + 0.5);
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(n) || ISNAN(k))
        return n + k;
#endif
    if (n < 0)
        ML_ERR_return_NAN; /* could be defined instead as
              (-1)^k (-n + k - 1 \\ k) {k>=0}*/
    if (k < 0 || n < k)
        return 0.;

    return floor(exp(lfastchoose(n, k)) + 0.5);
}
