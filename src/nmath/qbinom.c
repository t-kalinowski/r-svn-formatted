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
 *	The quantile function of the binomial distribution.
 *
 *  METHOD
 *
 *	Uses the Cornish-Fisher Expansion to include a skewness
 *	correction to a normal approximation.  This gives an
 *	initial value which never seems to be off by more than
 *	1 or 2.	 A search is then conducted of values close to
 *	this initial start point.
 */
#include "nmath.h"
#include "dpq.h"

double qbinom(double p, double n, double pr, int lower_tail, int log_p)
{
    double q, mu, sigma, gamma, z, y;

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(n) || ISNAN(pr))
        return p + n + pr;
#endif
    if (!R_FINITE(p) || !R_FINITE(n) || !R_FINITE(pr))
        ML_ERR_return_NAN;
    R_Q_P01_check(p);

    n = floor(n + 0.5);
    if (pr <= 0 || pr >= 1 || n <= 0)
        ML_ERR_return_NAN;

    if (p == R_DT_0)
        return 0.;
    if (p == R_DT_1)
        return n;

    q = 1 - pr;
    mu = n * pr;
    sigma = sqrt(n * pr * q);
    gamma = (q - pr) / sigma;

#ifdef DEBUG_qbinom
    REprintf("qbinom(p=%7g, n=%g, pr=%7g, l.t.=%d, log=%d): sigm=%g, gam=%g\n", p, n, pr, lower_tail, log_p, sigma,
             gamma);
#endif
    /* Note : "same" code in qpois.c, qbinom.c, qnbinom.c --
     * FIXME: This is far from optimal [cancellation for p ~= 1, etc]: */
    if (!lower_tail || log_p)
    {
        p = R_DT_qIv(p); /* need check again (cancellation!): */
        if (p == 0.)
            return 0.;
        if (p == 1.)
            return n;
    }
    /* temporary hack --- FIXME --- */
    if (p + 1.01 * DBL_EPSILON >= 1.)
        return n;

    /* y := approx.value (Cornish-Fisher expansion) :  */
    z = qnorm(p, 0., 1., /*lower_tail*/ LTRUE, /*log_p*/ LFALSE);
    y = floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5);
    if (y > n) /* way off */
        y = n;

#ifdef DEBUG_qbinom
    REprintf("  new (p,1-p)=(%7g,%7g), z=qnorm(..)=%7g, y=%5g\n", p, 1 - p, z, y);
#endif
    z = pbinom(y, n, pr, /*lower_tail*/ LTRUE, /*log_p*/ LFALSE);

    /* fuzz to ensure left continuity: */
    p *= 1 - 64 * DBL_EPSILON;

    if (y > n)
    {
        y = n;
        z = pbinom(y, n, pr, /*lower_tail*/ LTRUE, /*log_p*/ LFALSE);
    }
    else if (y < 0)
    {
        y = 0;
        z = pbinom(y, n, pr, /*lower_tail*/ LTRUE, /*log_p*/ LFALSE);
    }
#ifdef maybe_future
    if ((lower_tail && z >= p) || (!lower_tail && z <= p))
    {
#else
    if (z >= p)
    {
#endif
        /* search to the left */
#ifdef DEBUG_qbinom
        REprintf("\tnew z=%7g >= p = %7g  --> search to left (y--) ..\n", z, p);
#endif
        for (;;)
        {
            if (y == 0 || (z = pbinom(y - 1, n, pr, /*l._t.*/ LTRUE, /*log_p*/ LFALSE)) < p)
                return y;
            y = y - 1;
        }
    }
    else
    { /* search to the right */
#ifdef DEBUG_qbinom
        REprintf("\tnew z=%7g < p = %7g  --> search to right (y++) ..\n", z, p);
#endif
        for (;;)
        {
            y = y + 1;
            if (y == n || (z = pbinom(y, n, pr, /*l._t.*/ LTRUE, /*log_p*/ LFALSE)) >= p)
                return y;
        }
    }
}
