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
 *	The quantile function of the Poisson distribution.
 *
 *  METHOD
 *
 *	Uses the Cornish-Fisher Expansion to include a skewness
 *	correction to a normal approximation.  This gives an
 *	initial value which never seems to be off by more than
 *	1 or 2.	 A search is then conducted of values close to
 *	this initial start point.
 */

#include "Mathlib.h"
#include "dpq.h"

double qpois(double p, double lambda, int lower_tail, int log_p)
{
    double mu, sigma, gamma, z, y;
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(lambda))
        return p + lambda;
#endif
    if (!R_FINITE(lambda))
        ML_ERR_return_NAN;
    R_Q_P01_check(p);
    if (lambda <= 0)
        ML_ERR_return_NAN;

    if (p == R_DT_0)
        return 0;
#ifdef IEEE_754
    if (p == R_DT_1)
        return ML_POSINF;
#endif
    mu = lambda;
    sigma = sqrt(lambda);
    gamma = sigma;

    /* FIXME: This is far from optimal :
       -- "same" code in qpois.c, qbinom.c, qnbinom.c */
    if (!lower_tail || log_p)
        p = R_DT_qIv(p);

    z = qnorm(p, 0., 1., /*lower_tail*/ LTRUE, /*log_p*/ LFALSE);
    y = floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5);

    z = ppois(y, lambda, /*lower_tail*/ LTRUE, /*log_p*/ LFALSE);

    /* fuzz to ensure left continuity; 1 - 1e-7 may lose too much : */
    p *= 1 - 64 * DBL_EPSILON;
#ifdef maybe_future
    if ((lower_tail && z >= p) || (!lower_tail && z <= p))
    {
#else
    if (z >= p)
    {
#endif
        /* search to the left */
        for (;;)
        {
            if (y == 0 || (z = ppois(y - 1, lambda, /*l._t.*/ LTRUE, /*log_p*/ LFALSE)) < p)
                return y;
            y = y - 1;
        }
    }
    else
    { /* search to the right */
        for (;;)
        {
            y = y + 1;
            if ((z = ppois(y, lambda, /*l._t.*/ LTRUE, /*log_p*/ LFALSE)) >= p)
                return y;
        }
    }
}
