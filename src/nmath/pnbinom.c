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
 *	The distribution function of the negative binomial distribution.
 *
 *  NOTES
 *
 *	x = the number of failures before the n-th success
 */

#include "Mathlib.h"

double pnbinom(double x, double n, double p, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n) || ISNAN(p))
        return x + n + p;
    if (!R_FINITE(n) || !R_FINITE(p))
        ML_ERR_return_NAN;
#endif
    if (n <= 0 || p <= 0 || p >= 1)
        ML_ERR_return_NAN;

    x = floor(x + 1e-7);
    if (x < 0)
        return R_DT_0;
#ifdef IEEE_754
    if (!R_FINITE(x))
        return R_DT_1;
#endif
    return pbeta(p, n, x + 1, lower_tail, log_p);
}
