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
 *	The distribution function of the Cauchy distribution.
 */

#include "nmath.h"
#include "dpq.h"

double pcauchy(double x, double location, double scale, int lower_tail, int log_p)
{
    const double x_big = pow(5 * DBL_EPSILON, -0.25); /* = 5478.3 for IEEE */

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(location) || ISNAN(scale))
        return x + location + scale;
#endif
    if (scale <= 0)
        ML_ERR_return_NAN;

    x = (x - location) / scale;
    if (ISNAN(x))
        ML_ERR_return_NAN;
#ifdef IEEE_754
    if (!R_FINITE(x))
    {
        if (x < 0)
            return R_DT_0;
        else
            return R_DT_1;
    }
#endif
    /* for large (negative || "upper tail & positive")  x,
       the standard formula suffers from cancellation */
    if (x < 0 && lower_tail && (-x) > x_big)
    {
        /* P = 1/(pi*(-x)) *(1 - 1/(3* x^2)) */
        if (log_p)
            return -log(-M_PI * x) + log1p(-1. / (3 * x * x));
        /* else */
        return (1 - 1. / (3 * x * x)) / (-M_PI * x);
    }
    /* else */
    if (x > 0 && !lower_tail && x > x_big)
    {
        /* P = 1 - 1/(pi*x) *(1 - 1/(3* x^2))  */
        if (log_p)
            return -log(M_PI * x) + log1p(-1. / (3 * x * x));
        /* else */
        return (1 - 1. / (3 * x * x)) / (M_PI * x);
    }
    /* else */
    return R_DT_val(0.5 + atan(x) / M_PI);
}
