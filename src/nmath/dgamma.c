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
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  DESCRIPTION
 *
 *     Computes the density of the gamma distribution.
 */

#include "nmath.h"
#include "dpq.h"

double dgamma(double x, double shape, double scale, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(shape) || ISNAN(scale))
        return x + shape + scale;
#endif
    if (shape <= 0 || scale <= 0)
        ML_ERR_return_NAN;

    if (x < 0)
        return R_D__0;
    if (x == 0)
    {
        if (shape < 1)
            ML_ERR_return_NAN;
        if (shape > 1)
            return R_D__0;

        return give_log ? -log(scale) : 1 / scale;
    }
    x /= scale;
    return give_log ? ((shape - 1) * log(x) - lgammafn(shape) - x) - log(scale)
                    : exp((shape - 1) * log(x) - lgammafn(shape) - x) / scale;
}
