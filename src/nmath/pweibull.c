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
 *    The distribution function of the Weibull distribution.
 */

#include "Mathlib.h"
#include "dpq.h"

double pweibull(double x, double shape, double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(shape) || ISNAN(scale))
        return x + shape + scale;
#endif
    if (shape <= 0 || scale <= 0)
        ML_ERR_return_NAN;

    if (x <= 0)
        return R_DT_0;
    if (log_p && !lower_tail)
        return -pow(x / scale, shape);
    return R_DT_Cval(exp(-pow(x / scale, shape)));
}
