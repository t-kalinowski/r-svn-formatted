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
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  SYNOPSIS
 *
 *    #include "Mathlib.h"
 *    double qcauchy(double x, double location, double scale);
 *
 *  DESCRIPTION
 *
 *    The quantile function of the Cauchy distribution.
 */

#include "Mathlib.h"

double qcauchy(double x, double location, double scale)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(location) || ISNAN(scale))
        return x + location + scale;
    if (!FINITE(x) || !FINITE(location) || !FINITE(scale))
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }
#endif

    if (scale <= 0)
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }
    return location + scale * tan(M_PI * (x - 0.5));
}
