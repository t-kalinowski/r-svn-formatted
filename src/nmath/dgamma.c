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
 *    double dgamma(double x, double shape, double scale);
 *
 *  DESCRIPTION
 *
 *    Computes the density of the gamma distribution.
 */

#include "Mathlib.h"

double dgamma(double x, double shape, double scale)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(shape) || ISNAN(scale))
        return x + shape + scale;
#endif
    if (shape <= 0 || scale <= 0)
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }
    if (x < 0)
        return 0;
    if (x == 0)
    {
        if (shape < 1)
        {
            ML_ERROR(ME_RANGE);
            return ML_POSINF;
        }
        if (shape > 1)
        {
            return 0;
        }
        return 1 / scale;
    }
    x = x / scale;
    return exp((shape - 1) * log(x) - lgamma(shape) - x) / scale;
}
