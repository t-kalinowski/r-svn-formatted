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
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *  SYNOPSIS
 *
 *    #include "Mathlib.h"
 *    double pf(double x, double n1, double n2);
 *
 *  DESCRIPTION
 *
 *    The distribution function of the F distribution.
 */

#include "Mathlib.h"

double pf(double x, double n1, double n2)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n1) || ISNAN(n2))
        return x + n2 + n1;
#endif
    if (n1 <= 0.0 || n2 <= 0.0)
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }
    if (x <= 0.0)
        return 0.0;
    x = 1.0 - pbeta(n2 / (n2 + n1 * x), n2 / 2.0, n1 / 2.0);
    return ML_VALID(x) ? x : ML_NAN;
}
