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
 *    double runif(double a, double b);
 *
 *  DESCRIPTION
 *
 *    Random variates from the uniform distribution.
 */

#include "Mathlib.h"

double runif(double a, double b)
{
    if (
#ifdef IEEE_754
        !R_FINITE(a) || !R_FINITE(b) ||
#endif
        b < a)
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }
    if (a == b)
        return a;
    else
        return a + (b - a) * sunif();
}
