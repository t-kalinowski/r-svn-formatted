/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka and the R Development Core Team.
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
 *    double rgeom(double p);
 *
 *  DESCRIPTION
 *
 *    Random variates from the geometric distribution.
 *
 *  NOTES
 *
 *    We generate lambda as exponential with scale parameter
 *    p / (1 - p).  Return a Poisson deviate with mean lambda.
 *
 *  REFERENCE
 *
 *    Devroye, L. (1986).
 *    Non-Uniform Random Variate Generation.
 *    New York: Springer-Verlag.
 *    Page 480.
 */

#include "Mathlib.h"

double rgeom(double p)
{
    if (
#ifdef IEEE_754
        isnan(p) ||
#endif
        p <= 0 || p >= 1)
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }
    return rpois(sexp() * ((1 - p) / p));
}
