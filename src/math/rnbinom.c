/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
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
 */

/*
 *	Reference:
 *
 *	Devroye, L. (1980).
 *	Non-Uniform Random Variate Generation.
 *	New York:Springer-Verlag. Page 480.
 *
 *	Method:
 *
 *	Generate lambda as gamma with shape parameter n and scale
 *	parameter p/(1-p).  Return a Poisson deviate with mean lambda.
 */

#include "Mathlib.h"

double rnbinom(double n, double p)
{
    n = floor(n + 0.5);
    if (n <= 0.0 || p <= 0.0 || p >= 1.0)
        DOMAIN_ERROR;
    return rpois(rgamma(n, (1.0 - p) / p));
}
