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
 *  SYNOPSIS
 *
 *    #include "Mathlib.h"
 *    double rnorm(double mu, double sigma);
 *
 *  DESCRIPTION
 *
 *    Random variates from the normal distribution.
 *
 */

#include "Mathlib.h"

double rnorm(double mu, double sigma)
{
    if (
#ifdef IEEE_754
        !R_FINITE(mu) || !R_FINITE(sigma) ||
#endif
        sigma < 0.)
        ML_ERR_return_NAN;

    if (sigma == 0.)
        return mu;
    else
        return mu + sigma * snorm();
}
