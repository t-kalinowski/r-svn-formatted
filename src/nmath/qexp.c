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
 *    The quantile function of the exponential distribution.
 *
 */

#include "Mathlib.h"

double qexp(double p, double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(scale))
        return p + scale;
#endif
    R_Q_P01_check(p);
    if (scale < 0)
        ML_ERR_return_NAN;

    if (p == R_DT_0)
        return 0.;

    return -scale * R_DT_Clog(p);
}
