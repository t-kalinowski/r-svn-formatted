/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
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
 */

#include "Mathlib.h"

double qlogis(double p, double location, double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(location) || ISNAN(scale))
        return p + location + scale;
#endif
    R_Q_P01_check(p);
    if (scale < 0.)
        ML_ERR_return_NAN;
    if (scale == 0.)
        return location;

    if (R_DT_qIv(p) <= 0)
        return ML_NEGINF;
    if (p == R_DT_1)
        return ML_POSINF;

    /* p := logit(p) = log( p / (1-p) )	 : */
    if (log_p)
    {
        if (lower_tail)
            p = p - logrelerr(-exp(p));
        else
            p = logrelerr(-exp(p)) - p;
    }
    else
        p = log(lower_tail ? (p / (1. - p)) : ((1. - p) / p));

    return location + scale * p;
}
