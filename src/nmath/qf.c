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
 *    The quantile function of the F distribution.
 */

#include "nmath.h"
#include "dpq.h"

double qf(double p, double n1, double n2, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(n1) || ISNAN(n2))
        return p + n1 + n2;
#endif
    if (n1 <= 0. || n2 <= 0.)
        ML_ERR_return_NAN;

    R_Q_P01_check(p);
    if (p == R_DT_0)
        return 0;

    /* fudge the extreme DF cases -- qbeta doesn't do this well */

    if (n2 > 4e5)
        return qchisq(p, n1, lower_tail, log_p) / n1;

    if (n1 > 4e5)
        return 1 / qchisq(p, n2, !lower_tail, log_p) * n2;

    p = (1. / qbeta(R_DT_CIv(p), n2 / 2, n1 / 2, LTRUE, LFALSE) - 1.) * (n2 / n1);
    return ML_VALID(p) ? p : ML_NAN;
}
