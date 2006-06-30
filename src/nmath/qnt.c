/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2006 The R Development Core Team
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
 *  A copy of the GNU General Public License is available via WWW at
 *  http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
 *  writing to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "nmath.h"
#include "dpq.h"

double qnt(double p, double df, double delta, int lower_tail, int log_p)
{
    const static double accu = 1e-13;
    const static double Eps = 1e-11; /* must be > accu */

    double ux, lx, nx, pp;

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(df) || ISNAN(delta))
        return p + df + delta;
#endif
    if (!R_FINITE(df))
        ML_ERR_return_NAN;

    /* Was
     * df = floor(df + 0.5);
     * if (df < 1 || delta < 0) ML_ERR_return_NAN;
     */
    if (df < 0 || delta < 0)
        ML_ERR_return_NAN;

    R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);

    p = R_D_qIv(p);
    if (!lower_tail)
        p = 1 - p;

    /* Invert pnt(.) :
     * 1. finding an upper and lower bound */
    if (p > 1 - DBL_EPSILON)
        return ML_POSINF;
    pp = fmin2(1 - DBL_EPSILON, p * (1 + Eps));
    for (ux = delta; ux < DBL_MAX && pnt(ux, df, delta, TRUE, FALSE) < pp; ux *= 2)
        ;
    pp = p * (1 - Eps);
    for (lx = fmin2(-1, -delta); lx > -DBL_MAX && pnt(lx, df, delta, TRUE, FALSE) > pp; lx *= 2)
        ;

    /* 2. interval (lx,ux)  halving : */
    do
    {
        nx = 0.5 * (lx + ux);
        if (pnt(nx, df, delta, TRUE, FALSE) > p)
            ux = nx;
        else
            lx = nx;
    } while ((ux - lx) / fabs(nx) > accu);

    return 0.5 * (lx + ux);
}
