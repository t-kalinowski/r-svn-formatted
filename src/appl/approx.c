/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995-2001   Robert Gentleman, Ross Ihaka and the
 *			      R Development Core Team
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
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "R_ext/Arith.h"
#include "R_ext/Error.h"
#include "R_ext/Applic.h"

/* Linear and Step Function Interpolation */

/* Assumes that ordinates are in ascending order
 * The right interval is found by bisection
 * Linear/constant interpolation then takes place on that interval
 */

static double ylow;
static double yhigh;
static double f1;
static double f2;

static double approx1(double v, double *x, double *y, int n, int method)
{
    /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */
    int i, j, ij;

    i = 0;
    j = n - 1;

    /* handle out-of-domain points */

    if (v < x[i])
        return ylow;
    if (v > x[j])
        return yhigh;

    /* find the correct interval by bisection */

    while (i < j - 1)
    {
        ij = (i + j) / 2;
        if (v < x[ij])
            j = ij;
        else
            i = ij;
    }

    /* interpolation */

    if (v == x[i])
        return y[i];
    if (v == x[j])
        return y[j];
    if (x[j] == x[i])
        return y[i];

    if (method == 1)
    { /* linear */
        return y[i] + (y[j] - y[i]) * ((v - x[i]) / (x[j] - x[i]));
    }
    else
    { /* 2 : constant */
        return y[i] * f1 + y[j] * f2;
    }
}

/* R Frontend for Linear and Constant Interpolation */

void R_approx(double *x, double *y, int *nxy, double *xout, int *nout, int *method, double *yleft, double *yright,
              double *f)
{
    int i;

    /* check interpolation method */

    switch (*method)
    {
    case 1: /* linear */
        break;
    case 2: /* constant */
        if (!R_FINITE(*f) || *f < 0 || *f > 1)
            error("approx(): invalid f value");
        f2 = *f;
        f1 = 1 - *f;
        break;
    default:
        error("approx(): invalid interpolation method");
        break;
    }

    for (i = 0; i < *nxy; i++)
        if (ISNA(x[i]) || ISNA(y[i]))
            error("approx(): attempted to interpolate NA values");

    ylow = *yleft;
    yhigh = *yright;

    for (i = 0; i < *nout; i++)
        if (!ISNA(xout[i]))
            xout[i] = approx1(xout[i], x, y, *nxy, *method);
}
