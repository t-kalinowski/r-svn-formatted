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
 *    double fprec(double x, double digits);
 *
 *  DESCRIPTION
 *
 *    Returns the value of x rounded to "digits" significant
 *    decimal digits.
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *    Some modifications have been made so that the routines
 *    conform to the IEEE 754 standard.
 */

#include "Mathlib.h"

/* Improvements by Martin Maechler, May 1997,
/* further ones, Feb.2000:
   Replace  pow(x, (double)i) by  R_pow_di(x, i) {and use  int dig} */

/* DBL_DIG := digits of precision of a double, usually 15 : */
#define MAXPLACES DBL_DIG

double fprec(double x, double digits)
{
    double l10, pow10, sgn, p10, P10;
    int e10, e2, do_round, dig;
    /* Max.expon. of 10 (=308.2547) */
    const double max10e = DBL_MAX_EXP * M_LOG10_2;

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(digits))
        return x + digits;
    if (!R_FINITE(x))
        return x;
    if (!R_FINITE(digits))
    {
        if (digits > 0)
            return x;
        else
            return 0;
    }
#endif
    if (x == 0)
        return x;
    dig = (int)floor(digits + 0.5);
    if (dig > MAXPLACES)
    {
        return x;
    }
    else if (dig < 1)
        dig = 1;

    sgn = 1.0;
    if (x < 0.0)
    {
        sgn = -sgn;
        x = -x;
    }
    l10 = log10(x);
    e10 = (int)(dig - 1 - floor(l10));
    if (fabs(l10) < max10e - 2)
    {
        pow10 = R_pow_di(10., e10);
        return (sgn * floor(x * pow10 + 0.5) / pow10);
    }
    else
    { /* -- LARGE or small -- */
        do_round = max10e - l10 >= R_pow_di(10., -dig);
        e2 = (e10 > 0) ? 16 : -16;
        p10 = R_pow_di(10., e2);
        x *= p10;
        P10 = R_pow_di(10., e10 - e2);
        x *= P10;
        /*-- p10 * P10 = 10 ^ e10 */
        if (do_round)
            x += 0.5;
        x = floor(x) / p10;
        return (sgn * x / P10);
    }
}
