/*
 *  R : A Computer Language for Statistical Data Analysis
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

#include "Mathlib.h"
#include "Error.h"
#include "Arith.h"

/* bincode  cuts up the data using half open intervals defined as [a,b)
   bincode2 cuts up the data using half open intervals defined as (a,b]

 MM: Shouldn't we afford to collapse these two into one,
     There would have to be 'n' extra "if" evaluations only ... ?

*/
int bincode(double *x, int *pn, double *breaks, int *pnb, int *code, int *include_border, int *naok)
{
    int i, lo, hi;
    int n, nb1, new;

    n = *pn;
    nb1 = *pnb - 1;

    for (i = 0; i < n; i++)
        if (FINITE(x[i]))
        {
            lo = 0;
            hi = nb1;
            if (x[i] < breaks[lo] || breaks[hi] < x[i] || (x[i] == breaks[hi] && !*include_border))
                code[i] = NA_INTEGER;
            else
            {
                while (hi - lo >= 2)
                {
                    new = (hi + lo) / 2;
                    if (x[i] >= breaks[new])
                        lo = new;
                    else
                        hi = new;
                }
                code[i] = lo + 1;
            }
        }
        else if (!*naok)
            error("NA's in .C(\"bincode\",... NAOK=FALSE)");
    return 0;
}

int bincode2(double *x, int *pn, double *breaks, int *pnb, int *code, int *include_border, int *naok)
{
    int i, lo, hi;
    int n, nb1, new;

    n = *pn;
    nb1 = *pnb - 1;

    for (i = 0; i < n; i++)
        if (FINITE(x[i]))
        {
            lo = 0;
            hi = nb1;
            if (x[i] < breaks[lo] || breaks[hi] < x[i] || (x[i] == breaks[lo] && !*include_border))
                code[i] = NA_INTEGER;
            else
            {
                while (hi - lo >= 2)
                {
                    new = (hi + lo) / 2;
                    if (x[i] > breaks[new])
                        lo = new;
                    else
                        hi = new;
                }
                code[i] = lo + 1;
            }
        }
        else if (!*naok)
            error("NA's in .C(\"bincode\",... NAOK=FALSE)");
    return 0;
}

/* bincount is called by  hist(.)  [only]
 *
 * bincount *counts* like bincode2, i.e. half open intervals defined as (a,b]
 */
int bincount(double *x, int *pn, double *breaks, int *pnb, int *count, int *include_border, int *naok)
{
    int i, lo, hi;
    int n, nb1, new;

    n = *pn;
    nb1 = *pnb - 1;

    for (i = 0; i < nb1; i++)
        count[i] = 0;

    for (i = 0; i < n; i++)
        if (FINITE(x[i]))
        {
            lo = 0;
            hi = nb1;
            if (breaks[lo] <= x[i] && (x[i] < breaks[hi] || (x[i] == breaks[hi] && *include_border)))
            {
                while (hi - lo >= 2)
                {
                    new = (hi + lo) / 2;
                    if (x[i] > breaks[new])
                        lo = new;
                    else
                        hi = new;
                }
                count[lo] += 1;
            }
        }
        else if (!*naok)
            error("NA's in .C(\"bincode\",... NAOK=FALSE)");
    return 0;
}
