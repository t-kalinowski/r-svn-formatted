/*
 *  modreg/src/ksmooth.c by B. D. Ripley  Copyright (C) 1998
 */

#include <S.h>
#include <math.h>

#define NA_REAL R_NaReal
extern double R_NaReal;

static double dokern(double x, int kern)
{
    if (kern == 1)
        return (1.0);
    if (kern == 2)
        return (exp(-0.5 * x * x));
}

void BDRksmooth(double *x, double *y, int *n, double *xp, double *yp, long *np, long *kern, double *bandwidth)
{
    int i, imin, j;
    double cutoff, num, den, x0, w, bw = *bandwidth;

    /* bandwidth is in units of half inter-quartile range. */
    if (*kern == 1)
    {
        bw *= 0.5;
        cutoff = bw;
    }
    if (*kern == 2)
    {
        bw *= 0.3706506;
        cutoff = 4 * bw;
    }
    while (x[imin] < xp[0] - cutoff && imin < *n)
        imin++;
    for (j = 0; j < *np; j++)
    {
        num = den = 0.0;
        x0 = xp[j];
        for (i = imin; i < *n; i++)
        {
            if (x[i] < x0 - cutoff)
                imin = i;
            if (x[i] > x0 + cutoff)
                break;
            w = dokern(fabs(x[i] - x0) / bw, *kern);
            num += w * y[i];
            den += w;
        }
        if (den > 0)
            yp[j] = num / den;
        else
            yp[j] = NA_REAL;
    }
}
