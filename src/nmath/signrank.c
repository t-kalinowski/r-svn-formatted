/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1999 R Core Team
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
 *
 *  SYNOPSIS
 *
 *    #include "Mathlib.h"
 *    double dsignrank(double x, double n)
 *    double dsignrank(double x, double n)
 *    double dsignrank(double x, double n)
 *    double rsignrank(double n)
 *
 *  DESCRIPTION
 *
 *    dsignrank    The density of the Wilcoxon Signed Rank distribution.
 *    psignrank    The distribution function of the Wilcoxon Signed Rank
 *                 distribution.
 *    qsignrank    The quantile function of the Wilcoxon Signed Rank
 *                 distribution.
 *    rsignrank    Random variates from the Wilcoxon Signed Rank
 *                 distribution.
 */

#include "Mathlib.h"

static double **w;

static void w_free(int n)
{
    int i;

    n = imax2(n, SIGNRANK_MAX);
    for (i = n; i >= 0; i--)
    {
        free((void *)w[i]);
    }
    free((void *)w);
    w = 0;
}

static void w_free_maybe(int n)
{
    if (n > SIGNRANK_MAX)
        w_free(n);
}

static void w_init_maybe(int n)
{
    if (w && (n > SIGNRANK_MAX))
        w_free(SIGNRANK_MAX);

    if (!w)
    {
        n = imax2(n, SIGNRANK_MAX);
        w = (double **)calloc(n + 1, sizeof(double *));
        if (!w)
            error("signrank allocation error");
    }
}

static double csignrank(int k, int n)
{
    int c, u, i;

    u = n * (n + 1) / 2;
    c = (int)(u / 2);

    if ((k < 0) || (k > u))
        return (0);
    if (k > c)
        k = u - k;
    if (w[n] == 0)
    {
        w[n] = (double *)calloc(c + 1, sizeof(double));
        for (i = 0; i <= c; i++)
            w[n][i] = -1;
    }
    if (w[n][k] < 0)
    {
        if (n == 0)
            w[n][k] = (k == 0);
        else
            w[n][k] = csignrank(k - n, n - 1) + csignrank(k, n - 1);
    }
    return (w[n][k]);
}

double dsignrank(double x, double n)
{
    double d;

#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(x) || ISNAN(n))
        return x + n;
#endif
    n = floor(n + 0.5);
    if (n <= 0)
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }

    x = floor(x + 0.5);

    if ((x < 0) || (x > (n * (n + 1) / 2)))
        return 0;

    w_init_maybe(n);
    d = exp(log(csignrank(x, n)) - n * log(2));
    w_free_maybe(n);

    return (d);
}

double psignrank(double x, double n)
{
    int i;
    double f, p;

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n))
        return x + n;
    if (!FINITE(n))
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }
#endif
    n = floor(n + 0.5);
    if (n <= 0)
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }

    x = floor(x + 0.5);
    if (x < 0.0)
        return 0;
    if (x >= n * (n + 1) / 2)
        return 1;

    w_init_maybe(n);
    f = exp(-log(n) * 2);
    if (x <= (n * (n + 1) / 4))
    {
        p = 0;
        for (i = 0; i <= x; i++)
            p += csignrank(i, n) * f;
    }
    else
    {
        p = 1;
        for (i = 0; i <= x; i++)
            p -= csignrank(i, n) * f;
    }
    w_free_maybe(n);

    return (p);
}

double qsignrank(double x, double n)
{
    double f, p, q;

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n))
        return x + n;
    if (!FINITE(x) || !FINITE(n))
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }
#endif
    n = floor(n + 0.5);
    if (x < 0 || x > 1 || n <= 0)
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }

    if (x == 0)
        return (0.0);
    if (x == 1)
        return (n * (n + 1) / 2);

    w_init_maybe(n);
    f = exp(-log(n) * 2);
    q = 0;
    if (x <= 0.5)
    {
        p = 0;
        for (;;)
        {
            p += csignrank(q, n) * f;
            if (p >= x)
                break;
            q++;
        }
    }
    else
    {
        p = 1;
        for (;;)
        {
            p -= csignrank(q, n) * f;
            if (p < x)
            {
                q = n * (n + 1) / 2 - q;
                break;
            }
            q++;
        }
    }
    w_free_maybe(n);

    return (q);
}

double rsignrank(double n)
{
    int i, k;
    double r;

#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(n))
        return (n);
#endif
    n = floor(n + 0.5);
    if (n < 0)
    {
        ML_ERROR(ME_DOMAIN);
        return ML_NAN;
    }

    if (n == 0)
        return (0);
    r = 0.0;
    k = (int)n;
    for (i = 0; i < k;)
    {
        r += (++i) * floor(sunif() + 0.5);
    }
    return (r);
}
