/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1999-2000  The R Development Core Team
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
 *    #include "Rmath.h"
 *    double dsignrank(double x, double n, int give_log)
 *    double psignrank(double x, double n, int lower_tail, int log_p)
 *    double qsignrank(double x, double n, int lower_tail, int log_p)
 *    double rsignrank(double n)
 *
 *  DESCRIPTION
 *
 *    dsignrank	   The density of the Wilcoxon Signed Rank distribution.
 *    psignrank	   The distribution function of the Wilcoxon Signed Rank
 *		   distribution.
 *    qsignrank	   The quantile function of the Wilcoxon Signed Rank
 *		   distribution.
 *    rsignrank	   Random variates from the Wilcoxon Signed Rank
 *		   distribution.
 */

#include "nmath.h"
#include "dpq.h"

static double **w;

/* The idea is to allocate w of size SIGNRANK_MAX on the first small call, and
   to reallocate only for n > SIGNRANK_MAX, although for some reason
   realloc is not used */

static void w_free(int n)
{
    int i;

    if (!w)
        return;
    n = imax2(n, SIGNRANK_MAX);
    for (i = n; i >= 0; i--)
        if (w[i])
        {
            free((void *)w[i]);
            w[i] = 0;
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
            MATHLIB_ERROR("%s", "signrank allocation error");
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
        if (!w[n])
        {
            /* free up memory iff n is large */
            w_free_maybe(n);
            MATHLIB_ERROR("%s", "signrank allocation error");
        }
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

double dsignrank(double x, double n, int give_log)
{
    double d;

#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(x) || ISNAN(n))
        return (x + n);
#endif
    n = floor(n + 0.5);
    if (n <= 0)
        ML_ERR_return_NAN;

    if (fabs(x - floor(x + 0.5)) > 1e-7)
        return (R_D__0);
    x = floor(x + 0.5);
    if ((x < 0) || (x > (n * (n + 1) / 2)))
        return (R_D__0);

    w_init_maybe(n);
    d = R_D_exp(log(csignrank(x, n)) - n * M_LN2);
    w_free_maybe(n);

    return (d);
}

double psignrank(double x, double n, int lower_tail, int log_p)
{
    int i;
    double f, p;

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n))
        return (x + n);
#endif
    if (!R_FINITE(n))
        ML_ERR_return_NAN;
    n = floor(n + 0.5);
    if (n <= 0)
        ML_ERR_return_NAN;

    x = floor(x + 1e-7);
    if (x < 0.0)
        return (R_DT_0);
    if (x >= n * (n + 1) / 2)
        return (R_DT_1);

    w_init_maybe(n);
    f = exp(-n * M_LN2);
    p = 0;
    if (x <= (n * (n + 1) / 4))
    {
        for (i = 0; i <= x; i++)
            p += csignrank(i, n) * f;
    }
    else
    {
        x = n * (n + 1) / 2 - x;
        for (i = 0; i < x; i++)
            p += csignrank(i, n) * f;
        lower_tail = !lower_tail; /* p = 1 - p; */
    }
    w_free_maybe(n);

    return (R_DT_val(p));
} /* psignrank() */

double qsignrank(double x, double n, int lower_tail, int log_p)
{
    double f, p, q;

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n))
        return (x + n);
#endif
    if (!R_FINITE(x) || !R_FINITE(n))
        ML_ERR_return_NAN;
    R_Q_P01_check(x);

    n = floor(n + 0.5);
    if (n <= 0)
        ML_ERR_return_NAN;

    if (x == R_DT_0)
        return (0);
    if (x == R_DT_1)
        return (n * (n + 1) / 2);

    if (log_p || !lower_tail)
        x = R_DT_qIv(x); /* lower_tail,non-log "p" */

    w_init_maybe(n);
    f = exp(-n * M_LN2);
    p = 0;
    q = 0;
    if (x <= 0.5)
    {
        x = x - 10 * DBL_EPSILON;
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
        x = 1 - x + 10 * DBL_EPSILON;
        for (;;)
        {
            p += csignrank(q, n) * f;
            if (p > x)
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
        ML_ERR_return_NAN;

    if (n == 0)
        return (0);
    r = 0.0;
    k = (int)n;
    for (i = 0; i < k;)
    {
        r += (++i) * floor(unif_rand() + 0.5);
    }
    return (r);
}
