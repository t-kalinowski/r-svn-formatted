/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996, 1997  Robert Gentleman and Ross Ihaka
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

#include "Defn.h"
#include "FFDecl.h"

/* Fourier Transform for Univariate Spatial and Time Series */

SEXP do_fft(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP z, d;
    int i, inv, maxf, maxmaxf, maxmaxp, maxp, n, ndims, nseg, nspn;
    double *work;
    int *iwork;
    char *vmax;

    checkArity(op, args);

    z = CAR(args);

    switch (TYPEOF(z))
    {
    case INTSXP:
    case LGLSXP:
    case REALSXP:
        z = coerceVector(z, CPLXSXP);
        break;
    case CPLXSXP:
        if (NAMED(z))
            z = duplicate(z);
        break;
    default:
        errorcall(call, "non-numeric argument\n");
    }
    PROTECT(z);

    /* -2 for forward transform, complex values */
    /* +2 for backward transform, complex values */

    inv = asLogical(CADR(args));
    if (inv == NA_INTEGER || inv == 0)
        inv = -2;
    else
        inv = 2;

    if (LENGTH(z) > 1)
    {
        vmax = vmaxget();
        if (isNull(d = getAttrib(z, R_DimSymbol)))
        { /* temporal transform */
            n = length(z);
            fft_factor(n, &maxf, &maxp);
            if (maxf == 0)
                errorcall(call, "fft factorization error\n");
            work = (double *)R_alloc(4 * maxf, sizeof(double));
            iwork = (int *)R_alloc(maxp, sizeof(int));
            fft_work(&(COMPLEX(z)[0].r), &(COMPLEX(z)[0].i), 1, n, 1, inv, work, iwork);
        }
        else
        { /* spatial transform */
            maxmaxf = 1;
            maxmaxp = 1;
            ndims = LENGTH(d);
            for (i = 0; i < ndims; i++)
            {
                if (INTEGER(d)[i] > 1)
                {
                    fft_factor(INTEGER(d)[i], &maxf, &maxp);
                    if (maxf == 0)
                        errorcall(call, "fft factorization error\n");
                    if (maxf > maxmaxf)
                        maxmaxf = maxf;
                    if (maxp > maxmaxp)
                        maxmaxp = maxp;
                }
            }
            work = (double *)R_alloc(4 * maxmaxf, sizeof(double));
            iwork = (int *)R_alloc(maxmaxp, sizeof(int));
            nseg = LENGTH(z);
            n = 1;
            nspn = 1;
            for (i = 0; i < ndims; i++)
            {
                if (INTEGER(d)[i] > 1)
                {
                    nspn = nspn * n;
                    n = INTEGER(d)[i];
                    nseg = nseg / n;
                    fft_factor(n, &maxf, &maxp);
                    fft_work(&(COMPLEX(z)[0].r), &(COMPLEX(z)[0].i), nseg, n, nspn, inv, work, iwork);
                }
            }
        }
        vmaxset(vmax);
    }
    UNPROTECT(1);
    return z;
}

/* Fourier Transform for Vector-Valued Series */
/* Not to be confused with the spatial case. */

SEXP do_mvfft(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP z, d;
    int i, inv, maxf, maxp, n, p;
    double *work;
    int *iwork;
    char *vmax;

    checkArity(op, args);

    z = CAR(args);

    d = getAttrib(z, R_DimSymbol);
    if (length(d) > 2)
        errorcall(call, "vector-valued series required\n");
    n = INTEGER(d)[0];
    p = INTEGER(d)[1];

    switch (TYPEOF(z))
    {
    case INTSXP:
    case LGLSXP:
    case REALSXP:
        z = coerceVector(z, CPLXSXP);
        break;
    case CPLXSXP:
        if (NAMED(z))
            z = duplicate(z);
        break;
    default:
        errorcall(call, "non-numeric argument\n");
    }
    PROTECT(z);

    /* -2 for forward transform, complex values */
    /* +2 for backward transform, complex values */

    inv = asLogical(CADR(args));
    if (inv == NA_INTEGER || inv == 0)
        inv = -2;
    else
        inv = 2;

    if (n > 1)
    {
        vmax = vmaxget();
        fft_factor(n, &maxf, &maxp);
        if (maxf == 0)
            errorcall(call, "fft factorization error\n");
        work = (double *)R_alloc(4 * maxf, sizeof(double));
        iwork = (int *)R_alloc(maxp, sizeof(int));
        for (i = 0; i < p; i++)
            fft_work(&(COMPLEX(z)[i * n].r), &(COMPLEX(z)[i * n].i), 1, n, 1, inv, work, iwork);
        vmaxset(vmax);
    }
    UNPROTECT(1);
    return z;
}

static int ok(int n, int *f, int nf)
{
    int i;
    for (i = 0; i < nf; i++)
    {
        while (n % f[i] == 0)
        {
            if ((n = n / f[i]) == 1)
                return 1;
        }
    }
    return n == 1;
}

static int nextn(int n, int *f, int nf)
{
    while (!ok(n, f, nf))
        n++;
    return n;
}

SEXP do_nextn(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP n, f, ans;
    int i, nn, nf;
    checkArity(op, args);
    PROTECT(n = coerceVector(CAR(args), INTSXP));
    PROTECT(f = coerceVector(CADR(args), INTSXP));
    nn = LENGTH(n);
    nf = LENGTH(f);

    /* check the factors */

    if (nf == 0)
        errorcall(call, "no factors\n");
    for (i = 0; i < nf; i++)
        if (INTEGER(f)[i] == NA_INTEGER || INTEGER(f)[i] <= 1)
            errorcall(call, "invalid factors\n");

    ans = allocVector(INTSXP, nn);
    for (i = 0; i < nn; i++)
    {
        if (INTEGER(n)[i] == NA_INTEGER)
            INTEGER(ans)[i] = NA_INTEGER;
        else if (INTEGER(n)[i] <= 1)
            INTEGER(ans)[i] = 1;
        else
            INTEGER(ans)[i] = nextn(INTEGER(n)[i], INTEGER(f), nf);
    }
    UNPROTECT(2);
    return ans;
}
