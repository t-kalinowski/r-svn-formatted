/*
 *  R : A Computer Langage for Statistical Data Analysis
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

#include "Defn.h"
#include "Mathlib.h"

static int ZeroSD;

static void cov_pairwise1(int n, int ncx, double *x, double *ans, int cor)
{
    double sum, xmean, ymean, xsd, ysd, *xx, *yy;
    int i, j, k, nobs;

    for (i = 0; i < ncx; i++)
    {
        for (j = 0; j <= i; j++)
        {
            xx = &x[i * n];
            yy = &x[j * n];
            nobs = 0;
            xmean = 0;
            ymean = 0;
            for (k = 0; k < n; k++)
            {
                if (FINITE(xx[k]) && FINITE(yy[k]))
                {
                    nobs += 1;
                    xmean += xx[k];
                    ymean += yy[k];
                }
            }
            if (nobs >= 2)
            {
                xmean /= nobs;
                ymean /= nobs;
                xsd = 0.0;
                ysd = 0.0;
                sum = 0.0;
                for (k = 0; k < n; k++)
                {
                    if (FINITE(xx[k]) && FINITE(yy[k]))
                    {
                        sum += (xx[k] - xmean) * (yy[k] - ymean);
                        if (cor)
                        {
                            xsd += (xx[k] - xmean) * (xx[k] - xmean);
                            ysd += (yy[k] - ymean) * (yy[k] - ymean);
                        }
                    }
                }
                if (cor)
                {
                    xsd = sqrt(xsd / (nobs - 1));
                    ysd = sqrt(ysd / (nobs - 1));
                    if (xsd == 0.0 || ysd == 0.0)
                    {
                        ZeroSD = 1;
                        sum = NA_REAL;
                    }
                    else
                        sum = (sum / (nobs - 1)) / (xsd * ysd);
                }
                else
                    sum = sum / (nobs - 1);
                ans[j + i * ncx] = ans[i + j * ncx] = sum;
            }
            else
                ans[j + i * ncx] = ans[i + j * ncx] = NA_REAL;
        }
    }
}

static void cov_pairwise2(int n, int ncx, int ncy, double *x, double *y, double *ans, int cor)
{
    double sum, xmean, ymean, xsd, ysd, *xx, *yy;
    int i, j, k, nobs;

    for (i = 0; i < ncx; i++)
    {
        for (j = 0; j < ncy; j++)
        {
            xx = &x[i * n];
            yy = &y[j * n];
            nobs = 0;
            xmean = 0;
            ymean = 0;
            for (k = 0; k < n; k++)
            {
                if (FINITE(xx[k]) && FINITE(yy[k]))
                {
                    nobs += 1;
                    xmean += xx[k];
                    ymean += yy[k];
                }
            }
            if (nobs >= 2)
            {
                xmean /= nobs;
                ymean /= nobs;
                xsd = 0.0;
                ysd = 0.0;
                sum = 0.0;
                for (k = 0; k < n; k++)
                {
                    if (FINITE(xx[k]) && FINITE(yy[k]))
                    {
                        sum += (xx[k] - xmean) * (yy[k] - ymean);
                        if (cor)
                        {
                            xsd += (xx[k] - xmean) * (xx[k] - xmean);
                            ysd += (yy[k] - ymean) * (yy[k] - ymean);
                        }
                    }
                }
                if (cor)
                {
                    xsd = sqrt(xsd / (nobs - 1));
                    ysd = sqrt(ysd / (nobs - 1));
                    if (xsd == 0.0 || ysd == 0.0)
                    {
                        ZeroSD = 1;
                        sum = NA_REAL;
                    }
                    else
                        sum = (sum / (nobs - 1)) / (xsd * ysd);
                }
                else
                    sum = sum / (nobs - 1);
                ans[i + j * ncx] = sum;
            }
            else
                ans[i + j * ncx] = NA_REAL;
        }
    }
}

static int NaFail;

static void cov_complete1(int n, int ncx, double *x, double *xm, int *ind, double *ans, int cor)
{
    double sum, xxm, yym, *xx, *yy;
    int i, j, k, nobs;

    /* total number of complete observations */
    nobs = 0;
    for (k = 0; k < n; k++)
    {
        if (ind[k] != 0)
            nobs++;
    }
    if (nobs <= 1)
    {
        for (i = 0; i < ncx; i++)
            for (j = 0; j < ncx; j++)
                ans[i + j * ncx] = NA_REAL;
        return;
    }

    /* variable means */
    for (i = 0; i < ncx; i++)
    {
        xx = &x[i * n];
        sum = 0.0;
        for (k = 0; k < n; k++)
            if (ind[k] != 0)
                sum += xx[k];
        xm[i] = sum / nobs;
    }

    for (i = 0; i < ncx; i++)
    {
        xx = &x[i * n];
        xxm = xm[i];
        for (j = 0; j <= i; j++)
        {
            yy = &x[j * n];
            yym = xm[j];
            sum = 0.0;
            for (k = 0; k < n; k++)
                if (ind[k] != 0)
                    sum += (xx[k] - xxm) * (yy[k] - yym);
            ans[j + i * ncx] = ans[i + j * ncx] = sum / (nobs - 1);
        }
    }
    if (cor)
    {
        for (i = 0; i < ncx; i++)
            xm[i] = sqrt(ans[i + i * ncx]);
        for (i = 0; i < ncx; i++)
        {
            for (j = 0; j < i; j++)
            {
                if (xm[i] == 0 || xm[j] == 0)
                {
                    ZeroSD = 1;
                    ans[j + i * ncx] = ans[i + j * ncx] = NA_REAL;
                }
                else
                {
                    ans[j + i * ncx] = ans[i + j * ncx] = ans[i + j * ncx] / (xm[i] * xm[j]);
                }
            }
            ans[i + i * ncx] = 1.0;
        }
    }
}

static void cov_complete2(int n, int ncx, int ncy, double *x, double *y, double *xm, double *ym, int *ind, double *ans,
                          int cor)
{
    double sum, xxm, yym, *xx, *yy;
    int i, j, k, nobs;

    /* total number of complete observations */
    nobs = 0;
    for (k = 0; k < n; k++)
    {
        if (ind[k] != 0)
            nobs++;
    }
    if (nobs <= 1)
    {
        for (i = 0; i < ncx; i++)
            for (j = 0; j < ncy; j++)
                ans[i + j * ncx] = NA_REAL;
        return;
    }

    /* variable means */
    for (j = 0; j < ncx; j++)
    {
        xx = &x[j * n];
        sum = 0.0;
        for (i = 0; i < n; i++)
            if (ind[i] != 0)
                sum += xx[i];
        xm[j] = sum / nobs;
    }
    for (j = 0; j < ncy; j++)
    {
        yy = &y[j * n];
        sum = 0.0;
        for (i = 0; i < n; i++)
            if (ind[i] != 0)
                sum += yy[i];
        ym[j] = sum / nobs;
    }

    for (i = 0; i < ncx; i++)
    {
        xx = &x[i * n];
        xxm = xm[i];
        for (j = 0; j < ncy; j++)
        {
            yy = &y[j * n];
            yym = ym[j];
            sum = 0.0;
            for (k = 0; k < n; k++)
                if (ind[k] != 0)
                    sum += (xx[k] - xxm) * (yy[k] - yym);
            ans[i + j * ncx] = sum / (nobs - 1);
        }
    }

    if (cor)
    {
        for (i = 0; i < ncx; i++)
        {
            xx = &x[i * n];
            xxm = xm[i];
            sum = 0.0;
            for (k = 0; k < n; k++)
                if (ind[k] != 0)
                    sum += (xx[k] - xxm) * (xx[k] - xxm);
            xm[i] = sqrt(sum / (nobs - 1));
        }
        for (j = 0; j < ncy; j++)
        {
            yy = &y[j * n];
            yym = ym[j];
            sum = 0.0;
            for (k = 0; k < n; k++)
                if (ind[k] != 0)
                    sum += (yy[k] - yym) * (yy[k] - yym);
            ym[j] = sqrt(sum / (nobs - 1));
        }
        for (i = 0; i < ncx; i++)
        {
            for (j = 0; j < ncy; j++)
            {
                if (xm[i] == 0.0 || ym[j] == 0.0)
                {
                    ZeroSD = 1;
                    ans[i + j * ncx] = NA_REAL;
                }
                else
                    ans[i + j * ncx] = ans[i + j * ncx] / (xm[i] * ym[j]);
            }
        }
    }
}

/* This might look slightly inefficient, but it is designed to */
/* optimise paging in virtual memory systems ... */

static void complete1(int n, int ncx, double *x, int *ind)
{
    double *z;
    int i, j;
    for (i = 0; i < n; i++)
        ind[i] = 1;
    for (j = 0; j < ncx; j++)
    {
        z = &x[j * n];
        for (i = 0; i < n; i++)
            if (!FINITE(z[i]))
            {
                if (NaFail)
                    error("missing observations in cov/cor\n");
                else
                    ind[i] = 0;
            }
    }
}

static void complete2(int n, int ncx, int ncy, double *x, double *y, int *ind)
{
    double *z;
    int i, j;
    for (i = 0; i < n; i++)
        ind[i] = 1;
    for (j = 0; j < ncx; j++)
    {
        z = &x[j * n];
        for (i = 0; i < n; i++)
            if (!FINITE(z[i]))
            {
                if (NaFail)
                    error("missing observations in cov/cor\n");
                else
                    ind[i] = 0;
            }
    }
    for (j = 0; j < ncy; j++)
    {
        z = &y[j * n];
        for (i = 0; i < n; i++)
            if (!FINITE(z[i]))
            {
                if (NaFail)
                    error("missing observations in cov/cor\n");
                else
                    ind[i] = 0;
            }
    }
}

SEXP do_cov(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP x, y, ans, xm, ym, ind;
    int ansmat, cor, method, n, ncx, ncy, pair;

    checkArity(op, args);
    cor = PRIMVAL(op);

    /* Argument-1: x */

    x = CAR(args) = coerceVector(CAR(args), REALSXP);
    if ((ansmat = isMatrix(x)))
    {
        n = nrows(x);
        ncx = ncols(x);
    }
    else
    {
        n = length(x);
        ncx = 1;
    }
    args = CDR(args);

    /* Argument-2: y */

    if (isNull(CAR(args)))
    {
        y = R_NilValue;
        ncy = ncx;
    }
    else
    {
        y = CAR(args) = coerceVector(CAR(args), REALSXP);
        if (isMatrix(y))
        {
            if (nrows(y) != n)
                errorcall(call, "incompatible dimensions\n");
            ncy = ncols(y);
        }
        else
        {
            if (length(y) != n)
                errorcall(call, "incompatible dimensions\n");
            ncy = 1;
        }
        ansmat = (ansmat || isMatrix(y));
    }
    args = CDR(args);

    /* Argument-3: method */

    method = asInteger(CAR(args));
    switch (method)
    {
    case 1: /* none */
        NaFail = 1;
        pair = 0;
        break;
    case 2: /* complete */
        NaFail = 0;
        pair = 0;
        break;
    case 3: /* pairwise */
        NaFail = 0;
        pair = 1;
        break;
    default:
        errorcall(call, "invalid computational method\n");
    }
    if (pair == NA_INTEGER)
        pair = 0;

    if (ansmat)
        PROTECT(ans = allocMatrix(REALSXP, ncx, ncy));
    else
        PROTECT(ans = allocVector(REALSXP, ncx * ncy));

    ZeroSD = 0;

    if (isNull(y))
    {
        if (pair == 0)
        { /* complete */
            PROTECT(xm = allocVector(REALSXP, ncx));
            PROTECT(ind = allocVector(INTSXP, n));
            complete1(n, ncx, REAL(x), INTEGER(ind));
            cov_complete1(n, ncx, REAL(x), REAL(xm), INTEGER(ind), REAL(ans), cor);
            UNPROTECT(2);
        }
        else
        { /* pairwise */
            cov_pairwise1(n, ncx, REAL(x), REAL(ans), cor);
        }
    }
    else
    {
        if (pair == 0)
        { /* complete */
            PROTECT(xm = allocVector(REALSXP, ncx));
            PROTECT(ym = allocVector(REALSXP, ncy));
            PROTECT(ind = allocVector(INTSXP, n));
            complete2(n, ncx, ncy, REAL(x), REAL(y), INTEGER(ind));
            cov_complete2(n, ncx, ncy, REAL(x), REAL(y), REAL(xm), REAL(ym), INTEGER(ind), REAL(ans), cor);
            UNPROTECT(3);
        }
        else
        { /* pairwise */
            cov_pairwise2(n, ncx, ncy, REAL(x), REAL(y), REAL(ans), cor);
        }
    }
    if (ansmat)
    {
        if (isNull(y))
        {
            x = getAttrib(x, R_DimNamesSymbol);
            if (!isNull(x) && !isNull(CADR(x)))
            {
                PROTECT(ind = allocList(2));
                CAR(ind) = CADR(x);
                CADR(ind) = CADR(x);
                setAttrib(ans, R_DimNamesSymbol, ind);
                UNPROTECT(1);
            }
        }
        else
        {
            x = getAttrib(x, R_DimNamesSymbol);
            y = getAttrib(y, R_DimNamesSymbol);
            if ((!isNull(x) && !isNull(CADR(x))) || (!isNull(y) && !isNull(CADR(y))))
            {
                PROTECT(ind = allocList(2));
                CAR(ind) = CADR(x);
                CADR(ind) = CADR(y);
                setAttrib(ans, R_DimNamesSymbol, ind);
                UNPROTECT(1);
            }
        }
    }
    UNPROTECT(1);
    if (ZeroSD)
        warning("standard deviation equal to zero in cor\n");
    return ans;
}
