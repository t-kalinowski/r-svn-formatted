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

#include "Defn.h"

int OneIndex(SEXP x, SEXP s, int partial, SEXP *newname)
{
    SEXP names;
    int i, index, len, nx;

    if (length(s) > 1)
        error("attempt to select more than one element\n");
    if (length(s) < 1)
        error("attempt to select less than one element\n");

    index = -1;
    *newname = R_NilValue;
    switch (TYPEOF(s))
    {
    case LGLSXP:
    case INTSXP:
        index = INTEGER(s)[0] - 1;
        break;
    case REALSXP:
        index = REAL(s)[0] - 1;
        break;
    case STRSXP:
        nx = length(x);
        names = getAttrib(x, R_NamesSymbol);
        if (names != R_NilValue)
        {
            /* Try for exact match */
            for (i = 0; i < nx; i++)
                if (streql(CHAR(STRING(names)[i]), CHAR(STRING(s)[0])))
                {
                    index = i;
                    break;
                }
            /* Try for partial match */
            if (partial && index < 0)
            {
                len = strlen(CHAR(STRING(s)[0]));
                for (i = 0; i < nx; i++)
                {
                    if (!strncmp(CHAR(STRING(names)[i]), CHAR(STRING(s)[0]), len))
                    {
                        if (index == -1)
                            index = i;
                        else
                            index = -2;
                    }
                }
            }
        }
        if (index == -1)
            index = nx;
        *newname = STRING(s)[0];
        break;
    case SYMSXP:
        nx = length(x);
        names = getAttrib(x, R_NamesSymbol);
        if (names != R_NilValue)
        {
            for (i = 0; i < nx; i++)
                if (streql(CHAR(STRING(names)[i]), CHAR(PRINTNAME(s))))
                {
                    index = i;
                    break;
                }
        }
        if (index == -1)
            index = nx;
        *newname = STRING(s)[0];
        break;
    default:
        error("invalid subscript type\n");
    }
    return index;
}

/* Get a single index for the [[ operator. */
/* Check that only one index is being selected. */

int get1index(SEXP s, SEXP names, int pok)
{
    int index, i, len;

    if (length(s) > 1)
        error("attempt to select more than one element\n");
    if (length(s) < 1)
        error("attempt to select less than one element\n");

    index = -1;
    switch (TYPEOF(s))
    {
    case LGLSXP:
    case INTSXP:
        index = INTEGER(s)[0] - 1;
        break;
    case REALSXP:
        index = REAL(s)[0] - 1;
        break;
    case STRSXP:
        /* Try for exact match */
        for (i = 0; i < length(names); i++)
            if (streql(CHAR(STRING(names)[i]), CHAR(STRING(s)[0])))
            {
                index = i;
                break;
            }
        /* Try for partial match */
        if (pok && index < 0)
        {
            len = strlen(CHAR(STRING(s)[0]));
            for (i = 0; i < length(names); i++)
            {
                if (!strncmp(CHAR(STRING(names)[i]), CHAR(STRING(s)[0]), len))
                {
                    if (index == -1)
                        index = i;
                    else
                        index = -2;
                }
            }
        }
        break;
    case SYMSXP:
        for (i = 0; i < length(names); i++)
            if (streql(CHAR(STRING(names)[i]), CHAR(PRINTNAME(s))))
            {
                index = i;
                break;
            }
    default:
        error("invalid subscript type\n");
    }
    return index;
}

/* Special Matrix Subscripting: Handles the case x[i] where */
/* x is an n-way array and i is a matrix with n columns. */
/* This code returns a vector containing the integer subscripts */
/* to be extracted when x is regarded as unravelled. */

SEXP mat2indsub(SEXP dims, SEXP s)
{
    int tdim, j, i, nrs = nrows(s);
    SEXP rvec;

    PROTECT(rvec = allocVector(INTSXP, nrs));
    s = coerceVector(s, INTSXP);
    setIVector(INTEGER(rvec), nrs, 0);

    /* compute 0-based subscripts */
    for (i = 0; i < nrs; i++)
    {
        tdim = 1;
        for (j = 0; j < LENGTH(dims); j++)
        {
            if (INTEGER(s)[i + j * nrs] == NA_INTEGER)
            {
                INTEGER(rvec)[i] = NA_INTEGER;
                break;
            }
            if (INTEGER(s)[i + j * nrs] > INTEGER(dims)[j])
                error("subscript out of bounds\n");
            INTEGER(rvec)[i] += (INTEGER(s)[i + j * nrs] - 1) * tdim;
            tdim *= INTEGER(dims)[j];
        }
        /* transform to 1 based subscripting */
        if (INTEGER(rvec)[i] != NA_INTEGER)
            INTEGER(rvec)[i]++;
    }
    UNPROTECT(1);
    return (rvec);
}

static SEXP nullSubscript(int n)
{
    int i;
    SEXP index;
    index = allocVector(INTSXP, n);
    for (i = 0; i < n; i++)
        INTEGER(index)[i] = i + 1;
    return index;
}

static SEXP logicalSubscript(SEXP s, int ns, int nx)
{
    int count, i;
    SEXP index;
    if (ns > nx)
        error("subscript (\%d) out of bounds, should be at most %d\n", ns, nx);
    count = 0;
    for (i = 0; i < nx; i++)
        if (LOGICAL(s)[i % ns])
            count++;
    index = allocVector(INTSXP, count);
    count = 0;
    for (i = 0; i < nx; i++)
        if (LOGICAL(s)[i % ns])
        {
            if (LOGICAL(s)[i % ns] == NA_LOGICAL)
                INTEGER(index)[count++] = NA_INTEGER;
            else
                INTEGER(index)[count++] = i + 1;
        }
    return index;
}

static SEXP negativeSubscript(SEXP s, int ns, int nx)
{
    SEXP index;
    int i;
    PROTECT(index = allocVector(INTSXP, nx));
    for (i = 0; i < nx; i++)
        INTEGER(index)[i] = 1;
    for (i = 0; i < ns; i++)
        if (INTEGER(s)[i] != 0)
            INTEGER(index)[-INTEGER(s)[i] - 1] = 0;
    s = logicalSubscript(index, nx, nx);
    UNPROTECT(1);
    return s;
}

static SEXP positiveSubscript(SEXP s, int ns, int nx)
{
    SEXP index;
    int i, zct = 0;
    for (i = 0; i < ns; i++)
    {
        if (INTEGER(s)[i] == 0)
            zct++;
    }
    if (zct)
    {
        index = allocVector(INTSXP, (ns - zct));
        for (i = 0, zct = 0; i < ns; i++)
            if (INTEGER(s)[i] != 0)
                INTEGER(index)[zct++] = INTEGER(s)[i];
        return index;
    }
    else
        return s;
}

static SEXP integerSubscript(SEXP s, int ns, int nx, int *stretch)
{
    int i, ii, min, max, canstretch;
    canstretch = *stretch;
    *stretch = 0;
    min = 0;
    max = 0;
    for (i = 0; i < ns; i++)
    {
        ii = INTEGER(s)[i];
        if (ii != NA_INTEGER)
        {
            if (ii < min)
                min = ii;
            if (ii > max)
                max = ii;
        }
    }
    if (min < -nx)
        error("subscript out of bounds\n");
    if (max > nx)
    {
        if (canstretch)
            *stretch = max;
        else
            error("subscript out of bounds\n");
    }
    if (min < 0)
    {
        if (max == 0)
            return negativeSubscript(s, ns, nx);
        else
            error("only 0's may mix with negative subscripts\n");
    }
    else
        return positiveSubscript(s, ns, nx);
    return R_NilValue;
}

static SEXP stringSubscript(SEXP s, int ns, int nx, SEXP names, int *stretch)
{
    SEXP index, indexnames;
    int i, j, nnames, sub, extra;
    PROTECT(s);
    PROTECT(names);
    PROTECT(index = allocVector(INTSXP, ns));
    PROTECT(indexnames = allocVector(STRSXP, ns));
    nnames = nx;
    extra = nnames;

    /* Process each of the subscripts */
    /* First we compare with the names on the vector */
    /* and then (if there is no match) with each of */
    /* the previous subscripts. */

    for (i = 0; i < ns; i++)
    {
        sub = 0;
        if (names != R_NilValue)
        {
            for (j = 0; j < nnames; j++)
                if (NonNullStringMatch(STRING(s)[i], STRING(names)[j]))
                {
                    sub = j + 1;
                    STRING(indexnames)[i] = R_NilValue;
                    break;
                }
        }
        if (sub == 0)
        {
            for (j = 0; j < i; j++)
                if (NonNullStringMatch(STRING(s)[i], STRING(s)[j]))
                {
                    sub = INTEGER(index)[j];
                    STRING(indexnames)[i] = STRING(indexnames)[sub - 1];
                    break;
                }
        }
        if (sub == 0)
        {
            extra += 1;
            sub = extra;
            STRING(indexnames)[i] = STRING(s)[i];
        }
        INTEGER(index)[i] = sub;
    }
    /* Ghastly hack!  We attach the new names to the attribute */
    /* slot on the returned subscript vector. */
    if (extra != nnames)
    {
        *stretch = extra;
        ATTRIB(index) = indexnames;
    }
    UNPROTECT(4);
    return index;
}

/* Array Subscripts.  dim is the dimension (0 to k-1), s is */
/* the subscript list, x is the array to be subscripted. */

SEXP arraySubscript(int dim, SEXP s, SEXP x)
{
    int i, nd, ns, stretch = 0;
    SEXP dims, dnames;
    ns = length(s);
    dims = getAttrib(x, R_DimSymbol);
    nd = INTEGER(dims)[dim];

    switch (TYPEOF(s))
    {
    case NILSXP:
        return allocVector(INTSXP, 0);
    case LGLSXP:
        return logicalSubscript(s, ns, nd);
    case INTSXP:
        return integerSubscript(s, ns, nd, &stretch);
    case REALSXP:
        return integerSubscript(coerceVector(s, INTSXP), ns, nd, &stretch);
    case STRSXP:
        dnames = getAttrib(x, R_DimNamesSymbol);
        if (dnames == R_NilValue)
            error("no dimnames attribute for array\n");
#ifdef NEWLIST
        dnames = VECTOR(dnames)[dim];
#else
        for (i = 0; i < dim; i++)
            dnames = CDR(dnames);
        dnames = CAR(dnames);
#endif
        return stringSubscript(s, ns, nd, dnames, &stretch);
    case SYMSXP:
        if (s == R_MissingArg)
            return nullSubscript(nd);
    default:
        error("invalid subscript\n");
    }
    return R_NilValue;
}

/* Subscript creation.  The first thing we do is check to see */
/* if there are any user supplied NULL's, these result in */
/* returning a vector of length 0. */

SEXP makeSubscript(SEXP x, SEXP s, int *stretch)
{
    int nx, ns;
    SEXP ans, names, tmp;

    ans = R_NilValue;
    if (isVector(x) || isList(x) || isLanguage(x))
    {
        names = getAttrib(x, R_NamesSymbol);
        nx = length(x);
        ns = length(s);
        switch (TYPEOF(s))
        {
        case NILSXP:
            *stretch = 0;
            ans = allocVector(INTSXP, 0);
            break;
        case LGLSXP:
            *stretch = 0;
            ans = logicalSubscript(s, ns, nx);
            break;
        case INTSXP:
            ans = integerSubscript(s, ns, nx, stretch);
            break;
        case REALSXP:
            PROTECT(tmp = coerceVector(s, INTSXP));
            ans = integerSubscript(tmp, ns, nx, stretch);
            UNPROTECT(1);
            break;
        case STRSXP:
            *stretch = 0;
            ans = stringSubscript(s, ns, nx, names, stretch);
            break;
        case SYMSXP:
            *stretch = 0;
            if (s == R_MissingArg)
            {
                ans = nullSubscript(nx);
                break;
            }
        default:
            error("invalid subscript type\n");
        }
    }
    else
        error("subscripting on non-vector\n");
    return ans;
}
