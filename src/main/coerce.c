/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995-1998  Robert Gentleman and Ross Ihaka
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
#include "Print.h"

/* This section of code handles type conversion for elements */
/* of data vectors.  Type coersion throughout R should use these */
/* routines to ensure consistency. */

static char *truenames[] = {
    "T", "True", "TRUE", "true", (char *)0,
};

static char *falsenames[] = {
    "F", "False", "FALSE", "false", (char *)0,
};

#define WARN_NA 1
#define WARN_INACC 2
#define WARN_IMAG 4

void CoersionWarning(int warn)
{
    if (warn & WARN_NA)
        warning("NAs introduced by coersion\n");
    if (warn & WARN_INACC)
        warning("inaccurate integer conversion in coersion\n");
    if (warn & WARN_IMAG)
        warning("imaginary parts discarded in coersion\n");
}

int LogicalFromInteger(int x, int *warn)
{
    return (x == NA_INTEGER) ? NA_LOGICAL : (x != 0);
}

int LogicalFromReal(double x, int *warn)
{
    return ISNAN(x) ? NA_LOGICAL : (x != 0);
}

int LogicalFromComplex(complex x, int *warn)
{
    return (ISNAN(x.r) || ISNAN(x.i)) ? NA_LOGICAL : (x.r != 0 || x.i != 0);
}

int LogicalFromString(SEXP x, int *warn)
{
    if (x != R_NaString)
    {
        int i;
        for (i = 0; truenames[i]; i++)
            if (!strcmp(CHAR(x), truenames[i]))
                return 1;
        for (i = 0; falsenames[i]; i++)
            if (!strcmp(CHAR(x), falsenames[i]))
                return 0;
    }
    return NA_LOGICAL;
}

int IntegerFromLogical(int x, int *warn)
{
    return (x == NA_LOGICAL) ? NA_INTEGER : x;
}

int IntegerFromReal(double x, int *warn)
{
    if (ISNAN(x))
        return NA_INTEGER;
    else if (x > INT_MAX)
    {
        *warn |= WARN_INACC;
        return INT_MAX;
    }
    else if (x <= INT_MIN)
    {
        *warn |= WARN_INACC;
        return INT_MIN + 1;
    }
    return x;
}

int IntegerFromComplex(complex x, int *warn)
{
    if (ISNAN(x.r) || ISNAN(x.i))
        return NA_INTEGER;
    else if (x.r > INT_MAX)
    {
        *warn |= WARN_INACC;
        return INT_MAX;
    }
    else if (x.r <= INT_MIN)
    {
        *warn |= WARN_INACC;
        return INT_MIN + 1;
    }
    if (x.i != 0)
        *warn |= WARN_IMAG;
    return x.r;
}

int IntegerFromString(SEXP x, int *warn)
{
    double xdouble;
    char *endp;
    if (x != R_NaString)
    {
        xdouble = strtod(CHAR(x), &endp);
        if (*endp == '\0')
        {
            if (xdouble > INT_MAX)
            {
                *warn |= WARN_INACC;
                return INT_MAX;
            }
            else if (xdouble < INT_MIN + 1)
            {
                *warn |= WARN_INACC;
                return INT_MIN;
            }
            else
                return xdouble;
        }
        else
            *warn |= WARN_NA;
    }
    return NA_INTEGER;
}

double RealFromLogical(int x, int *warn)
{
    return (x == NA_LOGICAL) ? NA_REAL : x;
}

double RealFromInteger(int x, int *warn)
{
    if (x == NA_INTEGER)
        return NA_REAL;
    else
        return x;
}

double RealFromComplex(complex x, int *warn)
{
    if (ISNAN(x.r) || ISNAN(x.i))
        return NA_INTEGER;
    if (x.i != 0)
        *warn |= WARN_IMAG;
    return x.r;
}

double RealFromString(SEXP x, int *warn)
{
    double xdouble;
    char *endp;
    if (x != R_NaString)
    {
        xdouble = strtod(CHAR(x), &endp);
        if (*endp == '\0')
            return xdouble;
        else
            *warn |= WARN_NA;
    }
    return NA_INTEGER;
}

complex ComplexFromLogical(int x, int *warn)
{
    complex z;
    if (x == NA_LOGICAL)
    {
        z.r = NA_REAL;
        z.i = NA_REAL;
    }
    else
    {
        z.r = x;
        z.i = 0;
    }
    return z;
}

complex ComplexFromInteger(int x, int *warn)
{
    complex z;
    if (x == NA_INTEGER)
    {
        z.r = NA_REAL;
        z.i = NA_REAL;
    }
    else
    {
        z.r = x;
        z.i = 0;
    }
    return z;
}

complex ComplexFromReal(double x, int *warn)
{
    complex z;
    if (ISNAN(x))
    {
        z.r = NA_REAL;
        z.i = NA_REAL;
    }
    else
    {
        z.r = x;
        z.i = 0;
    }
    return z;
}

complex ComplexFromString(SEXP x, int *warn)
{
    double xr, xi;
    complex z;
    char *endp = CHAR(x);
    ;
    z.r = z.i = NA_REAL;
    if (x != R_NaString)
    {
        xr = strtod(endp, &endp);
        if (*endp == '\0')
        {
            z.r = xr;
            z.i = 0.0;
        }
        else if (*endp == '+' || *endp == '-')
        {
            xi = strtod(endp, &endp);
            if (endp[0] == 'i' && endp[1] == '\0')
            {
                z.r = xr;
                z.i = xi;
            }
            else
                *warn |= WARN_NA;
        }
        else
            *warn |= WARN_NA;
    }
    return z;
}

SEXP StringFromLogical(int x, int *warn)
{
    int w;
    formatLogical(&x, 1, &w);
    return mkChar(EncodeLogical(x, w));
}

SEXP StringFromInteger(int x, int *warn)
{
    int w;
    formatInteger(&x, 1, &w);
    return mkChar(EncodeInteger(x, w));
}

SEXP StringFromReal(double x, int *warn)
{
    int w, d, e;
    formatReal(&x, 1, &w, &d, &e);
    return mkChar(EncodeReal(x, w, d, e));
}

SEXP StringFromComplex(complex x, int *warn)
{
    int wr, dr, er, wi, di, ei;
    formatComplex(&x, 1, &wr, &dr, &er, &wi, &di, &ei);
    return mkChar(EncodeComplex(x, wr, dr, er, wi, di, ei));
}

/* Conversion between the two list types (LISTSXP and VECSXP). */

SEXP PairToVectorList(SEXP x)
{
    SEXP xptr, xnew, xnames, blank;
    int i, len = 0, named = 0;
    for (xptr = x; xptr != R_NilValue; xptr = CDR(xptr))
    {
        named = named | (TAG(xptr) != R_NilValue);
        len++;
    }
    PROTECT(x);
    PROTECT(xnew = allocVector(VECSXP, len));
    for (i = 0, xptr = x; i < len; i++, xptr = CDR(xptr))
        VECTOR(xnew)[i] = CAR(xptr);
    if (named)
    {
        PROTECT(xnames = allocVector(STRSXP, len));
        xptr = x;
        for (i = 0, xptr = x; i < len; i++, xptr = CDR(xptr))
        {
            if (TAG(xptr) == R_NilValue)
                STRING(xnames)[i] = R_BlankString;
            else
                STRING(xnames)[i] = PRINTNAME(TAG(xptr));
        }
        setAttrib(xnew, R_NamesSymbol, xnames);
        UNPROTECT(1);
    }
    copyMostAttrib(x, xnew);
    UNPROTECT(2);
    return xnew;
}

SEXP VectorToPairList(SEXP x)
{
    SEXP xptr, xnew, xnames;
    int i, len, named;
    len = length(x);
    PROTECT(x);
    PROTECT(xnew = allocList(len));
    PROTECT(xnames = getAttrib(x, R_NamesSymbol));
    named = (xnames != R_NilValue);
    xptr = xnew;
    for (i = 0; i < len; i++)
    {
        CAR(xptr) = VECTOR(x)[i];
        if (named && CHAR(STRING(xnames)[i])[0] != '\0')
            TAG(xptr) = install(CHAR(STRING(xnames)[i]));
        xptr = CDR(xptr);
    }
    copyMostAttrib(x, xnew);
    UNPROTECT(3);
    return xnew;
}

static SEXP coerceToSymbol(SEXP v)
{
    SEXP ans;
    int warn;
    if (length(v) <= 0)
        error("Invalid data of mode \"%s\" (too short)\n", type2str(TYPEOF(v)));
    PROTECT(v);
    switch (TYPEOF(v))
    {
    case LGLSXP:
        ans = StringFromLogical(LOGICAL(v)[0], &warn);
        break;
    case INTSXP:
        ans = StringFromInteger(INTEGER(v)[0], &warn);
        break;
    case REALSXP:
        ans = StringFromReal(REAL(v)[0], &warn);
        break;
    case CPLXSXP:
        ans = StringFromComplex(COMPLEX(v)[0], &warn);
        break;
    case STRSXP:
        ans = STRING(v)[0];
        break;
    }
    ans = install(CHAR(ans));
    UNPROTECT(1);
    return ans;
}

static SEXP coerceToLogical(SEXP v)
{
    SEXP ans;
    int i, n, warn = 0;
    PROTECT(ans = allocVector(LGLSXP, n = length(v)));
    ATTRIB(ans) = duplicate(ATTRIB(v));
    switch (TYPEOF(v))
    {
    case INTSXP:
        for (i = 0; i < n; i++)
            LOGICAL(ans)[i] = LogicalFromInteger(INTEGER(v)[i], &warn);
        break;
    case REALSXP:
        for (i = 0; i < n; i++)
            LOGICAL(ans)[i] = LogicalFromReal(REAL(v)[i], &warn);
        break;
    case CPLXSXP:
        for (i = 0; i < n; i++)
            LOGICAL(ans)[i] = LogicalFromComplex(COMPLEX(v)[i], &warn);
        break;
    case STRSXP:
        for (i = 0; i < n; i++)
            LOGICAL(ans)[i] = LogicalFromString(STRING(v)[i], &warn);
        break;
    }
    if (warn)
        CoersionWarning(warn);
    UNPROTECT(1);
    return ans;
}

static SEXP coerceToInteger(SEXP v)
{
    SEXP ans;
    int i, n, warn = 0;
    PROTECT(ans = allocVector(INTSXP, n = LENGTH(v)));
    ATTRIB(ans) = duplicate(ATTRIB(v));
    switch (TYPEOF(v))
    {
    case LGLSXP:
        for (i = 0; i < n; i++)
            INTEGER(ans)[i] = IntegerFromLogical(LOGICAL(v)[i], &warn);
        break;
    case REALSXP:
        for (i = 0; i < n; i++)
            INTEGER(ans)[i] = IntegerFromReal(REAL(v)[i], &warn);
        break;
    case CPLXSXP:
        for (i = 0; i < n; i++)
            INTEGER(ans)[i] = IntegerFromComplex(COMPLEX(v)[i], &warn);
        break;
    case STRSXP:
        for (i = 0; i < n; i++)
            INTEGER(ans)[i] = IntegerFromString(STRING(v)[i], &warn);
        break;
    }
    if (warn)
        CoersionWarning(warn);
    UNPROTECT(1);
    return ans;
}

static SEXP coerceToReal(SEXP v)
{
    SEXP ans;
    int i, n, warn = 0;
    PROTECT(ans = allocVector(REALSXP, n = LENGTH(v)));
    ATTRIB(ans) = duplicate(ATTRIB(v));
    switch (TYPEOF(v))
    {
    case LGLSXP:
        for (i = 0; i < n; i++)
            REAL(ans)[i] = RealFromLogical(LOGICAL(v)[i], &warn);
        break;
    case INTSXP:
        for (i = 0; i < n; i++)
            REAL(ans)[i] = RealFromInteger(INTEGER(v)[i], &warn);
        break;
    case CPLXSXP:
        for (i = 0; i < n; i++)
            REAL(ans)[i] = RealFromComplex(COMPLEX(v)[i], &warn);
        break;
    case STRSXP:
        for (i = 0; i < n; i++)
            REAL(ans)[i] = RealFromString(STRING(v)[i], &warn);
        break;
    }
    if (warn)
        CoersionWarning(warn);
    UNPROTECT(1);
    return ans;
}

static SEXP coerceToComplex(SEXP v)
{
    SEXP ans;
    int i, n, warn = 0;
    PROTECT(ans = allocVector(CPLXSXP, n = LENGTH(v)));
    ATTRIB(ans) = duplicate(ATTRIB(v));
    switch (TYPEOF(v))
    {
    case LGLSXP:
        for (i = 0; i < n; i++)
            COMPLEX(ans)[i] = ComplexFromLogical(LOGICAL(v)[i], &warn);
        break;
    case INTSXP:
        for (i = 0; i < n; i++)
            COMPLEX(ans)[i] = ComplexFromInteger(INTEGER(v)[i], &warn);
        break;
    case REALSXP:
        for (i = 0; i < n; i++)
            COMPLEX(ans)[i] = ComplexFromReal(REAL(v)[i], &warn);
        break;
    case STRSXP:
        for (i = 0; i < n; i++)
            COMPLEX(ans)[i] = ComplexFromString(STRING(v)[i], &warn);
        break;
    }
    if (warn)
        CoersionWarning(warn);
    UNPROTECT(1);
    return ans;
}

static SEXP coerceToString(SEXP v)
{
    SEXP ans;
    int i, n, savedigits, warn = 0;
    PROTECT(ans = allocVector(STRSXP, n = LENGTH(v)));
    ATTRIB(ans) = duplicate(ATTRIB(v));
    switch (TYPEOF(v))
    {
    case LGLSXP:
        for (i = 0; i < n; i++)
            STRING(ans)[i] = StringFromLogical(LOGICAL(v)[i], &warn);
        break;
    case INTSXP:
        for (i = 0; i < n; i++)
            STRING(ans)[i] = StringFromInteger(INTEGER(v)[i], &warn);
        break;
    case REALSXP:
        PrintDefaults(R_NilValue);
        savedigits = print_digits;
        print_digits = DBL_DIG; /*- MAXIMAL precision */
        for (i = 0; i < n; i++)
            STRING(ans)[i] = StringFromReal(REAL(v)[i], &warn);
        break;
        print_digits = savedigits;
    case CPLXSXP:
        PrintDefaults(R_NilValue);
        savedigits = print_digits;
        print_digits = DBL_DIG; /*- MAXIMAL precision */
        for (i = 0; i < n; i++)
            STRING(ans)[i] = StringFromComplex(COMPLEX(v)[i], &warn);
        break;
        print_digits = savedigits;
    }
    UNPROTECT(1);
    return (ans);
}

static SEXP coerceToExpression(SEXP v)
{
    SEXP ans, tmp;
    int i, n;
    if (isVectorObject(v))
    {
        n = LENGTH(v);
        PROTECT(ans = allocVector(EXPRSXP, n));
        switch (TYPEOF(v))
        {
        case LGLSXP:
            for (i = 0; i < n; i++)
                VECTOR(ans)[i] = ScalarLogical(LOGICAL(v)[i]);
            break;
        case INTSXP:
            for (i = 0; i < n; i++)
                VECTOR(ans)[i] = ScalarLogical(INTEGER(v)[i]);
            break;
        case REALSXP:
            for (i = 0; i < n; i++)
                VECTOR(ans)[i] = ScalarReal(REAL(v)[i]);
            break;
        case CPLXSXP:
            for (i = 0; i < n; i++)
                VECTOR(ans)[i] = ScalarComplex(COMPLEX(v)[i]);
            break;
        case STRSXP:
            for (i = 0; i < n; i++)
                VECTOR(ans)[i] = ScalarString(STRING(v)[i]);
            break;
        case VECSXP:
            for (i = 0; i < n; i++)
                VECTOR(ans)[i] = VECTOR(v)[i];
            break;
        }
    }
    else
    {
        PROTECT(ans = allocVector(EXPRSXP, 1));
        VECTOR(ans)[0] = duplicate(v);
    }
    UNPROTECT(1);
    return ans;
}

static SEXP coerceToVectorList(SEXP v)
{
    SEXP ans, tmp;
    int i, n;
    n = length(v);
    PROTECT(ans = allocVector(VECSXP, n));
    switch (TYPEOF(v))
    {
    case LGLSXP:
    case INTSXP:
        for (i = 0; i < n; i++)
            VECTOR(ans)[i] = ScalarInteger(INTEGER(v)[i]);
        break;
    case REALSXP:
        for (i = 0; i < n; i++)
            VECTOR(ans)[i] = ScalarReal(REAL(v)[i]);
        break;
    case CPLXSXP:
        for (i = 0; i < n; i++)
            VECTOR(ans)[i] = ScalarComplex(COMPLEX(v)[i]);
        break;
    case STRSXP:
        for (i = 0; i < n; i++)
            VECTOR(ans)[i] = ScalarString(STRING(v)[i]);
        break;
    case LISTSXP:
    case LANGSXP:
        tmp = v;
        for (i = 0; i < n; i++)
        {
            VECTOR(ans)[i] = CAR(tmp);
            tmp = CDR(tmp);
        }
        break;
    default:
        UNIMPLEMENTED("coerceToPairList");
    }
    tmp = getAttrib(v, R_NamesSymbol);
    if (tmp != R_NilValue)
        setAttrib(ans, R_NamesSymbol, tmp);
    UNPROTECT(1);
    return (ans);
}

static SEXP coerceToPairList(SEXP v)
{
    SEXP ans, ansp;
    int i, n;
    n = LENGTH(v);
    PROTECT(ansp = ans = allocList(n));
    for (i = 0; i < n; i++)
    {
        switch (TYPEOF(v))
        {
        case LGLSXP:
            CAR(ansp) = allocVector(LGLSXP, 1);
            INTEGER(CAR(ansp))[0] = INTEGER(v)[i];
            break;
        case INTSXP:
            CAR(ansp) = allocVector(INTSXP, 1);
            INTEGER(CAR(ansp))[0] = INTEGER(v)[i];
            break;
        case REALSXP:
            CAR(ansp) = allocVector(REALSXP, 1);
            REAL(CAR(ansp))[0] = REAL(v)[i];
            break;
        case CPLXSXP:
            CAR(ansp) = allocVector(CPLXSXP, 1);
            COMPLEX(CAR(ansp))[0] = COMPLEX(v)[i];
            break;
        case STRSXP:
            CAR(ansp) = allocVector(STRSXP, 1);
            STRING(CAR(ansp))[0] = STRING(v)[i];
            break;
        case VECSXP:
            CAR(ansp) = VECTOR(v)[i];
            break;
        case EXPRSXP:
            CAR(ansp) = VECTOR(v)[i];
            break;
        default:
            UNIMPLEMENTED("coerceToPairList");
        }
        ansp = CDR(ansp);
    }
    ansp = getAttrib(v, R_NamesSymbol);
    if (ansp != R_NilValue)
        setAttrib(ans, R_NamesSymbol, ansp);
    UNPROTECT(1);
    return (ans);
}

SEXP PairToVectorList(SEXP);

/* Coerce a list to the given type */
static SEXP coercePairList(SEXP v, SEXPTYPE type)
{
    int i, n;
    SEXP rval, vp, names;

    names = v;
    if (type == EXPRSXP)
    {
        PROTECT(rval = allocVector(type, 1));
        VECTOR(rval)[0] = v;
        UNPROTECT(1);
        return rval;
    }
    else if (type == STRSXP)
    {
        n = length(v);
        PROTECT(rval = allocVector(type, n));
        for (vp = v, i = 0; vp != R_NilValue; vp = CDR(vp), i++)
        {
            if (isString(CAR(vp)) && length(CAR(vp)) == 1)
                STRING(rval)[i] = STRING(CAR(vp))[0];
            else
                STRING(rval)[i] = STRING(deparse1(CAR(vp), 0))[0];
        }
    }
    else if (type == VECSXP)
    {
        rval = PairToVectorList(v);
        return rval;
    }
    else if (isVectorizable(v))
    {
        n = length(v);
        PROTECT(rval = allocVector(type, n));
        switch (type)
        {
        case LGLSXP:
            for (i = 0, vp = v; i < n; i++, vp = CDR(vp))
                LOGICAL(rval)[i] = asLogical(CAR(vp));
            break;
        case INTSXP:
            for (i = 0, vp = v; i < n; i++, vp = CDR(vp))
                INTEGER(rval)[i] = asInteger(CAR(vp));
            break;
        case REALSXP:
            for (i = 0, vp = v; i < n; i++, vp = CDR(vp))
                REAL(rval)[i] = asReal(CAR(vp));
            break;
        case CPLXSXP:
            for (i = 0, vp = v; i < n; i++, vp = CDR(vp))
                COMPLEX(rval)[i] = asComplex(CAR(vp));
            break;
        default:
            UNIMPLEMENTED("coerceList");
        }
    }
    else
        error("object cannot be coerced to vector type\n");

    /* If any tags are non-null then we */
    /* need to add a names attribute. */
    for (vp = v, i = 0; vp != R_NilValue; vp = CDR(vp))
        if (TAG(vp) != R_NilValue)
            i = 1;

    if (i)
    {
        i = 0;
        names = allocVector(STRSXP, n);
        for (vp = v; vp != R_NilValue; vp = CDR(vp), i++)
            if (TAG(vp) != R_NilValue)
                STRING(names)[i] = PRINTNAME(TAG(vp));
        setAttrib(rval, R_NamesSymbol, names);
    }
    UNPROTECT(1);
    return rval;
}

static SEXP coerceVectorList(SEXP v, SEXPTYPE type)
{
    int i, n;
    SEXP rval, names;

    names = v;
    if (type == EXPRSXP)
    {
        PROTECT(rval = allocVector(type, 1));
        VECTOR(rval)[0] = v;
        UNPROTECT(1);
        return rval;
    }
    else if (type == STRSXP)
    {
        n = length(v);
        PROTECT(rval = allocVector(type, n));
        for (i = 0; i < n; i++)
        {
            if (isString(VECTOR(v)[i]) && length(VECTOR(v)[i]) == 1)
                STRING(rval)[i] = STRING(VECTOR(v)[i])[0];
            else
                STRING(rval)[i] = STRING(deparse1(VECTOR(v)[i], 0))[0];
        }
    }
    else if (isVectorizable(v))
    {
        n = length(v);
        PROTECT(rval = allocVector(type, n));
        switch (type)
        {
        case LGLSXP:
            for (i = 0; i < n; i++)
                LOGICAL(rval)[i] = asLogical(VECTOR(v)[i]);
            break;
        case INTSXP:
            for (i = 0; i < n; i++)
                INTEGER(rval)[i] = asInteger(VECTOR(v)[i]);
            break;
        case REALSXP:
            for (i = 0; i < n; i++)
                REAL(rval)[i] = asReal(VECTOR(v)[i]);
            break;
        case CPLXSXP:
            for (i = 0; i < n; i++)
                COMPLEX(rval)[i] = asComplex(VECTOR(v)[i]);
            break;
        default:
            UNIMPLEMENTED("coerceList");
        }
    }
    else
        error("object cannot be coerced to vector type\n");

    names = getAttrib(v, R_NamesSymbol);
    if (names != R_NilValue)
        setAttrib(rval, R_NamesSymbol, names);
    UNPROTECT(1);
    return rval;
}

SEXP coerceVector(SEXP v, SEXPTYPE type)
{
    SEXP ans;
    if (TYPEOF(v) == type)
        return v;

    switch (TYPEOF(v))
    {
#ifdef NOTYET
    case NILSXP:
        ans = coerceNull(v, type);
        break;
    case SYMSXP:
        ans = coerceSymbol(v, type);
        break;
#endif
    case NILSXP:
    case LISTSXP:
    case LANGSXP:
        ans = coercePairList(v, type);
        break;
    case VECSXP:
    case EXPRSXP:
        ans = coerceVectorList(v, type);
        break;
    case ENVSXP:
        error("environments cannot be coerced to other types\n");
        break;
    case LGLSXP:
    case INTSXP:
    case REALSXP:
    case CPLXSXP:
    case STRSXP:
        switch (type)
        {
        case SYMSXP:
            ans = coerceToSymbol(v);
            break;
        case LGLSXP:
            ans = coerceToLogical(v);
            break;
        case INTSXP:
            ans = coerceToInteger(v);
            break;
        case REALSXP:
            ans = coerceToReal(v);
            break;
        case CPLXSXP:
            ans = coerceToComplex(v);
            break;
        case STRSXP:
            ans = coerceToString(v);
            break;
        case EXPRSXP:
            ans = coerceToExpression(v);
            break;
        case VECSXP:
            ans = coerceToVectorList(v);
            break;
        case LISTSXP:
            ans = coerceToPairList(v);
            break;
        }
    }
    return ans;
}

SEXP CreateTag(SEXP x)
{
    if (isNull(x) || isSymbol(x))
        return x;
    if (isString(x) && length(x) >= 1 && length(STRING(x)[0]) >= 1)
        x = install(CHAR(STRING(x)[0]));
    else
        x = install(CHAR(STRING(deparse1(x, 1))[0]));
    return x;
}

static SEXP asFunction(SEXP x)
{
    SEXP f, pf;
    int n;
    if (isFunction(x))
        return x;
    PROTECT(f = allocSExp(CLOSXP));
    CLOENV(f) = R_GlobalEnv;
    if (NAMED(x))
        PROTECT(x = duplicate(x));
    else
        PROTECT(x);

    if (isNull(x) || !isList(x))
    {
        FORMALS(f) = R_NilValue;
        BODY(f) = x;
    }
    else
    {
        n = length(x);
        pf = FORMALS(f) = allocList(n - 1);
        while (--n)
        {
            if (TAG(x) == R_NilValue)
            {
                TAG(pf) = CreateTag(CAR(x));
                CAR(pf) = R_MissingArg;
            }
            else
            {
                CAR(pf) = CAR(x);
                TAG(pf) = TAG(x);
            }
            pf = CDR(pf);
            x = CDR(x);
        }
        BODY(f) = CAR(x);
    }
    UNPROTECT(2);
    return f;
}

static SEXP ascommon(SEXP call, SEXP u, int type)
{
    SEXP v;
#ifdef OLD
    if (type == SYMSXP)
    {
        if (TYPEOF(u) == SYMSXP)
            return u;
        if (!isString(u) || LENGTH(u) < 0 || CHAR(STRING(u)[0])[0] == '\0')
            errorcall(call, "character argument required\n");
        return install(CHAR(STRING(u)[0]));
    }
    else
#endif
        if (type == CLOSXP)
    {
        return asFunction(u);
    }
    else if (isVector(u) || isList(u) || isLanguage(u))
    {
        if (NAMED(u))
            v = duplicate(u);
        else
            v = u;
        if (type != ANYSXP)
        {
            PROTECT(v);
            v = coerceVector(v, type);
            UNPROTECT(1);
        }
        if (type == LISTSXP &&
            !(TYPEOF(u) == LANGSXP || TYPEOF(u) == LISTSXP || TYPEOF(u) == EXPRSXP || TYPEOF(u) == VECSXP))
        {
            ATTRIB(v) = R_NilValue;
            OBJECT(v) = 0;
        }
        return v;
    }
    else if (isSymbol(u) && type == STRSXP)
    {
        v = allocVector(STRSXP, 1);
        STRING(v)[0] = PRINTNAME(u);
        return v;
    }
    else
        errorcall(call, "cannot coerce to vector\n");
    return u; /* -Wall */
}

/* as.logical */
/* as.integer */
/* as.real */
/* as.numeric*/
/* as.complex */
/* as.character */
/* as.list */
/* as.expression */
/* as.function */
/* as.name */

SEXP do_as(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    checkArity(op, args);
    return ascommon(call, CAR(args), PRIMVAL(op));
}

SEXP do_asvector(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP ans;
    int type;

    if (DispatchOrEval(call, op, args, rho, &ans, 1))
        return (ans);

    /* Method dispatch has failed, we now just */
    /* run the generic internal code */

    PROTECT(args = ans);
    checkArity(op, args);

    if (!isString(CADR(args)) || LENGTH(CADR(args)) < 1)
        errorcall(call, "invalid \"mode\" of argument\n");

    if (!strcmp("function", (CHAR(STRING(CADR(args))[0]))))
        type = CLOSXP;
    else
        type = str2type(CHAR(STRING(CADR(args))[0]));

    switch (type)
    {
    case SYMSXP:
    case LGLSXP:
    case INTSXP:
    case REALSXP:
    case CPLXSXP:
    case STRSXP:
    case EXPRSXP:
#ifdef NEWLIST
    case VECSXP:
#endif
    case LISTSXP:
    case CLOSXP:
    case ANYSXP:
        break;
    default:
        errorcall(call, "invalid mode\n");
    }
    ans = ascommon(call, CAR(args), type);
    ATTRIB(ans) = R_NilValue;
    UNPROTECT(1);
    return ans;
}

SEXP do_asfunction(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP arglist, envir, names, pargs;
    int i, n;

    checkArity(op, args);

    arglist = CAR(args);
    if (!isNewList(arglist))
        errorcall(call, "list argument expected\n");
    envir = CADR(args);
    if (!isNull(envir) && !isEnvironment(envir))
        errorcall(call, "invalid environment\n");

    n = length(arglist);
    if (n < 1)
        errorcall(call, "argument must have length at least 1\n");
    names = getAttrib(arglist, R_NamesSymbol);
    PROTECT(args = allocList(n - 1));
    for (i = 0, pargs = args; i < n - 1; i++, pargs = CDR(pargs))
    {
        CAR(pargs) = VECTOR(arglist)[i];
        if (names != R_NilValue)
            TAG(pargs) = install(CHAR(STRING(names)[i]));
    }
    CheckFormals(args);
    args = mkCLOSXP(args, VECTOR(arglist)[n - 1], envir);
    UNPROTECT(1);
    return args;
}

SEXP do_ascall(SEXP call, SEXP op, SEXP args, SEXP rho)
{
#ifdef NEWLIST
    SEXP ap, ans, names;
    int i, n;
    checkArity(op, args);
    args = CAR(args);
    switch (TYPEOF(args))
    {
    case LANGSXP:
        ans = args;
        break;
    case VECSXP:
    case EXPRSXP:
        names = getAttrib(args, R_NamesSymbol);
        n = length(args);
        PROTECT(ap = ans = allocList(n));
        for (i = 0; i < n; i++)
        {
            CAR(ap) = VECTOR(args)[i];
            if (names != R_NilValue && !StringBlank(STRING(names)[i]))
                TAG(ap) = install(CHAR(STRING(names)[i]));
            ap = CDR(ap);
        }
        UNPROTECT(1);
        break;
    case LISTSXP:
        ans = duplicate(args);
        break;
    default:
        errorcall(call, "invalid argument list\n");
        ans = R_NilValue;
    }
    TYPEOF(ans) = LANGSXP;
    TAG(ans) = R_NilValue;
    return ans;
#else
    SEXP s;
    checkArity(op, args);
    if (isLanguage(CAR(args)))
        return CAR(args);
    if (!isList(CAR(args)) || length(CAR(args)) < 1)
        errorcall(call, "invalid argument list\n");
    s = duplicate(args);
    TYPEOF(CAR(s)) = LANGSXP;
    TAG(CAR(s)) = R_NilValue;
    return CAR(s);
#endif
}

/* return the type of the SEXP */
SEXP do_typeof(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP ans;
    checkArity(op, args);
    PROTECT(ans = allocVector(STRSXP, 1));
    STRING(ans)[0] = type2str(TYPEOF(CAR(args)));
    UNPROTECT(1);
    return ans;
}

SEXP do_is(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP ans;
    checkArity(op, args);
    PROTECT(ans = allocVector(LGLSXP, 1));
    switch (PRIMVAL(op))
    {
    case NILSXP: /* is.null */
#ifdef NEWLIST
        LOGICAL(ans)[0] = isNull(CAR(args));
#else
        LOGICAL(ans)[0] = (CAR(args) == R_NilValue);
#endif
        break;
    case LGLSXP: /* is.logical */
        LOGICAL(ans)[0] = (TYPEOF(CAR(args)) == LGLSXP);
        break;
    case INTSXP: /* is.integer */
        LOGICAL(ans)[0] = (TYPEOF(CAR(args)) == INTSXP);
        break;
    case REALSXP: /* is.double */
        LOGICAL(ans)[0] = (TYPEOF(CAR(args)) == REALSXP);
        break;
    case CPLXSXP: /* is.complex */
        LOGICAL(ans)[0] = (TYPEOF(CAR(args)) == CPLXSXP);
        break;
    case STRSXP: /* is.character */
        LOGICAL(ans)[0] = (TYPEOF(CAR(args)) == STRSXP);
        break;
    case SYMSXP: /* is.name */
        LOGICAL(ans)[0] = (TYPEOF(CAR(args)) == SYMSXP);
        break;
    case ENVSXP: /* is.environment */
        LOGICAL(ans)[0] = (TYPEOF(CAR(args)) == ENVSXP);
        break;
#ifdef NEWLIST
    case LISTSXP: /* is.list */
        LOGICAL(ans)[0] = ((TYPEOF(CAR(args)) == VECSXP) || TYPEOF(CAR(args)) == NILSXP);
#else
    case LISTSXP: /* is.list */
        LOGICAL(ans)[0] = ((TYPEOF(CAR(args)) == LISTSXP) || TYPEOF(CAR(args)) == NILSXP);
#endif
        break;
    case EXPRSXP: /* is.expression */
        LOGICAL(ans)[0] = (TYPEOF(CAR(args)) == EXPRSXP);
        break;
    case 50: /* is.object */
        LOGICAL(ans)[0] = OBJECT(CAR(args));
        break;
    case 75: /* is.factor */
        LOGICAL(ans)[0] = isFactor(CAR(args));
        break;
    case 80:
        LOGICAL(ans)[0] = isFrame(CAR(args));
        break;

    case 100: /* is.numeric */
        LOGICAL(ans)[0] = (isNumeric(CAR(args)) && !isLogical(CAR(args)));
        break;
    case 101: /* is.matrix */
        LOGICAL(ans)[0] = isMatrix(CAR(args));
        break;
    case 102: /* is.array */
        LOGICAL(ans)[0] = isArray(CAR(args));
        break;
    case 103: /* is.ts */
        LOGICAL(ans)[0] = isTs(CAR(args));
        break;

    case 200: /* is.atomic */
        switch (TYPEOF(CAR(args)))
        {
        case NILSXP:
        case CHARSXP:
        case LGLSXP:
        case INTSXP:
        case REALSXP:
        case CPLXSXP:
        case STRSXP:
            LOGICAL(ans)[0] = 1;
            break;
        default:
            LOGICAL(ans)[0] = 0;
            break;
        }
        break;
    case 201: /* is.recursive */
        switch (TYPEOF(CAR(args)))
        {
        case LISTSXP:
        case CLOSXP:
        case ENVSXP:
        case PROMSXP:
        case LANGSXP:
        case SPECIALSXP:
        case BUILTINSXP:
        case DOTSXP:
        case ANYSXP:
            LOGICAL(ans)[0] = 1;
            break;
        default:
            LOGICAL(ans)[0] = 0;
            break;
        }
        break;

    case 300: /* is.call */
        LOGICAL(ans)[0] = TYPEOF(CAR(args)) == LANGSXP;
        break;
    case 301: /* is.language */
        LOGICAL(ans)
        [0] = (TYPEOF(CAR(args)) == SYMSXP) || (TYPEOF(CAR(args)) == LANGSXP) || (TYPEOF(CAR(args)) == EXPRSXP);
        break;
    case 302: /* is.function */
        LOGICAL(ans)
        [0] = (TYPEOF(CAR(args)) == CLOSXP) || (TYPEOF(CAR(args)) == SPECIALSXP) || (TYPEOF(CAR(args)) == BUILTINSXP);
        break;

    case 999: /* is.single */
        errorcall(call, "type \"single\" unimplemented in R\n");
    default:
        errorcall(call, "unimplemented predicate\n");
    }
    UNPROTECT(1);
    return (ans);
}

/* What should is.vector do ?
 * In S, if an object has no attributes it is a vector, otherwise it isn't.
 * It seems to make more sense to check for a dim attribute.
 */

SEXP do_isvector(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP ans, a;
    checkArity(op, args);
    if (!isString(CADR(args)) || LENGTH(CADR(args)) <= 0)
        errorcall(call, "invalid \"mode\" of argument\n");
    PROTECT(ans = allocVector(LGLSXP, 1));
    if (streql(CHAR(STRING(CADR(args))[0]), "any"))
    {
        LOGICAL(ans)[0] = isVector(CAR(args));
    }
    else if (streql(CHAR(STRING(CADR(args))[0]), CHAR(type2str(TYPEOF(CAR(args))))))
    {
        LOGICAL(ans)[0] = 1;
    }
    else
        LOGICAL(ans)[0] = 0;
    UNPROTECT(1);
    /* We allow a "names" attribute on any vector. */
    if (ATTRIB(CAR(args)) != R_NilValue)
    {
        a = ATTRIB(CAR(args));
        while (a != R_NilValue)
        {
            if (TAG(a) != R_NamesSymbol)
            {
                LOGICAL(ans)[0] = 0;
                return ans;
            }
            a = CDR(a);
        }
    }
    return (ans);
}

SEXP do_isna(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP ans, dims, names, x;
    int i, n;
    if (DispatchOrEval(call, op, args, rho, &ans, 1))
        return (ans);
    PROTECT(args = ans);
    checkArity(op, args);
#ifdef stringent_is
    if (!isList(CAR(args)) && !isVector(CAR(args)))
        errorcall(call, "is.na applies only to lists and vectors\n");
#endif
    ans = allocVector(LGLSXP, length(CAR(args)));
    x = CAR(args);
    n = length(x);
    if (isVector(x))
    {
        PROTECT(dims = getAttrib(x, R_DimSymbol));
        if (isArray(x))
            PROTECT(names = getAttrib(x, R_DimNamesSymbol));
        else
            PROTECT(names = getAttrib(x, R_NamesSymbol));
    }
    else
        dims = names = R_NilValue;
    switch (TYPEOF(x))
    {
    case LGLSXP:
    case INTSXP:
        for (i = 0; i < length(x); i++)
            LOGICAL(ans)[i] = (INTEGER(x)[i] == NA_INTEGER);
        break;
    case REALSXP:
        for (i = 0; i < length(x); i++)
            LOGICAL(ans)[i] = ISNAN(REAL(x)[i]);
        break;
    case CPLXSXP:
        for (i = 0; i < length(x); i++)
            LOGICAL(ans)[i] = (ISNAN(COMPLEX(x)[i].r) || ISNAN(COMPLEX(x)[i].i));
        break;
    case STRSXP:
        for (i = 0; i < length(x); i++)
            LOGICAL(ans)[i] = (STRING(x)[i] == NA_STRING);
        break;
    case LISTSXP:
        n = length(x);
        for (i = 0; i < n; i++)
        {
            if (!isVector(CAR(x)) || length(CAR(x)) > 1)
                LOGICAL(ans)[i] = 0;
            else
            {
                switch (TYPEOF(CAR(x)))
                {
                case LGLSXP:
                case INTSXP:
                    LOGICAL(ans)[i] = (INTEGER(CAR(x))[0] == NA_INTEGER);
                    break;
                case REALSXP:
                    LOGICAL(ans)[i] = ISNAN(REAL(CAR(x))[0]);
                    break;
                case STRSXP:
                    LOGICAL(ans)[i] = (STRING(CAR(x))[0] == NA_STRING);
                    break;
                case CPLXSXP:
                    LOGICAL(ans)[i] = (ISNAN(COMPLEX(CAR(x))[0].r) || ISNAN(COMPLEX(CAR(x))[0].i));
                    break;
                }
            }
            x = CDR(x);
        }
        break;
    default:
        warningcall(call, "is.na() applied to non-(list or vector)\n");
        for (i = 0; i < n; i++)
            LOGICAL(ans)[i] = 0;
    }
    if (dims != R_NilValue)
        setAttrib(ans, R_DimSymbol, dims);
    if (names != R_NilValue)
    {
        if (isArray(x))
            setAttrib(ans, R_DimNamesSymbol, names);
        else
            setAttrib(ans, R_NamesSymbol, names);
    }
    if (isVector(x))
        UNPROTECT(2);
    UNPROTECT(1);
    return ans;
}

SEXP do_isnan(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP ans, dims, names, x;
    int i, n;

    if (DispatchOrEval(call, op, args, rho, &ans, 1))
        return (ans);

    PROTECT(args = ans);
    checkArity(op, args);

#ifdef stringent_is
    if (!isList(CAR(args)) && !isVector(CAR(args)))
        errorcall(call, "is.nan applies only to lists and vectors\n");
#endif
    ans = allocVector(LGLSXP, length(CAR(args)));
    x = CAR(args);
    n = length(x);
    if (isVector(x))
    {
        PROTECT(dims = getAttrib(x, R_DimSymbol));
        if (isArray(x))
            PROTECT(names = getAttrib(x, R_DimNamesSymbol));
        else
            PROTECT(names = getAttrib(x, R_NamesSymbol));
    }
    else
        dims = names = R_NilValue;
    switch (TYPEOF(x))
    {
    case LGLSXP:
    case INTSXP:
    case STRSXP:
        for (i = 0; i < n; i++)
            LOGICAL(ans)[i] = 0;
        break;
    case REALSXP:
        for (i = 0; i < n; i++)
#ifdef IEEE_754
            LOGICAL(ans)[i] = R_IsNaN(REAL(x)[i]);
#else
            LOGICAL(ans)[i] = 0;
#endif
        break;
    case CPLXSXP:
        for (i = 0; i < n; i++)
#ifdef IEEE_754
            LOGICAL(ans)[i] = (R_IsNaN(COMPLEX(x)[i].r) || R_IsNaN(COMPLEX(x)[i].i));
#else
            LOGICAL(ans)[i] = 0;
#endif
        break;
    case LISTSXP:
        for (i = 0; i < n; i++)
        {
            if (!isVector(CAR(x)) || length(CAR(x)) > 1)
                LOGICAL(ans)[i] = 0;
            else
            {
                switch (TYPEOF(CAR(x)))
                {
                case LGLSXP:
                case INTSXP:
                case STRSXP:
                    LOGICAL(ans)[i] = 1;
                    break;
                case REALSXP:
#ifdef IEEE_754
                    LOGICAL(ans)[i] = R_IsNaN(REAL(CAR(x))[0]);
#else
                    LOGICAL(ans)[i] = 0;
#endif
                    break;
                case CPLXSXP:
#ifdef IEEE_754
                    LOGICAL(ans)[i] = (R_IsNaN(COMPLEX(CAR(x))[0].r) || R_IsNaN(COMPLEX(CAR(x))[0].i));
#else
                    LOGICAL(ans)[i] = 0;
#endif
                    break;
                }
            }
            x = CDR(x);
        }
        break;
    default:
        warningcall(call, "is.nan() applied to non-(list or vector)\n");
        for (i = 0; i < n; i++)
            LOGICAL(ans)[i] = 0;
    }
    if (dims != R_NilValue)
        setAttrib(ans, R_DimSymbol, dims);
    if (names != R_NilValue)
    {
        if (isArray(x))
            setAttrib(ans, R_DimNamesSymbol, names);
        else
            setAttrib(ans, R_NamesSymbol, names);
    }
    if (isVector(x))
        UNPROTECT(2);
    UNPROTECT(1);
    return ans;
}

SEXP do_isfinite(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP ans, x, names, dims;
    int i, n;
    checkArity(op, args);
#ifdef stringent_is
    if (!isList(CAR(args)) && !isVector(CAR(args)))
        errorcall(call, "is.finite applies only to lists and vectors\n");
#endif
    x = CAR(args);
    n = length(x);
    ans = allocVector(LGLSXP, n);
    if (isVector(x))
    {
        dims = getAttrib(x, R_DimSymbol);
        if (isArray(x))
            names = getAttrib(x, R_DimNamesSymbol);
        else
            names = getAttrib(x, R_NamesSymbol);
    }
    else
        dims = names = R_NilValue;
    switch (TYPEOF(x))
    {
    case LGLSXP:
    case INTSXP:
        for (i = 0; i < n; i++)
            LOGICAL(ans)[i] = (INTEGER(x)[i] != NA_INTEGER);
        break;
    case REALSXP:
        for (i = 0; i < n; i++)
            LOGICAL(ans)[i] = FINITE(REAL(x)[i]);
        break;
    case CPLXSXP:
        for (i = 0; i < n; i++)
            LOGICAL(ans)[i] = (FINITE(COMPLEX(x)[i].r) && FINITE(COMPLEX(x)[i].i));
        break;
    default:
        for (i = 0; i < n; i++)
            LOGICAL(ans)[i] = 0;
    }
    if (dims != R_NilValue)
        setAttrib(ans, R_DimSymbol, dims);
    if (names != R_NilValue)
    {
        if (isArray(x))
            setAttrib(ans, R_DimNamesSymbol, names);
        else
            setAttrib(ans, R_NamesSymbol, names);
    }
    return ans;
}

SEXP do_isinfinite(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP ans, x, names, dims;
    double xr, xi;
    int i, n;
    checkArity(op, args);
#ifdef stringent_is
    if (!isList(CAR(args)) && !isVector(CAR(args)))
        errorcall(call, "is.infinite applies only to list and vectors\n");
#endif
    x = CAR(args);
    n = length(x);
    ans = allocVector(LGLSXP, n);
    if (isVector(x))
    {
        dims = getAttrib(x, R_DimSymbol);
        if (isArray(x))
            names = getAttrib(x, R_DimNamesSymbol);
        else
            names = getAttrib(x, R_NamesSymbol);
    }
    else
        dims = names = R_NilValue;
#ifdef IEEE_754
    switch (TYPEOF(x))
    {
    case REALSXP:
        for (i = 0; i < n; i++)
        {
            xr = REAL(x)[i];
            if (xr != xr /*NaN*/ || FINITE(xr))
                LOGICAL(ans)[i] = 0;
            else
                LOGICAL(ans)[i] = 1;
        }
        break;
    case CPLXSXP:
        for (i = 0; i < n; i++)
        {
            xr = COMPLEX(x)[i].r;
            xi = COMPLEX(x)[i].i;
            if ((xr != xr || FINITE(xr)) && (xi != xi || FINITE(xi)))
                LOGICAL(ans)[i] = 0;
            else
                LOGICAL(ans)[i] = 1;
        }
        break;
    default:
        for (i = 0; i < n; i++)
            LOGICAL(ans)[i] = 0;
    }
#else
    for (i = 0; i < n; i++)
        LOGICAL(ans)[i] = 0;
#endif
    if (!isNull(dims))
        setAttrib(ans, R_DimSymbol, dims);
    if (!isNull(names))
    {
        if (isArray(x))
            setAttrib(ans, R_DimNamesSymbol, names);
        else
            setAttrib(ans, R_NamesSymbol, names);
    }
    return ans;
}

SEXP do_call(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP rest, evargs, rfun;

    PROTECT(rfun = eval(CAR(args), rho));
    if (!isString(rfun) || length(rfun) <= 0 || streql(CHAR(STRING(rfun)[0]), ""))
        errorcall(call, "first argument must be a character string\n");
    PROTECT(rfun = install(CHAR(STRING(rfun)[0])));
    PROTECT(evargs = duplicate(CDR(args)));
    for (rest = evargs; rest != R_NilValue; rest = CDR(rest))
        CAR(rest) = eval(CAR(rest), rho);
    rfun = LCONS(rfun, evargs);
    UNPROTECT(3);
    return (rfun);
}

SEXP do_docall(SEXP call, SEXP op, SEXP args, SEXP rho)
{
#ifdef NEWLIST
    SEXP c;
    int i, n;
#endif
    checkArity(op, args);
    if (!isString(CAR(args)) || length(CAR(args)) < 0 || streql(CHAR(STRING(CAR(args))[0]), ""))
        errorcall(call, "first argument must be a string\n");
#ifdef NEWLIST
    if (!isNull(CADR(args)) && !isNewList(CADR(args)))
        errorcall(call, "second argument must be a list\n");
    n = length(CADR(args));
    PROTECT(call = allocList(n + 1));
    CAR(call) = install(CHAR(STRING(CAR(args))[0]));
    args = CADR(args);
    c = CDR(call);
    for (i = 0; i < n; i++)
    {
        CAR(c) = VECTOR(args)[i];
        c = CDR(c);
    }
    UNPROTECT(1);
#else
    if (!isNull(CADR(args)) && !isList(CADR(args)))
        errorcall(call, "second argument must be a list\n");
    call = install(CHAR(STRING(CAR(args))[0]));
    PROTECT(call = LCONS(call, CADR(args)));
    call = eval(call, rho);
    UNPROTECT(1);
#endif
    return call;
}

/* do_substitute has two arguments, an expression and an environment */
/* (optional).  Symbols found in the expression are substituted with their */
/* values as found in the environment.  There is no inheritance so only */
/* the supplied environment is searched. If no environment is specified */
/* the environment in which substitute was called is used.  If the */
/* specified environment is R_NilValue then R_GlobalEnv is used. */
/* Arguments to do_substitute should not be evaluated. */

SEXP substituteList(SEXP, SEXP);

SEXP substitute(SEXP lang, SEXP rho)
{
    SEXP t;
    switch (TYPEOF(lang))
    {
    case PROMSXP:
        return substitute(PREXPR(lang), rho);
    case SYMSXP:
        t = findVarInFrame(FRAME(rho), lang);
        if (t != R_UnboundValue)
        {
            if (TYPEOF(t) == PROMSXP)
            {
                do
                {
                    t = PREXPR(t);
                } while (TYPEOF(t) == PROMSXP);
                return t;
#ifdef OLD
                return substitute(PREXPR(t), rho);
                return PREXPR(t);
#endif
            }
            else if (TYPEOF(t) == DOTSXP)
            {
                error("... used in an incorrect context\n");
            }
            if (rho != R_GlobalEnv)
                return t;
        }
        return (lang);
    case LANGSXP:
        return substituteList(lang, rho);
    default:
        return (lang);
    }
}

/* Work through a list doing substitute on the */
/* elements taking particular care to handle ... */

SEXP substituteList(SEXP el, SEXP rho)
{
    SEXP h, t;
    if (isNull(el))
        return el;
    if (CAR(el) == R_DotsSymbol)
    {
        h = findVar(CAR(el), rho);
        if (h == R_NilValue)
            return substituteList(CDR(el), rho);
        if (TYPEOF(h) != DOTSXP)
        {
            if (h == R_MissingArg)
                return substituteList(CDR(el), rho);
            error("... used in an incorrect context\n");
        }
        PROTECT(h = substituteList(h, rho));
        PROTECT(t = substituteList(CDR(el), rho));
        t = listAppend(h, t);
        UNPROTECT(2);
        return t;
    }
#ifdef OLD
    else if (CAR(el) == R_MissingArg)
    {
        return substituteList(CDR(el), rho);
    }
#endif
    else
    {
        PROTECT(h = substitute(CAR(el), rho));
        PROTECT(t = substituteList(CDR(el), rho));
        if (isLanguage(el))
            t = LCONS(h, t);
        else
            t = CONS(h, t);
        TAG(t) = TAG(el);
        UNPROTECT(2);
        return t;
    }
}

SEXP do_substitute(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP env, s, t;
    /* set up the environment for substitution */
    if (length(args) == 1)
        env = rho;
    else
        env = eval(CADR(args), rho);
    if (env == R_NilValue)
        env = R_GlobalEnv;
    if (TYPEOF(env) == LISTSXP)
    {
        PROTECT(s = duplicate(env));
        PROTECT(env = allocSExp(ENVSXP));
        FRAME(env) = s;
        UNPROTECT(2);
    }
    if (TYPEOF(env) != ENVSXP)
        errorcall(call, "invalid environment specified\n");

    PROTECT(env);
    PROTECT(t = duplicate(args));
    CDR(t) = R_NilValue;
    s = substituteList(t, env);
    UNPROTECT(2);
    return CAR(s);
}
