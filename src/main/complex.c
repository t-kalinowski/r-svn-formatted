/*
 *  R : A Computer Language for Statistical Data Analysis
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
#include "Mathlib.h"
#include "Fortran.h" /*- for POW_DI */

static int naflag;

SEXP complex_unary(int code, SEXP s1)
{
    int i, n;
    complex x;
    SEXP ans;

    switch (code)
    {
    case PLUSOP:
        return s1;
    case MINUSOP:
        ans = duplicate(s1);
        n = LENGTH(s1);
        for (i = 0; i < n; i++)
        {
            x = COMPLEX(s1)[i];
#ifdef IEEE_745
            COMPLEX(ans)[i].r = -x.r;
            COMPLEX(ans)[i].i = -x.i;
#else
            if (ISNAN(x.r) || ISNAN(x.i))
            {
                COMPLEX(ans)[i].r = NA_REAL;
                COMPLEX(ans)[i].i = NA_REAL;
            }
            else
            {
                COMPLEX(ans)[i].r = -x.r;
                COMPLEX(ans)[i].i = -x.i;
            }
#endif
        }
        return ans;
    default:
        error("illegal complex unary operator\n");
        return R_NilValue; /* -Wall*/
    }
}

static void complex_div(complex *c, complex *a, complex *b)
{
    double ratio, den;
    double abr, abi;

    if ((abr = b->r) < 0)
        abr = -abr;
    if ((abi = b->i) < 0)
        abi = -abi;
    if (abr <= abi)
    {
#ifndef IEEE_754
        if (abi == 0)
        {
            c->r = NA_REAL;
            c->i = NA_REAL;
            return;
        }
#endif
        ratio = b->r / b->i;
        den = b->i * (1 + ratio * ratio);
        c->r = (a->r * ratio + a->i) / den;
        c->i = (a->i * ratio - a->r) / den;
    }
    else
    {
        ratio = b->i / b->r;
        den = b->r * (1 + ratio * ratio);
        c->r = (a->r + a->i * ratio) / den;
        c->i = (a->i - a->r * ratio) / den;
    }
}

static void complex_pow(complex *r, complex *a, complex *b)
{
    double logr, logi, x, y;
    int ib;

    if (b->i == 0.)
    { /* be fast (and more accurate)*/
        if (b->r == 1.)
        { /* a^1 */
            r->r = a->r;
            r->i = a->i;
            return;
        }
        if (a->i == 0.)
        {
            r->r = pow(a->r, b->r);
            r->i = 0.;
            return;
        }
        if (a->r == 0. && b->r == (ib = (int)b->r))
        { /* (|a|*i)^b */
            x = POW_DI(&(a->i), &ib);
            if (ib % 2)
            { /* ib is odd ==> imaginary r */
                r->r = 0.;
                r->i = ((ib > 0 && ib % 4 == 3) || (ib < 0 && (-ib) % 4 == 1)) ? -x : x;
            }
            else
            { /* even exponent b : real r */
                r->r = (ib % 4) ? -x : x;
                r->i = 0.;
            }
            return;
        }
    }
    logr = log(hypot(a->r, a->i));
    logi = atan2(a->i, a->r);
    x = exp(logr * b->r - logi * b->i);
    y = logr * b->i + logi * b->r;
    r->r = x * cos(y);
    r->i = x * sin(y);
}

/* FIXME : Use the trick in arithmetic.c to eliminate "modulo" ops */

SEXP complex_binary(int code, SEXP s1, SEXP s2)
{
    int i, n, n1, n2;
    complex x1, x2;
    SEXP ans;
    n1 = LENGTH(s1);
    n2 = LENGTH(s2);
    n = (n1 > n2) ? n1 : n2;
    /* Note: "s1" and "s1" are protected in the calling code. */
    ans = allocVector(CPLXSXP, n);
    if (n1 < 1 || n2 < 1)
    {
        for (i = 0; i < n; i++)
        {
            COMPLEX(ans)[i].r = NA_REAL;
            COMPLEX(ans)[i].i = NA_REAL;
        }
        return ans;
    }
    switch (code)
    {
    case PLUSOP:
        for (i = 0; i < n; i++)
        {
            x1 = COMPLEX(s1)[i % n1];
            x2 = COMPLEX(s2)[i % n2];
#ifdef IEEE_754
            COMPLEX(ans)[i].r = x1.r + x2.r;
            COMPLEX(ans)[i].i = x1.i + x2.i;
#else
            if (ISNAN(x1.r) || ISNAN(x1.i) || ISNAN(x2.r) || ISNAN(x2.i))
            {
                COMPLEX(ans)[i].r = NA_REAL;
                COMPLEX(ans)[i].i = NA_REAL;
            }
            else
            {
                COMPLEX(ans)[i].r = MATH_CHECK(x1.r + x2.r);
                COMPLEX(ans)[i].i = MATH_CHECK(x1.i + x2.i);
            }
#endif
        }
        break;
    case MINUSOP:
        for (i = 0; i < n; i++)
        {
            x1 = COMPLEX(s1)[i % n1];
            x2 = COMPLEX(s2)[i % n2];
#ifdef IEEE_754
            COMPLEX(ans)[i].r = x1.r - x2.r;
            COMPLEX(ans)[i].i = x1.i - x2.i;
#else
            if (ISNAN(x1.r) || ISNAN(x1.i) || ISNAN(x2.r) || ISNAN(x2.i))
            {
                COMPLEX(ans)[i].r = NA_REAL;
                COMPLEX(ans)[i].i = NA_REAL;
            }
            else
            {
                COMPLEX(ans)[i].r = MATH_CHECK(x1.r - x2.r);
                COMPLEX(ans)[i].i = MATH_CHECK(x1.i - x2.i);
            }
#endif
        }
        break;
    case TIMESOP:
        for (i = 0; i < n; i++)
        {
            x1 = COMPLEX(s1)[i % n1];
            x2 = COMPLEX(s2)[i % n2];
#ifdef IEEE_754
            COMPLEX(ans)[i].r = x1.r * x2.r - x1.i * x2.i;
            COMPLEX(ans)[i].i = x1.r * x2.i + x1.i * x2.r;
#else
            if (ISNAN(x1.r) || ISNAN(x1.i) || ISNAN(x2.r) || ISNAN(x2.i))
            {
                COMPLEX(ans)[i].r = NA_REAL;
                COMPLEX(ans)[i].i = NA_REAL;
            }
            else
            {
                COMPLEX(ans)[i].r = MATH_CHECK(x1.r * x2.r - x1.i * x2.i);
                COMPLEX(ans)[i].i = MATH_CHECK(x1.r * x2.i + x1.i * x2.r);
            }
#endif
        }
        break;
    case DIVOP:
        for (i = 0; i < n; i++)
        {
            x1 = COMPLEX(s1)[i % n1];
            x2 = COMPLEX(s2)[i % n2];
#ifdef IEEE_754
            complex_div(&COMPLEX(ans)[i], &x1, &x2);
#else
            if (ISNAN(x1.r) || ISNAN(x1.i) || ISNAN(x2.r) || ISNAN(x2.i))
            {
                COMPLEX(ans)[i].r = NA_REAL;
                COMPLEX(ans)[i].i = NA_REAL;
            }
            else
                complex_div(&COMPLEX(ans)[i], &x1, &x2);
#endif
        }
        break;
    case POWOP:
        for (i = 0; i < n; i++)
        {
            x1 = COMPLEX(s1)[i % n1];
            x2 = COMPLEX(s2)[i % n2];
#ifdef IEEE_754
            complex_pow(&COMPLEX(ans)[i], &x1, &x2);
#else
            if (ISNAN(x1.r) || ISNAN(x1.i) || ISNAN(x2.r) || ISNAN(x2.i))
            {
                COMPLEX(ans)[i].r = NA_REAL;
                COMPLEX(ans)[i].i = NA_REAL;
            }
            else
                complex_pow(&COMPLEX(ans)[i], &x1, &x2);
#endif
        }
        break;
    default:
        error("unimplemented complex operation\n");
    }
    /* Copy attributes from longest argument. */
    if (n1 > n2)
        copyMostAttrib(s1, ans);
    else if (n1 == n2)
    {
        copyMostAttrib(s2, ans);
        copyMostAttrib(s1, ans);
    }
    else
        copyMostAttrib(s2, ans);
    return ans;
}

/* FIXME : Use the trick in arithmetic.c to eliminate "modulo" ops */

SEXP do_cmathfuns(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP x, y = R_NilValue; /* -Wall*/
    int i, n;

    checkArity(op, args);
    x = CAR(args);
    n = length(x);
    if (isComplex(x))
    {
        switch (PRIMVAL(op))
        {
        case 1: /* Re */
            y = allocVector(REALSXP, n);
            for (i = 0; i < n; i++)
                REAL(y)[i] = COMPLEX(x)[i].r;
            break;
        case 2: /* Im */
            y = allocVector(REALSXP, n);
            for (i = 0; i < n; i++)
                REAL(y)[i] = COMPLEX(x)[i].i;
            break;
        case 3: /* Mod */
            y = allocVector(REALSXP, n);
            for (i = 0; i < n; i++)
            {
#ifdef IEEE_754
                REAL(y)[i] = hypot(COMPLEX(x)[i].r, COMPLEX(x)[i].i);
#else
                if (ISNAN(COMPLEX(x)[i].r) || ISNAN(COMPLEX(x)[i].i))
                {
                    REAL(y)[i] = NA_REAL;
                }
                else
                {
                    REAL(y)[i] = hypot(COMPLEX(x)[i].r, COMPLEX(x)[i].i);
                }
#endif
            }
            break;
        case 4: /* Arg */
            y = allocVector(REALSXP, n);
            for (i = 0; i < n; i++)
            {
#ifdef IEEE_754
                REAL(y)[i] = atan2(COMPLEX(x)[i].i, COMPLEX(x)[i].r);
#else
                if (ISNAN(COMPLEX(x)[i].r) || ISNAN(COMPLEX(x)[i].i))
                {
                    REAL(y)[i] = NA_REAL;
                }
                else
                {
                    REAL(y)[i] = atan2(COMPLEX(x)[i].i, COMPLEX(x)[i].r);
                }
#endif
            }
            break;
        case 5: /* Conj */
            y = allocVector(CPLXSXP, n);
            for (i = 0; i < n; i++)
            {
#ifdef IEEE_754
                COMPLEX(y)[i].r = COMPLEX(x)[i].r;
                COMPLEX(y)[i].i = -COMPLEX(x)[i].i;
#else
                if (ISNAN(COMPLEX(x)[i].r) || ISNAN(COMPLEX(x)[i].i))
                {
                    COMPLEX(y)[i].r = NA_REAL;
                    COMPLEX(y)[i].i = NA_REAL;
                }
                else
                {
                    COMPLEX(y)[i].r = COMPLEX(x)[i].r;
                    COMPLEX(y)[i].i = -COMPLEX(x)[i].i;
                }
#endif
            }
            break;
        }
    }
    else if (isNumeric(x))
    {
        if (isReal(x))
            PROTECT(x);
        else
            PROTECT(x = coerceVector(x, REALSXP));
        switch (PRIMVAL(op))
        {
        case 1: /* Re */
        case 5: /* Conj */
            y = allocVector(REALSXP, n);
            for (i = 0; i < n; i++)
                REAL(y)[i] = REAL(x)[i];
            break;
        case 2: /* Im */
        case 4: /* Arg */
            y = allocVector(REALSXP, n);
            for (i = 0; i < n; i++)
                if (ISNAN(REAL(x)[i]))
                    REAL(y)[i] = REAL(x)[i];
                else
                    REAL(y)[i] = 0;
            break;
        case 3: /* Mod */
            y = allocVector(REALSXP, n);
            for (i = 0; i < n; i++)
            {
#ifdef IEEE_754
                REAL(y)[i] = fabs(REAL(x)[i]);
#else
                if (ISNAN(REAL(x)[i]))
                    REAL(y)[i] = REAL(x)[i];
                else
                    REAL(y)[i] = fabs(REAL(x)[i]);
#endif
            }
            break;
        }
        UNPROTECT(1);
    }
    else
        errorcall(call, "non-numeric argument to function\n");
    PROTECT(x);
    PROTECT(y);
    ATTRIB(y) = duplicate(ATTRIB(x));
    OBJECT(y) = OBJECT(x);
    UNPROTECT(2);
    return y;
}

static void z_rround(complex *r, complex *x, complex *p)
{
    r->r = rround(x->r, p->r);
    r->i = rround(x->i, p->r);
}

/* Question:  This treats real and imaginary parts */
/* separately.	Should it do them jointly? */

static void z_prec(complex *r, complex *x, complex *p)
{
    r->r = prec(x->r, p->r);
    r->i = prec(x->i, p->r);
}

static void z_log(complex *r, complex *z)
{
    r->i = atan2(z->i, z->r);
    r->r = log(hypot(z->r, z->i));
}

static void z_logbase(complex *r, complex *z, complex *base)
{
    complex t1, t2;
    z_log(&t1, z);
    z_log(&t2, base);
    complex_div(r, &t1, &t2);
}

static void z_exp(complex *r, complex *z)
{
    double expx;
    expx = exp(z->r);
    r->r = expx * cos(z->i);
    r->i = expx * sin(z->i);
}

static void z_sqrt(complex *r, complex *z)
{
    double mag;

    if ((mag = hypot(z->r, z->i)) == 0.0)
        r->r = r->i = 0.0;
    else if (z->r > 0)
    {
        r->r = sqrt(0.5 * (mag + z->r));
        r->i = z->i / r->r / 2;
    }
    else
    {
        r->i = sqrt(0.5 * (mag - z->r));
        if (z->i < 0)
            r->i = -r->i;
        r->r = z->i / r->i / 2;
    }
}

static void z_cos(complex *r, complex *z)
{
    r->r = cos(z->r) * cosh(z->i);
    r->i = -sin(z->r) * sinh(z->i);
}

static void z_sin(complex *r, complex *z)
{
    r->r = sin(z->r) * cosh(z->i);
    r->i = cos(z->r) * sinh(z->i);
}

static void z_tan(complex *r, complex *z)
{
    double x2, y2, den;
    x2 = 2.0 * z->r;
    y2 = 2.0 * z->i;
    den = cos(x2) + cosh(y2);
    r->r = sin(x2) / den;
    r->i = sinh(y2) / den;
}

/* Complex Arcsin and Arccos Functions */
/* Equation (4.4.37) Abramowitz and Stegun */

static void z_asin(complex *r, complex *z)
{
    double alpha, beta, t1, t2, x, y;
    x = z->r;
    y = z->i;
    t1 = 0.5 * sqrt((x + 1) * (x + 1) + y * y);
    t2 = 0.5 * sqrt((x - 1) * (x - 1) + y * y);
    alpha = t1 + t2;
    beta = t1 - t2;
    r->r = asin(beta);
    r->i = log(alpha + sqrt(alpha * alpha - 1));
}

static void z_acos(complex *r, complex *z)
{
    static double pi2 = 1.57079632679489661923;
    complex asin;
    z_asin(&asin, z);
    r->r = pi2 - asin.r;
    r->i = -asin.i;
}

/* Complex Arctangent Function */
/* Equation (4.4.39) Abramowitz and Stegun */

static void z_atan(complex *r, complex *z)
{
    double x, y;
    x = z->r;
    y = z->i;
    r->r = 0.5 * atan(2 * x / (1 - x * x - y * y));
    r->i = 0.25 * log((x * x + (y + 1) * (y + 1)) / (x * x + (y - 1) * (y - 1)));
}

static void z_atan2(complex *r, complex *csn, complex *ccs)
{
    static double pi = 3.14159265358979323846;
    complex tmp;
    if (ccs->r == 0 && ccs->i == 0)
    {
        if (csn->r == 0 && csn->r == 0)
        {
            r->r = NA_REAL;
            r->i = NA_REAL;
        }
        else
        {
            r->r = fsign(0.5 * pi, csn->r);
            r->i = 0;
        }
    }
    else
    {
        complex_div(&tmp, csn, ccs);
        z_atan(r, &tmp);
        if (ccs->r < 0)
            r->r += pi;
        if (r->r > pi)
            r->r -= 2 * pi;
    }
}

static void z_acosh(complex *r, complex *z)
{
    complex a;
    z_acos(&a, z);
    r->r = -a.i;
    r->i = a.r;
}

static void z_asinh(complex *r, complex *z)
{
    complex a, b;
    b.r = -z->i;
    b.i = z->r;
    z_asin(&a, &b);
    r->r = a.i;
    r->i = -a.r;
}

static void z_atanh(complex *r, complex *z)
{
    complex a, b;
    b.r = -z->i;
    b.i = z->r;
    z_atan(&a, &b);
    r->r = a.i;
    r->i = -a.r;
}

static void z_cosh(complex *r, complex *z)
{
    complex a;
    a.r = -z->i;
    a.i = z->r;
    z_cos(r, &a);
}

static void z_sinh(complex *r, complex *z)
{
    complex a, b;
    b.r = -z->i;
    b.i = z->r;
    z_sin(&a, &b);
    r->r = a.i;
    r->i = -a.r;
}

static void z_tanh(complex *r, complex *z)
{
    complex a, b;
    b.r = -z->i;
    b.i = z->r;
    z_tan(&a, &b);
    r->r = a.i;
    r->i = -a.r;
}

static void cmath1(void (*f)(), complex *x, complex *y, int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        if (ISNA(x[i].r) || ISNA(x[i].i))
        {
            y[i].r = NA_REAL;
            y[i].i = NA_REAL;
        }
        else
        {
            f(&y[i], &x[i]);
#ifndef IEEE_754
            if (ISNA(y[i].r) || ISNA(y[i].i))
            {
                y[i].r = NA_REAL;
                y[i].i = NA_REAL;
                naflag = 1;
            }
#endif
        }
    }
}

SEXP complex_math1(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP x, y;
    int n;
    x = CAR(args);
    n = length(x);
    y = allocVector(CPLXSXP, n);
    naflag = 0;

    switch (PRIMVAL(op))
    {
    case 10002:
        cmath1(z_atan, COMPLEX(x), COMPLEX(y), n);
        break;
    case 10003:
        cmath1(z_log, COMPLEX(x), COMPLEX(y), n);
        break;

    case 0:
        errorcall(call, "abs() unimplemented for complex; use Mod()\n");
        break;
    case 3:
        cmath1(z_sqrt, COMPLEX(x), COMPLEX(y), n);
        break;

    case 10:
        cmath1(z_exp, COMPLEX(x), COMPLEX(y), n);
        break;

    case 20:
        cmath1(z_cos, COMPLEX(x), COMPLEX(y), n);
        break;
    case 21:
        cmath1(z_sin, COMPLEX(x), COMPLEX(y), n);
        break;
    case 22:
        cmath1(z_tan, COMPLEX(x), COMPLEX(y), n);
        break;
    case 23:
        cmath1(z_acos, COMPLEX(x), COMPLEX(y), n);
        break;
    case 24:
        cmath1(z_asin, COMPLEX(x), COMPLEX(y), n);
        break;

    case 30:
        cmath1(z_cosh, COMPLEX(x), COMPLEX(y), n);
        break;
    case 31:
        cmath1(z_sinh, COMPLEX(x), COMPLEX(y), n);
        break;
    case 32:
        cmath1(z_tanh, COMPLEX(x), COMPLEX(y), n);
        break;
    case 33:
        cmath1(z_acosh, COMPLEX(x), COMPLEX(y), n);
        break;
    case 34:
        cmath1(z_asinh, COMPLEX(x), COMPLEX(y), n);
        break;
    case 35:
        cmath1(z_atanh, COMPLEX(x), COMPLEX(y), n);
        break;

#ifdef NOTYET
        MATH1(40, lgammafn);
        MATH1(41, gammafn);
#endif

    default:
        errorcall(call, "unimplemented complex function\n");
    }
    if (naflag)
        warning("NAs produced in function \"%s\"\n", PRIMNAME(op));
    return y;
}

/* FIXME : Use the trick in arithmetic.c to eliminate "modulo" ops */

static SEXP cmath2(SEXP op, SEXP sa, SEXP sb, void (*f)())
{
    int i, n, na, nb;
    complex ai, bi, *a, *b, *y;
    SEXP sy;

    na = length(sa);
    nb = length(sb);
    n = (na < nb) ? nb : na;
    PROTECT(sa = coerceVector(sa, CPLXSXP));
    PROTECT(sb = coerceVector(sb, CPLXSXP));
    PROTECT(sy = allocVector(CPLXSXP, n));
    a = COMPLEX(sa);
    b = COMPLEX(sb);
    y = COMPLEX(sy);
    if (na < 1 || na < 1)
    {
        for (i = 0; i < n; i++)
        {
            y[i].r = NA_REAL;
            y[i].i = NA_REAL;
        }
    }
    else
    {
        naflag = 0;
        for (i = 0; i < n; i++)
        {
            ai = a[i % na];
            bi = b[i % nb];
            if (ISNA(ai.r) && ISNA(ai.i) && ISNA(bi.r) && ISNA(bi.i))
            {
                y[i].r = NA_REAL;
                y[i].i = NA_REAL;
            }
            else
            {
                f(&y[i], &ai, &bi);
#ifndef IEEE_754
                if (ISNA(y[i].r) || ISNA(y[i].i))
                {
                    y[i].r = NA_REAL;
                    y[i].i = NA_REAL;
                    naflag = 1;
                }
#endif
            }
        }
    }
    if (naflag)
        warning("NAs produced in function \"%s\"\n", PRIMNAME(op));
    if (n == na)
    {
        ATTRIB(sy) = duplicate(ATTRIB(sa));
        OBJECT(sy) = OBJECT(sa);
    }
    else if (n == nb)
    {
        ATTRIB(sy) = duplicate(ATTRIB(sb));
        OBJECT(sy) = OBJECT(sb);
    }
    UNPROTECT(3);
    return sy;
}

/* Complex Functions of Two Arguments */

SEXP complex_math2(SEXP call, SEXP op, SEXP args, SEXP env)
{
    switch (PRIMVAL(op))
    {
    case 10001:
        return cmath2(op, CAR(args), CADR(args), z_rround);
    case 10002:
        return cmath2(op, CAR(args), CADR(args), z_atan2);
    case 10003:
        return cmath2(op, CAR(args), CADR(args), z_logbase);
    case 10004:
        return cmath2(op, CAR(args), CADR(args), z_prec);
    case 0:
        return cmath2(op, CAR(args), CADR(args), z_atan2);
    default:
        errorcall(call, "unimplemented complex function\n");
        return call; /* just for -Wall */
    }
}

SEXP do_complex(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    /* complex(length, real, imaginary) */
    SEXP ans, re, im;
    int i, na, nr, ni;
    na = asInteger(CAR(args));
    if (na == NA_INTEGER || na < 0)
        errorcall(call, "invalid length\n");
    PROTECT(re = coerceVector(CADR(args), REALSXP));
    PROTECT(im = coerceVector(CADDR(args), REALSXP));
    nr = length(re);
    ni = length(im);
    /* is always true: if (na >= 0) {*/
    na = (nr > na) ? nr : na;
    na = (ni > na) ? ni : na;
    /* }*/
    ans = allocVector(CPLXSXP, na);
    for (i = 0; i < na; i++)
    {
        COMPLEX(ans)[i].r = 0;
        COMPLEX(ans)[i].i = 0;
    }
    UNPROTECT(2);
    if (na > 0 && nr > 0)
    {
        for (i = 0; i < na; i++)
            COMPLEX(ans)[i].r = REAL(re)[i % nr];
    }
    if (na > 0 && ni > 0)
    {
        for (i = 0; i < na; i++)
            COMPLEX(ans)[i].i = REAL(im)[i % ni];
    }
    return ans;
}

int F77_SYMBOL(cpoly)(double *, double *, int *, double *, double *, int *);

SEXP do_polyroot(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP z, zr, zi, r, rr, ri;
    int degree, fail, i, n;
    checkArity(op, args);
    z = CAR(args);
    switch (TYPEOF(z))
    {
    case CPLXSXP:
        PROTECT(z);
        break;
    case REALSXP:
    case INTSXP:
    case LGLSXP:
        PROTECT(z = coerceVector(z, CPLXSXP));
        break;
    default:
        errorcall(call, "invalid argument type\n");
    }
    n = length(z);
    degree = n - 1;
    if (degree >= 1)
    {
        if (n > 49)
            errorcall(call, "polynomial degree too high (49 max)\n");
        /* <==>	 #define NMAX 50  in  ../appl/cpoly.c */

        if (COMPLEX(z)[n - 1].r == 0.0 && COMPLEX(z)[n - 1].i == 0.0)
            errorcall(call, "highest power has coefficient 0\n");

        PROTECT(rr = allocVector(REALSXP, n));
        PROTECT(ri = allocVector(REALSXP, n));
        PROTECT(zr = allocVector(REALSXP, n));
        PROTECT(zi = allocVector(REALSXP, n));

        for (i = 0; i < n; i++)
        {
            if (!FINITE(COMPLEX(z)[i].r) || !FINITE(COMPLEX(z)[i].i))
                errorcall(call, "invalid polynomial coefficient\n");
            REAL(zr)[degree - i] = COMPLEX(z)[i].r;
            REAL(zi)[degree - i] = COMPLEX(z)[i].i;
        }
        F77_SYMBOL(cpoly)(REAL(zr), REAL(zi), &degree, REAL(rr), REAL(ri), &fail);
        if (fail)
            errorcall(call, "root finding code failed\n");
        UNPROTECT(2);
        r = allocVector(CPLXSXP, degree);
        for (i = 0; i < n; i++)
        {
            COMPLEX(r)[i].r = REAL(rr)[i];
            COMPLEX(r)[i].i = REAL(ri)[i];
        }
        UNPROTECT(3);
    }
    else
    {
        UNPROTECT(1);
        r = allocVector(CPLXSXP, 0);
    }
    return r;
}
