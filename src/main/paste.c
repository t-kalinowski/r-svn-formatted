/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997--2005  The R Development Core Team
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
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 *  See ./printutils.c	 for general remarks on Printing
 *                       and the Encode.. utils.
 *
 *  See ./format.c	 for the  format_Foo_  functions.
 */

/* <UTF8> char here is handled as a whole */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Defn.h"
#include "Rmath.h" /* imax2 */
#include "Print.h"

/*  .Internal(paste(args, sep, collapse))
 *
 * do_paste uses two passes to paste the arguments (in CAR(args)) together.
 * The first pass calculates the width of the paste buffer,
 * then it is alloc-ed and the second pass stuffs the information in.
 */
SEXP do_paste(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans, collapse, sep, x, tmpchar;
    int i, j, k, maxlen, nx, pwidth, sepw;
    char *s, *buf;

    checkArity(op, args);

    /* We use formatting and so we */
    /* must initialize printing. */

    PrintDefaults(env);

    /* Check the arguments */

    x = CAR(args);
    if (!isVectorList(x))
        errorcall(call, _("invalid first argument"));

    sep = CADR(args);
    if (!isString(sep) || LENGTH(sep) <= 0)
        errorcall(call, _("invalid separator"));
    sep = STRING_ELT(sep, 0);
    sepw = strlen(CHAR(sep)); /* not LENGTH as might contain \0 */

    collapse = CADDR(args);
    if (!isNull(collapse))
        if (!isString(collapse) || LENGTH(collapse) <= 0)
            errorcall(call, _("invalid '%s' argument"), "collapse");

    /* Maximum argument length and */
    /* check for arguments of list type */

    nx = length(x);
    maxlen = 0;
    for (j = 0; j < nx; j++)
    {
        if (!isString(VECTOR_ELT(x, j)))
            error(_("non-string argument to Internal paste"));
        if (length(VECTOR_ELT(x, j)) > maxlen)
            maxlen = length(VECTOR_ELT(x, j));
    }
    if (maxlen == 0)
        return (!isNull(collapse)) ? mkString("") : allocVector(STRSXP, 0);

    PROTECT(ans = allocVector(STRSXP, maxlen));

    for (i = 0; i < maxlen; i++)
    {
        pwidth = 0;
        for (j = 0; j < nx; j++)
        {
            k = length(VECTOR_ELT(x, j));
            if (k > 0)
                pwidth += strlen(CHAR(STRING_ELT(VECTOR_ELT(x, j), i % k)));
        }
        pwidth += (nx - 1) * sepw;
        tmpchar = allocString(pwidth);
        buf = CHAR(tmpchar);
        for (j = 0; j < nx; j++)
        {
            k = length(VECTOR_ELT(x, j));
            if (k > 0)
            {
                s = CHAR(STRING_ELT(VECTOR_ELT(x, j), i % k));
                strcpy(buf, s);
                buf += strlen(s);
            }
            if (j != nx - 1 && sepw != 0)
            {
                strcpy(buf, CHAR(sep));
                buf += sepw;
            }
        }
        SET_STRING_ELT(ans, i, tmpchar);
    }

    /* Now collapse, if required. */

    if (collapse != R_NilValue && (nx = LENGTH(ans)) != 0)
    {
        sep = STRING_ELT(collapse, 0);
        sepw = strlen(CHAR(sep));
        pwidth = 0;
        for (i = 0; i < nx; i++)
            pwidth += strlen(CHAR(STRING_ELT(ans, i)));
        pwidth += (nx - 1) * sepw;
        tmpchar = allocString(pwidth);
        buf = CHAR(tmpchar);
        for (i = 0; i < nx; i++)
        {
            if (i > 0)
            {
                strcpy(buf, CHAR(sep));
                buf += sepw;
            }
            s = CHAR(STRING_ELT(ans, i));
            strcpy(buf, s);
            while (*buf)
                buf++;
        }
        PROTECT(tmpchar);
        ans = allocVector(STRSXP, 1);
        UNPROTECT(1);
        SET_STRING_ELT(ans, 0, tmpchar);
    }
    UNPROTECT(1);
    return ans;
}

/* format.default(x, trim, nsmall) : ../library/base/R/format.R
 * --------------   See "FIXME" in that file !
 */
SEXP do_format(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP l, x, y, swd;
    int i, il, n, digits, trim = 0, nsmall = 0, wd = 0, adj = -1;
    int w, d, e;
    int wi, di, ei;
    char *strp;

    checkArity(op, args);
    PrintDefaults(env);

    if (!isVector(x = CAR(args)))
        errorcall(call, _("first argument must be atomic"));
    args = CDR(args);

    trim = asLogical(CAR(args));
    if (trim == NA_INTEGER)
        errorcall(call, _("invalid '%s' argument"), "trim");
    args = CDR(args);

    if (!isNull(CAR(args)))
    {
        digits = asInteger(CAR(args));
        if (digits == NA_INTEGER || digits < 0)
            errorcall(call, _("invalid '%s' argument"), "digits");
        R_print.digits = digits;
    }
    args = CDR(args);

    nsmall = asInteger(CAR(args));
    if (nsmall == NA_INTEGER || nsmall < 0 || nsmall > 20)
        errorcall(call, _("invalid '%s' argument"), "nsmall");
    args = CDR(args);

    if (isNull(swd = CAR(args)))
        wd = 0;
    else
        wd = asInteger(swd);
    if (wd == NA_INTEGER)
        errorcall(call, _("invalid '%s' argument"), "width");
    args = CDR(args);

    adj = asInteger(CAR(args));
    if (adj == NA_INTEGER || adj < -1 || adj > 2)
        errorcall(call, _("invalid '%s' argument"), "justify");

    if ((n = LENGTH(x)) <= 0)
        return allocVector(STRSXP, 0);

    switch (TYPEOF(x))
    {

    case LGLSXP:
        PROTECT(y = allocVector(STRSXP, n));
        if (trim)
            w = 0;
        else
            formatLogical(LOGICAL(x), n, &w);
        w = imax2(w, wd);
        for (i = 0; i < n; i++)
        {
            strp = EncodeLogical(LOGICAL(x)[i], w);
            SET_STRING_ELT(y, i, mkChar(strp));
        }
        UNPROTECT(1);
        break;

    case INTSXP:
        PROTECT(y = allocVector(STRSXP, n));
        if (trim)
            w = 0;
        else
            formatInteger(INTEGER(x), n, &w);
        w = imax2(w, wd);
        for (i = 0; i < n; i++)
        {
            strp = EncodeInteger(INTEGER(x)[i], w);
            SET_STRING_ELT(y, i, mkChar(strp));
        }
        UNPROTECT(1);
        break;

    case REALSXP:
        formatReal(REAL(x), n, &w, &d, &e, nsmall);
        if (trim)
            w = 0;
        w = imax2(w, wd);
        PROTECT(y = allocVector(STRSXP, n));
        for (i = 0; i < n; i++)
        {
            strp = EncodeReal(REAL(x)[i], w, d, e, OutDec);
            SET_STRING_ELT(y, i, mkChar(strp));
        }
        UNPROTECT(1);
        break;

    case CPLXSXP:
        formatComplex(COMPLEX(x), n, &w, &d, &e, &wi, &di, &ei, nsmall);
        if (trim)
            wi = w = 0;
        w = imax2(w, wd);
        wi = imax2(wi, wd);
        PROTECT(y = allocVector(STRSXP, n));
        for (i = 0; i < n; i++)
        {
            strp = EncodeComplex(COMPLEX(x)[i], w, d, e, wi, di, ei, OutDec);
            SET_STRING_ELT(y, i, mkChar(strp));
        }
        UNPROTECT(1);
        break;

    case STRSXP: {
        /* this has to be different from formatString/EncodeString as
           we don't actually want to encode here */
        char *s, *buff, *q;
        int b, b0, cnt = 0, j;

        w = 0;
        if (adj >= 0)
            for (i = 0; i < n; i++)
                if (STRING_ELT(x, i) != NA_STRING)
                {
                    il = Rstrlen(STRING_ELT(x, i), 0);
                    if (il > w)
                        w = il;
                }
        w = imax2(w, wd);
        /* now calculate the buffer size needed, in bytes */
        for (i = 0; i < n; i++)
            if (STRING_ELT(x, i) != NA_STRING)
            {
                il = Rstrlen(STRING_ELT(x, i), 0);
                cnt = imax2(cnt, LENGTH(STRING_ELT(x, i)) + imax2(0, w - il));
            }
        buff = alloca(cnt + 1);
        PROTECT(y = allocVector(STRSXP, n));
        for (i = 0; i < n; i++)
        {
            if (STRING_ELT(x, i) == NA_STRING)
            {
                SET_STRING_ELT(y, i, NA_STRING);
            }
            else
            {
                q = buff;
                s = CHAR(STRING_ELT(x, i));
                cnt = LENGTH(STRING_ELT(x, i));
                il = Rstrlen(STRING_ELT(x, i), 0);
                b = w - il;
                if (b > 0 && adj != Rprt_adj_left)
                {
                    b0 = (adj == Rprt_adj_centre) ? b / 2 : b;
                    for (j = 0; j < b0; j++)
                        *q++ = ' ';
                    b -= b0;
                }
                for (j = 0; j < cnt; j++)
                    *q++ = *s++;
                if (b > 0 && adj != Rprt_adj_right)
                {
                    for (j = 0; j < b; j++)
                        *q++ = ' ';
                }
                *q = '\0';
                SET_STRING_ELT(y, i, mkChar(buff));
            }
        }
        UNPROTECT(1);
    }
    break;
    default:
        errorcall(call, _("Impossible mode ( x )"));
        y = R_NilValue; /* -Wall */
    }
    PROTECT(y);
    if ((l = getAttrib(x, R_DimSymbol)) != R_NilValue)
    {
        setAttrib(y, R_DimSymbol, l);
        if ((l = getAttrib(x, R_DimNamesSymbol)) != R_NilValue)
            setAttrib(y, R_DimNamesSymbol, l);
    }
    else if ((l = getAttrib(x, R_NamesSymbol)) != R_NilValue)
        setAttrib(y, R_NamesSymbol, l);

    UNPROTECT(1);
    return y;
}

/* format.info(obj)  --> 3 integers  (w,d,e) with the formatting information
 *			w = total width (#{chars}) per item
 *			d = #{digits} to RIGHT of "."
 *			e = {0:2}.   0: Fixpoint;
 *				   1,2: exponential with 2/3 digit expon.
 *
 * for complex : 2 x 3 integers for (Re, Im)
 */

SEXP do_formatinfo(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP x;
    int n, digits, nsmall, no = 1, w, d, e, wi, di, ei;

    checkArity(op, args);
    x = CAR(args);
    n = LENGTH(x);
    PrintDefaults(env);

    digits = asInteger(CADR(args));
    if (!isNull(CADR(args)))
    {
        digits = asInteger(CADR(args));
        if (digits == NA_INTEGER || digits < 0)
            errorcall(call, _("invalid '%s' argument"), "digits");
        R_print.digits = digits;
    }
    nsmall = asInteger(CADDR(args));
    if (nsmall == NA_INTEGER || nsmall < 0 || nsmall > 20)
        errorcall(call, _("invalid '%s' argument"), "nsmall");

    w = 0;
    d = 0;
    e = 0;
    switch (TYPEOF(x))
    {

    case RAWSXP:
        formatRaw(RAW(x), n, &w);
        break;

    case LGLSXP:
        formatLogical(LOGICAL(x), n, &w);
        break;

    case INTSXP:
        formatInteger(INTEGER(x), n, &w);
        break;

    case REALSXP:
        no = 3;
        formatReal(REAL(x), n, &w, &d, &e, nsmall);
        break;

    case CPLXSXP:
        no = 6;
        wi = di = ei = 0;
        formatComplex(COMPLEX(x), n, &w, &d, &e, &wi, &di, &ei, nsmall);
        break;

    case STRSXP: {
        int i, il;
        for (i = 0; i < n; i++)
            if (STRING_ELT(x, i) != NA_STRING)
            {
                il = Rstrlen(STRING_ELT(x, i), 0);
                if (il > w)
                    w = il;
            }
    }
    break;

    default:
        errorcall(call, _("atomic vector arguments only"));
    }
    x = allocVector(INTSXP, no);
    INTEGER(x)[0] = w;
    if (no > 1)
    {
        INTEGER(x)[1] = d;
        INTEGER(x)[2] = e;
    }
    if (no > 3)
    {
        INTEGER(x)[3] = wi;
        INTEGER(x)[4] = di;
        INTEGER(x)[5] = ei;
    }
    return x;
}
