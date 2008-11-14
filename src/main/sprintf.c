/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2002-8     the R Development Core Team
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
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 * Originally written by Jonathan Rougier
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Defn.h>
#include "RBufferUtils.h"

#define MAXLINE MAXELTSIZE

/*
   This is passed a format that started with % and may include other
   chars, e.g. '.2f abc'.  It's aim is to show that this is a valid
   format from one of the types given in pattern.
*/

static const char *findspec(const char *str)
{
    const char *p = str;

    if (*p != '%')
        return p;
    for (p++;; p++)
    {
        if (*p == '-' || *p == '+' || *p == ' ' || *p == '.')
            continue;
        if ((*p == '*' || (*p >= '0' && *p <= '9')) && *(p + 1) == '.')
        {
            p++;
            continue;
        }
        if (*p == '.' && (*(p + 1) == '*' || (*(p + 1) >= '0' && *(p + 1) <= '9')))
        {
            p++;
            continue;
        }
        if (*p == '*' || (*p >= '0' && *p <= '9'))
        {
            if (*(p + 1) != '.')
                continue;
            if (*(p + 2) == '*' || (*(p + 2) >= '0' && *(p + 2) <= '9'))
            {
                /* so have m.n pattern */
                p += 2; /* for loop will skip third byte */
                continue;
            }
        }
        break;
    }
    return p;
}

/*   FALSE is success, TRUE is an error. */
static Rboolean checkfmt(const char *fmt, const char *pattern)
{
    const char *p = fmt;

    if (*p != '%')
        return 1;
    p = findspec(fmt);
    return strcspn(p, pattern) != 0;
}

SEXP attribute_hidden do_sprintf(SEXP call, SEXP op, SEXP args, SEXP env)
{
    int i, nargs, cnt, v, thislen, nfmt, nprotect = 1;
    const char *formatString, *ss;
    char *starc;
    /* fmt2 is a copy of fmt with '*' expanded.
       bit will hold numeric formats and %<w>s, so be quite small. */
    char fmt[MAXLINE + 1], fmt2[MAXLINE + 10], *fmtp, bit[MAXLINE + 1], *outputString;
    size_t n, cur, chunk;

    SEXP format, ans, _this, a[100], tmp;
    int ns, maxlen, lens[100], nthis, has_star, star_arg = 0, nstar;
    static R_StringBuffer outbuff = {NULL, 0, MAXELTSIZE};
    Rboolean use_UTF8;

    outputString = R_AllocStringBuffer(0, &outbuff);

    /* grab the format string */
    nargs = length(args);
    format = CAR(args);
    if (!isString(format) || length(format) == 0)
        error(_("'fmt' is not a non-empty character vector"));
    args = CDR(args);
    nargs--;
    if (nargs >= 100)
        error(_("only 100 arguments are allowed"));

    /* record the args for possible coercion and later re-ordering */
    for (i = 0; i < nargs; i++, args = CDR(args))
        a[i] = CAR(args);

    maxlen = nfmt = length(format);
    for (i = 0; i < nargs; i++)
    {
        lens[i] = length(a[i]);
        if (lens[i] == 0)
            error(_("zero-length argument"));
        if (maxlen < lens[i])
            maxlen = lens[i];
    }
    if (maxlen % length(format))
        error(_("arguments cannot be recycled to the same length"));
    for (i = 0; i < nargs; i++)
        if (maxlen % lens[i])
            error(_("arguments cannot be recycled to the same length"));

    /* We do the format analysis a row at a time */
    PROTECT(ans = allocVector(STRSXP, maxlen));
    for (ns = 0; ns < maxlen; ns++)
    {
        outputString[0] = '\0';
        use_UTF8 = getCharCE(STRING_ELT(format, ns % nfmt)) == CE_UTF8;
        if (!use_UTF8)
        {
            for (i = 0; i < nargs; i++)
            {
                if (!isString(a[i]))
                    continue;
                if (getCharCE(STRING_ELT(a[i], ns % lens[i])) == CE_UTF8)
                {
                    use_UTF8 = TRUE;
                    break;
                }
            }
        }
        if (use_UTF8)
            formatString = translateCharUTF8(STRING_ELT(format, ns % nfmt));
        else
            formatString = translateChar(STRING_ELT(format, ns % nfmt));
        n = strlen(formatString);
        if (n > MAXLINE)
            error(_("'fmt' length exceeds maximal format length %d"), MAXLINE);
        /* process the format string */
        for (cur = 0, cnt = 0; cur < n; cur += chunk)
        {
            ss = NULL;
            if (formatString[cur] == '%')
            { /* handle special format command */

                if (cur < n - 1 && formatString[cur + 1] == '%')
                {
                    /* take care of %% in the format */
                    chunk = 2;
                    strcpy(bit, "%");
                }
                else
                {
                    /* recognise selected types from Table B-1 of K&R */
                    /* This is MBCS-OK, as we are in a format spec */
                    chunk = strcspn(formatString + cur + 1, "aAdisfeEgGxX%") + 2;
                    if (cur + chunk > n)
                        error(_("unrecognised format at end of string"));

                    strncpy(fmt, formatString + cur, chunk);
                    fmt[chunk] = '\0';

                    nthis = -1;
                    /* now look for %n$ or %nn$ form */
                    if (strlen(fmt) > 3 && fmt[1] >= '1' && fmt[1] <= '9')
                    {
                        v = fmt[1] - '0';
                        if (fmt[2] == '$')
                        {
                            if (v > nargs)
                                error(_("reference to non-existent argument %d"), v);
                            nthis = v - 1;
                            memmove(fmt + 1, fmt + 3, strlen(fmt) - 2);
                        }
                        else if (fmt[2] >= '1' && fmt[2] <= '9' && fmt[3] == '$')
                        {
                            v = 10 * v + fmt[2] - '0';
                            if (v > nargs)
                                error(_("reference to non-existent argument %d"), v);
                            nthis = v - 1;
                            memmove(fmt + 1, fmt + 4, strlen(fmt) - 3);
                        }
                    }

                    has_star = 0;
                    starc = Rf_strchr(fmt, '*');
                    if (starc)
                    { /* handle * format if present */
                        nstar = -1;
                        if (strlen(starc) > 3 && starc[1] >= '1' && starc[1] <= '9')
                        {
                            v = starc[1] - '0';
                            if (starc[2] == '$')
                            {
                                if (v > nargs)
                                    error(_("reference to non-existent argument %d"), v);
                                nstar = v - 1;
                                memmove(starc + 1, starc + 3, strlen(starc) - 2);
                            }
                            else if (starc[2] >= '1' && starc[2] <= '9' && starc[3] == '$')
                            {
                                v = 10 * v + starc[2] - '0';
                                if (v > nargs)
                                    error(_("reference to non-existent argument %d"), v);
                                nstar = v - 1;
                                memmove(starc + 1, starc + 4, strlen(starc) - 3);
                            }
                        }

                        if (nstar < 0)
                        {
                            if (cnt >= nargs)
                                error(_("too few arguments"));
                            nstar = cnt++;
                        }

                        if (Rf_strchr(starc + 1, '*'))
                            error(_("at most one asterisk '*' is supported in each conversion specification"));

                        _this = a[nstar];
                        if (ns == 0 && TYPEOF(_this) == REALSXP)
                        {
                            _this = coerceVector(_this, INTSXP);
                            PROTECT(a[nstar] = _this);
                            nprotect++;
                        }
                        if (TYPEOF(_this) != INTSXP || LENGTH(_this) < 1 ||
                            INTEGER(_this)[ns % LENGTH(_this)] == NA_INTEGER)
                            error(_("argument for '*' conversion specification must be a number"));
                        has_star = 1;
                        star_arg = INTEGER(_this)[ns % LENGTH(_this)];
                    }

                    if (fmt[strlen(fmt) - 1] == '%')
                    {
                        /* handle % with formatting options */
                        if (has_star)
                            sprintf(bit, fmt, star_arg);
                        else
                            sprintf(bit, fmt);
                    }
                    else
                    {
                        if (nthis < 0)
                        {
                            if (cnt >= nargs)
                                error(_("too few arguments"));
                            nthis = cnt++;
                        }
                        _this = a[nthis];
                        if (has_star)
                        {
                            char *p, *q = fmt2;
                            for (p = fmt; *p; p++)
                                if (*p == '*')
                                    q += sprintf(q, "%d", star_arg);
                                else
                                    *q++ = *p;
                            *q = '\0';
                            fmtp = fmt2;
                        }
                        else
                            fmtp = fmt;

                        if (ns == 0)
                        {
                            /* Now let us see if some minimal coercion
                               would be sensible, but only do so once. */
                            switch (*findspec(fmtp))
                            {
                            case 'd':
                            case 'i':
                            case 'x':
                            case 'X':
                                if (TYPEOF(_this) == REALSXP)
                                {
                                    double r = REAL(_this)[0];
                                    if ((double)((int)r) == r)
                                        _this = coerceVector(_this, INTSXP);
                                    PROTECT(a[nthis] = _this);
                                    nprotect++;
                                }
                                break;
                            case 'a':
                            case 'A':
                            case 'e':
                            case 'f':
                            case 'g':
                            case 'E':
                            case 'G':
                                if (TYPEOF(_this) != REALSXP)
                                {
                                    PROTECT(tmp = lang2(install("as.double"), _this));
                                    _this = eval(tmp, env);
                                    UNPROTECT(1);
                                    PROTECT(a[nthis] = _this);
                                    nprotect++;
                                }
                                break;
                            case 's':
                                if (TYPEOF(_this) != STRSXP)
                                {
                                    PROTECT(tmp = lang2(install("as.character"), _this));
                                    _this = eval(tmp, env);
                                    UNPROTECT(1);
                                    PROTECT(a[nthis] = _this);
                                    nprotect++;
                                }
                                break;
                            default:
                                break;
                            }
                        }

                        PROTECT(_this);
                        thislen = length(_this);
                        if (thislen == 0)
                            error(_("coercion has changed vector length to 0"));

                        switch (TYPEOF(_this))
                        {
                        case LGLSXP: {
                            int x = LOGICAL(_this)[ns % thislen];
                            if (checkfmt(fmtp, "di"))
                                error("%s", _("use format %d or %i for logical objects"));
                            if (x == NA_LOGICAL)
                            {
                                fmtp[strlen(fmtp) - 1] = 's';
                                sprintf(bit, fmtp, "NA");
                            }
                            else
                            {
                                sprintf(bit, fmtp, x);
                            }
                            break;
                        }
                        case INTSXP: {
                            int x = INTEGER(_this)[ns % thislen];
                            if (checkfmt(fmtp, "dixX"))
                                error("%s", _("use format %d, %i, %x or %X for integer objects"));
                            if (x == NA_INTEGER)
                            {
                                fmtp[strlen(fmtp) - 1] = 's';
                                sprintf(bit, fmtp, "NA");
                            }
                            else
                            {
                                sprintf(bit, fmtp, x);
                            }
                            break;
                        }
                        case REALSXP: {
                            double x = REAL(_this)[ns % thislen];
                            if (checkfmt(fmtp, "aAfeEgG"))
                                error("%s", _("use format %f, %e, %g or %a for numeric objects"));
                            if (R_FINITE(x))
                            {
                                sprintf(bit, fmtp, x);
                            }
                            else
                            {
                                char *p = Rf_strchr(fmtp, '.');
                                if (p)
                                {
                                    *p++ = 's';
                                    *p = '\0';
                                }
                                else
                                    fmtp[strlen(fmtp) - 1] = 's';
                                if (ISNA(x))
                                {
                                    if (strcspn(fmtp, " ") < strlen(fmtp))
                                        sprintf(bit, fmtp, " NA");
                                    else
                                        sprintf(bit, fmtp, "NA");
                                }
                                else if (ISNAN(x))
                                {
                                    if (strcspn(fmtp, " ") < strlen(fmtp))
                                        sprintf(bit, fmtp, " NaN");
                                    else
                                        sprintf(bit, fmtp, "NaN");
                                }
                                else if (x == R_PosInf)
                                {
                                    if (strcspn(fmtp, "+") < strlen(fmtp))
                                        sprintf(bit, fmtp, "+Inf");
                                    else if (strcspn(fmtp, " ") < strlen(fmtp))
                                        sprintf(bit, fmtp, " Inf");
                                    else
                                        sprintf(bit, fmtp, "Inf");
                                }
                                else if (x == R_NegInf)
                                    sprintf(bit, fmtp, "-Inf");
                            }
                            break;
                        }
                        case STRSXP:
                            /* NA_STRING will be printed as 'NA' */
                            if (checkfmt(fmtp, "s"))
                                error("%s", _("use format %s for character objects"));
                            if (use_UTF8)
                                ss = translateCharUTF8(STRING_ELT(_this, ns % thislen));
                            else
                                ss = translateChar(STRING_ELT(_this, ns % thislen));
                            if (fmtp[1] != 's')
                            {
                                if (strlen(ss) > MAXLINE)
                                    warning(_("likely truncation of character string to %d characters"), MAXLINE - 1);
                                snprintf(bit, MAXLINE, fmtp, ss);
                                bit[MAXLINE] = '\0';
                                ss = NULL;
                            }
                            break;

                        default:
                            error(_("unsupported type"));
                            break;
                        }

                        UNPROTECT(1);
                    }
                }
            }
            else
            {                                                  /* not '%' : handle string part */
                char *ch = Rf_strchr(formatString + cur, '%'); /* MBCS-aware
                                          version used */
                if (ch)
                    chunk = ch - formatString - cur;
                else
                    chunk = strlen(formatString + cur);
                strncpy(bit, formatString + cur, chunk);
                bit[chunk] = '\0';
            }
            if (ss)
            {
                outputString = R_AllocStringBuffer(strlen(outputString) + strlen(ss) + 1, &outbuff);
                strcat(outputString, ss);
            }
            else
            {
                outputString = R_AllocStringBuffer(strlen(outputString) + strlen(bit) + 1, &outbuff);
                strcat(outputString, bit);
            }
        }
        SET_STRING_ELT(ans, ns, mkCharCE(outputString, use_UTF8 ? CE_UTF8 : CE_NATIVE));
    }

    UNPROTECT(nprotect);
    R_FreeStringBufferL(&outbuff);
    return ans;
}

/* Local Variables: */
/* indent-tabs-mode: t */
/* c-basic-offset: 4 */
/* End: */
