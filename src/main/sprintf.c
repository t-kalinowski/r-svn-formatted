/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2002     the R Development Core Team
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
 * Written by Jonathan Rougier, email J.C.Rougier@durham.ac.uk
 */

/* <UTF8> char here is handled as a whole string */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Defn.h>

#define MAXLINE MAXELTSIZE

/* Simple wrapper for C sprintf function: now (1.6.0) checks the
   types and handles the R specials.

   We make no attempt to allow fmt to be MBCS.
*/
SEXP do_sprintf(SEXP call, SEXP op, SEXP args, SEXP env)
{
    int nargs;
    char *formatString;
    char fmt[MAXLINE + 1], bit[MAXLINE + 1], outputString[MAXLINE + 1] = "";
    size_t n, cur, chunk;

    SEXP format, ans;

    /* grab the format string */

    nargs = length(args);
    format = CAR(args);
    if (!isString(format) || LENGTH(format) != 1)
        errorcall(call, _("'fmt' is not a character string of length 1"));
    formatString = CHAR(STRING_ELT(format, 0));
    n = strlen(formatString);
    if (n > MAXLINE)
        errorcall(call, _("string length exceeds maximal buffer length %d"), MAXLINE);

    /* process the format string */
    for (cur = 0; cur < n; cur += chunk)
    {

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

                chunk = strcspn(formatString + cur, "disfeEgG") + 1;
                if (cur + chunk > n)
                    errorcall(call, _("unrecognised format at end of string"));

                strncpy(fmt, formatString + cur, chunk);
                fmt[chunk] = '\0';

                if (--nargs > 0)
                    args = CDR(args);
                else
                    errorcall(call, _("too few arguments"));

                if (LENGTH(CAR(args)) < 1)
                    error(_("zero-length argument"));
                switch (TYPEOF(CAR(args)))
                {
                case LGLSXP: {
                    int x = LOGICAL(CAR(args))[0];
                    if (strcspn(fmt, "di") >= strlen(fmt))
                        error("%s", _("use format %d or %i for logical objects"));
                    if (x == NA_LOGICAL)
                    {
                        fmt[chunk - 1] = 's';
                        sprintf(bit, fmt, "NA");
                    }
                    else
                        sprintf(bit, fmt, x);
                    break;
                }
                case INTSXP: {
                    int x = INTEGER(CAR(args))[0];
                    if (strcspn(fmt, "di") >= strlen(fmt))
                        error("%s", _("use format %d or %i for integer objects"));
                    if (x == NA_INTEGER)
                    {
                        fmt[chunk - 1] = 's';
                        sprintf(bit, fmt, "NA");
                    }
                    else
                        sprintf(bit, fmt, x);
                    break;
                }
                case REALSXP: {
                    double x = REAL(CAR(args))[0];
                    if (strcspn(fmt, "feEgG") >= strlen(fmt))
                        error("%s", _("use format %f, %e or %g for numeric objects"));
                    if (R_FINITE(x))
                    {
                        sprintf(bit, fmt, x);
                    }
                    else
                    {
                        char *p = strchr(fmt, '.');
                        if (p)
                        {
                            *p++ = 's';
                            *p = '\0';
                        }
                        else
                            fmt[chunk - 1] = 's';
                        if (ISNA(x))
                        {
                            if (strcspn(fmt, " ") < strlen(fmt))
                                sprintf(bit, fmt, " NA");
                            else
                                sprintf(bit, fmt, "NA");
                        }
                        else if (ISNAN(x))
                        {
                            if (strcspn(fmt, " ") < strlen(fmt))
                                sprintf(bit, fmt, " NaN");
                            else
                                sprintf(bit, fmt, "NaN");
                        }
                        else if (x == R_PosInf)
                        {
                            if (strcspn(fmt, "+") < strlen(fmt))
                                sprintf(bit, fmt, "+Inf");
                            else if (strcspn(fmt, " ") < strlen(fmt))
                                sprintf(bit, fmt, " Inf");
                            else
                                sprintf(bit, fmt, "Inf");
                        }
                        else if (x == R_NegInf)
                            sprintf(bit, fmt, "-Inf");
                    }
                    break;
                }
                case STRSXP:
                    /* NA_STRING will be printed as `NA' */
                    if (strcspn(fmt, "s") >= strlen(fmt))
                        error("%s", _("use format %s for character objects"));
                    if (strlen(CHAR(STRING_ELT(CAR(args), 0))) > MAXLINE)
                        warning(_("Likely truncation of character string"));
                    snprintf(bit, MAXLINE, fmt, CHAR(STRING_ELT(CAR(args), 0)));
                    bit[MAXLINE] = '\0';
                    break;
                default:
                    errorcall(call, _("unsupported type"));
                    break;
                }
            }
        }
        else
        { /* not '%' : handle string part */

            chunk = strcspn(formatString + cur, "%");
            strncpy(bit, formatString + cur, chunk);
            bit[chunk] = '\0';
        }

        if (strlen(outputString) + strlen(bit) > MAXLINE)
            errorcall(call, _("String length exceeds buffer size of %d"), MAXLINE);
        strcat(outputString, bit);
    }

    /* return outputString as SEXP */

    if (nargs > 1)
        warning(_("Unused arguments"));

    PROTECT(ans = allocVector(STRSXP, 1));
    SET_STRING_ELT(ans, 0, mkChar(outputString));
    UNPROTECT(1);
    return ans;
}
