/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995-1997, 1998  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998-2000  the R Development Core Team.
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
 *  EXPORTS	printVector()
 *		printNamedVector()
 *		printRealVector()
 *		printIntegerVector()
 *		printComplexVector()
 *
 *  See ./printutils.c	 for remarks on Printing and the Encoding utils.
 *  See ./format.c	 for the formatXXXX functions used below.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Print.h"

#define DO_first_lab                                                                                                   \
    if (indx)                                                                                                          \
    {                                                                                                                  \
        labwidth = IndexWidth(n) + 2;                                                                                  \
        /* labwidth may well be                                                                                        \
           one more than desired ..*/                                                                                  \
        VectorIndex(1, labwidth);                                                                                      \
        width = labwidth;                                                                                              \
    }                                                                                                                  \
    else                                                                                                               \
        width = 0

#define DO_newline                                                                                                     \
    Rprintf("\n");                                                                                                     \
    if (indx)                                                                                                          \
    {                                                                                                                  \
        VectorIndex(i + 1, labwidth);                                                                                  \
        width = labwidth;                                                                                              \
    }                                                                                                                  \
    else                                                                                                               \
        width = 0

void printLogicalVector(int *x, int n, int indx)
{
    int i, w, labwidth = 0, width;

    DO_first_lab;
    formatLogical(x, n, &w);
    w += R_print.gap;

    for (i = 0; i < n; i++)
    {
        if (width + w > R_print.width)
        {
            DO_newline;
        }
        Rprintf("%s", EncodeLogical(x[i], w));
        width += w;
    }
    Rprintf("\n");
}

void printIntegerVector(int *x, int n, int indx)
{
    int i, w, labwidth = 0, width;

    DO_first_lab;
    formatInteger(x, n, &w);
    w += R_print.gap;

    for (i = 0; i < n; i++)
    {
        if (width + w > R_print.width)
        {
            DO_newline;
        }
        Rprintf("%s", EncodeInteger(x[i], w));
        width += w;
    }
    Rprintf("\n");
}

void printRealVector(double *x, int n, int indx)
{
    int i, w, d, e, labwidth = 0, width;

    DO_first_lab;
    formatReal(x, n, &w, &d, &e);
    w += R_print.gap;

    for (i = 0; i < n; i++)
    {
        if (width + w > R_print.width)
        {
            DO_newline;
        }
        Rprintf("%s", EncodeReal(x[i], w, d, e));
        width += w;
    }
    Rprintf("\n");
}

void printComplexVector(Rcomplex *x, int n, int indx)
{
    int i, w, wr, dr, er, wi, di, ei, labwidth = 0, width;

    DO_first_lab;
    formatComplex(x, n, &wr, &dr, &er, &wi, &di, &ei);

    w = wr + wi + 2; /* +2 for "+" and "i" */
    w += R_print.gap;

    for (i = 0; i < n; i++)
    {
        if (width + w > R_print.width)
        {
            DO_newline;
        }
        if (ISNA(x[i].r) || ISNA(x[i].i))
            Rprintf("%s", EncodeReal(NA_REAL, w, 0, 0));
        else
            Rprintf("%s", EncodeComplex(x[i], wr + R_print.gap, dr, er, wi, di, ei));
        width += w;
    }
    Rprintf("\n");
}

static void printStringVector(SEXP *x, int n, int quote, int indx)
{
    int i, w, labwidth = 0, width;

    DO_first_lab;
    formatString(x, n, &w, quote);

    for (i = 0; i < n; i++)
    {
        if (i > 0 && width + w + R_print.gap > R_print.width)
        {
            DO_newline;
        }
        Rprintf("%*s%s", R_print.gap, "", EncodeString(CHAR(x[i]), w, quote, Rprt_adj_left));
        width += w + R_print.gap;
    }
    Rprintf("\n");
}

void printVector(SEXP x, int indx, int quote)
{
    /* print R vector x[];	if(indx) print indices; if(quote) quote strings */
    int n;

    if ((n = LENGTH(x)) != 0)
        switch (TYPEOF(x))
        {
        case LGLSXP:
            printLogicalVector(LOGICAL(x), n, indx);
            break;
        case INTSXP:
            printIntegerVector(INTEGER(x), n, indx);
            break;
        case REALSXP:
            printRealVector(REAL(x), n, indx);
            break;
        case STRSXP:
            if (quote)
                printStringVector(STRING_PTR(x), n, '"', indx);
            else
                printStringVector(STRING_PTR(x), n, 0, indx);
            break;
        case CPLXSXP:
            printComplexVector(COMPLEX(x), n, indx);
            break;
        }
    else
#define PRINT_V_0                                                                                                      \
    switch (TYPEOF(x))                                                                                                 \
    {                                                                                                                  \
    case LGLSXP:                                                                                                       \
        Rprintf("logical(0)\n");                                                                                       \
        break;                                                                                                         \
    case INTSXP:                                                                                                       \
    case REALSXP:                                                                                                      \
        Rprintf("numeric(0)\n");                                                                                       \
        break;                                                                                                         \
    case CPLXSXP:                                                                                                      \
        Rprintf("complex(0)\n");                                                                                       \
        break;                                                                                                         \
    case STRSXP:                                                                                                       \
        Rprintf("character(0)\n");                                                                                     \
        break;                                                                                                         \
    }
        PRINT_V_0;
}

#undef DO_first_lab
#undef DO_newline

/* The following code prints vectors which have every element named.

 * Primitives for each type of vector are presented first, followed
 * by the main (dispatching) function.
 * 1) These primitives are almost identical... ==> use PRINT_N_VECTOR macro
 * 2) S prints a _space_ in the first column for named vectors; we dont.
 */

#define PRINT_N_VECTOR(INI_FORMAT, PRINT_1)                                                                            \
    {                                                                                                                  \
        int i, j, k, nlines, nperline, w, wn;                                                                          \
        INI_FORMAT;                                                                                                    \
                                                                                                                       \
        formatString(names, n, &wn, 0);                                                                                \
        if (w < wn)                                                                                                    \
            w = wn;                                                                                                    \
        nperline = R_print.width / (w + R_print.gap);                                                                  \
        if (nperline <= 0)                                                                                             \
            nperline = 1;                                                                                              \
        nlines = n / nperline;                                                                                         \
        if (n % nperline)                                                                                              \
            nlines += 1;                                                                                               \
                                                                                                                       \
        for (i = 0; i < nlines; i++)                                                                                   \
        {                                                                                                              \
            if (i)                                                                                                     \
                Rprintf("\n");                                                                                         \
            for (j = 0; j < nperline && (k = i * nperline + j) < n; j++)                                               \
                Rprintf("%s%*s", EncodeString(CHAR(names[k]), w, 0, Rprt_adj_right), R_print.gap, "");                 \
            Rprintf("\n");                                                                                             \
            for (j = 0; j < nperline && (k = i * nperline + j) < n; j++)                                               \
                PRINT_1;                                                                                               \
        }                                                                                                              \
        Rprintf("\n");                                                                                                 \
    }

static void printNamedLogicalVector(int *x, int n, SEXP *names)
    PRINT_N_VECTOR(formatLogical(x, n, &w), Rprintf("%s%*s", EncodeLogical(x[k], w), R_print.gap, ""))

        static void printNamedIntegerVector(int *x, int n, SEXP *names)
            PRINT_N_VECTOR(formatInteger(x, n, &w), Rprintf("%s%*s", EncodeInteger(x[k], w), R_print.gap, ""))

#undef INI_F_REAL
#define INI_F_REAL                                                                                                     \
    int d, e;                                                                                                          \
    formatReal(x, n, &w, &d, &e)

                static void printNamedRealVector(double *x, int n, SEXP *names)
                    PRINT_N_VECTOR(INI_F_REAL, Rprintf("%s%*s", EncodeReal(x[k], w, d, e), R_print.gap, ""))

#undef INI_F_CPLX
#define INI_F_CPLX                                                                                                     \
    int wr, dr, er, wi, di, ei;                                                                                        \
    formatComplex(x, n, &wr, &dr, &er, &wi, &di, &ei);                                                                 \
    w = wr + wi + 2

#undef P_IMAG_NA
#ifdef IEEE_754
#define P_IMAG_NA                                                                                                      \
    if (ISNAN(x[k].i))                                                                                                 \
        Rprintf("+%si", "NaN");                                                                                        \
    else
#else
#define P_IMAG_NA
#endif

                        static void printNamedComplexVector(Rcomplex *x, int n, SEXP *names)
                            PRINT_N_VECTOR(INI_F_CPLX,
                                           { /* PRINT_1 */
                                             if (j)
                                                 Rprintf("%*s", R_print.gap, "");
                                             if (ISNA(x[j].r) || ISNA(x[j].i))
                                             {
                                                 Rprintf("%s", EncodeReal(NA_REAL, w, 0, 0));
                                             }
                                             else
                                             {
                                                 Rprintf("%s", EncodeReal(x[k].r, wr, dr, er));
                                                 P_IMAG_NA
                                                 if (x[k].i >= 0)
                                                     Rprintf("+%si", EncodeReal(x[k].i, wi, di, ei));
                                                 else
                                                     Rprintf("-%si", EncodeReal(-x[k].i, wi, di, ei));
                                             }
                                           })

                                static void printNamedStringVector(SEXP *x, int n, int quote, SEXP *names)
                                    PRINT_N_VECTOR(formatString(x, n, &w, quote),
                                                   Rprintf("%s%*s", EncodeString(CHAR(x[k]), w, quote, Rprt_adj_right),
                                                           R_print.gap, ""))

                                        void printNamedVector(SEXP x, SEXP names, int quote, char *title)
{
    int n;

    if (title != NULL)
        Rprintf("%s\n", title);

    if ((n = LENGTH(x)) != 0)
        switch (TYPEOF(x))
        {
        case LGLSXP:
            printNamedLogicalVector(LOGICAL(x), n, STRING_PTR(names));
            break;
        case INTSXP:
            printNamedIntegerVector(INTEGER(x), n, STRING_PTR(names));
            break;
        case REALSXP:
            printNamedRealVector(REAL(x), n, STRING_PTR(names));
            break;
        case CPLXSXP:
            printNamedComplexVector(COMPLEX(x), n, STRING_PTR(names));
            break;
        case STRSXP:
            if (quote)
                quote = '"';
            printNamedStringVector(STRING_PTR(x), n, quote, STRING_PTR(names));
            break;
        }
    else
    {
        Rprintf("named ");
        PRINT_V_0;
    }
}
