/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996	Robert Gentleman and Ross Ihaka
 *  Copyright (C) 2000--2006	The R Development Core Team.
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
 *  Foundation, Inc., 51 Franklin Street Fifth Floor, Boston, MA 02110-1301  USA
 *
 *
 *  EXPORTS	printMatrix()
 *		printArray()
 *
 *  See ./printutils.c	 for general remarks on Printing
 *			 and the Encode.. utils.
 *
 *  See ./format.c	 for the  format_FOO_  functions used below.
 */

/* <UTF8> char here is handled as a whole,
   but lengths were used as display widths */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Defn.h"
#include "Print.h"

/* We need display width of a string */
int Rstrwid(char *str, int slen, int quote); /* from printutils.c */
#define strwidth(x) Rstrwid(x, strlen(x), 0)

/* This is the first (of 6)  print<TYPE>Matrix()  functions.
 * We define macros that will be re-used in the other functions,
 * and comment the common code here (only):
 */
static void printLogicalMatrix(SEXP sx, int offset, int r_pr, int r, int c, SEXP rl, SEXP cl, char *rn, char *cn)
{
    int *x;

/* initialization; particularly of row labels, rl= dimnames(.)[[1]] and
 * rn = names(dimnames(.))[1] : */
#define _PRINT_INIT_rl_rn                                                                                              \
    SEXP sw;                                                                                                           \
    int *w;                                                                                                            \
    int width, rlabw = -1, clabw = -1; /* -Wall */                                                                     \
    int i, j, jmin = 0, jmax = 0, lbloff = 0;                                                                          \
                                                                                                                       \
    if (!isNull(rl))                                                                                                   \
        formatString(STRING_PTR(rl), r, &rlabw, 0);                                                                    \
    else                                                                                                               \
        rlabw = IndexWidth(r + 1) + 3;                                                                                 \
                                                                                                                       \
    if (rn)                                                                                                            \
    {                                                                                                                  \
        int rnw = strwidth(rn);                                                                                        \
        if (rnw < rlabw + R_MIN_LBLOFF)                                                                                \
            lbloff = R_MIN_LBLOFF;                                                                                     \
        else                                                                                                           \
            lbloff = rnw - rlabw;                                                                                      \
                                                                                                                       \
        rlabw += lbloff;                                                                                               \
    }

    _PRINT_INIT_rl_rn;

    sw = allocVector(INTSXP, c);
    x = LOGICAL(sx) + offset;
    w = INTEGER(sw);
    /* compute w[j] = column-width of j(+1)-th column : */
    for (j = 0; j < c; j++)
    {
        formatLogical(&x[j * r], r, &w[j]);

#define _PRINT_SET_clabw                                                                                               \
                                                                                                                       \
    if (!isNull(cl))                                                                                                   \
    {                                                                                                                  \
        if (STRING_ELT(cl, j) == NA_STRING)                                                                            \
            clabw = R_print.na_width_noquote;                                                                          \
        else                                                                                                           \
            clabw = strwidth(CHAR(STRING_ELT(cl, j)));                                                                 \
    }                                                                                                                  \
    else                                                                                                               \
        clabw = IndexWidth(j + 1) + 3;

        _PRINT_SET_clabw;

        if (w[j] < clabw)
            w[j] = clabw;
        w[j] += R_print.gap;
    }

#define _PRINT_DEAL_c_eq_0                                                                                             \
                                                                                                                       \
    if (c == 0)                                                                                                        \
    {                                                                                                                  \
        for (i = 0; i < r; i++)                                                                                        \
            MatrixRowLabel(rl, i, rlabw, lbloff);                                                                      \
        Rprintf("\n");                                                                                                 \
        return;                                                                                                        \
    }
    _PRINT_DEAL_c_eq_0;

    while (jmin < c)
    {
        /* print columns  jmin:(jmax-1)	 where jmax has to be determined first */

        width = rlabw;
        /* initially, jmax = jmin */
        do
        {
            width += w[jmax];
            jmax++;

        } while (jmax < c && width + w[jmax] < R_print.width);

#define _PRINT_ROW_LAB                                                                                                 \
                                                                                                                       \
    if (cn != NULL)                                                                                                    \
        Rprintf("%*s%s\n", rlabw, "", cn);                                                                             \
                                                                                                                       \
    if (rn != NULL)                                                                                                    \
        Rprintf("%*s", -rlabw, rn);                                                                                    \
    else                                                                                                               \
        Rprintf("%*s", rlabw, "");

        _PRINT_ROW_LAB;

        for (j = jmin; j < jmax; j++)
            MatrixColumnLabel(cl, j, w[j]);
        for (i = 0; i < r_pr; i++)
        {
            MatrixRowLabel(rl, i, rlabw, lbloff); /* starting with an "\n" */
            for (j = jmin; j < jmax; j++)
            {
                Rprintf("%s", EncodeLogical(x[i + j * r], w[j]));
            }
        }
        Rprintf("\n");
        jmin = jmax;
    }
}

static void printIntegerMatrix(SEXP sx, int offset, int r_pr, int r, int c, SEXP rl, SEXP cl, char *rn, char *cn)
{
    int *x;

    _PRINT_INIT_rl_rn;

    sw = allocVector(INTSXP, c);
    x = INTEGER(sx) + offset;
    w = INTEGER(sw);
    for (j = 0; j < c; j++)
    {
        formatInteger(&x[j * r], r, &w[j]);
        _PRINT_SET_clabw;
        if (w[j] < clabw)
            w[j] = clabw;
        w[j] += R_print.gap;
    }
    _PRINT_DEAL_c_eq_0;
    while (jmin < c)
    {
        width = rlabw;
        do
        {
            width += w[jmax];
            jmax++;
        } while (jmax < c && width + w[jmax] < R_print.width);

        _PRINT_ROW_LAB;

        for (j = jmin; j < jmax; j++)
            MatrixColumnLabel(cl, j, w[j]);
        for (i = 0; i < r_pr; i++)
        {
            MatrixRowLabel(rl, i, rlabw, lbloff);
            for (j = jmin; j < jmax; j++)
            {
                Rprintf("%s", EncodeInteger(x[i + j * r], w[j]));
            }
        }
        Rprintf("\n");
        jmin = jmax;
    }
}

static void printRealMatrix(SEXP sx, int offset, int r_pr, int r, int c, SEXP rl, SEXP cl, char *rn, char *cn)
{
    SEXP sd, se;
    double *x;
    int *d, *e;
    _PRINT_INIT_rl_rn;

    PROTECT(sd = allocVector(INTSXP, c));
    PROTECT(se = allocVector(INTSXP, c));
    sw = allocVector(INTSXP, c);
    UNPROTECT(2);
    x = REAL(sx) + offset;
    d = INTEGER(sd);
    e = INTEGER(se);
    w = INTEGER(sw);

    for (j = 0; j < c; j++)
    {
        formatReal(&x[j * r], r, &w[j], &d[j], &e[j], 0);
        _PRINT_SET_clabw;
        if (w[j] < clabw)
            w[j] = clabw;
        w[j] += R_print.gap;
    }
    _PRINT_DEAL_c_eq_0;
    while (jmin < c)
    {
        width = rlabw;
        do
        {
            width += w[jmax];
            jmax++;
        } while (jmax < c && width + w[jmax] < R_print.width);

        _PRINT_ROW_LAB;

        for (j = jmin; j < jmax; j++)
            MatrixColumnLabel(cl, j, w[j]);
        for (i = 0; i < r_pr; i++)
        {
            MatrixRowLabel(rl, i, rlabw, lbloff);
            for (j = jmin; j < jmax; j++)
            {
                Rprintf("%s", EncodeReal(x[i + j * r], w[j], d[j], e[j], OutDec));
            }
        }
        Rprintf("\n");
        jmin = jmax;
    }
}

static void printComplexMatrix(SEXP sx, int offset, int r_pr, int r, int c, SEXP rl, SEXP cl, char *rn, char *cn)
{
    SEXP sdr, ser, swr, sdi, sei, swi;
    Rcomplex *x;
    int *dr, *er, *wr, *di, *ei, *wi;
    _PRINT_INIT_rl_rn;

    PROTECT(sdr = allocVector(INTSXP, c));
    PROTECT(ser = allocVector(INTSXP, c));
    PROTECT(swr = allocVector(INTSXP, c));
    PROTECT(sdi = allocVector(INTSXP, c));
    PROTECT(sei = allocVector(INTSXP, c));
    PROTECT(swi = allocVector(INTSXP, c));
    PROTECT(sw = allocVector(INTSXP, c));
    UNPROTECT(7);
    x = COMPLEX(sx) + offset;
    dr = INTEGER(sdr);
    er = INTEGER(ser);
    wr = INTEGER(swr);
    di = INTEGER(sdi);
    ei = INTEGER(sei);
    wi = INTEGER(swi);
    w = INTEGER(sw);

    /* Determine the column widths */

    for (j = 0; j < c; j++)
    {
        formatComplex(&x[j * r], r, &wr[j], &dr[j], &er[j], &wi[j], &di[j], &ei[j], 0);
        _PRINT_SET_clabw;
        w[j] = wr[j] + wi[j] + 2;
        if (w[j] < clabw)
            w[j] = clabw;
        w[j] += R_print.gap;
    }

    _PRINT_DEAL_c_eq_0;
    while (jmin < c)
    {
        width = rlabw;
        do
        {
            width += w[jmax];
            jmax++;
        } while (jmax < c && width + w[jmax] < R_print.width);

        _PRINT_ROW_LAB;

        for (j = jmin; j < jmax; j++)
            MatrixColumnLabel(cl, j, w[j]);
        for (i = 0; i < r_pr; i++)
        {
            MatrixRowLabel(rl, i, rlabw, lbloff);
            for (j = jmin; j < jmax; j++)
            {
                if (ISNA(x[i + j * r].r) || ISNA(x[i + j * r].i))
                    Rprintf("%s", EncodeReal(NA_REAL, w[j], 0, 0, OutDec));
                else
                    Rprintf("%s", EncodeComplex(x[i + j * r], wr[j] + R_print.gap, dr[j], er[j], wi[j], di[j], ei[j],
                                                OutDec));
            }
        }
        Rprintf("\n");
        jmin = jmax;
    }
}

static void printStringMatrix(SEXP sx, int offset, int r_pr, int r, int c, int quote, int right, SEXP rl, SEXP cl,
                              char *rn, char *cn)
{
    SEXP *x;
    _PRINT_INIT_rl_rn;

    sw = allocVector(INTSXP, c);
    x = STRING_PTR(sx) + offset;
    w = INTEGER(sw);
    for (j = 0; j < c; j++)
    {
        formatString(&x[j * r], r, &w[j], quote);
        _PRINT_SET_clabw;
        if (w[j] < clabw)
            w[j] = clabw;
    }
    _PRINT_DEAL_c_eq_0;
    while (jmin < c)
    {
        width = rlabw;
        do
        {
            width += w[jmax] + R_print.gap;
            jmax++;
        } while (jmax < c && width + w[jmax] + R_print.gap < R_print.width);

        _PRINT_ROW_LAB;

        if (right)
        {
            for (j = jmin; j < jmax; j++)
                RightMatrixColumnLabel(cl, j, w[j]);
        }
        else
        {
            for (j = jmin; j < jmax; j++)
                LeftMatrixColumnLabel(cl, j, w[j]);
        }
        for (i = 0; i < r_pr; i++)
        {
            MatrixRowLabel(rl, i, rlabw, lbloff);
            for (j = jmin; j < jmax; j++)
            {
                Rprintf("%*s%s", R_print.gap, "", EncodeString(x[i + j * r], w[j], quote, right));
            }
        }
        Rprintf("\n");
        jmin = jmax;
    }
}

static void printRawMatrix(SEXP sx, int offset, int r_pr, int r, int c, SEXP rl, SEXP cl, char *rn, char *cn)
{
    Rbyte *x;
    _PRINT_INIT_rl_rn;

    sw = allocVector(INTSXP, c);
    x = RAW(sx) + offset;
    w = INTEGER(sw);
    for (j = 0; j < c; j++)
    {
        formatRaw(&x[j * r], r, &w[j]);
        _PRINT_SET_clabw;
        if (w[j] < clabw)
            w[j] = clabw;
        w[j] += R_print.gap;
    }
    _PRINT_DEAL_c_eq_0;
    while (jmin < c)
    {
        width = rlabw;
        do
        {
            width += w[jmax];
            jmax++;
        } while (jmax < c && width + w[jmax] < R_print.width);

        _PRINT_ROW_LAB;

        for (j = jmin; j < jmax; j++)
            MatrixColumnLabel(cl, j, w[j]);
        for (i = 0; i < r_pr; i++)
        {
            MatrixRowLabel(rl, i, rlabw, lbloff);
            for (j = jmin; j < jmax; j++)
                Rprintf("%*s%s", w[j] - 2, "", EncodeRaw(x[i + j * r]));
        }
        Rprintf("\n");
        jmin = jmax;
    }
}

void printMatrix(SEXP x, int offset, SEXP dim, int quote, int right, SEXP rl, SEXP cl, char *rn, char *cn)
{
    /* 'rl' and 'cl' are dimnames(.)[[1]] and dimnames(.)[[2]]  whereas
     * 'rn' and 'cn' are the  names(dimnames(.))
     */
    int r, c, r_pr;

    r = INTEGER(dim)[0];
    c = INTEGER(dim)[1];
    /* PR#850 */
    if ((rl != R_NilValue) && (r > length(rl)))
        error(_("too few row labels"));
    if ((cl != R_NilValue) && (c > length(cl)))
        error(_("too few column labels"));
    if (r == 0 && c == 0)
    {
        Rprintf("<0 x 0 matrix>\n");
        return;
    }
    r_pr = r;
    if (c > 0 && R_print.max / c < r) /* avoid integer overflow */
        r_pr = R_print.max / c;
    switch (TYPEOF(x))
    {
    case LGLSXP:
        printLogicalMatrix(x, offset, r_pr, r, c, rl, cl, rn, cn);
        break;
    case INTSXP:
        printIntegerMatrix(x, offset, r_pr, r, c, rl, cl, rn, cn);
        break;
    case REALSXP:
        printRealMatrix(x, offset, r_pr, r, c, rl, cl, rn, cn);
        break;
    case CPLXSXP:
        printComplexMatrix(x, offset, r_pr, r, c, rl, cl, rn, cn);
        break;
    case STRSXP:
        if (quote)
            quote = '"';
        printStringMatrix(x, offset, r_pr, r, c, quote, right, rl, cl, rn, cn);
        break;
    case RAWSXP:
        printRawMatrix(x, offset, r_pr, r, c, rl, cl, rn, cn);
        break;
    default:
        UNIMPLEMENTED_TYPE("printMatrix", x);
    }
    if (r_pr < r)
        Rprintf(" [[ reached getOption(\"max.print\") -- omitted %d more rows ]]\n", r - r_pr);
}

static void printArrayGeneral(SEXP x, SEXP dim, int quote, int right, SEXP dimnames)
{
    /* == printArray(.) */

    int ndim = LENGTH(dim);
    char *rn = NULL, *cn = NULL;

    if (ndim == 1)
        printVector(x, 1, quote);
    else if (ndim == 2)
    {
        SEXP rl, cl;
        GetMatrixDimnames(x, &rl, &cl, &rn, &cn);
        printMatrix(x, 0, dim, quote, 0, rl, cl, rn, cn);
    }
    else
    { /* ndim >= 3 */
        SEXP dn, dnn, dn0, dn1;
        int i, j, has_dimnames, has_dnn, nb;
        int nr = INTEGER(dim)[0];
        int nc = INTEGER(dim)[1];
        int b = nr * nc;

        if (dimnames == R_NilValue)
        {
            has_dimnames = 0;
            has_dnn = 0;
            dn0 = R_NilValue;
            dn1 = R_NilValue;
            dnn = R_NilValue; /* -Wall */
        }
        else
        {
            dn0 = VECTOR_ELT(dimnames, 0);
            dn1 = VECTOR_ELT(dimnames, 1);
            has_dimnames = 1;
            dnn = getAttrib(dimnames, R_NamesSymbol);
            has_dnn = !isNull(dnn);
            if (has_dnn)
            {
                rn = CHAR(STRING_ELT(dnn, 0));
                cn = CHAR(STRING_ELT(dnn, 1));
            }
        }
        /* nb := #{entries} in a slice such as x[1,1,..] or equivalently,
         *       the number of matrix slices   x[ , , *, ..]  which
         *       are printed as matrices -- if options("max.print") allows */
        for (i = 2, nb = 1; i < ndim; i++)
            nb *= INTEGER(dim)[i];
        if (b > 0 && R_print.max / b < nb)
            nb = R_print.max / b;

        for (i = 0; i < nb; i++)
        {
            int k = 1;
            Rprintf(", ");
            for (j = 2; j < ndim; j++)
            {
                int l = (i / k) % INTEGER(dim)[j] + 1;
                if (has_dimnames && ((dn = VECTOR_ELT(dimnames, j)) != R_NilValue))
                {
                    if (has_dnn)
                        Rprintf(", %s = %s", CHAR(STRING_ELT(dnn, j)), CHAR(STRING_ELT(dn, l - 1)));
                    else
                        Rprintf(", %s", CHAR(STRING_ELT(dn, l - 1)));
                }
                else
                    Rprintf(", %d", l);
                k = k * INTEGER(dim)[j];
            }
            Rprintf("\n\n");
            switch (TYPEOF(x))
            {
            case LGLSXP:
                printLogicalMatrix(x, i * b, nr /*FIXME*/, nr, nc, dn0, dn1, rn, cn);
                break;
            case INTSXP:
                printIntegerMatrix(x, i * b, nr /*FIXME*/, nr, nc, dn0, dn1, rn, cn);
                break;
            case REALSXP:
                printRealMatrix(x, i * b, nr /*FIXME*/, nr, nc, dn0, dn1, rn, cn);
                break;
            case CPLXSXP:
                printComplexMatrix(x, i * b, nr /*FIXME*/, nr, nc, dn0, dn1, rn, cn);
                break;
            case STRSXP:
                if (quote)
                    quote = '"';
                printStringMatrix(x, i * b, nr /*FIXME*/, nr, nc, quote, right, dn0, dn1, rn, cn);
                break;
            case RAWSXP:
                printRawMatrix(x, i * b, nr /*FIXME*/, nr, nc, dn0, dn1, rn, cn);
                break;
            }
            Rprintf("\n");
        }
    }
}

void printArray(SEXP x, SEXP dim, int quote, int right, SEXP dimnames)
{
    printArrayGeneral(x, dim, quote, right, dimnames);
}
