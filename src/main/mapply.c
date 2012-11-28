/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2003-12   The R Core Team
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
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Defn.h>

SEXP attribute_hidden do_mapply(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    checkArity(op, args);

    SEXP f = CAR(args), varyingArgs = CADR(args), constantArgs = CADDR(args);
    int i, j, m, named, zero = 0;
    R_xlen_t *lengths, *counters, longest = 0;
    SEXP vnames, fcall = R_NilValue, mindex, nindex, tmp1, tmp2, ans;
    static SEXP length_op = NULL;

    /* Store the .Primitive for 'length' for DispatchOrEval to use. */
    if (length_op == NULL)
    {
        SEXP R_lengthSymbol = install("length");
        length_op = eval(R_lengthSymbol, R_BaseEnv);
        if (TYPEOF(length_op) != BUILTINSXP)
        {
            length_op = NULL;
            error("'length' is not a BUILTIN");
        }
        R_PreserveObject(length_op);
    }

    m = length(varyingArgs);
    vnames = PROTECT(getAttrib(varyingArgs, R_NamesSymbol));
    named = vnames != R_NilValue;

    lengths = (R_xlen_t *)R_alloc(m, sizeof(R_xlen_t));
    for (i = 0; i < m; i++)
    {
        int dispatch_ok = 0;
        tmp1 = VECTOR_ELT(varyingArgs, i);
        if (isObject(tmp1))
        {                                // possibly dispatch on  length():
            PROTECT(tmp2 = list1(tmp1)); // DispatchOrEval() needs 'args', a pairlist
            dispatch_ok = DispatchOrEval(call, length_op, "length", tmp2, rho, &ans, 0, 1);
            UNPROTECT(1);
        }
        lengths[i] = dispatch_ok ? asInteger(ans) : xlength(tmp1);
        if (lengths[i] == 0)
            zero++;
        if (lengths[i] > longest)
            longest = lengths[i];
    }
    if (zero && longest)
        error(_("zero-length inputs cannot be mixed with those of non-zero length"));

    counters = (R_xlen_t *)R_alloc(m, sizeof(R_xlen_t));
    memset(counters, 0, m * sizeof(R_xlen_t));
    // for(i = 0; i < m; counters[i++] = 0);

    mindex = PROTECT(allocVector(VECSXP, m));
    nindex = PROTECT(allocVector(VECSXP, m));

    /* build a call like
       f(dots[[1]][[4]], dots[[2]][[4]], dots[[3]][[4]], d=7)
    */

    if (constantArgs == R_NilValue)
        PROTECT(fcall = R_NilValue);
    else if (isVectorList(constantArgs))
        PROTECT(fcall = VectorToPairList(constantArgs));
    else
        error(_("argument 'MoreArgs' of 'mapply' is not a list"));

    Rboolean realIndx = longest > INT_MAX;
    for (j = m - 1; j >= 0; j--)
    {
        SET_VECTOR_ELT(mindex, j, ScalarInteger(j + 1));
        SET_VECTOR_ELT(nindex, j, allocVector(realIndx ? REALSXP : INTSXP, 1));
        PROTECT(tmp1 = lang3(R_Bracket2Symbol, install("dots"), VECTOR_ELT(mindex, j)));
        PROTECT(tmp2 = lang3(R_Bracket2Symbol, tmp1, VECTOR_ELT(nindex, j)));
        UNPROTECT(3);
        PROTECT(fcall = LCONS(tmp2, fcall));
        if (named && CHAR(STRING_ELT(vnames, j))[0] != '\0')
            SET_TAG(fcall, install(translateChar(STRING_ELT(vnames, j))));
    }

    UNPROTECT(1);
    PROTECT(fcall = LCONS(f, fcall));

    PROTECT(ans = allocVector(VECSXP, longest));

    for (i = 0; i < longest; i++)
    {
        for (j = 0; j < m; j++)
        {
            counters[j] = (++counters[j] > lengths[j]) ? 1 : counters[j];
            if (realIndx)
                REAL(VECTOR_ELT(nindex, j))[0] = (double)counters[j];
            else
                INTEGER(VECTOR_ELT(nindex, j))[0] = (int)counters[j];
        }
        SET_VECTOR_ELT(ans, i, eval(fcall, rho));
    }

    for (j = 0; j < m; j++)
        if (counters[j] != lengths[j])
            warning(_("longer argument not a multiple of length of shorter"));

    UNPROTECT(5);
    return (ans);
}
