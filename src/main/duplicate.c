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

/*  duplicate  -  object duplication  */

/*  Because we try to maintain the illusion of call by
 *  value, we often need to duplicate entire data
 *  objects.  There are a couple of points to note.
 *  First, duplication of list-like objects is done
 *  iteratively to prevent growth of the pointer
 *  protection stack, and second, the duplication of
 *  promises requires that the promises be forced and
 *  the value duplicated.  */

SEXP duplicate(SEXP s)
{
    SEXP h, t, sp;

    switch (TYPEOF(s))
    {
    case NILSXP:
    case SYMSXP:
    case ENVSXP:
    case SPECIALSXP:
    case BUILTINSXP:
        return s;
    case CLOSXP:
        PROTECT(s);
        t = allocSExp(CLOSXP);
        FORMALS(t) = FORMALS(s);
        BODY(t) = BODY(s);
        CLOENV(t) = CLOENV(s);
        ATTRIB(t) = duplicate(ATTRIB(s));
        UNPROTECT(1);
        break;
    case LISTSXP:
        PROTECT(sp = s);
        PROTECT(h = t = CONS(R_NilValue, R_NilValue));
        while (sp != R_NilValue)
        {
            CDR(t) = CONS(duplicate(CAR(sp)), R_NilValue);
            t = CDR(t);
            TAG(t) = TAG(sp);
            ATTRIB(t) = duplicate(ATTRIB(sp));
            sp = CDR(sp);
        }
        t = CDR(h);
        UNPROTECT(2);
        break;
    case LANGSXP:
        PROTECT(sp = s);
        PROTECT(h = t = CONS(R_NilValue, R_NilValue));
        while (sp != R_NilValue)
        {
            CDR(t) = CONS(duplicate(CAR(sp)), R_NilValue);
            t = CDR(t);
            TAG(t) = TAG(sp);
            ATTRIB(t) = duplicate(ATTRIB(sp));
            sp = CDR(sp);
        }
        t = CDR(h);
        TYPEOF(t) = LANGSXP;
        ATTRIB(t) = duplicate(ATTRIB(s));
        UNPROTECT(2);
        break;
    case CHARSXP:
        PROTECT(s);
        PROTECT(t = allocString(strlen(CHAR(s))));
        strcpy(CHAR(t), CHAR(s));
        ATTRIB(t) = duplicate(ATTRIB(s));
        UNPROTECT(2);
        break;
    case STRSXP:
    case LGLSXP:
    case FACTSXP:
    case ORDSXP:
    case INTSXP:
    case REALSXP:
    case CPLXSXP:
    case EXPRSXP:
        PROTECT(s);
        t = allocVector(TYPEOF(s), LENGTH(s));
        copyVector(t, s);
        LEVELS(t) = LEVELS(s);
        PROTECT(t);
        ATTRIB(t) = duplicate(ATTRIB(s));
        UNPROTECT(2);
        break;
    case PROMSXP: /* duplication requires that we evaluate the promise */
        if (PRVALUE(s) == R_UnboundValue)
        {
            t = eval(PREXPR(s), PRENV(s));
            PRVALUE(s) = t;
        }
        t = duplicate(PRVALUE(s));
        break;
    default:
        UNIMPLEMENTED("duplicate");
    }
    if (TYPEOF(t) == TYPEOF(s)) /* surely it only makes sense in this case*/
        OBJECT(t) = OBJECT(s);
    return t;
}

void copyVector(SEXP s, SEXP t)
{
    int i, ns, nt;

    nt = LENGTH(t);
    ns = LENGTH(s);
    switch (TYPEOF(s))
    {
    case STRSXP:
    case EXPRSXP:
        for (i = 0; i < ns; i++)
            VECTOR(s)[i] = VECTOR(t)[i % nt];
        break;
    case LGLSXP:
        for (i = 0; i < ns; i++)
            LOGICAL(s)[i] = LOGICAL(t)[i % nt];
        break;
    case FACTSXP:
    case ORDSXP:
        for (i = 0; i < ns; i++)
            FACTOR(s)[i] = FACTOR(t)[i % nt];
        break;
    case INTSXP:
        for (i = 0; i < ns; i++)
            INTEGER(s)[i] = INTEGER(t)[i % nt];
        break;
    case REALSXP:
        for (i = 0; i < ns; i++)
            REAL(s)[i] = REAL(t)[i % nt];
        break;
#ifdef COMPLEX_DATA
    case CPLXSXP:
        for (i = 0; i < ns; i++)
            COMPLEX(s)[i] = COMPLEX(t)[i % nt];
        break;
#endif
    default:
        UNIMPLEMENTED("copyVector");
    }
}

void copyListMatrix(SEXP s, SEXP t, int byrow)
{
    SEXP pt, tmp;
    int i, j, nr, nc, ns;

    nr = nrows(s);
    nc = ncols(s);
    ns = nr * nc;
    pt = t;
    if (byrow)
    {
        PROTECT(tmp = allocVector(STRSXP, nr * nc));
        for (i = 0; i < nr; i++)
            for (j = 0; j < nc; j++)
            {
                STRING(tmp)[i + j * nr] = duplicate(CAR(pt));
                pt = CDR(pt);
                if (pt == R_NilValue)
                    pt = t;
            }
        for (i = 0; i < ns; i++)
        {
            CAR(s) = STRING(tmp)[i++];
            s = CDR(s);
        }
        UNPROTECT(1);
    }
    else
    {
        for (i = 0; i < ns; i++)
        {
            CAR(s) = duplicate(CAR(pt));
            s = CDR(s);
            pt = CDR(pt);
            if (pt == R_NilValue)
                pt = t;
        }
    }
}

void copyMatrix(SEXP s, SEXP t, int byrow)
{
    int i, j, k, nr, nc, nt;

    nr = nrows(s);
    nc = ncols(s);
    nt = LENGTH(t);
    k = 0;

    if (byrow)
    {
        switch (TYPEOF(s))
        {
        case STRSXP:
            for (i = 0; i < nr; i++)
                for (j = 0; j < nc; j++)
                    STRING(s)[i + j * nr] = STRING(t)[k++ % nt];
            break;
        case LGLSXP:
            for (i = 0; i < nr; i++)
                for (j = 0; j < nc; j++)
                    LOGICAL(s)[i + j * nr] = LOGICAL(t)[k++ % nt];
            break;
        case FACTSXP:
        case ORDSXP:
            for (i = 0; i < nr; i++)
                for (j = 0; j < nc; j++)
                    FACTOR(s)[i + j * nr] = FACTOR(t)[k++ % nt];
            break;
        case INTSXP:
            for (i = 0; i < nr; i++)
                for (j = 0; j < nc; j++)
                    INTEGER(s)[i + j * nr] = INTEGER(t)[k++ % nt];
            break;
        case REALSXP:
            for (i = 0; i < nr; i++)
                for (j = 0; j < nc; j++)
                    REAL(s)[i + j * nr] = REAL(t)[k++ % nt];
            break;
#ifdef COMPLEX_DATA
        case CPLXSXP:
            for (i = 0; i < nr; i++)
                for (j = 0; j < nc; j++)
                    COMPLEX(s)[i + j * nr] = COMPLEX(t)[k++ % nt];
            break;
#endif
        default:
            UNIMPLEMENTED("copyMatrix");
        }
    }
    else
        copyVector(s, t);
}
