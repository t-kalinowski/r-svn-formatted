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
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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

/* This macro pulls out the common code in copying an atomic vector.
   The special handling of the scalar case (__n__ == 1) seems to make
   a small but measurable difference, at least for some cases. */
#define DUPLICATE_ATOMIC_VECTOR(type, fun, to, from)                                                                   \
    do                                                                                                                 \
    {                                                                                                                  \
        int __n__ = LENGTH(from);                                                                                      \
        PROTECT(from);                                                                                                 \
        PROTECT(to = allocVector(TYPEOF(from), __n__));                                                                \
        if (__n__ == 1)                                                                                                \
            fun(to)[0] = fun(from)[0];                                                                                 \
        else                                                                                                           \
        {                                                                                                              \
            int __i__;                                                                                                 \
            type *__fp__ = fun(from), *__tp__ = fun(to);                                                               \
            for (__i__ = 0; __i__ < __n__; __i__++)                                                                    \
                __tp__[__i__] = __fp__[__i__];                                                                         \
        }                                                                                                              \
        DUPLICATE_ATTRIB(to, from);                                                                                    \
        SET_TRUELENGTH(to, TRUELENGTH(from));                                                                          \
        UNPROTECT(2);                                                                                                  \
    } while (0)

/* The following macros avoid the cost of going through calls to the
   assignment functions (and duplicate in the case of ATTRIB) when the
   ATTRIB or TAG value to be stored is R_NilValue, the value the field
   will have been set to by the allocation function */
#define DUPLICATE_ATTRIB(to, from)                                                                                     \
    do                                                                                                                 \
    {                                                                                                                  \
        SEXP __a__ = ATTRIB(from);                                                                                     \
        if (__a__ != R_NilValue)                                                                                       \
            SET_ATTRIB(to, duplicate(__a__));                                                                          \
    } while (0)

#define COPY_TAG(to, from)                                                                                             \
    do                                                                                                                 \
    {                                                                                                                  \
        SEXP __tag__ = TAG(from);                                                                                      \
        if (__tag__ != R_NilValue)                                                                                     \
            SET_TAG(to, __tag__);                                                                                      \
    } while (0)

SEXP duplicate(SEXP s)
{
    SEXP h, t, sp;
    int i, n;

    switch (TYPEOF(s))
    {
    case NILSXP:
    case SYMSXP:
    case ENVSXP:
    case SPECIALSXP:
    case BUILTINSXP:
    case EXTPTRSXP:
        return s;
    case CLOSXP:
        PROTECT(s);
        PROTECT(t = allocSExp(CLOSXP));
        SET_FORMALS(t, FORMALS(s));
        SET_BODY(t, BODY(s));
        SET_CLOENV(t, CLOENV(s));
        DUPLICATE_ATTRIB(t, s);
        UNPROTECT(2);
        break;
    case LISTSXP:
        PROTECT(sp = s);
        PROTECT(h = t = CONS(R_NilValue, R_NilValue));
        while (sp != R_NilValue)
        {
            SETCDR(t, CONS(duplicate(CAR(sp)), R_NilValue));
            t = CDR(t);
            COPY_TAG(t, sp);
            DUPLICATE_ATTRIB(t, sp);
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
            SETCDR(t, CONS(duplicate(CAR(sp)), R_NilValue));
            t = CDR(t);
            COPY_TAG(t, sp);
            DUPLICATE_ATTRIB(t, sp);
            sp = CDR(sp);
        }
        t = CDR(h);
        SET_TYPEOF(t, LANGSXP);
        DUPLICATE_ATTRIB(t, s);
        UNPROTECT(2);
        break;
    case CHARSXP:
        PROTECT(s);
        PROTECT(t = allocString(strlen(CHAR(s))));
        strcpy(CHAR(t), CHAR(s));
        DUPLICATE_ATTRIB(t, s);
        UNPROTECT(2);
        break;
    case EXPRSXP:
    case VECSXP:
        n = LENGTH(s);
        PROTECT(s);
        PROTECT(t = allocVector(TYPEOF(s), n));
        for (i = 0; i < n; i++)
            SET_VECTOR_ELT(t, i, duplicate(VECTOR_ELT(s, i)));
        DUPLICATE_ATTRIB(t, s);
        SET_TRUELENGTH(t, TRUELENGTH(s));
        UNPROTECT(2);
        break;
    case LGLSXP:
        DUPLICATE_ATOMIC_VECTOR(int, LOGICAL, t, s);
        break;
    case INTSXP:
        DUPLICATE_ATOMIC_VECTOR(int, INTEGER, t, s);
        break;
    case REALSXP:
        DUPLICATE_ATOMIC_VECTOR(double, REAL, t, s);
        break;
    case CPLXSXP:
        DUPLICATE_ATOMIC_VECTOR(Rcomplex, COMPLEX, t, s);
        break;
    case STRSXP:
        /* direct copying and bypassing the write barrier is OK since
           t was just allocated and so it cannot be older than any of
           the elements in s.  LT */
        DUPLICATE_ATOMIC_VECTOR(SEXP, STRING_PTR, t, s);
        break;
    case PROMSXP: /* duplication requires that we evaluate the promise */
#ifdef OLD
        if (PRVALUE(s) == R_UnboundValue)
        {
            t = eval(PREXPR(s), PRENV(s));
            PRVALUE(s) = t;
        }
        t = duplicate(PRVALUE(s));
#endif
        return s;
        break;
    default:
        UNIMPLEMENTED("duplicate");
        t = s; /* for -Wall */
    }
    if (TYPEOF(t) == TYPEOF(s)) /* surely it only makes sense in this case*/
        SET_OBJECT(t, OBJECT(s));
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
            SET_VECTOR_ELT(s, i, VECTOR_ELT(t, i % nt));
        break;
    case LGLSXP:
        for (i = 0; i < ns; i++)
            LOGICAL(s)[i] = LOGICAL(t)[i % nt];
        break;
    case INTSXP:
        for (i = 0; i < ns; i++)
            INTEGER(s)[i] = INTEGER(t)[i % nt];
        break;
    case REALSXP:
        for (i = 0; i < ns; i++)
            REAL(s)[i] = REAL(t)[i % nt];
        break;
    case CPLXSXP:
        for (i = 0; i < ns; i++)
            COMPLEX(s)[i] = COMPLEX(t)[i % nt];
        break;
    default:
        UNIMPLEMENTED("copyVector");
    }
}

void copyListMatrix(SEXP s, SEXP t, Rboolean byrow)
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
                SET_STRING_ELT(tmp, i + j * nr, duplicate(CAR(pt)));
                pt = CDR(pt);
                if (pt == R_NilValue)
                    pt = t;
            }
        for (i = 0; i < ns; i++)
        {
            SETCAR(s, STRING_ELT(tmp, i++));
            s = CDR(s);
        }
        UNPROTECT(1);
    }
    else
    {
        for (i = 0; i < ns; i++)
        {
            SETCAR(s, duplicate(CAR(pt)));
            s = CDR(s);
            pt = CDR(pt);
            if (pt == R_NilValue)
                pt = t;
        }
    }
}

void copyMatrix(SEXP s, SEXP t, Rboolean byrow)
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
                    SET_STRING_ELT(s, i + j * nr, STRING_ELT(t, k++ % nt));
            break;
        case LGLSXP:
            for (i = 0; i < nr; i++)
                for (j = 0; j < nc; j++)
                    LOGICAL(s)[i + j * nr] = LOGICAL(t)[k++ % nt];
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
        case CPLXSXP:
            for (i = 0; i < nr; i++)
                for (j = 0; j < nc; j++)
                    COMPLEX(s)[i + j * nr] = COMPLEX(t)[k++ % nt];
            break;
        default:
            UNIMPLEMENTED("copyMatrix");
        }
    }
    else
        copyVector(s, t);
}
