/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997--1999  Robert Gentleman, Ross Ihaka and the
 *                            R Development Core Team
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

/* Code to handle list / vector switch */

#ifdef HAVE_CONFIG_H
#include <Rconfig.h>
#endif

#include "Defn.h"
#include "Mathlib.h" /* imax2 */

#define LIST_ASSIGN(x)                                                                                                 \
    {                                                                                                                  \
        VECTOR(ans_ptr)[ans_length] = x;                                                                               \
        ans_length++;                                                                                                  \
    }

static SEXP cbind(SEXP, SEXP, SEXPTYPE);
static SEXP rbind(SEXP, SEXP, SEXPTYPE);

/* The following code establishes the return type for the */
/* functions  unlist, c, cbind, and rbind and also determines */
/* whether the returned object is to have a names attribute. */

static int ans_flags;
static SEXP ans_ptr;
static int ans_length;
static SEXP ans_names;
static int ans_nnames;

static int HasNames(SEXP x)
{
    if (isVector(x))
    {
        if (!isNull(getAttrib(x, R_NamesSymbol)))
            return 1;
    }
    else if (isList(x))
    {
        while (!isNull(x))
        {
            if (!isNull(TAG(x)))
                return 1;
            x = CDR(x);
        }
    }
    return 0;
}

static void AnswerType(SEXP x, int recurse, int usenames)
{
    switch (TYPEOF(x))
    {
    case NILSXP:
        break;
    case LGLSXP:
        ans_flags |= 1;
        ans_length += LENGTH(x);
        break;
    case INTSXP:
        ans_flags |= 8;
        ans_length += LENGTH(x);
        break;
    case REALSXP:
        ans_flags |= 16;
        ans_length += LENGTH(x);
        break;
    case CPLXSXP:
        ans_flags |= 32;
        ans_length += LENGTH(x);
        break;
    case STRSXP:
        ans_flags |= 64;
        ans_length += LENGTH(x);
        break;
    case VECSXP:
    case EXPRSXP:
        if (recurse)
        {
            int i, n;
            n = length(x);
            if (usenames && !ans_names && !isNull(getAttrib(x, R_NamesSymbol)))
                ans_nnames = 1;
            for (i = 0; i < n; i++)
            {
                if (usenames && !ans_nnames)
                    ans_nnames = HasNames(VECTOR(x)[i]);
                AnswerType(VECTOR(x)[i], recurse, usenames);
            }
        }
        else
        {
            if (TYPEOF(x) == EXPRSXP)
                ans_flags |= 256;
            else
                ans_flags |= 128;
            ans_length += length(x);
        }
        break;
    case LISTSXP:
        if (recurse)
        {
            while (x != R_NilValue)
            {
                if (usenames && !ans_nnames)
                {
                    if (!isNull(TAG(x)))
                        ans_nnames = 1;
                    else
                        ans_nnames = HasNames(CAR(x));
                }
                AnswerType(CAR(x), recurse, usenames);
                x = CDR(x);
            }
        }
        else
        {
            ans_flags |= 128;
            ans_length += length(x);
        }
        break;
    default:
        ans_flags |= 128;
        ans_length += 1;
        break;
    }
}

/* The following functions are used to coerce arguments to */
/* the appropriate type for inclusion in the returned value. */

static void ListAnswer(SEXP x, int recurse)
{
    int i;

    switch (TYPEOF(x))
    {
    case NILSXP:
        break;
    case LGLSXP:
        for (i = 0; i < LENGTH(x); i++)
            LIST_ASSIGN(ScalarLogical(LOGICAL(x)[i]));
        break;
    case INTSXP:
        for (i = 0; i < LENGTH(x); i++)
            LIST_ASSIGN(ScalarInteger(INTEGER(x)[i]));
        break;
    case REALSXP:
        for (i = 0; i < LENGTH(x); i++)
            LIST_ASSIGN(ScalarReal(REAL(x)[i]));
        break;
    case CPLXSXP:
        for (i = 0; i < LENGTH(x); i++)
            LIST_ASSIGN(ScalarComplex(COMPLEX(x)[i]));
        break;
    case STRSXP:
        for (i = 0; i < LENGTH(x); i++)
            LIST_ASSIGN(ScalarString(STRING(x)[i]));
        break;
    case VECSXP:
    case EXPRSXP:
        if (recurse)
        {
            for (i = 0; i < LENGTH(x); i++)
                ListAnswer(VECTOR(x)[i], recurse);
        }
        else
        {
            for (i = 0; i < LENGTH(x); i++)
                LIST_ASSIGN(duplicate(VECTOR(x)[i]));
        }
        break;
    case LISTSXP:
        if (recurse)
        {
            while (x != R_NilValue)
            {
                ListAnswer(CAR(x), recurse);
                x = CDR(x);
            }
        }
        else
            while (x != R_NilValue)
            {
                LIST_ASSIGN(duplicate(CAR(x)));
                x = CDR(x);
            }
        break;
    default:
        LIST_ASSIGN(duplicate(x));
        break;
    }
}

static void StringAnswer(SEXP x)
{
    int i, n;
    switch (TYPEOF(x))
    {
    case NILSXP:
        break;
    case LISTSXP:
        while (x != R_NilValue)
        {
            StringAnswer(CAR(x));
            x = CDR(x);
        }
        break;
    case VECSXP:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            StringAnswer(VECTOR(x)[i]);
        break;
    default:
        PROTECT(x = coerceVector(x, STRSXP));
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            STRING(ans_ptr)[ans_length++] = STRING(x)[i];
        UNPROTECT(1);
        break;
    }
}

static void IntegerAnswer(SEXP x)
{
    int i, n;
    switch (TYPEOF(x))
    {
    case NILSXP:
        break;
    case LISTSXP:
        while (x != R_NilValue)
        {
            IntegerAnswer(CAR(x));
            x = CDR(x);
        }
        break;
    case VECSXP:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            IntegerAnswer(VECTOR(x)[i]);
        break;
    default:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            INTEGER(ans_ptr)[ans_length++] = INTEGER(x)[i];
        break;
    }
}

static void RealAnswer(SEXP x)
{
    int i, n, xi;
    switch (TYPEOF(x))
    {
    case NILSXP:
        break;
    case LISTSXP:
        while (x != R_NilValue)
        {
            RealAnswer(CAR(x));
            x = CDR(x);
        }
        break;
    case VECSXP:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            RealAnswer(VECTOR(x)[i]);
        break;
    case REALSXP:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            REAL(ans_ptr)[ans_length++] = REAL(x)[i];
        break;
    default:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
        {
            xi = INTEGER(x)[i];
            if (xi == NA_INTEGER)
                REAL(ans_ptr)[ans_length++] = NA_REAL;
            else
                REAL(ans_ptr)[ans_length++] = xi;
        }
        break;
    }
}

static void ComplexAnswer(SEXP x)
{
    int i, n, xi;
    switch (TYPEOF(x))
    {
    case NILSXP:
        break;
    case LISTSXP:
        while (x != R_NilValue)
        {
            ComplexAnswer(CAR(x));
            x = CDR(x);
        }
        break;
    case VECSXP:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            ComplexAnswer(VECTOR(x)[i]);
        break;
    case REALSXP:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
        {
            COMPLEX(ans_ptr)[ans_length].r = REAL(x)[i];
            COMPLEX(ans_ptr)[ans_length].i = 0.0;
            ans_length++;
        }
        break;
    case CPLXSXP:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            COMPLEX(ans_ptr)[ans_length++] = COMPLEX(x)[i];
        break;
    default:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
        {
            xi = INTEGER(x)[i];
            if (xi == NA_INTEGER)
                REAL(ans_ptr)[ans_length++] = NA_REAL;
            else
                REAL(ans_ptr)[ans_length++] = xi;
        }
        break;
    }
}

static SEXP NewBase(SEXP base, SEXP tag)
{
    SEXP ans;
    base = EnsureString(base);
    tag = EnsureString(tag);
    if (*CHAR(base) && *CHAR(tag))
    {
        ans = allocString(strlen(CHAR(tag)) + strlen(CHAR(base)) + 2);
        sprintf(CHAR(ans), "%s.%s", CHAR(base), CHAR(tag));
    }
    else if (*CHAR(tag))
    {
        ans = tag;
    }
    else if (*CHAR(base))
    {
        ans = base;
    }
    else
        ans = R_BlankString;
    return ans;
}

SEXP NewName(SEXP base, SEXP tag, int i, int n, int seqno)
{
    /* Construct a new Name/Tag, using
     *	base.tag
     *	base<seqno>	or
     *	tag
     *
     * NOTE: i,n   are NOT used currently */

    SEXP ans;
    base = EnsureString(base);
    tag = EnsureString(tag);
    if (*CHAR(base) && *CHAR(tag))
    {
        ans = allocString(strlen(CHAR(base)) + strlen(CHAR(tag)) + 1);
        sprintf(CHAR(ans), "%s.%s", CHAR(base), CHAR(tag));
    }
    else if (*CHAR(base))
    {
        ans = allocString(strlen(CHAR(base)) + IndexWidth(seqno));
        sprintf(CHAR(ans), "%s%d", CHAR(base), seqno);
    }
    else if (*CHAR(tag))
    {
        ans = allocString(strlen(CHAR(tag)));
        sprintf(CHAR(ans), "%s", CHAR(tag));
    }
    else
        ans = R_BlankString;
    return ans;
}

SEXP ItemName(SEXP names, int i)
{
    /* return  names[i]  if it is a character (>= 1 cgar), or NULL otherwise */
    if (names != R_NilValue && STRING(names)[i] != R_NilValue && CHAR(STRING(names)[i])[0] != '\0')
        return STRING(names)[i];
    else
        return R_NilValue;
}

/* NewExtractNames(v, base, tag, recurse):  For c() and	 unlist().
 * On entry, "base" is the naming component we have acquired by
 * recursing down from above.
 *	If we have a list and we are recursing, we append a new tag component
 * to the base tag (either by using the list tags, or their offsets),
 * and then we do the recursion.
 *	If we have a vector, we just create the tags for each element. */

static int count, seqno;
static int firstpos;

static void NewExtractNames(SEXP v, SEXP base, SEXP tag, int recurse)
{
    SEXP names, namei;
    int i, n, savecount = 0, saveseqno, savefirstpos = 0;

    /* If we beneath a new tag, we reset the index */
    /* sequence and create the new basename string. */

    if (tag != R_NilValue)
    {
        PROTECT(base = NewBase(base, tag));
        savefirstpos = firstpos;
        saveseqno = seqno;
        savecount = count;
        count = 0;
        seqno = 0;
        firstpos = -1;
    }
    else
        saveseqno = 0;

    n = length(v);
    names = getAttrib(v, R_NamesSymbol);

    switch (TYPEOF(v))
    {
    case NILSXP:
        break;
    case LISTSXP:
        for (i = 0; i < n; i++)
        {
            namei = ItemName(names, i);
            if (recurse)
            {
                NewExtractNames(CAR(v), base, namei, recurse);
            }
            else
            {
                if (namei == R_NilValue && count == 0)
                    firstpos = ans_nnames;
                count++;
                namei = NewName(base, namei, i, n, ++seqno);
                STRING(ans_names)[ans_nnames++] = namei;
            }
            v = CDR(v);
        }
        break;
    case VECSXP:
    case EXPRSXP:
        for (i = 0; i < n; i++)
        {
            namei = ItemName(names, i);
            if (recurse)
            {
                NewExtractNames(VECTOR(v)[i], base, namei, recurse);
            }
            else
            {
                if (namei == R_NilValue && count == 0)
                    firstpos = ans_nnames;
                count++;
                namei = NewName(base, namei, i, n, ++seqno);
                STRING(ans_names)[ans_nnames++] = namei;
            }
        }
        break;
    case LGLSXP:
    case INTSXP:
    case REALSXP:
    case CPLXSXP:
    case STRSXP:
        for (i = 0; i < n; i++)
        {
            namei = ItemName(names, i);
            if (namei == R_NilValue && count == 0)
                firstpos = ans_nnames;
            count++;
            namei = NewName(base, namei, i, n, ++seqno);
            STRING(ans_names)[ans_nnames++] = namei;
        }
        break;
    default:
        if (count == 0)
            firstpos = ans_nnames;
        count++;
        namei = NewName(base, R_NilValue, 0, 1, ++seqno);
        STRING(ans_names)[ans_nnames++] = base;
    }
    if (tag != R_NilValue)
    {
        if (firstpos >= 0 && count == 1)
            STRING(ans_names)[firstpos] = base;
        firstpos = savefirstpos;
        count = savecount;
        UNPROTECT(1);
    }
    seqno = seqno + saveseqno;
}

/* Code to extract the optional arguments to c().  We do it this */
/* way, rather than having an interpreted font-end do the job, */
/* because we want to avoid duplication at the top level. */
/* FIXME : is there another possibility? */

static SEXP ExtractOptionals(SEXP ans, int *recurse, int *usenames)
{
    SEXP a, n, r, u;
    int v;
    PROTECT(a = ans = CONS(R_NilValue, ans));
    r = install("recursive");
    u = install("use.names");
    while (a != R_NilValue && CDR(a) != R_NilValue)
    {
        n = TAG(CDR(a));
        if (n != R_NilValue && pmatch(r, n, 1))
        {
            if ((v = asLogical(CADR(a))) != NA_INTEGER)
            {
                *recurse = v;
            }
            CDR(a) = CDDR(a);
        }
        else if (n != R_NilValue && pmatch(u, n, 1))
        {
            if ((v = asLogical(CADR(a))) != NA_INTEGER)
            {
                *usenames = v;
            }
            CDR(a) = CDDR(a);
        }
        a = CDR(a);
    }
    UNPROTECT(1);
    return CDR(ans);
}

/* The change to lists based on dotted pairs has meant that it was */
/* necessary to separate the internal code for "c" and "unlist". */
/* Although the functions are quite similar, they operate on very */
/* different data structures. */

/* The major difference between the two functions is that the value of */
/* the "recursive" argument is FALSE by default for "c" and TRUE for */
/* "unlist".  In addition, "list" takes ... while "unlist" takes a single */
/* argument, and unlist has two optional arguments, while list has none. */

SEXP do_c(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans, t;
    int mode, recurse, usenames;

    checkArity(op, args);

    /* Attempt method dispatch. */

    if (DispatchOrEval(call, op, args, env, &ans, 1))
    {
        R_Visible = 1;
        return (ans);
    }
    R_Visible = 1;

    /* Method dispatch has failed; run the default code. */
    /* By default we do not recurse, but this can be over-ridden */
    /* by an optional "recursive" argument. */

    usenames = 1;
    recurse = 0;
    if (length(args) > 1)
        PROTECT(args = ExtractOptionals(ans, &recurse, &usenames));
    else
        PROTECT(args = ans);

    /* Determine the type of the returned value. */
    /* The strategy here is appropriate because the */
    /* object being operated on is a pair based list. */

    ans_flags = 0;
    ans_length = 0;
    ans_nnames = 0;

    for (t = args; t != R_NilValue; t = CDR(t))
    {
        if (usenames && !ans_nnames)
        {
            if (!isNull(TAG(t)))
                ans_nnames = 1;
            else
                ans_nnames = HasNames(CAR(t));
        }
        AnswerType(CAR(t), recurse, usenames);
    }

    /* If a non-vector argument was encountered (perhaps a list if */
    /* recursive is FALSE) then we must return a list.	Otherwise, */
    /* we use the natural coercion for vector types. */

    mode = NILSXP;
    if (ans_flags & 256)
        mode = EXPRSXP;
    else if (ans_flags & 128)
        mode = VECSXP;
    else if (ans_flags & 64)
        mode = STRSXP;
    else if (ans_flags & 32)
        mode = CPLXSXP;
    else if (ans_flags & 16)
        mode = REALSXP;
    else if (ans_flags & 8)
        mode = INTSXP;
    else if (ans_flags & 1)
        mode = LGLSXP;

    /* Allocate the return value and set up to pass through */
    /* the arguments filling in values of the returned object. */

    PROTECT(ans = allocVector(mode, ans_length));
    ans_ptr = ans;
    ans_length = 0;
    t = args;

    if (mode == VECSXP || mode == EXPRSXP)
    {
        if (!recurse)
        {
            while (args != R_NilValue)
            {
                ListAnswer(CAR(args), 0);
                args = CDR(args);
            }
        }
        else
            ListAnswer(args, recurse);
        ans_length = length(ans);
    }
    else if (mode == STRSXP)
        StringAnswer(args);
    else if (mode == CPLXSXP)
        ComplexAnswer(args);
    else if (mode == REALSXP)
        RealAnswer(args);
    else
        IntegerAnswer(args);
    args = t;

    /* Build and attach the names attribute for the returned object. */

    if (ans_nnames && ans_length > 0)
    {
#ifdef OLD
        PROTECT(ans_names = allocVector(STRSXP, ans_length));
        ans_nnames = 0;
        ExtractNames(args, recurse, 1, R_NilValue);
        setAttrib(ans, R_NamesSymbol, ans_names);
        UNPROTECT(1);
#else
        PROTECT(ans_names = allocVector(STRSXP, ans_length));
        ans_nnames = 0;
#ifdef EXPT
        if (!recurse)
        {
#endif
            while (args != R_NilValue)
            {
                seqno = 0;
                firstpos = 0;
                count = 0;
                NewExtractNames(CAR(args), R_NilValue, TAG(args), recurse);
                args = CDR(args);
            }
#ifdef EXPT
        }
        else
        {
            seqno = 0;
            firstpos = 0;
            count = 0;
            NewExtractNames(args, R_NilValue, TAG(args), recurse);
        }
#endif /*EXPT*/
        setAttrib(ans, R_NamesSymbol, ans_names);
        UNPROTECT(1);
#endif /* (not) OLD */
    }
    UNPROTECT(2);
    return ans;
} /* do_c */

SEXP do_unlist(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans, t;
    int mode, recurse, usenames;
    int i, n;

    checkArity(op, args);

    /* Attempt method dispatch. */

    if (DispatchOrEval(call, op, args, env, &ans, 1))
    {
        R_Visible = 1;
        return (ans);
    }
    R_Visible = 1;

    /* Method dispatch has failed; run the default code. */
    /* By default we recurse, but this can be over-ridden */
    /* by an optional "recursive" argument. */

    PROTECT(args = CAR(ans));
    recurse = asLogical(CADR(ans));
    usenames = asLogical(CADDR(ans));

    /* Determine the type of the returned value. */
    /* The strategy here is appropriate because the */
    /* object being operated on is a generic vector. */

    ans_flags = 0;
    ans_length = 0;
    ans_nnames = 0;

    n = 0; /* -Wall */
    if (isNewList(args))
    {
        n = length(args);
        if (usenames && getAttrib(args, R_NamesSymbol) != R_NilValue)
            ans_nnames = 1;
        for (i = 0; i < n; i++)
        {
            if (usenames && !ans_nnames)
                ans_nnames = HasNames(VECTOR(args)[i]);
            AnswerType(VECTOR(args)[i], recurse, usenames);
        }
    }
    else if (isList(args))
    {
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            if (usenames && !ans_nnames)
            {
                if (!isNull(TAG(t)))
                    ans_nnames = 1;
                else
                    ans_nnames = HasNames(CAR(t));
            }
            AnswerType(CAR(t), recurse, usenames);
        }
    }
    else
    {
        UNPROTECT(1);
        if (isVector(args))
            return args;
        else
            errorcall(call, "argument not a list\n");
    }

    /* If a non-vector argument was encountered (perhaps a list if */
    /* recursive = F) then we must return a list.  Otherwise, we use */
    /* the natural coercion for vector types. */

    mode = NILSXP;
    if (ans_flags & 128)
        mode = VECSXP;
    else if (ans_flags & 64)
        mode = STRSXP;
    else if (ans_flags & 32)
        mode = CPLXSXP;
    else if (ans_flags & 16)
        mode = REALSXP;
    else if (ans_flags & 8)
        mode = INTSXP;
    else if (ans_flags & 1)
        mode = LGLSXP;

    /* Allocate the return value and set up to pass through */
    /* the arguments filling in values of the returned object. */

    PROTECT(ans = allocVector(mode, ans_length));
    ans_ptr = ans;
    ans_length = 0;
    t = args;

    /* FIXME : The following assumes one of pair or vector */
    /* based lists applies.  It needs to handle both */

#ifdef OLD
    /* This is here only for historical interest */
    if (mode == LISTSXP)
    {
        if (!recurse)
        {
            while (args != R_NilValue)
            {
                ListAnswer(CAR(args), 0);
                args = CDR(args);
            }
        }
        else
            ListAnswer(args, recurse);
        ans_length = length(ans);
    }
#else
    if (mode == VECSXP)
    {
        if (!recurse)
        {
            for (i = 0; i < n; i++)
                ListAnswer(VECTOR(args)[i], 0);
        }
        else
            ListAnswer(args, recurse);
        ans_length = length(ans);
    }
#endif
    else if (mode == STRSXP)
        StringAnswer(args);
    else if (mode == CPLXSXP)
        ComplexAnswer(args);
    else if (mode == REALSXP)
        RealAnswer(args);
    else
        IntegerAnswer(args);
    args = t;

    /* Build and attach the names attribute for the returned object. */

    if (ans_nnames && ans_length > 0)
    {
        PROTECT(ans_names = allocVector(STRSXP, ans_length));
        if (!recurse)
        {
            if (TYPEOF(args) == VECSXP)
            {
                SEXP names = getAttrib(args, R_NamesSymbol);
                ans_nnames = 0;
                seqno = 0;
                firstpos = 0;
                count = 0;
                for (i = 0; i < n; i++)
                {
                    NewExtractNames(VECTOR(args)[i], R_NilValue, ItemName(names, i), recurse);
                }
            }
            else if (TYPEOF(args) == LISTSXP)
            {
                ans_nnames = 0;
                seqno = 0;
                firstpos = 0;
                count = 0;
                while (args != R_NilValue)
                {
                    NewExtractNames(CAR(args), R_NilValue, TAG(args), recurse);
                    args = CDR(args);
                }
            }
        }
        else
        {
            ans_nnames = 0;
            seqno = 0;
            firstpos = 0;
            count = 0;
            NewExtractNames(args, R_NilValue, R_NilValue, recurse);
        }
        setAttrib(ans, R_NamesSymbol, ans_names);
        UNPROTECT(1);
    }
    UNPROTECT(2);
    return ans;
} /* do_unlist */

#define LNAMBUF 100

static SEXP rho;

SEXP FetchMethod(char *generic, char *classname, SEXP env)
{
    char buf[LNAMBUF];
    SEXP method;
    if (strlen(generic) + strlen(classname) + 2 > LNAMBUF)
        error("class name too long in %s\n", generic);
    sprintf(buf, "%s.%s", generic, classname);
    method = findVar(install(buf), env);
    if (TYPEOF(method) != CLOSXP)
        method = R_NilValue;
    return method;
}

SEXP do_bind(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP a, t, obj, class, classlist, classname, method, classmethod;
    char *generic;
    int mode = ANYSXP;

    /* Lazy evaluation and method dispatch based on argument types are
     * fundamentally incompatible notions.  The results here are
     * ghastly.
     *
     * We build promises to evaluate the arguments and then force the
     * promises so that if we despatch to a closure below, the closure
     * is still in a position to use "substitute" to get the actual
     * expressions which generated the argument (for naming purposes).
     *
     * The dispatch rule here is as follows:
     *
     * 1) For each argument we get the list of possible class
     *    memberships from the class attribute.
     *
     * 2) We inspect each class in turn to see if there is an
     *    an applicable method.
     *
     * 3) If we find an applicable method we make sure that it is
     *    identical to any method determined for prior arguments.
     *    If it is identical, we proceed, otherwise we immediately
     *    drop through to the default code.
     */

    PROTECT(args = promiseArgs(args, env));
    mode = 0;
    generic = (PRIMVAL(op) == 1) ? "cbind" : "rbind";
    class = R_NilValue;
    method = R_NilValue;
    for (a = args; a != R_NilValue; a = CDR(a))
    {
        obj = eval(CAR(a), env);
        if (isObject(obj))
        {
            int i;
            classlist = getAttrib(obj, R_ClassSymbol);
            for (i = 0; i < length(classlist); i++)
            {
                classname = STRING(classlist)[i];
                classmethod = FetchMethod(generic, CHAR(classname), env);
                if (classmethod != R_NilValue)
                {
                    if (class == R_NilValue)
                    {
                        /* There is no previous class */
                        /* We use this method. */
                        class = classname;
                        method = classmethod;
                    }
                    else
                    {
                        /* Check compatibility with the */
                        /* previous class.  If the two are not */
                        /* compatible we drop through to the */
                        /* default method. */
                        if (strcmp(CHAR(class), CHAR(classname)))
                        {
                            method = R_NilValue;
                            break;
                        }
                    }
                }
            }
        }
    }
    if (method != R_NilValue)
    {
        PROTECT(method);
        args = applyClosure(call, method, args, env, R_NilValue);
        UNPROTECT(2);
        return args;
    }

    /* Despatch based on class membership has failed. */
    /* The default code for rbind/cbind.default follows */
    /* First, extract the evaluated arguments. */
    for (a = args; a != R_NilValue; a = CDR(a))
        CAR(a) = PRVALUE(CAR(a));

    rho = env; /* GLOBAL */
    ans_flags = 0;
    ans_length = 0;
    ans_nnames = 0;

    for (t = args; t != R_NilValue; t = CDR(t))
        AnswerType(CAR(t), 0, 0);

    /* zero-extent matrices shouldn't give NULL, anymore :
       if (ans_length == 0)
       return R_NilValue;
    */
    if (ans_flags >= 128)
    {
        if (ans_flags & 128)
            mode = LISTSXP;
    }
    else if (ans_flags >= 64)
    {
        if (ans_flags & 64)
            mode = STRSXP;
    }
    else
    {
        if (ans_flags & 1)
            mode = LGLSXP;
        if (ans_flags & 8)
            mode = INTSXP;
        if (ans_flags & 16)
            mode = REALSXP;
        if (ans_flags & 32)
            mode = CPLXSXP;
    }

    switch (mode)
    {
    case NILSXP:
    case LGLSXP:
    case INTSXP:
    case REALSXP:
    case CPLXSXP:
    case STRSXP:
        break;
    default:
        errorcall(call, "cannot create a matrix from these types\n");
    }

    if (PRIMVAL(op) == 1)
        a = cbind(call, args, mode);
    else
        a = rbind(call, args, mode);
    UNPROTECT(1);
    return a;
}

static void SetRowNames(SEXP dimnames, SEXP x)
{
    if (TYPEOF(dimnames) == VECSXP)
        VECTOR(dimnames)[0] = x;
    else if (TYPEOF(dimnames) == LISTSXP)
        CAR(dimnames) = x;
}

static void SetColNames(SEXP dimnames, SEXP x)
{
    if (TYPEOF(dimnames) == VECSXP)
        VECTOR(dimnames)[1] = x;
    else if (TYPEOF(dimnames) == LISTSXP)
        CADR(dimnames) = x;
}

static SEXP cbind(SEXP call, SEXP args, SEXPTYPE mode)
{
    int i, j, k, idx, n;
    int have_rnames, have_cnames;
    int nnames, mnames;
    int rows, cols, mrows;
    int warned;
    SEXP dn, t, u, result, dims;

    have_rnames = 0;
    have_cnames = 0;
    nnames = 0;
    mnames = 0;
    rows = 0;
    cols = 0;
    /* mrows = 0;*/
    mrows = -1;

    /* check conformability of matrix arguments */

    n = 0;
    for (t = args; t != R_NilValue; t = CDR(t))
    {
        if (length(CAR(t)) >= 0)
        {
            dims = getAttrib(CAR(t), R_DimSymbol);
            if (length(dims) == 2)
            {
                if (mrows == -1)
                    mrows = INTEGER(dims)[0];
                else if (mrows != INTEGER(dims)[0])
                    errorcall(call, "number of rows of matrices must match (see arg %d)\n", n + 1);
                cols += INTEGER(dims)[1];
            }
            else if (length(CAR(t)) > 0)
            {
                rows = imax2(rows, length(CAR(t)));
                cols += 1;
            }
        }
        n++;
    }
    if (mrows != -1)
        rows = mrows;

    /* Check conformability of vector arguments. -- Look for dimnames. */

    n = 0;
    warned = 0;
    for (t = args; t != R_NilValue; t = CDR(t))
    {
        n++;
        if (length(CAR(t)) >= 0)
        {
            dims = getAttrib(CAR(t), R_DimSymbol);
            if (length(dims) == 2)
            {
                dn = getAttrib(CAR(t), R_DimNamesSymbol);
                if (length(dn) == 2)
                {
                    if (VECTOR(dn)[1] != R_NilValue)
                        have_cnames = 1;
                    if (VECTOR(dn)[0] != R_NilValue)
                        mnames = mrows;
                }
            }
            else
            {
                k = length(CAR(t));
                if (!warned && k > 0 && (k > rows || rows % k))
                {
                    warned = 1;
                    PROTECT(call = substituteList(call, rho));
                    warningcall(call, "number of rows of result\n\tis not a multiple of vector length (arg %d)", n);
                    UNPROTECT(1);
                }
                dn = getAttrib(CAR(t), R_NamesSymbol);
                if (TAG(t) != R_NilValue)
                    have_cnames = 1;
                nnames = imax2(nnames, length(dn));
            }
        }
    }
    if (mnames || nnames == rows)
        have_rnames = 1;

    PROTECT(result = allocMatrix(mode, rows, cols));
    n = 0;

    if (mode == STRSXP)
    {
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            if (length(CAR(t)) > 0)
            { /* rbind() has ">=" here */
                u = CAR(t) = coerceVector(CAR(t), STRSXP);
                k = LENGTH(u);
                idx = (!isMatrix(CAR(t))) ? rows : k;
                for (i = 0; i < idx; i++)
                    STRING(result)[n++] = STRING(u)[i % k];
            }
        }
    }
    else if (mode == CPLXSXP)
    {
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            if (length(CAR(t)) > 0)
            {
                u = CAR(t) = coerceVector(CAR(t), CPLXSXP);
                k = LENGTH(u);
                idx = (!isMatrix(CAR(t))) ? rows : k;
                for (i = 0; i < idx; i++)
                    COMPLEX(result)[n++] = COMPLEX(u)[i % k];
            }
        }
    }
    else
    {
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            if (length(CAR(t)) > 0)
            {
                u = CAR(t);
                k = LENGTH(u);
                idx = (!isMatrix(CAR(t))) ? rows : k;
                if (TYPEOF(u) <= INTSXP)
                {
                    if (mode <= INTSXP)
                    {
                        for (i = 0; i < idx; i++)
                            INTEGER(result)[n++] = INTEGER(u)[i % k];
                    }
                    else
                    {
                        for (i = 0; i < idx; i++)
                            REAL(result)[n++] = (INTEGER(u)[i % k]) == NA_INTEGER ? NA_REAL : INTEGER(u)[i % k];
                    }
                }
                else
                {
                    for (i = 0; i < idx; i++)
                        REAL(result)[n++] = REAL(u)[i % k];
                }
            }
        }
    }

    /* Adjustment of dimnames attributes. */
    if (have_cnames | have_rnames)
    {
        SEXP nam, tnam;
        PROTECT(dn = allocVector(VECSXP, 2));
        if (have_cnames)
            nam = VECTOR(dn)[1] = allocVector(STRSXP, cols);
        else
            nam = R_NilValue; /* -Wall */
        j = 0;
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            if (length(CAR(t)) >= 0)
            {

                if (isMatrix(CAR(t)))
                {

                    u = getAttrib(CAR(t), R_DimNamesSymbol);
                    tnam = GetColNames(u);

                    if (have_rnames && GetRowNames(dn) == R_NilValue && GetRowNames(u) != R_NilValue)
                        SetRowNames(dn, duplicate(GetRowNames(u)));

                    /* rbind() does this only  if(have_?names) .. : */
                    if (tnam != R_NilValue)
                    {
                        for (i = 0; i < length(tnam); i++)
                            STRING(nam)[j++] = STRING(tnam)[i];
                    }
                    else if (have_cnames)
                    {
                        for (i = 0; i < ncols(CAR(t)); i++)
                            STRING(nam)[j++] = R_BlankString;
                    }
                }
                else if (length(CAR(t)) > 0)
                {

                    u = getAttrib(CAR(t), R_NamesSymbol);
                    tnam = u;

                    if (have_rnames && GetRowNames(dn) == R_NilValue && u != R_NilValue && length(u) == rows)
                        SetRowNames(dn, duplicate(u));

                    if (TAG(t) != R_NilValue)
                    {
                        STRING(nam)[j++] = PRINTNAME(TAG(t));
                    }
                    else if (have_cnames)
                    {
                        STRING(nam)[j++] = R_BlankString;
                    }
                }
            }
        }
        setAttrib(result, R_DimNamesSymbol, dn);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return result;
} /* cbind */

static SEXP rbind(SEXP call, SEXP args, SEXPTYPE mode)
{
    int i, j, k, n;
    int have_rnames, have_cnames;
    int nnames, mnames;
    int rows, cols, mcols, mrows;
    int warned;
    SEXP dn, t, u, result, dims;

    have_rnames = 0;
    have_cnames = 0;
    nnames = 0;
    mnames = 0;
    rows = 0;
    cols = 0;
    mcols = 0;

    /* check conformability of matrix arguments */

    n = 0;
    for (t = args; t != R_NilValue; t = CDR(t))
    {
        if (length(CAR(t)) >= 0)
        {
            dims = getAttrib(CAR(t), R_DimSymbol);
            if (length(dims) == 2)
            {
                if (mcols == 0)
                    mcols = INTEGER(dims)[1];
                else if (mcols != INTEGER(dims)[1])
                    errorcall(call, "number of columns of matrices must match (see arg %d)\n", n + 1);
                rows += INTEGER(dims)[0];
            }
            else if (length(CAR(t)) > 0)
            {
                cols = imax2(cols, length(CAR(t)));
                rows += 1;
            }
        }
        n++;
    }
    if (mcols != 0)
        cols = mcols;

    /* Check conformability of vector arguments. -- Look for dimnames. */

    n = 0;
    warned = 0;
    for (t = args; t != R_NilValue; t = CDR(t))
    {
        n++;
        if (length(CAR(t)) >= 0)
        {
            dims = getAttrib(CAR(t), R_DimSymbol);
            if (length(dims) == 2)
            {
                dn = getAttrib(CAR(t), R_DimNamesSymbol);
                if (length(dn) == 2)
                {
                    if (VECTOR(dn)[0] != R_NilValue)
                        have_rnames = 1;
                    if (VECTOR(dn)[1] != R_NilValue)
                        mnames = mcols;
                }
            }
            else
            {
                k = length(CAR(t));
                if (!warned && k > 0 && (k > cols || cols % k))
                {
                    warned = 1;
                    PROTECT(call = substituteList(call, rho));
                    warningcall(call, "number of columns of result\n\tnot a multiple of vector length (arg %d)", n);
                    UNPROTECT(1);
                }
                dn = getAttrib(CAR(t), R_NamesSymbol);
                if (TAG(t) != R_NilValue)
                    have_rnames = 1;
                nnames = imax2(nnames, length(dn));
            }
        }
    }
    if (mnames || nnames == cols)
        have_cnames = 1;

    PROTECT(result = allocMatrix(mode, rows, cols));
    n = 0;

    if (mode == STRSXP)
    {
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            if (length(CAR(t)) >= 0)
            { /* cbind() has ">" here */
                u = CAR(t) = coerceVector(CAR(t), STRSXP);
                k = LENGTH(u);
                mrows = (isMatrix(u)) ? nrows(u) : 1;
                if (k == 0)
                    mrows = 0;
                for (i = 0; i < mrows; i++)
                    for (j = 0; j < cols; j++)
                        STRING(result)[i + n + (j * rows)] = STRING(u)[(i + j * mrows) % k];
                n += mrows;
            }
        }
        /* 0.63.3 erroneously:	return result; */
    }
    else if (mode == CPLXSXP)
    {
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            if (length(CAR(t)) > 0)
            {
                u = CAR(t) = coerceVector(CAR(t), CPLXSXP);
                k = LENGTH(u);
                mrows = (isMatrix(u)) ? nrows(u) : 1;
                for (i = 0; i < mrows; i++)
                    for (j = 0; j < cols; j++)
                        COMPLEX(result)[i + n + (j * rows)] = COMPLEX(u)[(i + j * mrows) % k];
                n += mrows;
            }
        }
        /* 0.63.3 erroneously:	return result; */
    }
    else
    {
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            if (length(CAR(t)) >= 0)
            {
                u = CAR(t);
                k = LENGTH(u);
                /* mrows = (isMatrix(u)) ? nrows(u) : 1; */
                if (isMatrix(u))
                {
                    mrows = nrows(u);
                }
                else
                {
                    if (length(u) == 0)
                        mrows = 0;
                    else
                        mrows = 1;
                }
                if (TYPEOF(u) <= INTSXP)
                {
                    if (mode <= INTSXP)
                    {
                        for (i = 0; i < mrows; i++)
                            for (j = 0; j < cols; j++)
                                INTEGER(result)[i + n + (j * rows)] = INTEGER(u)[(i + j * mrows) % k];
                        n += mrows;
                    }
                    else
                    {
                        for (i = 0; i < mrows; i++)
                            for (j = 0; j < cols; j++)
                                REAL(result)
                                [i + n + (j * rows)] =
                                    (INTEGER(u)[(i + j * mrows) % k]) == NA_INTEGER ? NA_REAL
                                                                                    : INTEGER(u)[(i + j * mrows) % k];
                        n += mrows;
                    }
                }
                else
                {
                    for (i = 0; i < mrows; i++)
                        for (j = 0; j < cols; j++)
                            REAL(result)[i + n + (j * rows)] = REAL(u)[(i + j * mrows) % k];
                    n += mrows;
                }
            }
        }
    }
    /* Adjustment of dimnames attributes. */
    if (have_rnames | have_cnames)
    {
        SEXP nam, tnam;
        PROTECT(dn = allocVector(VECSXP, 2));
        if (have_rnames)
            nam = VECTOR(dn)[0] = allocVector(STRSXP, rows);
        else
            nam = R_NilValue; /* -Wall */
        j = 0;
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            if (length(CAR(t)) >= 0)
            {

                if (isMatrix(CAR(t)))
                {

                    u = getAttrib(CAR(t), R_DimNamesSymbol);
                    tnam = GetRowNames(u);

                    if (have_cnames && GetColNames(dn) == R_NilValue && GetColNames(u) != R_NilValue)
                        SetColNames(dn, duplicate(GetColNames(u)));

                    /* cbind() doesn't test have_?names BEFORE tnam!=Nil..:*/
                    if (have_rnames)
                    {
                        if (tnam != R_NilValue)
                        {
                            for (i = 0; i < length(tnam); i++)
                                STRING(nam)[j++] = STRING(tnam)[i];
                        }
                        else
                        {
                            for (i = 0; i < nrows(CAR(t)); i++)
                                STRING(nam)[j++] = R_BlankString;
                        }
                    }
                }
                else if (length(CAR(t)) > 0)
                {

                    u = getAttrib(CAR(t), R_NamesSymbol);
                    tnam = u;

                    if (have_cnames && GetColNames(dn) == R_NilValue && u != R_NilValue && length(u) == cols)
                        SetColNames(dn, duplicate(u));

                    if (TAG(t) != R_NilValue)
                    {
                        STRING(nam)[j++] = PRINTNAME(TAG(t));
                    }
                    else if (have_rnames)
                    {
                        STRING(nam)[j++] = R_BlankString;
                    }
                }
            }
        }
        setAttrib(result, R_DimNamesSymbol, dn);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return result;
} /* rbind */
