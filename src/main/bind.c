/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997--2002  Robert Gentleman, Ross Ihaka and the
 *			      R Development Core Team
 *  Copyright (C) 2002-3      The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available via WWW at
 *  http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
 *  writing to the Free Software Foundation, Inc., 59 Temple Place,
 *  Suite 330, Boston, MA  02111-1307  USA.
 */

/* Code to handle list / vector switch */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Defn.h>
#include <Rmath.h> /* imax2 */

#define LIST_ASSIGN(x)                                                                                                 \
    {                                                                                                                  \
        SET_VECTOR_ELT(data->ans_ptr, data->ans_length, x);                                                            \
        data->ans_length++;                                                                                            \
    }

static SEXP cbind(SEXP, SEXP, SEXPTYPE, SEXP);
static SEXP rbind(SEXP, SEXP, SEXPTYPE, SEXP);

/* The following code establishes the return type for the */
/* functions  unlist, c, cbind, and rbind and also determines */
/* whether the returned object is to have a names attribute. */

struct BindData
{
    int ans_flags;
    SEXP ans_ptr;
    int ans_length;
    SEXP ans_names;
    int ans_nnames;
    int deparse_level; /* Initialize to 1. */
};

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

static void AnswerType(SEXP x, int recurse, int usenames, struct BindData *data)
{
    switch (TYPEOF(x))
    {
    case NILSXP:
        break;
    case LGLSXP:
        data->ans_flags |= 1;
        data->ans_length += LENGTH(x);
        break;
    case INTSXP:
        data->ans_flags |= 8;
        data->ans_length += LENGTH(x);
        break;
    case REALSXP:
        data->ans_flags |= 16;
        data->ans_length += LENGTH(x);
        break;
    case CPLXSXP:
        data->ans_flags |= 32;
        data->ans_length += LENGTH(x);
        break;
    case STRSXP:
        data->ans_flags |= 64;
        data->ans_length += LENGTH(x);
        break;
    case VECSXP:
    case EXPRSXP:
        if (recurse)
        {
            int i, n;
            n = length(x);
            if (usenames && !data->ans_names && !isNull(getAttrib(x, R_NamesSymbol)))
                data->ans_nnames = 1;
            for (i = 0; i < n; i++)
            {
                if (usenames && !data->ans_nnames)
                    data->ans_nnames = HasNames(VECTOR_ELT(x, i));
                AnswerType(VECTOR_ELT(x, i), recurse, usenames, data);
            }
        }
        else
        {
            if (TYPEOF(x) == EXPRSXP)
                data->ans_flags |= 256;
            else
                data->ans_flags |= 128;
            data->ans_length += length(x);
        }
        break;
    case LISTSXP:
        if (recurse)
        {
            while (x != R_NilValue)
            {
                if (usenames && !data->ans_nnames)
                {
                    if (!isNull(TAG(x)))
                        data->ans_nnames = 1;
                    else
                        data->ans_nnames = HasNames(CAR(x));
                }
                AnswerType(CAR(x), recurse, usenames, data);
                x = CDR(x);
            }
        }
        else
        {
            data->ans_flags |= 128;
            data->ans_length += length(x);
        }
        break;
    default:
        data->ans_flags |= 128;
        data->ans_length += 1;
        break;
    }
}

/* The following functions are used to coerce arguments to */
/* the appropriate type for inclusion in the returned value. */

static void ListAnswer(SEXP x, int recurse, struct BindData *data)
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
            LIST_ASSIGN(ScalarString(STRING_ELT(x, i)));
        break;
    case VECSXP:
    case EXPRSXP:
        if (recurse)
        {
            for (i = 0; i < LENGTH(x); i++)
                ListAnswer(VECTOR_ELT(x, i), recurse, data);
        }
        else
        {
            for (i = 0; i < LENGTH(x); i++)
                LIST_ASSIGN(duplicate(VECTOR_ELT(x, i)));
        }
        break;
    case LISTSXP:
        if (recurse)
        {
            while (x != R_NilValue)
            {
                ListAnswer(CAR(x), recurse, data);
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

static void StringAnswer(SEXP x, struct BindData *data)
{
    int i, n;
    switch (TYPEOF(x))
    {
    case NILSXP:
        break;
    case LISTSXP:
        while (x != R_NilValue)
        {
            StringAnswer(CAR(x), data);
            x = CDR(x);
        }
        break;
    case VECSXP:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            StringAnswer(VECTOR_ELT(x, i), data);
        break;
    default:
        PROTECT(x = coerceVector(x, STRSXP));
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            SET_STRING_ELT(data->ans_ptr, data->ans_length++, STRING_ELT(x, i));
        UNPROTECT(1);
        break;
    }
}

static void IntegerAnswer(SEXP x, struct BindData *data)
{
    int i, n;
    switch (TYPEOF(x))
    {
    case NILSXP:
        break;
    case LISTSXP:
        while (x != R_NilValue)
        {
            IntegerAnswer(CAR(x), data);
            x = CDR(x);
        }
        break;
    case VECSXP:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            IntegerAnswer(VECTOR_ELT(x, i), data);
        break;
    default:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            INTEGER(data->ans_ptr)[data->ans_length++] = INTEGER(x)[i];
        break;
    }
}

static void RealAnswer(SEXP x, struct BindData *data)
{
    int i, n, xi;
    switch (TYPEOF(x))
    {
    case NILSXP:
        break;
    case LISTSXP:
        while (x != R_NilValue)
        {
            RealAnswer(CAR(x), data);
            x = CDR(x);
        }
        break;
    case VECSXP:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            RealAnswer(VECTOR_ELT(x, i), data);
        break;
    case REALSXP:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            REAL(data->ans_ptr)[data->ans_length++] = REAL(x)[i];
        break;
    default:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
        {
            xi = INTEGER(x)[i];
            if (xi == NA_INTEGER)
                REAL(data->ans_ptr)[data->ans_length++] = NA_REAL;
            else
                REAL(data->ans_ptr)[data->ans_length++] = xi;
        }
        break;
    }
}

static void ComplexAnswer(SEXP x, struct BindData *data)
{
    int i, n, xi;
    switch (TYPEOF(x))
    {
    case NILSXP:
        break;
    case LISTSXP:
        while (x != R_NilValue)
        {
            ComplexAnswer(CAR(x), data);
            x = CDR(x);
        }
        break;
    case VECSXP:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            ComplexAnswer(VECTOR_ELT(x, i), data);
        break;
    case REALSXP:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
        {
            COMPLEX(data->ans_ptr)[data->ans_length].r = REAL(x)[i];
            COMPLEX(data->ans_ptr)[data->ans_length].i = 0.0;
            data->ans_length++;
        }
        break;
    case CPLXSXP:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
            COMPLEX(data->ans_ptr)[data->ans_length++] = COMPLEX(x)[i];
        break;
    default:
        n = LENGTH(x);
        for (i = 0; i < n; i++)
        {
            xi = INTEGER(x)[i];
            if (xi == NA_INTEGER)
            {
                COMPLEX(data->ans_ptr)[data->ans_length].r = NA_REAL;
                COMPLEX(data->ans_ptr)[data->ans_length].i = NA_REAL;
            }
            else
            {
                COMPLEX(data->ans_ptr)[data->ans_length].r = xi;
                COMPLEX(data->ans_ptr)[data->ans_length].i = 0.0;
            }
            data->ans_length++;
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
        ans = allocString(strlen(CHAR(tag)) + strlen(CHAR(base)) + 1);
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

static SEXP NewName(SEXP base, SEXP tag, int i, int n, int seqno)
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
        if (tag == NA_STRING)
            ans = NA_STRING;
        else
        {
            ans = allocString(strlen(CHAR(tag)));
            sprintf(CHAR(ans), "%s", CHAR(tag));
        }
    }
    else
        ans = R_BlankString;
    return ans;
}

SEXP ItemName(SEXP names, int i)
{
    /* return  names[i]  if it is a character (>= 1 char), or NULL otherwise */
    if (names != R_NilValue && STRING_ELT(names, i) != R_NilValue && CHAR(STRING_ELT(names, i))[0] != '\0')
        return STRING_ELT(names, i);
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

struct NameData
{
    int count;
    int seqno;
    int firstpos;
};

static void NewExtractNames(SEXP v, SEXP base, SEXP tag, int recurse, struct BindData *data, struct NameData *nameData)
{
    SEXP names, namei;
    int i, n, savecount = 0, saveseqno, savefirstpos = 0;

    /* If we beneath a new tag, we reset the index */
    /* sequence and create the new basename string. */

    if (tag != R_NilValue)
    {
        PROTECT(base = NewBase(base, tag));
        savefirstpos = nameData->firstpos;
        saveseqno = nameData->seqno;
        savecount = nameData->count;
        nameData->count = 0;
        nameData->seqno = 0;
        nameData->firstpos = -1;
    }
    else
        saveseqno = 0;

    n = length(v);
    PROTECT(names = getAttrib(v, R_NamesSymbol));

    switch (TYPEOF(v))
    {
    case NILSXP:
        break;
    case LISTSXP:
        for (i = 0; i < n; i++)
        {
            PROTECT(namei = ItemName(names, i));
            if (recurse)
            {
                NewExtractNames(CAR(v), base, namei, recurse, data, nameData);
            }
            else
            {
                if (namei == R_NilValue && nameData->count == 0)
                    nameData->firstpos = data->ans_nnames;
                nameData->count++;
                namei = NewName(base, namei, i, n, ++(nameData->seqno));
                SET_STRING_ELT(data->ans_names, (data->ans_nnames)++, namei);
            }
            v = CDR(v);
            UNPROTECT(1); /*namei*/
        }
        break;
    case VECSXP:
    case EXPRSXP:
        for (i = 0; i < n; i++)
        {
            namei = ItemName(names, i);
            if (recurse)
            {
                NewExtractNames(VECTOR_ELT(v, i), base, namei, recurse, data, nameData);
            }
            else
            {
                if (namei == R_NilValue && nameData->count == 0)
                    nameData->firstpos = data->ans_nnames;
                nameData->count++;
                namei = NewName(base, namei, i, n, ++(nameData->seqno));
                SET_STRING_ELT(data->ans_names, (data->ans_nnames)++, namei);
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
            if (namei == R_NilValue && nameData->count == 0)
                nameData->firstpos = data->ans_nnames;
            nameData->count++;
            namei = NewName(base, namei, i, n, ++(nameData->seqno));
            SET_STRING_ELT(data->ans_names, (data->ans_nnames)++, namei);
        }
        break;
    default:
        if (nameData->count == 0)
            nameData->firstpos = data->ans_nnames;
        nameData->count++;
        namei = NewName(base, R_NilValue, 0, 1, ++(nameData->seqno));
        SET_STRING_ELT(data->ans_names, (data->ans_nnames)++, namei);
    }
    if (tag != R_NilValue)
    {
        if (nameData->firstpos >= 0 && nameData->count == 1)
            SET_STRING_ELT(data->ans_names, nameData->firstpos, base);
        nameData->firstpos = savefirstpos;
        nameData->count = savecount;
        UNPROTECT(1);
    }
    UNPROTECT(1); /*names*/
    nameData->seqno = nameData->seqno + saveseqno;
}

/* Code to extract the optional arguments to c().  We do it this */
/* way, rather than having an interpreted front-end do the job, */
/* because we want to avoid duplication at the top level. */
/* FIXME : is there another possibility? */

static SEXP ExtractOptionals(SEXP ans, int *recurse, int *usenames)
{
    SEXP a, n, last = NULL, next = NULL;
    int v, n_recurse = 0, n_usenames = 0;

    for (a = ans; a != R_NilValue; a = next)
    {
        n = TAG(a);
        next = CDR(a);
        if (n != R_NilValue && pmatch(R_RecursiveSymbol, n, 1))
        {
            if (n_recurse++ == 1)
                error("repeated formal argument 'recursive'");
            if ((v = asLogical(CAR(a))) != NA_INTEGER)
            {
                *recurse = v;
            }
            if (last == NULL)
                last = ans = next;
            else
                SETCDR(last, next);
        }
        else if (n != R_NilValue && pmatch(R_UseNamesSymbol, n, 1))
        {
            if (n_usenames++ == 1)
                error("repeated formal argument 'use.names'");
            if ((v = asLogical(CAR(a))) != NA_INTEGER)
            {
                *usenames = v;
            }
            if (last == NULL)
                last = ans = next;
            else
                SETCDR(last, next);
        }
        else
            last = a;
    }
    return ans;
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
    SEXP ans;

    checkArity(op, args);

    /* Attempt method dispatch. */

    if (DispatchOrEval(call, op, "c", args, env, &ans, 1, 0))
    {
        R_Visible = 1;
        return (ans);
    }
    return do_c_dflt(call, op, ans, env);
}

SEXP do_c_dflt(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans, t;
    int mode, recurse, usenames;
    struct BindData data;
    struct NameData nameData;

    data.deparse_level = 1; /* Initialize this early. */
    R_Visible = 1;

    /* Method dispatch has failed; run the default code. */
    /* By default we do not recurse, but this can be over-ridden */
    /* by an optional "recursive" argument. */

    usenames = 1;
    recurse = 0;
    /* this was only done for length(args) > 1 prior to 1.5.0,
       _but_ `recursive' might be the only argument */
    PROTECT(args = ExtractOptionals(args, &recurse, &usenames));

    /* Determine the type of the returned value. */
    /* The strategy here is appropriate because the */
    /* object being operated on is a pair based list. */

    data.ans_flags = 0;
    data.ans_length = 0;
    data.ans_nnames = 0;

    for (t = args; t != R_NilValue; t = CDR(t))
    {
        if (usenames && !data.ans_nnames)
        {
            if (!isNull(TAG(t)))
                data.ans_nnames = 1;
            else
                data.ans_nnames = HasNames(CAR(t));
        }
        AnswerType(CAR(t), recurse, usenames, &data);
    }

    /* If a non-vector argument was encountered (perhaps a list if */
    /* recursive is FALSE) then we must return a list.	Otherwise, */
    /* we use the natural coercion for vector types. */

    mode = NILSXP;
    if (data.ans_flags & 256)
        mode = EXPRSXP;
    else if (data.ans_flags & 128)
        mode = VECSXP;
    else if (data.ans_flags & 64)
        mode = STRSXP;
    else if (data.ans_flags & 32)
        mode = CPLXSXP;
    else if (data.ans_flags & 16)
        mode = REALSXP;
    else if (data.ans_flags & 8)
        mode = INTSXP;
    else if (data.ans_flags & 1)
        mode = LGLSXP;

    /* Allocate the return value and set up to pass through */
    /* the arguments filling in values of the returned object. */

    PROTECT(ans = allocVector(mode, data.ans_length));
    data.ans_ptr = ans;
    data.ans_length = 0;
    t = args;

    if (mode == VECSXP || mode == EXPRSXP)
    {
        if (!recurse)
        {
            while (args != R_NilValue)
            {
                ListAnswer(CAR(args), 0, &data);
                args = CDR(args);
            }
        }
        else
            ListAnswer(args, recurse, &data);
        data.ans_length = length(ans);
    }
    else if (mode == STRSXP)
        StringAnswer(args, &data);
    else if (mode == CPLXSXP)
        ComplexAnswer(args, &data);
    else if (mode == REALSXP)
        RealAnswer(args, &data);
    else
        IntegerAnswer(args, &data);
    args = t;

    /* Build and attach the names attribute for the returned object. */

    if (data.ans_nnames && data.ans_length > 0)
    {
        PROTECT(data.ans_names = allocVector(STRSXP, data.ans_length));
        data.ans_nnames = 0;
#ifdef EXPT
        if (!recurse)
        {
#endif
            while (args != R_NilValue)
            {
                nameData.seqno = 0;
                nameData.firstpos = 0;
                nameData.count = 0;
                NewExtractNames(CAR(args), R_NilValue, TAG(args), recurse, &data, &nameData);
                args = CDR(args);
            }
#ifdef EXPT
        }
        else
        {
            nameData.seqno = 0;
            nameData.firstpos = 0;
            nameData.count = 0;
            NewExtractNames(args, R_NilValue, TAG(args), recurse);
        }
#endif /*EXPT*/
        setAttrib(ans, R_NamesSymbol, data.ans_names);
        UNPROTECT(1);
    }
    UNPROTECT(2);
    return ans;
} /* do_c */

SEXP do_unlist(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans, t;
    int mode, recurse, usenames;
    int i, n;
    struct BindData data;
    struct NameData nameData;
    data.deparse_level = 1;
    checkArity(op, args);

    /* Attempt method dispatch. */

    if (DispatchOrEval(call, op, "unlist", args, env, &ans, 1, 0))
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

    data.ans_flags = 0;
    data.ans_length = 0;
    data.ans_nnames = 0;

    n = 0; /* -Wall */
    if (isNewList(args))
    {
        n = length(args);
        if (usenames && getAttrib(args, R_NamesSymbol) != R_NilValue)
            data.ans_nnames = 1;
        for (i = 0; i < n; i++)
        {
            if (usenames && !data.ans_nnames)
                data.ans_nnames = HasNames(VECTOR_ELT(args, i));
            AnswerType(VECTOR_ELT(args, i), recurse, usenames, &data);
        }
    }
    else if (isList(args))
    {
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            if (usenames && !data.ans_nnames)
            {
                if (!isNull(TAG(t)))
                    data.ans_nnames = 1;
                else
                    data.ans_nnames = HasNames(CAR(t));
            }
            AnswerType(CAR(t), recurse, usenames, &data);
        }
    }
    else
    {
        UNPROTECT(1);
        if (isVector(args))
            return args;
        else
            errorcall(call, "argument not a list");
    }

    /* If a non-vector argument was encountered (perhaps a list if */
    /* recursive = F) then we must return a list.  Otherwise, we use */
    /* the natural coercion for vector types. */

    mode = NILSXP;
    if (data.ans_flags & 128)
        mode = VECSXP;
    else if (data.ans_flags & 64)
        mode = STRSXP;
    else if (data.ans_flags & 32)
        mode = CPLXSXP;
    else if (data.ans_flags & 16)
        mode = REALSXP;
    else if (data.ans_flags & 8)
        mode = INTSXP;
    else if (data.ans_flags & 1)
        mode = LGLSXP;

    /* Allocate the return value and set up to pass through */
    /* the arguments filling in values of the returned object. */

    PROTECT(ans = allocVector(mode, data.ans_length));
    data.ans_ptr = ans;
    data.ans_length = 0;
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
        data.ans_length = length(ans);
    }
#else
    if (mode == VECSXP)
    {
        if (!recurse)
        {
            for (i = 0; i < n; i++)
                ListAnswer(VECTOR_ELT(args, i), 0, &data);
        }
        else
            ListAnswer(args, recurse, &data);
        data.ans_length = length(ans);
    }
#endif
    else if (mode == STRSXP)
        StringAnswer(args, &data);
    else if (mode == CPLXSXP)
        ComplexAnswer(args, &data);
    else if (mode == REALSXP)
        RealAnswer(args, &data);
    else
        IntegerAnswer(args, &data);
    args = t;

    /* Build and attach the names attribute for the returned object. */

    if (data.ans_nnames && data.ans_length > 0)
    {
        PROTECT(data.ans_names = allocVector(STRSXP, data.ans_length));
        if (!recurse)
        {
            if (TYPEOF(args) == VECSXP)
            {
                SEXP names = getAttrib(args, R_NamesSymbol);
                data.ans_nnames = 0;
                nameData.seqno = 0;
                nameData.firstpos = 0;
                nameData.count = 0;
                for (i = 0; i < n; i++)
                {
                    NewExtractNames(VECTOR_ELT(args, i), R_NilValue, ItemName(names, i), recurse, &data, &nameData);
                }
            }
            else if (TYPEOF(args) == LISTSXP)
            {
                data.ans_nnames = 0;
                nameData.seqno = 0;
                nameData.firstpos = 0;
                nameData.count = 0;
                while (args != R_NilValue)
                {
                    NewExtractNames(CAR(args), R_NilValue, TAG(args), recurse, &data, &nameData);
                    args = CDR(args);
                }
            }
        }
        else
        {
            data.ans_nnames = 0;
            nameData.seqno = 0;
            nameData.firstpos = 0;
            nameData.count = 0;
            NewExtractNames(args, R_NilValue, R_NilValue, recurse, &data, &nameData);
        }
        setAttrib(ans, R_NamesSymbol, data.ans_names);
        UNPROTECT(1);
    }
    UNPROTECT(2);
    return ans;
} /* do_unlist */

#define LNAMBUF 100

SEXP FetchMethod(char *generic, char *classname, SEXP env)
{
    char buf[LNAMBUF];
    SEXP method;
    if (strlen(generic) + strlen(classname) + 2 > LNAMBUF)
        error("class name too long in %s", generic);
    sprintf(buf, "%s.%s", generic, classname);
    method = findVar(install(buf), env);
    if (TYPEOF(method) == PROMSXP)
        method = eval(method, env);
    if (TYPEOF(method) != CLOSXP)
        method = R_NilValue;
    return method;
}

/* cbind(...) and rbind(...) : */
SEXP do_bind(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP a, t, obj, class, classlist, classname, method, classmethod, rho;
    char *generic;
    int mode;
    struct BindData data;
    data.deparse_level = 1;

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
     *	  memberships from the class attribute.
     *
     * 2) We inspect each class in turn to see if there is an
     *	  an applicable method.
     *
     * 3) If we find an applicable method we make sure that it is
     *	  identical to any method determined for prior arguments.
     *	  If it is identical, we proceed, otherwise we immediately
     *	  drop through to the default code.
     */

    PROTECT(args = promiseArgs(args, env));
    generic = (PRIMVAL(op) == 1) ? "cbind" : "rbind";
    class = R_NilValue;
    method = R_NilValue;
    for (a = args; a != R_NilValue; a = CDR(a))
    {
        PROTECT(obj = eval(CAR(a), env));
        if (isObject(obj))
        {
            int i;
            classlist = getAttrib(obj, R_ClassSymbol);
            for (i = 0; i < length(classlist); i++)
            {
                classname = STRING_ELT(classlist, i);
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
        UNPROTECT(1);
    }
    if (method != R_NilValue)
    {
        PROTECT(method);
        args = applyClosure(call, method, args, env, R_NilValue);
        UNPROTECT(2);
        return args;
    }

    /* Dispatch based on class membership has failed. */
    /* The default code for rbind/cbind.default follows */
    /* First, extract the evaluated arguments. */
#ifdef OLD
    for (a = args; a != R_NilValue; a = CDR(a))
        CAR(a) = PRVALUE(CAR(a));
#endif

    rho = env;
    data.ans_flags = 0;
    data.ans_length = 0;
    data.ans_nnames = 0;
    for (t = args; t != R_NilValue; t = CDR(t))
        AnswerType(PRVALUE(CAR(t)), 0, 0, &data);

    /* zero-extent matrices shouldn't give NULL, but cbind(NULL) should: */
    if (!data.ans_flags && !data.ans_length)
    {
        UNPROTECT(1);
        return R_NilValue;
    }

    mode = NILSXP;
    if (data.ans_flags >= 128)
    {
        if (data.ans_flags & 128)
            mode = LISTSXP;
    }
    else if (data.ans_flags >= 64)
    {
        if (data.ans_flags & 64)
            mode = STRSXP;
    }
    else
    {
        if (data.ans_flags & 1)
            mode = LGLSXP;
        if (data.ans_flags & 8)
            mode = INTSXP;
        if (data.ans_flags & 16)
            mode = REALSXP;
        if (data.ans_flags & 32)
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
        errorcall(call, "cannot create a matrix from these types");
    }

    if (PRIMVAL(op) == 1)
        a = cbind(call, args, mode, rho);
    else
        a = rbind(call, args, mode, rho);
    UNPROTECT(1);
    R_Visible = 1; /* assignment in arguments would set this to zero */
    return a;
}

static void SetRowNames(SEXP dimnames, SEXP x)
{
    if (TYPEOF(dimnames) == VECSXP)
        SET_VECTOR_ELT(dimnames, 0, x);
    else if (TYPEOF(dimnames) == LISTSXP)
        SETCAR(dimnames, x);
}

static void SetColNames(SEXP dimnames, SEXP x)
{
    if (TYPEOF(dimnames) == VECSXP)
        SET_VECTOR_ELT(dimnames, 1, x);
    else if (TYPEOF(dimnames) == LISTSXP)
        SETCADR(dimnames, x);
}

/*
 * Apparently i % 0 could occur here (PR#2541).  But it should not,
 * as zero-length vectors are ignored and
 * zero-length matrices must have zero columns,
 * unless the result has zero rows, hence is of length zero and no
 * copying will be done.
 */
static SEXP cbind(SEXP call, SEXP args, SEXPTYPE mode, SEXP rho)
{
    int i, j, k, idx, n;
    int have_rnames, have_cnames;
    int nnames, mnames;
    int rows, cols, mrows, lenmin = 0;
    int warned;
    SEXP dn, t, u, result, dims, expr;
    struct BindData data;
    data.deparse_level = 1;

    have_rnames = 0;
    have_cnames = 0;
    nnames = 0;
    mnames = 0;
    rows = 0;
    cols = 0;
    mrows = -1;

    /* check if we are in the zero-row case */

    for (t = args; t != R_NilValue; t = CDR(t))
    {
        u = PRVALUE(CAR(t));
        if ((isMatrix(u) ? nrows(u) : length(u)) > 0)
        {
            lenmin = 1;
            break;
        }
    }

    /* check conformability of matrix arguments */

    n = 0;
    for (t = args; t != R_NilValue; t = CDR(t))
    {
        u = PRVALUE(CAR(t));
        dims = getAttrib(u, R_DimSymbol);
        if (length(dims) == 2)
        {
            if (mrows == -1)
                mrows = INTEGER(dims)[0];
            else if (mrows != INTEGER(dims)[0])
                errorcall(call, "number of rows of matrices must match (see arg %d)", n + 1);
            cols += INTEGER(dims)[1];
        }
        else if (length(u) >= lenmin)
        {
            rows = imax2(rows, length(u));
            cols += 1;
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
        u = PRVALUE(CAR(t));
        n++;
        dims = getAttrib(u, R_DimSymbol);
        if (length(dims) == 2)
        {
            dn = getAttrib(u, R_DimNamesSymbol);
            if (length(dn) == 2)
            {
                if (VECTOR_ELT(dn, 1) != R_NilValue)
                    have_cnames = 1;
                if (VECTOR_ELT(dn, 0) != R_NilValue)
                    mnames = mrows;
            }
        }
        else
        {
            k = length(u);
            if (!warned && k > 0 && (k > rows || rows % k))
            {
                warned = 1;
                PROTECT(call = substituteList(call, rho));
                warningcall(call, "number of rows of result\n\tis not a multiple of vector length (arg %d)", n);
                UNPROTECT(1);
            }
            dn = getAttrib(u, R_NamesSymbol);
            if (k >= lenmin &&
                (TAG(t) != R_NilValue || ((data.deparse_level == 1) && isSymbol(substitute(CAR(t), R_NilValue)))))
                have_cnames = 1;
            nnames = imax2(nnames, length(dn));
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
            u = PRVALUE(CAR(t));
            if (isMatrix(u) || length(u) >= lenmin)
            {
                u = coerceVector(u, STRSXP);
                k = LENGTH(u);
                idx = (!isMatrix(u)) ? rows : k;
                for (i = 0; i < idx; i++)
                    SET_STRING_ELT(result, n++, STRING_ELT(u, i % k));
            }
        }
    }
    else if (mode == CPLXSXP)
    {
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            u = PRVALUE(CAR(t));
            if (isMatrix(u) || length(u) >= lenmin)
            {
                u = coerceVector(u, CPLXSXP);
                k = LENGTH(u);
                idx = (!isMatrix(u)) ? rows : k;
                for (i = 0; i < idx; i++)
                    COMPLEX(result)[n++] = COMPLEX(u)[i % k];
            }
        }
    }
    else
    {
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            u = PRVALUE(CAR(t));
            if (isMatrix(u) || length(u) >= lenmin)
            {
                k = LENGTH(u);
                idx = (!isMatrix(u)) ? rows : k;
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
        SEXP nam, tnam, v;
        PROTECT(dn = allocVector(VECSXP, 2));
        if (have_cnames)
            nam = SET_VECTOR_ELT(dn, 1, allocVector(STRSXP, cols));
        else
            nam = R_NilValue; /* -Wall */
        j = 0;
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            u = PRVALUE(CAR(t));
            if (isMatrix(u))
            {

                v = getAttrib(u, R_DimNamesSymbol);
                tnam = GetColNames(v);

                if (have_rnames && GetRowNames(dn) == R_NilValue && GetRowNames(v) != R_NilValue)
                    SetRowNames(dn, duplicate(GetRowNames(v)));

                /* rbind() does this only  if(have_?names) .. : */
                if (tnam != R_NilValue)
                {
                    for (i = 0; i < length(tnam); i++)
                        SET_STRING_ELT(nam, j++, STRING_ELT(tnam, i));
                }
                else if (have_cnames)
                {
                    for (i = 0; i < ncols(u); i++)
                        SET_STRING_ELT(nam, j++, R_BlankString);
                }
            }
            else if (length(u) >= lenmin)
            {

                u = getAttrib(u, R_NamesSymbol);
                tnam = u;

                if (have_rnames && GetRowNames(dn) == R_NilValue && u != R_NilValue && length(u) == rows)
                    SetRowNames(dn, duplicate(u));

                if (TAG(t) != R_NilValue)
                    SET_STRING_ELT(nam, j++, PRINTNAME(TAG(t)));
                else if ((data.deparse_level == 1) && isSymbol(expr = substitute(CAR(t), R_NilValue)))
                    SET_STRING_ELT(nam, j++, PRINTNAME(expr));
                else if (have_cnames)
                    SET_STRING_ELT(nam, j++, R_BlankString);
            }
        }
        setAttrib(result, R_DimNamesSymbol, dn);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return result;
} /* cbind */

static SEXP rbind(SEXP call, SEXP args, SEXPTYPE mode, SEXP rho)
{
    int i, j, k, n;
    int have_rnames, have_cnames;
    int nnames, mnames;
    int rows, cols, mcols, mrows, have_mcols, lenmin = 0;
    int warned;
    SEXP dn, t, u, result, dims, expr;
    struct BindData data;
    data.deparse_level = 1;

    have_rnames = 0;
    have_cnames = 0;
    nnames = 0;
    mnames = 0;
    rows = 0;
    cols = 0;
    mcols = 0;
    have_mcols = 0;

    /* check if we are in the zero-cols case */

    for (t = args; t != R_NilValue; t = CDR(t))
    {
        u = PRVALUE(CAR(t));
        if ((isMatrix(u) ? ncols(u) : length(u)) > 0)
        {
            lenmin = 1;
            break;
        }
    }

    /* check conformability of matrix arguments */

    n = 0;
    for (t = args; t != R_NilValue; t = CDR(t))
    {
        u = PRVALUE(CAR(t));
        dims = getAttrib(u, R_DimSymbol);
        if (length(dims) == 2)
        {
            if (!have_mcols)
            {
                mcols = INTEGER(dims)[1];
                have_mcols = 1;
            }
            else if (mcols != INTEGER(dims)[1])
                errorcall(call, "number of columns of matrices must match (see arg %d)", n + 1);
            rows += INTEGER(dims)[0];
        }
        else if (length(u) >= lenmin)
        {
            cols = imax2(cols, length(u));
            rows += 1;
        }
        n++;
    }
    if (have_mcols)
        cols = mcols;

    /* Check conformability of vector arguments. -- Look for dimnames. */

    n = 0;
    warned = 0;
    for (t = args; t != R_NilValue; t = CDR(t))
    {
        u = PRVALUE(CAR(t));
        n++;
        dims = getAttrib(u, R_DimSymbol);
        if (length(dims) == 2)
        {
            dn = getAttrib(u, R_DimNamesSymbol);
            if (length(dn) == 2)
            {
                if (VECTOR_ELT(dn, 0) != R_NilValue)
                    have_rnames = 1;
                if (VECTOR_ELT(dn, 1) != R_NilValue)
                    mnames = mcols;
            }
        }
        else
        {
            k = length(u);
            if (!warned && k > 0 && (k > cols || cols % k))
            {
                warned = 1;
                PROTECT(call = substituteList(call, rho));
                warningcall(call, "number of columns of result\n\tnot a multiple of vector length (arg %d)", n);
                UNPROTECT(1);
            }
            dn = getAttrib(u, R_NamesSymbol);
            if (k >= lenmin &&
                (TAG(t) != R_NilValue || ((data.deparse_level == 1) && isSymbol(substitute(CAR(t), R_NilValue)))))
                have_rnames = 1;
            nnames = imax2(nnames, length(dn));
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
            u = PRVALUE(CAR(t));
            if (isMatrix(u) || length(u) >= lenmin)
            {
                u = coerceVector(u, STRSXP);
                k = LENGTH(u);
                mrows = (isMatrix(u)) ? nrows(u) : (k > 0);
                for (i = 0; i < mrows; i++)
                    for (j = 0; j < cols; j++)
                        SET_STRING_ELT(result, i + n + (j * rows), STRING_ELT(u, (i + j * mrows) % k));
                n += mrows;
            }
        }
    }
    else if (mode == CPLXSXP)
    {
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            u = PRVALUE(CAR(t));
            if (isMatrix(u) || length(u) >= lenmin)
            {
                u = coerceVector(u, CPLXSXP);
                k = LENGTH(u);
                mrows = (isMatrix(u)) ? nrows(u) : (k > 0);
                for (i = 0; i < mrows; i++)
                    for (j = 0; j < cols; j++)
                        COMPLEX(result)[i + n + (j * rows)] = COMPLEX(u)[(i + j * mrows) % k];
                n += mrows;
            }
        }
    }
    else
    {
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            u = PRVALUE(CAR(t));
            if (isMatrix(u) || length(u) >= lenmin)
            {
                k = LENGTH(u);
                mrows = (isMatrix(u)) ? nrows(u) : (k > 0);
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
        SEXP nam, tnam, v;
        PROTECT(dn = allocVector(VECSXP, 2));
        if (have_rnames)
            nam = SET_VECTOR_ELT(dn, 0, allocVector(STRSXP, rows));
        else
            nam = R_NilValue; /* -Wall */
        j = 0;
        for (t = args; t != R_NilValue; t = CDR(t))
        {
            u = PRVALUE(CAR(t));
            if (isMatrix(u))
            {

                v = getAttrib(u, R_DimNamesSymbol);
                tnam = GetRowNames(v);

                if (have_cnames && GetColNames(dn) == R_NilValue && GetColNames(v) != R_NilValue)
                    SetColNames(dn, duplicate(GetColNames(v)));

                /* cbind() doesn't test have_?names BEFORE tnam!=Nil..:*/
                if (have_rnames)
                {
                    if (tnam != R_NilValue)
                    {
                        for (i = 0; i < length(tnam); i++)
                            SET_STRING_ELT(nam, j++, STRING_ELT(tnam, i));
                    }
                    else
                    {
                        for (i = 0; i < nrows(u); i++)
                            SET_STRING_ELT(nam, j++, R_BlankString);
                    }
                }
            }
            else if (length(u) >= lenmin)
            {

                u = getAttrib(u, R_NamesSymbol);
                tnam = u;

                if (have_cnames && GetColNames(dn) == R_NilValue && u != R_NilValue && length(u) == cols)
                    SetColNames(dn, duplicate(u));

                if (TAG(t) != R_NilValue)
                    SET_STRING_ELT(nam, j++, PRINTNAME(TAG(t)));
                else if ((data.deparse_level == 1) && isSymbol(expr = substitute(CAR(t), R_NilValue)))
                    SET_STRING_ELT(nam, j++, PRINTNAME(expr));
                else if (have_rnames)
                    SET_STRING_ELT(nam, j++, R_BlankString);
            }
        }
        setAttrib(result, R_DimNamesSymbol, dn);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return result;
} /* rbind */
