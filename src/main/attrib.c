/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997--1998  Robert Gentleman, Ross Ihaka and the R Core team.
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
#include "Mathlib.h"

static void checkNames(SEXP, SEXP);
static SEXP installAttrib(SEXP, SEXP, SEXP);
static SEXP removeAttrib(SEXP, SEXP);

SEXP comment(SEXP);
SEXP commentgets(SEXP, SEXP);

static SEXP stripAttrib(SEXP tag, SEXP lst)
{
    if (lst == R_NilValue)
        return lst;
    if (tag == TAG(lst))
        return stripAttrib(tag, CDR(lst));
    CDR(lst) = stripAttrib(tag, CDR(lst));
    return lst;
}

SEXP getAttrib(SEXP vec, SEXP name)
{
    SEXP s, blank;
    int len, i, any;

    if (isString(name))
        name = install(CHAR(STRING(name)[0]));

    if (name == R_NamesSymbol)
    {
        if (isVector(vec) || isList(vec) || isLanguage(vec))
        {
            s = getAttrib(vec, R_DimSymbol);
            if (TYPEOF(s) == INTSXP && length(s) == 1)
            {
                s = getAttrib(vec, R_DimNamesSymbol);
                if (!isNull(s))
                    return CAR(s);
            }
        }
        if (isList(vec) || isLanguage(vec))
        {
            len = length(vec);
            PROTECT(s = allocVector(STRSXP, len));
            blank = mkChar("");
            i = 0;
            any = 0;
            for (; vec != R_NilValue; vec = CDR(vec), i++)
            {
                if (TAG(vec) == R_NilValue)
                    STRING(s)[i] = blank;
                else if (isSymbol(TAG(vec)))
                {
                    any = 1;
                    STRING(s)[i] = PRINTNAME(TAG(vec));
                }
                else
                    error("getAttrib: invalid type for TAG\n");
            }
            UNPROTECT(1);
            if (any)
                return (s);
            return R_NilValue;
        }
    }
    for (s = ATTRIB(vec); s != R_NilValue; s = CDR(s))
        if (TAG(s) == name)
        {
            NAMED(CAR(s)) = NAMED(vec);
            return CAR(s);
        }
    return R_NilValue;
}

SEXP setAttrib(SEXP vec, SEXP name, SEXP val)
{
    if (isString(name))
        name = install(CHAR(STRING(name)[0]));
    if (val == R_NilValue)
        return removeAttrib(vec, name);

    if (vec == R_NilValue)
        error("attempt to set an attribute on NULL\n");

    PROTECT(vec);
    PROTECT(name);
    val = duplicate(val);
    UNPROTECT(2);

    if (name == R_NamesSymbol)
        return namesgets(vec, val);
    else if (name == R_DimSymbol)
        return dimgets(vec, val);
    else if (name == R_DimNamesSymbol)
        return dimnamesgets(vec, val);
    else if (name == R_ClassSymbol)
        return classgets(vec, val);
    else if (name == R_TspSymbol)
        return tspgets(vec, val);
    else if (name == R_CommentSymbol)
        return commentgets(vec, val);
    else
        return installAttrib(vec, name, val);
}

/* This is called in the case of binary operations to copy */
/* most attributes from (one of) the input arguments to */
/* the output.	Note that the Dim and Names attributes */
/* should have been assigned elsewhere. */

void copyMostAttrib(SEXP inp, SEXP ans)
{
    SEXP s;
    PROTECT(ans);
    PROTECT(inp);
    for (s = ATTRIB(inp); s != R_NilValue; s = CDR(s))
    {
        if ((TAG(s) != R_NamesSymbol) && (TAG(s) != R_DimSymbol) && (TAG(s) != R_DimNamesSymbol))
        {
            installAttrib(ans, TAG(s), CAR(s));
        }
    }
    UNPROTECT(2);
}

static SEXP installAttrib(SEXP vec, SEXP name, SEXP val)
{
    SEXP s, t;

    PROTECT(vec);
    PROTECT(name);
    PROTECT(val);

    for (s = ATTRIB(vec); s != R_NilValue; s = CDR(s))
    {
        if (TAG(s) == name)
        {
            CAR(s) = val;
            UNPROTECT(3);
            return val;
        }
    }
    s = allocList(1);
    CAR(s) = val;
    TAG(s) = name;
    if (ATTRIB(vec) == R_NilValue)
        ATTRIB(vec) = s;
    else
    {
        t = nthcdr(ATTRIB(vec), length(ATTRIB(vec)) - 1);
        SETCDR(t, s);
    }
    UNPROTECT(3);
    return val;
}

static SEXP removeAttrib(SEXP vec, SEXP name)
{
    SEXP t;

    if (name == R_NamesSymbol && isList(vec))
    {
        for (t = vec; t != R_NilValue; t = CDR(t))
            TAG(t) = R_NilValue;
        return R_NilValue;
    }
    else
    {
        if (name == R_DimSymbol)
            ATTRIB(vec) = stripAttrib(R_DimNamesSymbol, ATTRIB(vec));
        ATTRIB(vec) = stripAttrib(name, ATTRIB(vec));
        if (name == R_ClassSymbol)
            OBJECT(vec) = 0;
    }
    return R_NilValue;
}

static void checkNames(SEXP x, SEXP s)
{
    if (isVector(x) || isList(x) || isLanguage(x))
    {
        if (!isVector(s) && !isList(s))
            error("invalid type for names: must be vector\n");
        if (length(x) != length(s))
            error("names attribute must be the same length as the vector\n");
    }
    else
        error("names applied to non-vector\n");
}

/* Time Series Parameters */

static void badtsp()
{
    error("invalid time series parameters specified\n");
}

SEXP tspgets(SEXP vec, SEXP val)
{
    double start, end, frequency;
    int n;

    if (!isNumeric(val) || length(val) != 3)
        error("tsp attribute must be numeric of length three\n");

    if (isReal(val))
    {
        start = REAL(val)[0];
        end = REAL(val)[1];
        frequency = REAL(val)[2];
    }
    else
    {
        start = (INTEGER(val)[0] == NA_INTEGER) ? NA_REAL : INTEGER(val)[0];
        end = (INTEGER(val)[1] == NA_INTEGER) ? NA_REAL : INTEGER(val)[1];
        frequency = (INTEGER(val)[2] == NA_INTEGER) ? NA_REAL : INTEGER(val)[2];
    }
    if (frequency <= 0)
        badtsp();
    n = nrows(vec);
    if (fabs(end - start - (n - 1) / frequency) > 1.e-5)
        badtsp();

    PROTECT(vec);
    val = allocVector(REALSXP, 3);
    PROTECT(val);
    REAL(val)[0] = start;
    REAL(val)[1] = end;
    REAL(val)[2] = frequency;
    installAttrib(vec, R_TspSymbol, val);
    UNPROTECT(2);
    return vec;
}

SEXP commentgets(SEXP vec, SEXP comment)
{
    if (isNull(comment) || isString(comment))
    {
        if (length(comment) <= 0)
        {
            ATTRIB(vec) = stripAttrib(R_CommentSymbol, vec);
        }
        else
        {
            installAttrib(vec, R_CommentSymbol, comment);
        }
        return R_NilValue;
    }
    error("attempt to set invalid comment attribute\n");
    return R_NilValue; /*- just for -Wall */
}

SEXP do_commentgets(SEXP call, SEXP op, SEXP args, SEXP env)
{
    checkArity(op, args);
    if (NAMED(CAR(args)) == 2)
        CAR(args) = duplicate(CAR(args));
    if (length(CADR(args)) == 0)
        CADR(args) = R_NilValue;
    setAttrib(CAR(args), R_CommentSymbol, CADR(args));
    return CAR(args);
}

SEXP do_comment(SEXP call, SEXP op, SEXP args, SEXP env)
{
    checkArity(op, args);
    return getAttrib(CAR(args), R_CommentSymbol);
}

SEXP classgets(SEXP vec, SEXP class)
{
    if (isNull(class) || isString(class))
    {
        if (length(class) <= 0)
        {
            ATTRIB(vec) = stripAttrib(R_ClassSymbol, vec);
            OBJECT(vec) = 0;
        }
        else
        {
            if (streql(CHAR(STRING(class)[0]), "data.frame") && !isList(vec))
                error("attempt to make non-list a data frame\n");
            installAttrib(vec, R_ClassSymbol, class);
            OBJECT(vec) = 1;
        }
        return R_NilValue;
    }
    error("attempt to set invalid class attribute\n");
    return R_NilValue; /*- just for -Wall */
}

SEXP do_classgets(SEXP call, SEXP op, SEXP args, SEXP env)
{
    checkArity(op, args);
    if (NAMED(CAR(args)) == 2)
        CAR(args) = duplicate(CAR(args));
    if (length(CADR(args)) == 0)
        CADR(args) = R_NilValue;
    setAttrib(CAR(args), R_ClassSymbol, CADR(args));
    return CAR(args);
}

SEXP do_class(SEXP call, SEXP op, SEXP args, SEXP env)
{
    checkArity(op, args);
    return getAttrib(CAR(args), R_ClassSymbol);
}

SEXP do_namesgets(SEXP call, SEXP op, SEXP args, SEXP env)
{
    checkArity(op, args);
    if (NAMED(CAR(args)) == 2)
        CAR(args) = duplicate(CAR(args));
    setAttrib(CAR(args), R_NamesSymbol, CADR(args));
    return CAR(args);
}

SEXP namesgets(SEXP vec, SEXP val)
{
    int i;
    SEXP s, rval;

    PROTECT(vec);
    PROTECT(val);

    /* Ensure that the labels are indeed */
    /* a vector of character strings */

    if (isList(val))
    {
        if (!isVectorizable(val))
            error("incompatible names argument\n");
        else
        {
            rval = allocVector(STRSXP, length(vec));
            PROTECT(rval);
            for (i = 0; i < length(vec); i++)
            {
                s = coerceVector(CAR(val), STRSXP);
                STRING(rval)[i] = STRING(s)[0];
            }
            UNPROTECT(1);
            val = rval;
        }
    }
    else
        val = coerceVector(val, STRSXP);
    UNPROTECT(1);
    PROTECT(val);

    /* Check that the lengths and types are compatible */

    checkNames(vec, val);

    /* Special treatment for one dimensional arrays */

    if (isVector(vec) || isList(vec) || isLanguage(vec))
    {
        s = getAttrib(vec, R_DimSymbol);
        if (TYPEOF(s) == INTSXP && length(s) == 1)
        {
            PROTECT(val = CONS(val, R_NilValue));
            setAttrib(vec, R_DimNamesSymbol, val);
            UNPROTECT(3);
            return vec;
        }
    }

    /* Cons-cell based objects */

    if (isList(vec) || isLanguage(vec))
    {
        i = 0;
        for (s = vec; s != R_NilValue; s = CDR(s), i++)
            if (STRING(val)[i] != R_NilValue && STRING(val)[i] != R_NaString && *CHAR(STRING(val)[i]) != 0)
                TAG(s) = install(CHAR(STRING(val)[i]));
            else
                TAG(s) = R_NilValue;
    }
    else if (isVector(vec))
        installAttrib(vec, R_NamesSymbol, val);
    else
        error("invalid type to set names attribute\n");
    UNPROTECT(2);
    return vec;
}

SEXP do_names(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP s;
    checkArity(op, args);
    s = CAR(args);
    if (isVector(s) || isList(s) || isLanguage(s))
        return getAttrib(s, R_NamesSymbol);
    return R_NilValue;
}

SEXP do_dimnamesgets(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans;
    if (DispatchOrEval(call, op, args, env, &ans, 0))
        return (ans);
    PROTECT(args = ans);
    checkArity(op, args);
    if (NAMED(CAR(args)) > 2)
        CAR(args) = duplicate(CAR(args));
    setAttrib(CAR(args), R_DimNamesSymbol, CADR(args));
    UNPROTECT(1);
    return CAR(args);
}

SEXP dimnamesgets(SEXP vec, SEXP val)
{
    SEXP dims, top;
    int i, k;

    PROTECT(vec);
    PROTECT(val);

    if (!isArray(vec) && !isList(vec))
        error("dimnames applied to non-array\n");
    if (!isList(val))
        error("invalid type for dimnames: must be a list\n");
    dims = getAttrib(vec, R_DimSymbol);
    if ((k = LENGTH(dims)) != length(val))
        error("dimnames: number of dimensions must equal number of names\n");
    top = val;
    for (i = 0; i < k; i++)
    {
        if (CAR(val) != R_NilValue)
        {
            if (!isVector(CAR(val)))
                error("invalid type for dim name must be a vector\n");
            if (INTEGER(dims)[i] != LENGTH(CAR(val)) && LENGTH(CAR(val)) != 0)
                error("length of namelist must equal dims\n");
            if (LENGTH(CAR(val)) == 0)
            {
                CAR(val) = R_NilValue;
            }
            else if (!isString(CAR(val)))
            {
                CAR(val) = coerceVector(CAR(val), STRSXP);
            }
        }
        val = CDR(val);
    }
    installAttrib(vec, R_DimNamesSymbol, top);
    if (isList(vec) && k == 1)
    {
        top = CAR(top);
        i = 0;
        for (val = vec; !isNull(val); val = CDR(val))
            TAG(val) = install(CHAR(STRING(top)[i++]));
    }
    UNPROTECT(2);
    return (vec);
}

SEXP do_dimnames(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans;
    if (DispatchOrEval(call, op, args, env, &ans, 0))
        return (ans);
    PROTECT(args = ans);
    checkArity(op, args);
    ans = getAttrib(CAR(args), R_DimNamesSymbol);
    UNPROTECT(1);
    return ans;
}

SEXP do_dim(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans;
    if (DispatchOrEval(call, op, args, env, &ans, 0))
        return (ans);
    PROTECT(args = ans);
    checkArity(op, args);
    ans = getAttrib(CAR(args), R_DimSymbol);
    UNPROTECT(1);
    return ans;
}

SEXP do_dimgets(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans;
    if (DispatchOrEval(call, op, args, env, &ans, 0))
        return (ans);
    PROTECT(args = ans);
    checkArity(op, args);
    if (NAMED(CAR(args)) > 1)
        CAR(args) = duplicate(CAR(args));
    setAttrib(CAR(args), R_DimSymbol, CADR(args));
    setAttrib(CAR(args), R_NamesSymbol, R_NilValue);
    UNPROTECT(1);
    return CAR(args);
}

SEXP dimgets(SEXP vec, SEXP val)
{
    int len, ndim, i, total;

    PROTECT(vec);
    PROTECT(val);

    if (!isVector(vec) && !isList(vec))
        error("dim<- : invalid first argument\n");

    if (!isVector(val) && !isList(val))
        error("dim<- : invalid second argument\n");
    val = coerceVector(val, INTSXP);
    UNPROTECT(1);
    PROTECT(val);

    len = length(vec);
    ndim = length(val);
    if (ndim == 0)
        error("dim: Invalid dimension vector\n");
    total = 1;
    for (i = 0; i < ndim; i++)
        total *= INTEGER(val)[i];
    if (total != len)
        error("dim<- length of dims do not match the length of object\n");
    removeAttrib(vec, R_DimNamesSymbol);
    installAttrib(vec, R_DimSymbol, val);
    UNPROTECT(2);
    return vec;
}

#ifdef NEWLIST
SEXP do_attributes(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP attrs, blank, names, namesattr, value;
    int nvalues;
    namesattr = R_NilValue;
    attrs = ATTRIB(CAR(args));
    nvalues = length(attrs);
    if (isList(CAR(args)))
    {
        namesattr = getAttrib(CAR(args), R_NamesSymbol);
        if (namesattr != R_NilValue)
            nvalues++;
    }
    PROTECT(value = allocVector(VECSXP, nvalues));
    PROTECT(names = allocVector(STRSXP, nvalues));
    PROTECT(blank = mkChar(""));
    nvalues = 0;
    if (namesattr != R_NilValue)
    {
        VECTOR(value)[nvalues] = namesattr;
        STRING(names)[nvalues] = PRINTNAME(R_NamesSymbol);
        nvalues++;
    }
    while (attrs != R_NilValue)
    {
        VECTOR(value)[nvalues] = CAR(attrs);
        if (TAG(attrs) == R_NilValue)
            STRING(names)[nvalues] = blank;
        else
            STRING(names)[nvalues] = PRINTNAME(TAG(attrs));
        attrs = CDR(attrs);
        nvalues++;
    }
    setAttrib(value, R_NamesSymbol, names);
    NAMED(value) = NAMED(CAR(args));
    UNPROTECT(3);
    return value;
}
#else
SEXP do_attributes(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP s;

    s = R_NilValue;
    if (isList(CAR(args)))
        s = getAttrib(CAR(args), R_NamesSymbol);
    PROTECT(s);
    if (s != R_NilValue)
    {
        s = CONS(s, ATTRIB(CAR(args)));
        TAG(s) = R_NamesSymbol;
    }
    else
        s = ATTRIB(CAR(args));
    UNPROTECT(1);
    NAMED(s) = NAMED(CAR(args));
    return s;
}
#endif

/* NOTE: The following code ensures that when an attribute list */
/* is attached to an object, that the "dim" attibute is always */
/* brought to the front of the list.  This ensures that when both */
/* "dim" and "dimnames" are set that the "dim" is attached first. */

#ifdef NEWLIST
SEXP do_attributesgets(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP object, attrs, names;
    int i, nattrs;

    /* If there are multiple references to the */
    /* object being mutated, we must duplicate. */

    if (NAMED(CAR(args)) == 2)
        CAR(args) = duplicate(CAR(args));

    /* Extract the arguments from the argument list */

    object = CAR(args);
    if (object == R_NilValue)
    {
        warning("attempt to set attributes on NULL\n");
        return (R_NilValue);
    }
    attrs = CADR(args);
    if (!isNewList(attrs))
        errorcall(call, "attributes must be in a list\n");
    names = getAttrib(attrs, R_NamesSymbol);
    if (names == R_NilValue)
        errorcall(call, "attributes must be named\n");

    /* Empty the existing attribute list */
    /* FIXME: the code below treats pair-based */
    /* structures in a special way.  This can */
    /* probably be dropped down the road. */
    /* Users should never encounter pair-based lists. */

    if (isList(object))
        setAttrib(object, R_NamesSymbol, R_NilValue);
    ATTRIB(object) = R_NilValue;
    OBJECT(object) = 0;

    /* We do two passes through the attributes */
    /* the first finding and transferring "dims" */
    /* and the second transferring the rest. */
    /* This is to ensure that "dim" occurs */
    /* in the attribute list before "dimnames". */

    nattrs = length(attrs);
    for (i = 0; i < nattrs; i++)
    {
        if (STRING(names)[i] == R_NilValue && CHAR(STRING(names)[i])[0] == '\0')
        {
            error("all attributes must have names\n");
        }
        if (!strcmp(CHAR(STRING(names)[i]), "dim"))
            setAttrib(object, R_DimSymbol, VECTOR(attrs)[i]);
    }
    for (i = 0; i < nattrs; i++)
    {
        if (strcmp(CHAR(STRING(attrs)[i]), "dim"))
            setAttrib(object, install(CHAR(STRING(names)[i])), VECTOR(attrs)[i]);
    }
    return object;
}
#else
static SEXP dimptr;

static SEXP TrimDim(SEXP l)
{
    if (l != R_NilValue)
    {
        if (TAG(l) == R_DimSymbol)
        {
            dimptr = l;
            return CDR(l);
        }
        else
        {
            CDR(l) = TrimDim(CDR(l));
            return l;
        }
    }
    return R_NilValue;
}

SEXP do_attributesgets(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP s, t;

    if (CAR(args) == R_NilValue)
    {
        warning("attempt to set attributes on NULL\n");
        return (R_NilValue);
    }
    if (NAMED(CAR(args)) == 2)
        CAR(args) = duplicate(CAR(args));
    s = CAR(args);
    t = CADR(args);
    if (isList(s))
        setAttrib(s, R_NamesSymbol, R_NilValue);
    ATTRIB(s) = R_NilValue;

    OBJECT(s) = 0;

    if (!isList(t))
        errorcall(call, "attributes must be in a list\n");

    /* Ghastly hack to ensure that "dim" */
    /* is always the first attribute */

    dimptr = R_NilValue;
    t = TrimDim(t);
    if (dimptr != R_NilValue)
    {
        CDR(dimptr) = t;
        t = dimptr;
    }

    for (; t != R_NilValue; t = CDR(t))
    {
        if (TAG(t) == R_NilValue)
            error("all attributes must have names\n");
        setAttrib(s, TAG(t), CAR(t));
    }
    return s;
}
#endif

SEXP do_attr(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP s, t;

    s = CAR(args);
    t = CADR(args);

    if (!isString(t) || length(t) == 0)
        error("attribute name must be of mode character\n");

    return getAttrib(s, install(CHAR(STRING(t)[0])));
}

SEXP do_attrgets(SEXP call, SEXP op, SEXP args, SEXP env)
{
    /*  attr(obj, "<name>")  <-  value  */
    SEXP obj, name, value;

    obj = eval(CAR(args), env);
    if (NAMED(obj) == 2)
        PROTECT(duplicate(obj));
    else
        PROTECT(obj);

    PROTECT(name = eval(CADR(args), env));
    if (!isString(name))
        error("attr<- : name must be of mode character\n");

    /* no eval(.), RHS is already evaluated: */
    PROTECT(value = CADDR(args));
    setAttrib(obj, name, value);
    UNPROTECT(3);
    return obj;
}
