/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998-2000   The R Development Core Team.
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
#include "Print.h"

/* The global var. R_Expressions is in Defn.h */
#define R_MIN_EXPRESSIONS_OPT 25
#define R_MAX_EXPRESSIONS_OPT 100000

/* Interface to the (polymorphous!)  options(...)  command.
 *
 * We have two kind of options:
 *   1) those used exclusively from R code,
 *	typically initialized in Rprofile.

 *	Their names need not appear here, but may, when we want
 *	to make sure that they are assigned `valid' values only.
 *
 *   2) Those used (and sometimes set) from C code;
 *	Either accessing and/or setting a global C variable,
 *	or just accessed by e.g.  GetOption(install("pager"), ..)
 *
 * A (complete?!) list of these (2):
 *
 *	"prompt"
 *	"continue"
 *	"editor"
 *	"expressions"
 *	"width"
 *	"digits"
 *	"contrasts"
 *	"echo"
 *	"verbose"
 *	"keep.source"
 *	"keep.source.pkgs"

 *	"de.cellwidth"		../unix/X11/ & ../gnuwin32/dataentry.c
 *	"device"
 *	"pager"
 *	"paper.size"		./devPS.c

 *	"timeout"		./connections.c

 *	"check.bounds"
 *	"error"
 *	"error.messages"
 *	"show.error.messages"
 *	"warn"
 *	"warning.length"
 *	"warning.expression"

 *
 * S additionally/instead has (and one might think about some)
 * "free",	"keep"
 * "length",	"memory"
 * "object.size"
 * "reference", "show"
 * "scrap"
 */

static SEXP Options(void)
{
    return install(".Options");
}

static SEXP FindTaggedItem(SEXP lst, SEXP tag)
{
    for (; lst != R_NilValue; lst = CDR(lst))
    {
        if (TAG(lst) == tag)
            return lst;
    }
    return R_NilValue;
}

static SEXP makeErrorCall(SEXP fun)
{
    SEXP call;
    PROTECT(call = allocList(1));
    SET_TYPEOF(call, LANGSXP);
    SETCAR(call, fun);
    UNPROTECT(1);
    return call;
}

SEXP GetOption(SEXP tag, SEXP rho)
{
    SEXP opt = findVar(Options(), R_NilValue);
    if (!isList(opt))
        error("corrupted options list");
    opt = FindTaggedItem(opt, tag);
    return CAR(opt);
}

int GetOptionWidth(SEXP rho)
{
    int w;
    w = asInteger(GetOption(install("width"), rho));
    if (w < R_MIN_WIDTH_OPT || w > R_MAX_WIDTH_OPT)
    {
        warning("invalid printing width, used 80");
        return 80;
    }
    return w;
}

int GetOptionDigits(SEXP rho)
{
    int d;
    d = asInteger(GetOption(install("digits"), rho));
    if (d < R_MIN_DIGITS_OPT || d > R_MAX_DIGITS_OPT)
    {
        warning("invalid printing digits, used 7");
        return 7;
    }
    return d;
}

/* Change the value of an option or add a new option or, */
/* if called with value R_NilValue, remove that option. */

static SEXP SetOption(SEXP tag, SEXP value)
{
    SEXP opt, old, t;
    t = opt = SYMVALUE(Options());
    if (!isList(opt))
        error("corrupted options list");
    opt = FindTaggedItem(opt, tag);

    /* The option is being removed. */
    if (value == R_NilValue)
    {
        for (; t != R_NilValue; t = CDR(t))
            if (TAG(CDR(t)) == tag)
            {
                old = CAR(t);
                SETCDR(t, CDDR(t));
                return old;
            }
        return R_NilValue;
    }
    /* If the option is new, a new slot */
    /* is added to the end of .Options */
    if (opt == R_NilValue)
    {
        while (CDR(t) != R_NilValue)
            t = CDR(t);
        PROTECT(value);
        SETCDR(t, allocList(1));
        UNPROTECT(1);
        opt = CDR(t);
        SET_TAG(opt, tag);
    }
    old = CAR(opt);
    SETCAR(opt, value);
    return old;
}

/* Set the width of lines for printing i.e. like options(width=...) */
/* Returns the previous value for the options. */

int R_SetOptionWidth(int w)
{
    SEXP t, v;
    if (w < R_MIN_WIDTH_OPT)
        w = R_MIN_WIDTH_OPT;
    if (w > R_MAX_WIDTH_OPT)
        w = R_MAX_WIDTH_OPT;
    PROTECT(t = install("width"));
    PROTECT(v = ScalarInteger(w));
    v = SetOption(t, v);
    UNPROTECT(2);
    return INTEGER(v)[0];
}

int R_SetOptionWarn(int w)
{
    SEXP t, v;

    t = install("warn");
    PROTECT(v = ScalarInteger(w));
    v = SetOption(t, v);
    UNPROTECT(1);
    return INTEGER(v)[0];
}

/* Note that options are stored as a dotted pair list */
/* This is barely historical, but is also useful. */

void InitOptions(void)
{
    SEXP t, val, v;
    char *p;

    PROTECT(v = val = allocList(13));

    SET_TAG(v, install("prompt"));
    SETCAR(v, mkString("> "));
    v = CDR(v);

    SET_TAG(v, install("continue"));
    SETCAR(v, mkString("+ "));
    v = CDR(v);

    SET_TAG(v, install("editor"));
    SETCAR(v, mkString("vi"));
    v = CDR(v);

    SET_TAG(v, install("expressions"));
    SETCAR(v, ScalarInteger(R_Expressions));
    v = CDR(v);

    SET_TAG(v, install("width"));
    SETCAR(v, ScalarInteger(80));
    v = CDR(v);

    SET_TAG(v, install("digits"));
    SETCAR(v, ScalarInteger(7));
    v = CDR(v);

    SET_TAG(v, install("contrasts"));
    SETCAR(v, allocVector(STRSXP, 2));
    SET_STRING_ELT(CAR(v), 0, mkChar("contr.treatment"));
    SET_STRING_ELT(CAR(v), 1, mkChar("contr.poly"));
    PROTECT(t = allocVector(STRSXP, 2));
    SET_STRING_ELT(t, 0, mkChar("unordered"));
    SET_STRING_ELT(t, 1, mkChar("ordered"));
    namesgets(CAR(v), t);
    v = CDR(v);

    SET_TAG(v, install("echo"));
    SETCAR(v, allocVector(LGLSXP, 1));
    LOGICAL(CAR(v))[0] = !R_Slave;
    v = CDR(v);

    SET_TAG(v, install("verbose"));
    SETCAR(v, allocVector(LGLSXP, 1));
    LOGICAL(CAR(v))[0] = R_Verbose;
    v = CDR(v);

    SET_TAG(v, install("check.bounds"));
    SETCAR(v, allocVector(LGLSXP, 1));
    LOGICAL(CAR(v))[0] = 0; /* no checking */
    v = CDR(v);

    p = getenv("R_KEEP_PKG_SOURCE");
    R_KeepSource = (p && (strcmp(p, "yes") == 0)) ? 1 : 0;

    SET_TAG(v, install("keep.source"));
    SETCAR(v, allocVector(LGLSXP, 1));
    LOGICAL(CAR(v))[0] = R_KeepSource;
    v = CDR(v);

    SET_TAG(v, install("keep.source.pkgs"));
    SETCAR(v, allocVector(LGLSXP, 1));
    LOGICAL(CAR(v))[0] = R_KeepSource;
    v = CDR(v);

    SET_TAG(v, install("error.messages"));
    SETCAR(v, allocVector(LGLSXP, 1));
    LOGICAL(CAR(v))[0] = 1;

    SET_TAG(v, install("warnings.length"));
    SETCAR(v, allocVector(INTSXP, 1));
    INTEGER(CAR(v))[0] = 1000;

    SET_SYMVALUE(install(".Options"), val);
    UNPROTECT(2);
}

SEXP do_options(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP argi = R_NilValue, argnames = R_NilValue, namei = R_NilValue, names, options, s, tag,
         value; /* = R_Nil..: -Wall */
    int i, k, n;

    /* Locate the options values in the symbol table.
       This will need to change if options are to live in the session
       frame.
       */

    options = SYMVALUE(Options());

    if (args == R_NilValue)
    {
        /* This is the zero argument case.
           We alloc up a real list and write the system values into it.
        */
        n = length(options);
        PROTECT(value = allocVector(VECSXP, n));
        PROTECT(names = allocVector(STRSXP, n));
        i = 0;
        while (options != R_NilValue)
        {
            SET_VECTOR_ELT(names, i, PRINTNAME(TAG(options)));
            SET_VECTOR_ELT(value, i, duplicate(CAR(options)));
            options = CDR(options);
            i++;
        }
        setAttrib(value, R_NamesSymbol, names);
        UNPROTECT(2);
        return value;
    }

    /* The arguments to "options" can either be a sequence of
       name = value form, or can be a single list.
       This means that we must code so that both forms will work.
       [ Vomits quietly onto shoes ... ]
       */

    n = length(args);
    if (n == 1 && (isPairList(CAR(args)) || isVectorList(CAR(args))) && TAG(args) == R_NilValue)
    {
        args = CAR(args);
        n = length(args);
    }
    PROTECT(value = allocVector(VECSXP, n));
    PROTECT(names = allocVector(STRSXP, n));

    switch (TYPEOF(args))
    {
    case NILSXP:
    case LISTSXP:
        argnames = R_NilValue;
        break;
    case VECSXP:
        argnames = getAttrib(args, R_NamesSymbol);
        if (LENGTH(argnames) != n)
            errorcall(call, "list argument has no names or invalid ones");
        break;
    }

    R_Visible = 0;
    for (i = 0; i < n; i++)
    { /* i-th argument */

        switch (TYPEOF(args))
        {
        case LISTSXP:
            argi = CAR(args);
            namei = EnsureString(TAG(args));
            args = CDR(args);
            break;
        case VECSXP:
            argi = VECTOR_ELT(args, i);
            namei = EnsureString(STRING_ELT(argnames, i));
            break;
        }

        if (*CHAR(namei))
        { /* name = value  ---> assignment */
            tag = install(CHAR(namei));
            if (streql(CHAR(namei), "width"))
            {
                k = asInteger(argi);
                if (k < R_MIN_WIDTH_OPT || k > R_MAX_WIDTH_OPT)
                    errorcall(call, "invalid width parameter, allowed %d...%d", R_MIN_WIDTH_OPT, R_MAX_WIDTH_OPT);
                SET_VECTOR_ELT(value, i, SetOption(tag, ScalarInteger(k)));
            }
            else if (streql(CHAR(namei), "digits"))
            {
                k = asInteger(argi);
                if (k < R_MIN_DIGITS_OPT || k > R_MAX_DIGITS_OPT)
                    errorcall(call, "invalid digits parameter, allowed %d...%d", R_MIN_DIGITS_OPT, R_MAX_DIGITS_OPT);
                SET_VECTOR_ELT(value, i, SetOption(tag, ScalarInteger(k)));
            }
            else if (streql(CHAR(namei), "expressions"))
            {
                k = asInteger(argi);
                if (k < R_MIN_EXPRESSIONS_OPT || k > R_MAX_EXPRESSIONS_OPT)
                    errorcall(call, "expressions parameter invalid, allowed %d...%d", R_MIN_EXPRESSIONS_OPT,
                              R_MAX_EXPRESSIONS_OPT);
                R_Expressions = k;
                SET_VECTOR_ELT(value, i, SetOption(tag, ScalarInteger(k)));
            }
            else if (streql(CHAR(namei), "keep.source"))
            {
                if (TYPEOF(argi) != LGLSXP || LENGTH(argi) != 1)
                    errorcall(call, "keep.source parameter invalid");
                k = asInteger(argi);
                R_KeepSource = k;
                SET_VECTOR_ELT(value, i, SetOption(tag, ScalarLogical(k)));
            }
            else if (streql(CHAR(namei), "editor"))
            {
                s = asChar(argi);
                if (s == NA_STRING || length(s) == 0)
                    errorcall(call, "invalid editor parameter");
                SET_VECTOR_ELT(value, i, SetOption(tag, ScalarString(s)));
            }
            else if (streql(CHAR(namei), "continue"))
            {
                s = asChar(argi);
                if (s == NA_STRING || length(s) == 0)
                    errorcall(call, "invalid continue parameter");
                SET_VECTOR_ELT(value, i, SetOption(tag, ScalarString(s)));
            }
            else if (streql(CHAR(namei), "prompt"))
            {
                s = asChar(argi);
                if (s == NA_STRING || length(s) == 0)
                    errorcall(call, "prompt parameter invalid");
                SET_VECTOR_ELT(value, i, SetOption(tag, ScalarString(s)));
            }
            else if (streql(CHAR(namei), "contrasts"))
            {
                if (TYPEOF(argi) != STRSXP || LENGTH(argi) != 2)
                    errorcall(call, "contrasts parameter invalid");
                SET_VECTOR_ELT(value, i, SetOption(tag, argi));
            }
            else if (streql(CHAR(namei), "check.bounds"))
            {
                if (TYPEOF(argi) != LGLSXP || LENGTH(argi) != 1)
                    errorcall(call, "check.bounds parameter invalid");
                k = asInteger(argi);
                /* R_CheckBounds = k; */
                SET_VECTOR_ELT(value, i, SetOption(tag, ScalarLogical(k)));
            }
            else if (streql(CHAR(namei), "warn"))
            {
                if (!isNumeric(argi) || length(argi) != 1)
                    errorcall(call, "warn parameter invalid");
                SET_VECTOR_ELT(value, i, SetOption(tag, argi));
            }
            else if (streql(CHAR(namei), "warning.length"))
            {
                k = asInteger(argi);
                if (k < 100 || k > 8192)
                    errorcall(call, "warning.length parameter invalid");
                R_WarnLength = k;
                SET_VECTOR_ELT(value, i, SetOption(tag, argi));
            }
            else if (streql(CHAR(namei), "warning.expression"))
            {
                if (!isLanguage(argi) && !isExpression(argi))
                    errorcall(call, "warning.expression parameter invalid");
                SET_VECTOR_ELT(value, i, SetOption(tag, argi));
            }
            else if (streql(CHAR(namei), "error"))
            {
                if (isFunction(argi))
                    argi = makeErrorCall(argi);
                else if (!isLanguage(argi) && !isExpression(argi))
                    errorcall(call, "error parameter invalid");
                SET_VECTOR_ELT(value, i, SetOption(tag, argi));
            }
            /* handle this here to avoid GetOption during error handling */
            else if (streql(CHAR(namei), "show.error.messages"))
            {
                if (!isLogical(argi) && length(argi) != 1)
                    errorcall(call, "show.error.messages parameter invalid");
                SET_VECTOR_ELT(value, i, SetOption(tag, argi));
                R_ShowErrorMessages = LOGICAL(argi)[0];
            }
            else if (streql(CHAR(namei), "echo"))
            {
                if (TYPEOF(argi) != LGLSXP || LENGTH(argi) != 1)
                    errorcall(call, "echo parameter invalid");
                k = asInteger(argi);
                /* Should be quicker than checking options(echo)
                   every time R prompts for input:
                   */
                R_Slave = !k;
                SET_VECTOR_ELT(value, i, SetOption(tag, ScalarLogical(k)));
            }
            else
            {
                SET_VECTOR_ELT(value, i, SetOption(tag, duplicate(argi)));
            }
            SET_STRING_ELT(names, i, namei);
        }
        else
        { /* querying arg */
            if (!isString(argi) || LENGTH(argi) <= 0)
                errorcall(call, R_MSG_IA);
            SET_VECTOR_ELT(value, i, duplicate(CAR(FindTaggedItem(options, install(CHAR(STRING_ELT(argi, 0)))))));
            SET_STRING_ELT(names, i, STRING_ELT(argi, 0));
            R_Visible = 1;
        }
    } /* for() */
    setAttrib(value, R_NamesSymbol, names);
    UNPROTECT(2);
    return value;
}
