/*
 *  R : A Computer Language for Statistical Data Analysis
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

#define ARGUSED(x) LEVELS(x)

#include "Defn.h"

SEXP do_browser(SEXP, SEXP, SEXP, SEXP);

#ifdef Macintosh

/* The universe will end if the Stack on the Mac grows til it hits the heap. */
/* This code places a limit on the depth to which eval can recurse. */

#define EVAL_LIMIT 500
void isintrpt()
{
}

#endif

#ifdef Win32
extern void ProcessEvents();
#endif

/* NEEDED: A fixup is needed in browser, because it */
/* can trap errors, and currently does not reset the */
/* limit to the right value. */

/* These are not ifdef'ed to keep main etc happy. */
extern int R_EvalDepth;
extern int R_EvalCount;

/* Return value of "e" evaluated in "rho". */

SEXP eval(SEXP e, SEXP rho)
{
    SEXP op, tmp, val;

    /* The use of depthsave below is necessary because of the possibility */
    /* of non-local returns from evaluation.  Without this an "expression */
    /* too complex error" is quite likely. */

#ifdef EVAL_LIMIT
    int depthsave = R_EvalDepth++;
    if (R_EvalDepth > EVAL_LIMIT)
        error("expression too complex for evaluator\n");
#endif
#ifdef Macintosh
    /* check for a user abort */
    if ((R_EvalCount++ % 100) == 0)
    {
        isintrpt();
    }
#endif
#ifdef Win32
    if ((R_EvalCount++ % 100) == 0)
    {
        ProcessEvents();
        /* Is this safe? R_EvalCount is not used in other part and
         * don't want to overflow it
         */
        R_EvalCount = 0;
    }
#endif

    tmp = R_NilValue; /* -Wall */

    R_Visible = 1;
    switch (TYPEOF(e))
    {
    case NILSXP:
    case LISTSXP:
    case LGLSXP:
    case INTSXP:
    case REALSXP:
    case STRSXP:
    case CPLXSXP:
    case SPECIALSXP:
    case BUILTINSXP:
    case ENVSXP:
    case CLOSXP:
    case VECSXP:
#ifndef OLD
    case EXPRSXP:
#endif
        tmp = e;
        break;
    case SYMSXP:
        R_Visible = 1;
        if (e == R_DotsSymbol)
            error("... used in an incorrect context\n");
        if (DDVAL(e))
            tmp = ddfindVar(e, rho);
        else
            tmp = findVar(e, rho);
        if (tmp == R_UnboundValue)
            error("Object \"%s\" not found\n", CHAR(PRINTNAME(e)));
        /* if ..d is missing then ddfindVar will signal */
        else if (tmp == R_MissingArg && !DDVAL(e))
        {
            char *n = CHAR(PRINTNAME(e));
            if (*n)
                error("Argument \"%s\" is missing, with no default\n", CHAR(PRINTNAME(e)));
            else
                error("Argument is missing, with no default\n");
        }
        else if (TYPEOF(tmp) == PROMSXP)
        {
            PROTECT(tmp);
            tmp = eval(tmp, rho);
#ifdef old
            if (NAMED(tmp) == 1)
                NAMED(tmp) = 2;
            else
                NAMED(tmp) = 1;
#else
            NAMED(tmp) = 2;
#endif
            UNPROTECT(1);
        }
#ifdef OLD
        else if (!isNull(tmp))
            NAMED(tmp) = 1;
#else
        else if (!isNull(tmp) && NAMED(tmp) < 1)
            NAMED(tmp) = 1;
#endif
        break;
    case PROMSXP:
        if (PRVALUE(e) == R_UnboundValue)
        {
            if (PRSEEN(e))
                errorcall(R_GlobalContext->call, "recursive default argument reference\n");
            PRSEEN(e) = 1;
            val = eval(PREXPR(e), PRENV(e));
            PRSEEN(e) = 0;
            PRVALUE(e) = val;
        }
        tmp = PRVALUE(e);
        break;
#ifdef OLD
    case EXPRSXP: {
        int i, n;
        n = LENGTH(e);
        for (i = 0; i < n; i++)
            tmp = eval(VECTOR(e)[i], rho);
    }
    break;
#endif
    case LANGSXP:
        if (TYPEOF(CAR(e)) == SYMSXP)
            PROTECT(op = findFun(CAR(e), rho));
        else
            PROTECT(op = eval(CAR(e), rho));
        if (TRACE(op))
        {
            Rprintf("trace: ");
            PrintValue(e);
        }
        if (TYPEOF(op) == SPECIALSXP)
        {
            int save = R_PPStackTop;
            PROTECT(CDR(e));
            R_Visible = 1 - PRIMPRINT(op);
            tmp = PRIMFUN(op)(e, op, CDR(e), rho);
            UNPROTECT(1);
            if (save != R_PPStackTop)
            {
                Rprintf("stack imbalance in %s, %d then %d\n", PRIMNAME(op), save, R_PPStackTop);
            }
        }
        else if (TYPEOF(op) == BUILTINSXP)
        {
            int save = R_PPStackTop;
            PROTECT(tmp = evalList(CDR(e), rho));
            R_Visible = 1 - PRIMPRINT(op);
            tmp = PRIMFUN(op)(e, op, tmp, rho);
            UNPROTECT(1);
            if (save != R_PPStackTop)
            {
                Rprintf("stack imbalance in %s, %d then %d\n", PRIMNAME(op), save, R_PPStackTop);
            }
        }
        else if (TYPEOF(op) == CLOSXP)
        {
            PROTECT(tmp = promiseArgs(CDR(e), rho));
            tmp = applyClosure(e, op, tmp, rho, R_NilValue);
            UNPROTECT(1);
        }
        else
            error("attempt to apply non-function\n");
        UNPROTECT(1);
        break;
    case DOTSXP:
        error("... used in an incorrect context\n");
    default:
        UNIMPLEMENTED("eval");
    }
#ifdef EVAL_LIMIT
    R_EvalDepth = depthsave;
#endif
    return (tmp);
}

/* Apply SEXP op of type CLOSXP to actuals */

SEXP applyClosure(SEXP call, SEXP op, SEXP arglist, SEXP rho, SEXP suppliedenv)
{
    int nargs;
    SEXP body, formals, actuals, savedrho;
    volatile SEXP newrho;
    SEXP f, a, tmp;
    RCNTXT cntxt;

    /* formals = list of formal parameters */
    /* actuals = values to be bound to formals */
    /* arglist = the tagged list of arguments */

    formals = FORMALS(op);
    body = BODY(op);
    savedrho = CLOENV(op);
    nargs = length(arglist);

    /*set up a context with the call in it so error has access to it */
    begincontext(&cntxt, CTXT_RETURN, call, newrho, rho, arglist);

    PROTECT(actuals = matchArgs(formals, arglist));
    PROTECT(newrho = extendEnv(savedrho, formals, actuals));

    /* Use default code for unbound formals */
    f = formals;
    a = actuals;
    while (f != R_NilValue)
    {
        if (CAR(a) == R_MissingArg && CAR(f) != R_MissingArg)
        {
            CAR(a) = mkPROMISE(CAR(f), newrho);
            MISSING(a) = 2;
        }
        f = CDR(f);
        a = CDR(a);
    }

    /* Fix up any extras that were supplied by usemethod. */
    if (suppliedenv != R_NilValue)
    {
        for (tmp = FRAME(suppliedenv); tmp != R_NilValue; tmp = CDR(tmp))
        {
            for (a = actuals; a != R_NilValue; a = CDR(a))
                if (TAG(a) == TAG(tmp))
                    break;
            if (a == R_NilValue)
            {
                FRAME(newrho) = CONS(CAR(tmp), FRAME(newrho));
                TAG(FRAME(newrho)) = TAG(tmp);
            }
        }
    }
    NARGS(newrho) = nargs;

    /* Terminate the previous context and start */
    /* a new one with the correct environment. */

    endcontext(&cntxt);

    /* If we have a generic function we need to use the */
    /* sysparent of the generic as the sysparent of the */
    /* method because the method is a straight substitution */
    /* of the generic. */

    if (R_GlobalContext->callflag == CTXT_GENERIC)
        begincontext(&cntxt, CTXT_RETURN, call, newrho, R_GlobalContext->sysparent, arglist);
    else
        begincontext(&cntxt, CTXT_RETURN, call, newrho, rho, arglist);

    /* The default return value */

    tmp = R_NilValue;

    /* Debugging */

    DEBUG(newrho) = DEBUG(op);
    if (DEBUG(op))
    {
        Rprintf("debugging in: ");
        PrintValueRec(call, rho);
        /* now find out if the body is function with only one statement */
        if (isSymbol(CAR(body)))
            tmp = findFun(CAR(body), rho);
        else
            tmp = eval(CAR(body), rho);
        if ((TYPEOF(tmp) == BUILTINSXP || TYPEOF(tmp) == SPECIALSXP) && !strcmp(PRIMNAME(tmp), "for") &&
            !strcmp(PRIMNAME(tmp), "{") && !strcmp(PRIMNAME(tmp), "repeat") && !strcmp(PRIMNAME(tmp), "while"))
            goto regdb;
        Rprintf("debug: ");
        PrintValue(body);
        do_browser(call, op, arglist, newrho);
    }
regdb:
    /* Set a longjmp target which will catch any */
    /* explicit returns from the function body. */

    if (SETJMP(cntxt.cjmpbuf))
    {
        PROTECT(tmp = R_ReturnedValue);
    }
    else
    {
        PROTECT(tmp = eval(body, newrho));
    }

    endcontext(&cntxt);

    if (DEBUG(op))
    {
        Rprintf("exiting from: ");
        PrintValueRec(call, rho);
    }
    UNPROTECT(3);
    return (tmp);
}

static SEXP EnsureLocal(SEXP symbol, SEXP rho)
{
    SEXP vl;

    if ((vl = findVarInFrame(FRAME(rho), symbol)) != R_UnboundValue)
    {
        vl = eval(symbol, rho); /* for promises */
        if (NAMED(vl) == 2)
        {
            PROTECT(vl = duplicate(vl));
            defineVar(symbol, vl, rho);
            UNPROTECT(1);
        }
        return vl;
    }

    vl = eval(symbol, ENCLOS(rho));
    if (vl == R_UnboundValue)
        error("Object \"%s\" not found\n", CHAR(PRINTNAME(symbol)));

    PROTECT(vl = duplicate(vl));
    defineVar(symbol, vl, rho);
    UNPROTECT(1);
    NAMED(vl) = 1;
    return vl;
}

/* Note: If val is a language object it must be protected */
/* to prevent evaluation.  As an example consider */
/* e <- quote(f(x=1,y=2); names(e) <- c("","a","b") */

static SEXP replaceCall(SEXP fun, SEXP val, SEXP args, SEXP rhs)
{
    SEXP tmp, ptmp;
    PROTECT(fun);
    PROTECT(args);
    PROTECT(rhs);
    PROTECT(val);
    ptmp = tmp = allocList(length(args) + 3);
    UNPROTECT(4);
    CAR(ptmp) = fun;
    ptmp = CDR(ptmp);
    CAR(ptmp) = val;
    ptmp = CDR(ptmp);
    while (args != R_NilValue)
    {
        CAR(ptmp) = CAR(args);
        ptmp = CDR(ptmp);
        args = CDR(args);
    }
    CAR(ptmp) = rhs;
    TAG(ptmp) = install("value");
    TYPEOF(tmp) = LANGSXP;
    return tmp;
}

static SEXP assignCall(SEXP op, SEXP symbol, SEXP fun, SEXP val, SEXP args, SEXP rhs)
{
    PROTECT(op);
    PROTECT(symbol);
    val = replaceCall(fun, val, args, rhs);
    UNPROTECT(2);
    return lang3(op, symbol, val);
}

SEXP do_if(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    int cond = asLogical(eval(CAR(args), rho));
    if (cond == NA_LOGICAL)
        errorcall(call, "missing value where logical needed\n");
    else if (cond)
        return (eval(CAR(CDR(args)), rho));
    else if (length(args) > 2)
        return (eval(CAR(CDR(CDR(args))), rho));
    R_Visible = 0;
    return R_NilValue;
}

SEXP do_for(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    int tmp, dbg;
    volatile int i, n, bgn;
    SEXP sym, body;
    volatile SEXP ans, v, val;
    RCNTXT cntxt;

    sym = CAR(args);
    val = CADR(args);
    body = CADDR(args);

    PROTECT(args);
    PROTECT(rho);
    PROTECT(val = eval(val, rho));
    defineVar(sym, R_NilValue, rho);
    if (isList(val) || isNull(val))
    {
        n = length(val);
        PROTECT(v = R_NilValue);
    }
    else
    {
        n = LENGTH(val);
        PROTECT(v = allocVector(TYPEOF(val), 1));
    }
    ans = R_NilValue;

    dbg = DEBUG(rho);
    if (isLanguage(body) && isSymbol(CAR(body)) && strcmp(CHAR(PRINTNAME(CAR(body))), "{"))
        bgn = 1;
    else
        bgn = 0;

    for (i = 0; i < n; i++)
    {
        if (DEBUG(rho) && bgn)
        {
            Rprintf("debug: ");
            PrintValue(body);
            do_browser(call, op, args, rho);
        }
        begincontext(&cntxt, CTXT_LOOP, R_NilValue, R_NilValue, R_NilValue, R_NilValue);
        if ((tmp = SETJMP(cntxt.cjmpbuf)))
        {
            if (tmp == CTXT_BREAK)
                break; /* break */
            else
                continue; /* next	 */
        }
        else
        {
            if (isVector(v))
            {
                UNPROTECT(1);
                PROTECT(v = allocVector(TYPEOF(val), 1));
            }
            switch (TYPEOF(val))
            {
            case LGLSXP:
                LOGICAL(v)[0] = LOGICAL(val)[i];
                setVar(sym, v, rho);
                ans = eval(body, rho);
                break;
            case INTSXP:
                INTEGER(v)[0] = INTEGER(val)[i];
                setVar(sym, v, rho);
                ans = eval(body, rho);
                break;
            case REALSXP:
                REAL(v)[0] = REAL(val)[i];
                setVar(sym, v, rho);
                ans = eval(body, rho);
                break;
            case CPLXSXP:
                COMPLEX(v)[0] = COMPLEX(val)[i];
                setVar(sym, v, rho);
                ans = eval(body, rho);
                break;
            case STRSXP:
                STRING(v)[0] = STRING(val)[i];
                setVar(sym, v, rho);
                ans = eval(body, rho);
                break;
            case EXPRSXP:
            case VECSXP:
                setVar(sym, VECTOR(val)[i], rho);
                ans = eval(body, rho);
                break;
            case LISTSXP:
                setVar(sym, CAR(val), rho);
                ans = eval(body, rho);
                val = CDR(val);
            }
            endcontext(&cntxt);
        }
    }
    UNPROTECT(4);
    R_Visible = 0;
    DEBUG(rho) = dbg;
    return ans;
}

SEXP do_while(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    int cond, dbg;
    volatile int bgn;
    volatile SEXP s, t;
    RCNTXT cntxt;

    checkArity(op, args);
    s = eval(CAR(args), rho); /* ??? */

    dbg = DEBUG(rho);
    t = CAR(CADR(args));
    if (isSymbol(t) && strcmp(CHAR(PRINTNAME(t)), "{"))
        bgn = 1;
    t = R_NilValue;
    for (;;)
    {
        if ((cond = asLogical(s)) == NA_LOGICAL)
            errorcall(call, "missing value where logical needed\n");
        else if (!cond)
            break;
        if (bgn && DEBUG(rho))
        {
            Rprintf("debug: ");
            PrintValue(CAR(args));
            do_browser(call, op, args, rho);
        }
        begincontext(&cntxt, CTXT_LOOP, R_NilValue, R_NilValue, R_NilValue, R_NilValue);
        if ((cond = SETJMP(cntxt.cjmpbuf)))
        {
            if (cond == CTXT_BREAK)
                break; /* break */
            else
                continue; /* next	 */
        }
        else
        {
            PROTECT(t = eval(CAR(CDR(args)), rho));
            s = eval(CAR(args), rho);
            UNPROTECT(1);
            endcontext(&cntxt);
        }
    }
    R_Visible = 0;
    DEBUG(rho) = dbg;
    return t;
}

SEXP do_repeat(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    int cond, dbg;
    volatile int bgn;
    volatile SEXP t;
    RCNTXT cntxt;

    checkArity(op, args);

    dbg = DEBUG(rho);
    if (isSymbol(CAR(args)) && strcmp(CHAR(PRINTNAME(CAR(args))), "{"))
        bgn = 1;

    t = R_NilValue;
    for (;;)
    {
        if (DEBUG(rho) && bgn)
        {
            Rprintf("debug: ");
            PrintValue(CAR(args));
            do_browser(call, op, args, rho);
        }
        begincontext(&cntxt, CTXT_LOOP, R_NilValue, R_NilValue, R_NilValue, R_NilValue);
        if ((cond = SETJMP(cntxt.cjmpbuf)))
        {
            if (cond == CTXT_BREAK)
                break; /*break */
            else
                continue; /* next */
        }
        else
        {
            t = eval(CAR(args), rho);
            endcontext(&cntxt);
        }
    }
    R_Visible = 0;
    DEBUG(rho) = dbg;
    return t;
}

SEXP do_break(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    findcontext(PRIMVAL(op), rho, R_NilValue);
    return R_NilValue;
}

SEXP do_paren(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    checkArity(op, args);
    return CAR(args);
}

SEXP do_begin(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP s;
    if (args == R_NilValue)
    {
        s = R_NilValue;
    }
    else
    {
        while (args != R_NilValue)
        {
            if (DEBUG(rho))
            {
                Rprintf("debug: ");
                PrintValue(CAR(args));
                do_browser(call, op, args, rho);
            }
            s = eval(CAR(args), rho);
            args = CDR(args);
        }
    }
    return s;
}

SEXP do_return(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP a, v, vals;
    int nv = 0;
    /* We do the evaluation here so that we */
    /* can tag any untagged return values if */
    /* they are specified by symbols */
    PROTECT(vals = evalList(args, rho));
    a = args;
    v = vals;
    while (!isNull(a))
    {
        nv += 1;
        if (isNull(TAG(a)) && isSymbol(CAR(a)))
            TAG(v) = CAR(a);
        if (NAMED(CAR(v)) > 1)
            CAR(v) = duplicate(CAR(v));
        a = CDR(a);
        v = CDR(v);
    }
    UNPROTECT(1);
    switch (nv)
    {
    case 0:
        v = R_NilValue;
        break;
    case 1:
        v = CAR(vals);
        break;
    default:
        v = PairToVectorList(vals);
        break;
    }
    if (R_BrowseLevel)
        findcontext(CTXT_BROWSER | CTXT_RETURN, rho, v);
    else
        findcontext(CTXT_RETURN, rho, v);
    return R_NilValue; /*NOTREACHED*/
}

SEXP do_function(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    if (length(args) < 2)
        WrongArgCount("lambda");
    CheckFormals(CAR(args));
    return mkCLOSXP(CAR(args), CADR(args), rho);
}

/*
 *  Assignments for complex LVAL specifications. This is the stuff that
 *  nightmares are made of ...	Note that "evalseq" preprocesses the LHS
 *  of an assignment.  Given an expression, it builds a list of partial
 *  values for the exression.  For example, the assignment x$a[3] <- 10
 *  with LHS x$a[3] yields the (improper) list:
 *
 *	 (eval(x$a[3])	eval(x$a)  eval(x)  .  x)
 *
 *  (Note the terminating symbol).  The partial evaluations are carried
 *  out efficiently using previously computed components.
 */

static SEXP evalseq(SEXP expr, SEXP rho, int forcelocal, SEXP tmploc)
{
    SEXP val, nval, nexpr;
    if (isNull(expr))
        error("invalid (NULL) left side of assignment\n");
    if (isSymbol(expr))
    {
        PROTECT(expr);
        if (forcelocal)
        {
            nval = EnsureLocal(expr, rho);
        }
        else
        {
            nval = eval(expr, rho);
        }
        UNPROTECT(1);
        return CONS(nval, expr);
    }
    else if (isLanguage(expr))
    {
        PROTECT(expr);
        PROTECT(val = evalseq(CADR(expr), rho, forcelocal, tmploc));
        CAR(tmploc) = CAR(val);
        PROTECT(nexpr = LCONS(TAG(tmploc), CDDR(expr)));
        PROTECT(nexpr = LCONS(CAR(expr), nexpr));
        nval = eval(nexpr, rho);
        UNPROTECT(4);
        return CONS(nval, val);
    }
    else
        error("Target of assignment expands to non-language object\n");
    return R_NilValue; /*NOTREACHED*/
}

/* Main entry point for complex assignments */
/* We have checked to see that CAR(args) is a LANGSXP */

static char *asym[] = {":=", "<-", "<<-"};

SEXP applydefine(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP expr, lhs, rhs, saverhs, tmp, tmploc, tmpsym;
    char buf[32];

    expr = CAR(args);

    /* it's important that the rhs get evaluated first */
    /* because assignment is right associative i.e. */
    /* a <- b <- c is parsed as a <- (b <- c) */

    PROTECT(saverhs = rhs = eval(CADR(args), rho));
    tmpsym = install("*tmp*");
    defineVar(tmpsym, R_NilValue, rho);
    tmploc = FRAME(rho);
    while (tmploc != R_NilValue && TAG(tmploc) != tmpsym)
        tmploc = CDR(tmploc);

    /* do a partial evaluation down through the lhs */
    lhs = evalseq(CADR(expr), rho, PRIMVAL(op) == 1, tmploc);

    PROTECT(lhs);
    PROTECT(rhs); /*just to get the loop right */

    while (isLanguage(CADR(expr)))
    {
        sprintf(buf, "%s<-", CHAR(PRINTNAME(CAR(expr))));
        tmp = install(buf);
        UNPROTECT(1);
        CAR(tmploc) = CAR(lhs);
        PROTECT(rhs = replaceCall(tmp, TAG(tmploc), CDDR(expr), rhs));
        rhs = eval(rhs, rho);
        UNPROTECT(1);
        PROTECT(rhs);
        lhs = CDR(lhs);
        expr = CADR(expr);
    }
    sprintf(buf, "%s<-", CHAR(PRINTNAME(CAR(expr))));
    CAR(tmploc) = CAR(lhs);
    PROTECT(expr = assignCall(install(asym[PRIMVAL(op)]), CDR(lhs), install(buf), TAG(tmploc), CDDR(expr), rhs));
    expr = eval(expr, rho);
    UNPROTECT(4);
    unbindVar(tmpsym, rho);
    return duplicate(saverhs);
}

SEXP do_alias(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    checkArity(op, args);
    NAMED(CAR(args)) = 0;
    return CAR(args);
}

/*  Assignment in its various forms  */

SEXP do_set(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP s, t;
    if (length(args) != 2)
        WrongArgCount(asym[PRIMVAL(op)]);
    if (isString(CAR(args)))
        CAR(args) = install(CHAR(STRING(CAR(args))[0]));

    switch (PRIMVAL(op))
    {
    case 1: /* <- */
        if (isSymbol(CAR(args)))
        {
            s = eval(CADR(args), rho);
            if (NAMED(s))
            {
                PROTECT(s);
                t = duplicate(s);
                UNPROTECT(1);
                s = t;
            }
            PROTECT(s);
            R_Visible = 0;
            defineVar(CAR(args), s, rho);
            UNPROTECT(1);
            NAMED(s) = 1;
            return (s);
        }
        else if (isLanguage(CAR(args)))
        {
            R_Visible = 0;
            return applydefine(call, op, args, rho);
        }
        else
            errorcall(call, "invalid (do_set) left-hand side to assignment\n");
    case 2: /* <<- */
        if (isSymbol(CAR(args)))
        {
            s = eval(CADR(args), rho);
            if (NAMED(s))
                s = duplicate(s);
            PROTECT(s);
            R_Visible = 0;
            setVar(CAR(args), s, ENCLOS(rho));
            UNPROTECT(1);
            NAMED(s) = 1;
            return s;
        }
        else if (isLanguage(CAR(args)))
            return applydefine(call, op, args, rho);
        else
            error("invalid assignment lhs\n");

    default:
        UNIMPLEMENTED("do_set");
    }
    return R_NilValue; /*NOTREACHED*/
}

/* Evaluate each expression in "el" in the environment "rho".  This is */
/* a naturally recursive algorithm, but we use the iterative form below */
/* because it is does not cause growth of the pointer protection stack, */
/* and because it is a little more efficient. */

SEXP evalList(SEXP el, SEXP rho)
{
    SEXP ans, h, tail;

    PROTECT(ans = tail = CONS(R_NilValue, R_NilValue));

    while (el != R_NilValue)
    {

        /* If we have a ... symbol, we look to see what it is */
        /* bound to.  If its binding is Null (i.e. zero length) */
        /* we just ignore it and return the cdr with all its */
        /* expressions evaluated; if it is bound to a ... list */
        /* of promises, we force all the promises and then splice */
        /* the list of resulting values into the return value.	*/
        /* Anything else bound to a ... symbol is an error */

        if (CAR(el) == R_DotsSymbol)
        {
            h = findVar(CAR(el), rho);
            if (TYPEOF(h) == DOTSXP || h == R_NilValue)
            {
                while (h != R_NilValue)
                {
                    CDR(tail) = CONS(eval(CAR(h), rho), R_NilValue);
                    TAG(CDR(tail)) = CreateTag(TAG(h));
                    tail = CDR(tail);
                    h = CDR(h);
                }
            }
            else if (h != R_MissingArg)
                error("... used in an incorrect context\n");
        }
        else if (CAR(el) != R_MissingArg)
        {
            CDR(tail) = CONS(eval(CAR(el), rho), R_NilValue);
            tail = CDR(tail);
            TAG(tail) = CreateTag(TAG(el));
        }
        el = CDR(el);
    }
    UNPROTECT(1);
    return CDR(ans);
}

/* A slight variation of evaluating each expression in "el" in "rho". */
/* This is a naturally recursive algorithm, but we use the iterative */
/* form below because it is does not cause growth of the pointer */
/* protection stack, and because it is a little more efficient. */

SEXP evalListKeepMissing(SEXP el, SEXP rho)
{
    SEXP ans, h, tail;

    PROTECT(ans = tail = CONS(R_NilValue, R_NilValue));

    while (el != R_NilValue)
    {

        /* If we have a ... symbol, we look to see what it is */
        /* bound to.  If its binding is Null (i.e. zero length) */
        /* we just ignore it and return the cdr with all its */
        /* expressions evaluated; if it is bound to a ... list */
        /* of promises, we force all the promises and then splice */
        /* the list of resulting values into the return value.	*/
        /* Anything else bound to a ... symbol is an error */

        if (CAR(el) == R_DotsSymbol)
        {
            h = findVar(CAR(el), rho);
            if (TYPEOF(h) == DOTSXP || h == R_NilValue)
            {
                while (h != R_NilValue)
                {
                    if (CAR(h) == R_MissingArg)
                        CDR(tail) = CONS(R_MissingArg, R_NilValue);
                    else
                        CDR(tail) = CONS(eval(CAR(h), rho), R_NilValue);
                    TAG(CDR(tail)) = CreateTag(TAG(h));
                    tail = CDR(tail);
                    h = CDR(h);
                }
            }
            else if (h != R_MissingArg)
                error("... used in an incorrect context\n");
        }
        else if (CAR(el) == R_MissingArg)
        {
            CDR(tail) = CONS(R_MissingArg, R_NilValue);
            tail = CDR(tail);
            TAG(tail) = CreateTag(TAG(el));
        }
        else
        {
            CDR(tail) = CONS(eval(CAR(el), rho), R_NilValue);
            tail = CDR(tail);
            TAG(tail) = CreateTag(TAG(el));
        }
        el = CDR(el);
    }
    UNPROTECT(1);
    return CDR(ans);
}

/* Create a promise to evaluate each argument.	Although this is most */
/* naturally attacked with a recursive algorithm, we use the iterative */
/* form below because it is does not cause growth of the pointer */
/* protection stack, and because it is a little more efficient. */

SEXP promiseArgs(SEXP el, SEXP rho)
{
    SEXP ans, h, tail;

    PROTECT(ans = tail = CONS(R_NilValue, R_NilValue));

    while (el != R_NilValue)
    {

        /* If we have a ... symbol, we look to see what it */
        /* is bound to.	 If its binding is Null (i.e. zero length) */
        /* we just ignore it and return the cdr with all its */
        /* expressions promised; if it is bound to a ... list */
        /* of promises, we repromise all the promises and then splice */
        /* the list of resulting values into the return value. */
        /* Anything else bound to a ... symbol is an error */

        /* Is this double promise mechanism really needed? */

        if (CAR(el) == R_DotsSymbol)
        {
            h = findVar(CAR(el), rho);
            if (TYPEOF(h) == DOTSXP || h == R_NilValue)
            {
                while (h != R_NilValue)
                {
                    CDR(tail) = CONS(mkPROMISE(CAR(h), rho), R_NilValue);
                    TAG(CDR(tail)) = CreateTag(TAG(h));
                    tail = CDR(tail);
                    h = CDR(h);
                }
            }
            else if (h != R_MissingArg)
                error("... used in an incorrect context\n");
        }
        else if (CAR(el) == R_MissingArg)
        {
            CDR(tail) = CONS(R_MissingArg, R_NilValue);
            tail = CDR(tail);
            TAG(tail) = CreateTag(TAG(el));
        }
        else
        {
            CDR(tail) = CONS(mkPROMISE(CAR(el), rho), R_NilValue);
            tail = CDR(tail);
            TAG(tail) = CreateTag(TAG(el));
        }
        el = CDR(el);
    }
    UNPROTECT(1);
    return CDR(ans);
}

/* Check that each formal is a symbol */

void CheckFormals(SEXP ls)
{
    if (isList(ls))
        for (; ls != R_NilValue; ls = CDR(ls))
            if (TYPEOF(TAG(ls)) != SYMSXP)
                error("invalid formal argument list for \"function\"\n");
}

/* "eval" and "eval.with.vis" : Evaluate the first argument */
/* in the environment specified by the second argument. */

SEXP do_eval(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP expr, env, encl, tmp;
    int nback;
    RCNTXT cntxt;

    checkArity(op, args);
    expr = CAR(args);
    env = CADR(args);
    encl = CADDR(args);
    if (!isNull(encl) && !isEnvironment(encl))
        errorcall(call, "invalid 3rd argument\n");
    switch (TYPEOF(env))
    {
    case NILSXP:
    case ENVSXP:
        PROTECT(env); /* so we can unprotect 2 at the end */
        break;
    case LISTSXP:
        PROTECT(env = allocSExp(ENVSXP));
        FRAME(env) = duplicate(CADR(args));
        ENCLOS(env) = encl;
        break;
    case VECSXP:
        PROTECT(env = allocSExp(ENVSXP));
        FRAME(env) = VectorToPairList(CADR(args));
        ENCLOS(env) = encl;
        break;
    case INTSXP:
    case REALSXP:
        nback = asInteger(env);
        if (nback == NA_INTEGER)
            errorcall(call, "invalid environment\n");
        if (nback > 0)
            nback -= framedepth(R_GlobalContext);
        nback = -nback;
        PROTECT(env = R_sysframe(nback, R_GlobalContext));
        break;
    default:
        errorcall(call, "invalid second argument\n");
    }
    if (isLanguage(expr) || isExpression(expr) || isSymbol(expr))
    {
        PROTECT(expr);
        begincontext(&cntxt, CTXT_RETURN, call, env, rho, args);
        if (!SETJMP(cntxt.cjmpbuf))
            expr = eval(expr, env);
        endcontext(&cntxt);
        UNPROTECT(1);
    }
    if (isExpression(expr))
    {
        int i, n;
        PROTECT(expr);
        n = LENGTH(expr);
        begincontext(&cntxt, CTXT_RETURN, call, env, rho, args);
        if (!SETJMP(cntxt.cjmpbuf))
            for (i = 0; i < n; i++)
                tmp = eval(VECTOR(expr)[i], env);
        endcontext(&cntxt);
        UNPROTECT(1);
        expr = tmp;
    }
    if (PRIMVAL(op))
    {
        PROTECT(expr);
        PROTECT(env = allocVector(VECSXP, 2));
        PROTECT(encl = allocVector(STRSXP, 2));
        STRING(encl)[0] = mkChar("value");
        STRING(encl)[1] = mkChar("visible");
        VECTOR(env)[0] = expr;
        VECTOR(env)[1] = ScalarLogical(R_Visible);
        setAttrib(env, R_NamesSymbol, encl);
        expr = env;
        UNPROTECT(3);
    }
    UNPROTECT(1);
    return expr;
}

SEXP do_recall(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    RCNTXT *cptr;
    SEXP t, s, ans;
    cptr = R_GlobalContext;
    /* get the args supplied */
    while (cptr != NULL)
    {
        if (cptr->callflag == CTXT_RETURN && cptr->cloenv == rho)
            break;
        cptr = cptr->nextcontext;
    }
    args = cptr->promargs;
    /* get the env recall was called from */
    s = R_GlobalContext->sysparent;
    while (cptr != NULL)
    {
        if (cptr->callflag == CTXT_RETURN && cptr->cloenv == s)
            break;
        cptr = cptr->nextcontext;
    }
    if (cptr == NULL)
        error("Recall called from outside a closure\n");
    if (TYPEOF(CAR(cptr->call)) == SYMSXP)
        PROTECT(s = findFun(CAR(cptr->call), cptr->sysparent));
    else
        PROTECT(s = eval(CAR(cptr->call), cptr->sysparent));
    ans = applyClosure(cptr->call, s, args, cptr->sysparent, R_NilValue);
    UNPROTECT(1);
    return ans;
}

SEXP EvalArgs(SEXP el, SEXP rho, int dropmissing)
{
    if (dropmissing)
        return evalList(el, rho);
    else
        return evalListKeepMissing(el, rho);
}

/* DispatchOrEval is used in internal functions which dispatch to
 * object methods (e.g. "[" or "[[").  The code either builds promises
 * and dispatches to the appropriate method, or it evaluates the
 * (unevaluated) arguments it comes in with and returns them so that
 * the generic built-in C code can continue.

 * To call this an ugly hack would be to insult all existing ugly hacks
 * at large in the world.
 */
int DispatchOrEval(SEXP call, SEXP op, SEXP args, SEXP rho, SEXP *ans, int dropmissing)
{
    SEXP x;
    char *pt, buf[128];
    RCNTXT cntxt;

    /* NEW */
    PROTECT(args = promiseArgs(args, rho));
    PROTECT(x = eval(CAR(args), rho));
    pt = strrchr(CHAR(PRINTNAME(CAR(call))), '.');

    if (isObject(x) && (pt == NULL || strcmp(pt, ".default")))
    {
        /* PROTECT(args = promiseArgs(args, rho)); */
        PRVALUE(CAR(args)) = x;
        sprintf(buf, "%s", CHAR(PRINTNAME(CAR(call))));
        begincontext(&cntxt, CTXT_RETURN, call, rho, rho, args);
        if (usemethod(buf, x, call, args, rho, ans))
        {
            endcontext(&cntxt);
            UNPROTECT(2);
            return 1;
        }
        endcontext(&cntxt);
    }
    /* else PROTECT(args); */
    *ans = CONS(x, EvalArgs(CDR(args), rho, dropmissing));
    TAG(*ans) = CreateTag(TAG(args));
    UNPROTECT(2);
    return 0;
}

/* Compare two class attributes and decide whether one dominates the */
/* other.  For now this is a pretty simple version but it could get */
/* more complicated comparing factors (ordered vs not) etc. */

static SEXP dominates(SEXP cclass, SEXP newobj)
{
    SEXP t;
    t = getAttrib(newobj, R_ClassSymbol);
    if (cclass == R_NilValue)
        return t;
    return cclass;
}

int DispatchGroup(char *group, SEXP call, SEXP op, SEXP args, SEXP rho, SEXP *ans)
{
    int i, j, nargs, whichclass, set;
    SEXP class, s, t, m, meth, sxp, gr, newrho;
    char buf[512], generic[128], *pt;
    /* check whether we are processing the default method */
    sprintf(buf, "%s", PRIMNAME(op));
    pt = strtok(buf, ".");
    pt = strtok(NULL, ".");

    if (pt != NULL && !strcmp(pt, "default"))
        return 0;

    if (!strcmp(group, "Ops"))
        nargs = length(args);
    else
        nargs = 1;

    if (nargs == 1 && !isObject(CAR(args)))
        return 0;

    if (!isObject(CAR(args)) && !isObject(CADR(args)))
        return 0;

    class = getAttrib(CAR(args), R_ClassSymbol);
    if (nargs == 2)
        class = dominates(class, CADR(args));

    sprintf(generic, "%s", PRIMNAME(op));

    j = length(class);
    sxp = R_NilValue;
    /* -Wall: */ gr = R_NilValue;
    meth = R_NilValue;
    /* Need to interleave looking for group and generic methods */
    /* eg if class(x) is "foo" "bar" then x>3 should invoke */
    /* "Ops.foo" rather than ">.bar" */
    for (whichclass = 0; whichclass < j; whichclass++)
    {
        sprintf(buf, "%s.%s", generic, CHAR(STRING(class)[whichclass]));
        meth = install(buf);
        sxp = findVar(meth, rho);
        if (isFunction(sxp))
        {
            PROTECT(gr = mkString(""));
            break;
        }
        sprintf(buf, "%s.%s", group, CHAR(STRING(class)[whichclass]));
        meth = install(buf);
        sxp = findVar(meth, rho);
        if (isFunction(sxp))
        {
            PROTECT(gr = mkString(group));
            break;
        }
    }
    if (!isFunction(sxp)) /* no generic or group method so use default*/
        return 0;

    /* we either have a group method or a class method */

    PROTECT(newrho = allocSExp(ENVSXP));
    PROTECT(m = allocVector(STRSXP, nargs));
    s = args;
    for (i = 0; i < nargs; i++)
    {
        t = getAttrib(CAR(s), R_ClassSymbol);
        set = 0;
        if (isString(t))
        {
            for (j = 0; j < length(t); j++)
            {
                if (!strcmp(CHAR(STRING(t)[j]), CHAR(STRING(class)[whichclass])))
                {
                    STRING(m)[i] = mkChar(buf);
                    set = 1;
                    break;
                }
            }
        }
        if (!set)
            STRING(m)[i] = R_BlankString;
        s = CDR(s);
    }

    defineVar(install(".Method"), m, newrho);
    UNPROTECT(1);
    PROTECT(t = mkString(generic));
    defineVar(install(".Generic"), t, newrho);
    UNPROTECT(1);
    defineVar(install(".Group"), gr, newrho);
    set = length(class) - whichclass;
    PROTECT(t = allocVector(STRSXP, set));
    for (j = 0; j < set; j++)
        STRING(t)[j] = duplicate(STRING(class)[whichclass++]);
    defineVar(install(".Class"), t, newrho);
    UNPROTECT(1);

    PROTECT(t = LCONS(meth, CDR(call)));

    /* the arguments have been evaluated; since we are passing them */
    /* out to a closure we need to wrap them in promises so that */
    /* they get duplicated and things like missing/substitute work. */

    PROTECT(s = promiseArgs(CDR(call), rho));
    if (length(s) != length(args))
        errorcall(call, "dispatch error\n");
    for (m = s; m != R_NilValue; m = CDR(m), args = CDR(args))
        PRVALUE(CAR(m)) = CAR(args);

    *ans = applyClosure(t, sxp, s, rho, newrho);
    UNPROTECT(4);
    return 1;
}
