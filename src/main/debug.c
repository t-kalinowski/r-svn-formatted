/*
 *  R : A Computer Langage for Statistical Data Analysis
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
 *  Foundation, Inc., 51 Franklin Street Fifth Floor, Boston, MA 02110-1301  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Defn.h"

SEXP attribute_hidden do_debug(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    checkArity(op, args);
#define find_char_fun                                                                                                  \
    if (isValidString(CAR(args)))                                                                                      \
    {                                                                                                                  \
        SEXP s;                                                                                                        \
        PROTECT(s = install(CHAR(STRING_ELT(CAR(args), 0))));                                                          \
        SETCAR(args, findFun(s, rho));                                                                                 \
        UNPROTECT(1);                                                                                                  \
    }
    find_char_fun

        if (TYPEOF(CAR(args)) != CLOSXP) errorcall(call, "argument must be a function");
    switch (PRIMVAL(op))
    {
    case 0:
        SET_DEBUG(CAR(args), 1);
        break;
    case 1:
        if (DEBUG(CAR(args)) != 1)
            warningcall(call, "argument is not being debugged");
        SET_DEBUG(CAR(args), 0);
        break;
    }
    return R_NilValue;
}

SEXP attribute_hidden do_trace(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    checkArity(op, args);

    find_char_fun

        if (TYPEOF(CAR(args)) != CLOSXP && TYPEOF(CAR(args)) != BUILTINSXP && TYPEOF(CAR(args)) != SPECIALSXP)
            errorcall(call, "argument must be a function");

    switch (PRIMVAL(op))
    {
    case 0:
        SET_TRACE(CAR(args), 1);
        break;
    case 1:
        SET_TRACE(CAR(args), 0);
        break;
    }
    return R_NilValue;
}

/* maintain global trace state */

static Rboolean tracing_state = TRUE;
#define GET_TRACE_STATE tracing_state
#define SET_TRACE_STATE(value) tracing_state = value

SEXP R_traceOnOff(SEXP onOff)
{
    SEXP value;
    Rboolean prev = GET_TRACE_STATE;
    if (length(onOff) > 0)
    {
        Rboolean new = asLogical(onOff);
        if (new == TRUE || new == FALSE)
            SET_TRACE_STATE(new);
        else
            error("Value for tracingState must be TRUE or FALSE");
    }
    value = allocVector(LGLSXP, 1);
    LOGICAL(value)[0] = prev;
    return value;
}

Rboolean R_current_trace_state()
{
    return GET_TRACE_STATE;
}

/* memory tracing */
/* report when a traced object is duplicated */

SEXP attribute_hidden do_memtrace(SEXP call, SEXP op, SEXP args, SEXP rho)
{
#ifdef R_MEMORY_PROFILING
    SEXP object;
    char buffer[20];

    checkArity(op, args);

    object = CAR(args);
    if (TYPEOF(object) == CLOSXP || TYPEOF(object) == BUILTINSXP || TYPEOF(object) == SPECIALSXP)
        errorcall(call, "argument must not be a function");

    if (object == R_NilValue)
        errorcall(call, "cannot trace NULL");

    if (TYPEOF(object) == ENVSXP || TYPEOF(object) == PROMSXP)
        errorcall(call, "memtrace is not useful for promise and environment objects");
    if (TYPEOF(object) == EXTPTRSXP || TYPEOF(object) == WEAKREFSXP)
        errorcall(call, "memtrace is not useful for weak reference or pointer objects");

    SET_TRACE(object, 1);
    sprintf(buffer, "<%p>", object);
    return mkString(buffer);
#else
    errorcall(call, "R not compiled with memory profiling");
    return R_NilValue;
#endif
}

SEXP attribute_hidden do_memuntrace(SEXP call, SEXP op, SEXP args, SEXP rho)
{
#ifdef R_MEMORY_PROFILING
    SEXP object;

    checkArity(op, args);

    object = CAR(args);
    if (TYPEOF(object) == CLOSXP || TYPEOF(object) == BUILTINSXP || TYPEOF(object) == SPECIALSXP)
        errorcall(call, "argument must not be a function");

    if (TRACE(object))
        SET_TRACE(object, 0);
#else
    error(call, "R not compiled with memory profiling");
#endif
    return R_NilValue;
}

#ifndef R_MEMORY_PROFILING
void memtrace_report(SEXP object)
{
    return;
}
#else
void memtrace_report(SEXP old, SEXP new)
{
    RCNTXT *cptr;

    if (!R_current_trace_state())
        return;
    Rprintf("memtrace[%p->%p]: ", old, new);
    for (cptr = R_GlobalContext; cptr; cptr = cptr->nextcontext)
    {
        if ((cptr->callflag & (CTXT_FUNCTION | CTXT_BUILTIN)) && TYPEOF(cptr->call) == LANGSXP)
        {
            SEXP fun = CAR(cptr->call);
            Rprintf("%s ", TYPEOF(fun) == SYMSXP ? CHAR(PRINTNAME(fun)) : "<Anonymous>");
        }
    }
    Rprintf("\n");
}

#endif /* R_MEMORY_PROFILING */
