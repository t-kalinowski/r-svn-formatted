/*
 *  Routines used in calculating least squares solutions in a
 *  nonlinear model in nls library for R.
 *
 *  Copyright 1999-2001 Douglas M. Bates <bates@stat.wisc.edu>,
 *                      Saikat DebRoy <saikat@stat.wisc.edu>
 *
 *  Copyright 2005 The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied
 *  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the GNU General Public License for more
 *  details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA.
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include "nls.h"

#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

/*
 * get the list element named str. names is the name attribute of list
 */

static SEXP getListElement(SEXP list, SEXP names, char *str)
{
    SEXP elmt = (SEXP)NULL;
    char *tempChar;
    int i;

    for (i = 0; i < length(list); i++)
    {
        tempChar = CHAR(STRING_ELT(names, i));
        if (strcmp(tempChar, str) == 0)
        {
            elmt = VECTOR_ELT(list, i);
            break;
        }
    }
    return elmt;
}

/*
 *  call to nls_iter from R -
 *  .Call("nls_iter", m, control, doTrace)
 *  where m and control are nlsModel and nlsControl objects.
 *  doTrace is a logical value.  The returned value is m.
 */

SEXP nls_iter(SEXP m, SEXP control, SEXP doTraceArg)
{

    double dev, fac, minFac, tolerance, newDev, convNew;
    int i, j, maxIter, hasConverged, nPars, doTrace, evaltotCnt, warnOnly, printEval;
    SEXP tmp, conv, incr, deviance, setPars, getPars, pars, newPars, trace;

    doTrace = asLogical(doTraceArg);

    if (!isNewList(control))
        error(_("'control' must be a list"));
    if (!isNewList(m))
        error(_("'m' must be a list"));

    PROTECT(tmp = getAttrib(control, R_NamesSymbol));

    conv = getListElement(control, tmp, "maxiter");
    if (conv == NULL || !isNumeric(conv))
        error(_("'%s' absent"), "control$maxiter");
    maxIter = asInteger(conv);

    conv = getListElement(control, tmp, "tol");
    if (conv == NULL || !isNumeric(conv))
        error(_("'%s' absent"), "control$tol");
    tolerance = asReal(conv);

    conv = getListElement(control, tmp, "minFactor");
    if (conv == NULL || !isNumeric(conv))
        error(_("'%s' absent"), "control$minFactor");
    minFac = asReal(conv);

    conv = getListElement(control, tmp, "warnOnly");
    if (conv == NULL || !isLogical(conv))
        error(_("'%s' absent"), "control$warnOnly");
    warnOnly = conv;

    conv = getListElement(control, tmp, "printEval");
    if (conv == NULL || !isLogical(conv))
        error(_("'%s' absent"), "control$printEval");
    printEval = conv;

    UNPROTECT(1);

    PROTECT(tmp = getAttrib(m, R_NamesSymbol));

    conv = getListElement(m, tmp, "conv");
    if (conv == NULL || !isFunction(conv))
        error(_("'%s' absent"), "m$conv()");
    PROTECT(conv = lang1(conv));

    incr = getListElement(m, tmp, "incr");
    if (incr == NULL || !isFunction(incr))
        error(_("'%s' absent"), "m$incr()");
    PROTECT(incr = lang1(incr));

    deviance = getListElement(m, tmp, "deviance");
    if (deviance == NULL || !isFunction(deviance))
        error(_("'%s' absent"), "m$deviance()");
    PROTECT(deviance = lang1(deviance));

    trace = getListElement(m, tmp, "trace");
    if (trace == NULL || !isFunction(trace))
        error(_("'%s' absent"), "m$trace()");
    PROTECT(trace = lang1(trace));

    setPars = getListElement(m, tmp, "setPars");
    if (setPars == NULL || !isFunction(setPars))
        error(_("'%s' absent"), "m$setPars()");
    PROTECT(setPars);

    getPars = getListElement(m, tmp, "getPars");
    if (getPars == NULL || !isFunction(getPars))
        error(_("'%s' absent"), "m$getPars()");
    PROTECT(getPars = lang1(getPars));

    PROTECT(pars = eval(getPars, R_GlobalEnv));
    nPars = LENGTH(pars);

    dev = asReal(eval(deviance, R_GlobalEnv));
    if (doTrace)
        eval(trace, R_GlobalEnv);

    fac = 1.0;
    hasConverged = FALSE;

    PROTECT(newPars = allocVector(REALSXP, nPars));
    if (printEval)
        evaltotCnt = 1;
    for (i = 0; i < maxIter; i++)
    {
        SEXP newIncr;
        int evalCnt;
        if ((convNew = asReal(eval(conv, R_GlobalEnv))) < tolerance)
        {
            hasConverged = TRUE;
            break;
        }
        PROTECT(newIncr = eval(incr, R_GlobalEnv));

        if (printEval)
            evalCnt = 1;

        while (fac >= minFac)
        {
            if (printEval)
            {
                Rprintf("\nIteration: %d   Evaluation: %d   Total eval.: %d \n", i + 1, evalCnt, evaltotCnt);
                evalCnt++;
                evaltotCnt++;
            }
            for (j = 0; j < nPars; j++)
                REAL(newPars)[j] = REAL(pars)[j] + fac * REAL(newIncr)[j];

            PROTECT(tmp = lang2(setPars, newPars));
            if (asLogical(eval(tmp, R_GlobalEnv)))
            { /* singular gradient */
                UNPROTECT(11);
                if (warnOnly)
                {
                    warning(_("singular gradient"));
                    return m;
                }
                else
                    error(_("singular gradient"));
            }
            UNPROTECT(1);

            newDev = asReal(eval(deviance, R_GlobalEnv));

            if (newDev <= dev)
            {
                dev = newDev;
                fac = MIN(2 * fac, 1);
                tmp = newPars;
                newPars = pars;
                pars = tmp;
                break;
            }
            fac /= 2.;
        }
        UNPROTECT(1);
        if (fac < minFac)
        {
            UNPROTECT(9);
            if (warnOnly)
            {
                warning(_("step factor %g reduced below 'minFactor' of %g"), fac, minFac);
                return m;
            }
            else
                error(_("step factor %g reduced below 'minFactor' of %g"), fac, minFac);
        }
        if (doTrace)
            eval(trace, R_GlobalEnv);
    }

    UNPROTECT(9);
    if (!hasConverged)
    {
        if (warnOnly)
            warning(_("number of iterations exceeded maximum of %d"), maxIter);
        else
            error(_("number of iterations exceeded maximum of %d"), maxIter);
    }
    return m;
}

/*
 *  call to numeric_deriv from R -
 *  .Call("numeric_deriv", expr, theta, rho)
 *  Returns: ans
 */

SEXP numeric_deriv(SEXP expr, SEXP theta, SEXP rho, SEXP dir)
{
    SEXP ans, gradient, pars;
    double eps = sqrt(DOUBLE_EPS), *rDir;
    int start, i, j, k, lengthTheta = 0;

    if (!isString(theta))
        error(_("'theta' should be of type character"));
    if (isNull(rho))
    {
        error(_("use of NULL environment is defunct"));
        rho = R_BaseEnv;
    }
    else if (!isEnvironment(rho))
        error(_("'rho' should be an environment"));
    if (TYPEOF(dir) != REALSXP || LENGTH(dir) != LENGTH(theta))
        error(_("'dir' is not a numeric vector of the correct length"));
    rDir = REAL(dir);

    PROTECT(pars = allocVector(VECSXP, LENGTH(theta)));

    if (TYPEOF(expr) == SYMSXP)
        PROTECT(ans = duplicate(eval(expr, rho)));
    else
        PROTECT(ans = eval(expr, rho));

    if (!isReal(ans))
    {
        SEXP temp = coerceVector(ans, REALSXP);
        UNPROTECT(1);
        PROTECT(ans = temp);
    }
    for (i = 0; i < LENGTH(ans); i++)
    {
        if (!R_FINITE(REAL(ans)[i]))
            error(_("Missing value or an infinity produced when evaluating the model"));
    }
    for (i = 0; i < LENGTH(theta); i++)
    {
        char *name = CHAR(STRING_ELT(theta, i));
        SEXP temp = findVar(install(name), rho);
        if (isInteger(temp))
            error(_("variable '%s' is integer, not numeric"), name);
        if (!isReal(temp))
            error(_("variable '%s' is not numeric"), name);
        SET_VECTOR_ELT(pars, i, temp);
        lengthTheta += LENGTH(VECTOR_ELT(pars, i));
    }
    PROTECT(gradient = allocMatrix(REALSXP, LENGTH(ans), lengthTheta));

    for (i = 0, start = 0; i < LENGTH(theta); i++)
    {
        for (j = 0; j < LENGTH(VECTOR_ELT(pars, i)); j++, start += LENGTH(ans))
        {
            SEXP ans_del;
            double origPar, xx, delta;

            origPar = REAL(VECTOR_ELT(pars, i))[j];
            xx = fabs(origPar);
            delta = (xx == 0) ? eps : xx * eps;
            REAL(VECTOR_ELT(pars, i))[j] += rDir[i] * delta;
            PROTECT(ans_del = eval(expr, rho));
            if (!isReal(ans_del))
                ans_del = coerceVector(ans_del, REALSXP);
            UNPROTECT(1);
            for (k = 0; k < LENGTH(ans); k++)
            {
                if (!R_FINITE(REAL(ans_del)[k]))
                    error(_("Missing value or an infinity produced when evaluating the model"));
                REAL(gradient)[start + k] = rDir[i] * (REAL(ans_del)[k] - REAL(ans)[k]) / delta;
            }
            REAL(VECTOR_ELT(pars, i))[j] = origPar;
        }
    }
    setAttrib(ans, install("gradient"), gradient);
    UNPROTECT(3);
    return ans;
}
