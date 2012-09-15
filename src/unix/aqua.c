/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999-2012 The R Core Team
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
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Defn.h>
#include <Internal.h>

#include "Runix.h"
#include <sys/types.h>
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif

#ifdef HAVE_AQUA

/* tell QuartzDevice.h to insert definitions for us (to maintain consistency) */
#define IN_AQUA_C 1

#include <R_ext/GraphicsEngine.h>
#include <R_ext/Rdynload.h>
#include <R_ext/QuartzDevice.h>

extern Rboolean useaqua; /* from src/unix/system.c */

/* These are in no header.  Their definitions are in
   Mac-GUI/REngine/Rinit.m, which sets them to functions in
   Mac-GUI/REngine/Rcallbacks.m

   So this is a essentially a private hook arrangement for R.app

   There's another one in src/main/systutils.c, ptr_CocoaSystem .

extern SEXP (*ptr_do_packagemanger)(SEXP, SEXP, SEXP, SEXP);
extern SEXP (*ptr_do_datamanger)(SEXP, SEXP, SEXP, SEXP);
extern SEXP (*ptr_do_browsepkgs)(SEXP, SEXP, SEXP, SEXP);
extern SEXP (*ptr_do_wsbrowser)(SEXP, SEXP, SEXP, SEXP);
extern SEXP (*ptr_do_hsbrowser)(SEXP, SEXP, SEXP, SEXP);
*/

DL_FUNC ptr_do_wsbrowser, ptr_GetQuartzParameters, ptr_do_browsepkgs, ptr_do_datamanger, ptr_do_packagemanger,
    ptr_do_hsbrowser;

int (*ptr_Raqua_CustomPrint)(const char *, SEXP);

/* called from Mac-GUI/RController.m */
static QuartzFunctions_t *qfn; // could this not be internal to function?
QuartzFunctions_t *getQuartzFunctions(void)
{
    if (qfn)
        return qfn;
    {
        QuartzFunctions_t *(*fn)(void);
        fn = (QuartzFunctions_t * (*)(void)) R_FindSymbol("getQuartzAPI", "grDevices", NULL);
        if (!fn)
        {
            /* we need to load grDevices - not sure if this is the best way, though ... */
            SEXP call = lang2(install("loadNamespace"), install("grDevices"));
            PROTECT(call);
            eval(call, R_GlobalEnv);
            UNPROTECT(1);
            fn = (QuartzFunctions_t * (*)(void)) R_FindSymbol("getQuartzAPI", "grDevices", NULL);
            if (!fn)
                error(_("unable to load Quartz"));
        }
        return fn();
    }
}

SEXP do_wsbrowser(SEXP call, SEXP op, SEXP args, SEXP env)
{
    return ptr_do_wsbrowser(call, op, args, env);
}

SEXP do_browsepkgs(SEXP call, SEXP op, SEXP args, SEXP env)
{
    return ptr_do_browsepkgs(call, op, args, env);
}

SEXP do_datamanger(SEXP call, SEXP op, SEXP args, SEXP env)
{
    return ptr_do_datamanger(call, op, args, env);
}

SEXP do_hsbrowser(SEXP call, SEXP op, SEXP args, SEXP env)
{
    return ptr_do_hsbrowser(call, op, args, env);
}

SEXP do_packagemanger(SEXP call, SEXP op, SEXP args, SEXP env)
{
    return ptr_do_packagemanger(call, op, args, env);
}

SEXP do_aqua_custom_print(SEXP call, SEXP op, SEXP args, SEXP env)
{
    const void *vm;
    const char *ct;
    int cpr;
    SEXP rv, objType, obj;

    if (!ptr_Raqua_CustomPrint)
        return R_NilValue;

    checkArity(op, args);

    vm = vmaxget();

    objType = CAR(args);
    args = CDR(args);
    obj = CAR(args);

    if (!isString(objType) || LENGTH(objType) < 1)
        errorcall(call, "invalid arguments");
    ct = CHAR(STRING_ELT(objType, 0));
    cpr = ptr_Raqua_CustomPrint(ct, obj);

    /* FIXME: trying to store a pointer in an integer is wrong */
    PROTECT(rv = allocVector(INTSXP, 1));
    INTEGER(rv)[0] = cpr;

    vmaxset(vm);
    UNPROTECT(1);

    return rv;
}
#endif
