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

/* This file provides support for R.app, the OS X front end */

/*
   There's another one in src/main/systutils.c, ptr_CocoaSystem .
*/

#include <Defn.h>
/*#include <Internal.h>

#include "Runix.h"
#include <sys/types.h>
#ifdef HAVE_SYS_STAT_H
# include <sys/stat.h>
#endif
*/

/* tell QuartzDevice.h to insert definitions for us (to maintain consistency) */
#define IN_AQUA_C 1

#include <R_ext/GraphicsEngine.h>
#include <R_ext/Rdynload.h>
#include <R_ext/QuartzDevice.h>

DL_FUNC ptr_GetQuartzParameters;

/* called from Mac-GUI/RController.m, before packages are loaded.
   If this fails, it hangs R.app */

/* FIXME: this should not be allowed: we were requiring symbols in
   grDevices.dll */
QuartzFunctions_t *getQuartzFunctions(void)
{
    /* presumably this was intended to cache the result.
       But,
       - it was never used
       - it would be unsafe as the namespace could be unloaded/reloaded
    static QuartzFunctions_t* qfn;
    if (qfn) return qfn;
    */

    QuartzFunctions_t *(*fn)(void);
    fn = (QuartzFunctions_t * (*)(void)) R_FindSymbol("getQuartzAPI", "grDevices", NULL);
    if (!fn)
    {
        SEXP call = lang2(install("loadNamespace"), mkString("grDevices"));
        PROTECT(call);
        eval(call, R_GlobalEnv);
        UNPROTECT(1);
        fn = (QuartzFunctions_t * (*)(void)) R_FindSymbol("getQuartzAPI", "grDevices", NULL);
        if (!fn)
            error("unable to get QuartzAPI");
    }
    return fn();
}
