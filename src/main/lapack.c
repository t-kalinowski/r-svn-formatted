/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001 The R Development Core Team
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

#include <Defn.h>
#include "R_ext/Rdynpriv.h"

typedef SEXP (*sDL_FUNC)();
static sDL_FUNC ptr_svd, ptr_rs, ptr_rg;

/*
SEXP La_svd(SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, SEXP v)
SEXP La_rs(SEXP x, SEXP only_values)
SEXP La_rg(SEXP x, SEXP only_values)
*/

static int initialized = 0;

static void La_Init(void)
{
    int res = moduleCdynload("lapack", 1, 1);
    initialized = -1;
    if (!res)
        return;

    if (!(ptr_svd = (sDL_FUNC)R_FindSymbol("La_svd", "lapack")))
        return;
    if (!(ptr_rs = (sDL_FUNC)R_FindSymbol("La_rs", "lapack")))
        return;
    if (!(ptr_rg = (sDL_FUNC)R_FindSymbol("La_rg", "lapack")))
        return;
    initialized = 1;
    return;
}

SEXP La_svd(SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, SEXP v)
{
    if (!initialized)
        La_Init();
    if (initialized > 0)
        return (*ptr_svd)(jobu, jobv, x, s, u, v);
    else
    {
        error("lapack routines cannot be loaded");
        return R_NilValue;
    }
}

SEXP La_rs(SEXP x, SEXP only_values)
{
    if (!initialized)
        La_Init();
    if (initialized > 0)
        return (*ptr_rs)(x, only_values);
    else
    {
        error("lapack routines cannot be loaded");
        return R_NilValue;
    }
}

SEXP La_rg(SEXP x, SEXP only_values)
{
    if (!initialized)
        La_Init();
    if (initialized > 0)
        return (*ptr_rg)(x, only_values);
    else
    {
        error("lapack routines cannot be loaded");
        return R_NilValue;
    }
}
