/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001-3   The R Development Core Team.
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

#include <R.h>
#include "ctest.h"
#include "modreg.h"
#include <R_ext/Rdynload.h>

#include <Rinternals.h>

static R_NativePrimitiveArgType chisqsim_t[11] = {INTSXP,  INTSXP, INTSXP,  INTSXP, INTSXP, INTSXP,
                                                  REALSXP, INTSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType d2_t[5] = {INTSXP, REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType dansari_t[4] = {INTSXP, REALSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType pansari_t[4] = {INTSXP, REALSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType qansari_t[4] = {INTSXP, REALSXP, INTSXP, INTSXP};

static R_NativePrimitiveArgType fexact_t[10] = {INTSXP,  INTSXP,  REALSXP, INTSXP,  REALSXP,
                                                REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};

static R_NativePrimitiveArgType pkendall_t[3] = {INTSXP, REALSXP, INTSXP};

static R_NativePrimitiveArgType pkstwo_t[3] = {INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType prho_t[5] = {INTSXP, REALSXP, REALSXP, INTSXP, LGLSXP};
static R_NativePrimitiveArgType psmirnov2x_t[3] = {REALSXP, INTSXP, INTSXP};

static R_NativePrimitiveArgType swilk_t[9] = {LGLSXP,    SINGLESXP, INTSXP,  INTSXP, INTSXP,
                                              SINGLESXP, REALSXP,   REALSXP, INTSXP};

static R_NativePrimitiveArgType Srunmed_t[6] = {REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, LGLSXP};
static R_NativePrimitiveArgType Trunmed_t[9] = {INTSXP, INTSXP,  REALSXP, REALSXP, INTSXP,
                                                INTSXP, REALSXP, INTSXP,  INTSXP};

static const R_CMethodDef cMethods[] = {{"chisqsim", (DL_FUNC)&chisqsim, 11, chisqsim_t},
                                        {"d2x2xk", (DL_FUNC)&d2x2xk, 5, d2_t},
                                        {"dansari", (DL_FUNC)&dansari, 4, dansari_t},
                                        {"fexact", (DL_FUNC)&fexact, 10, fexact_t},
                                        {"pansari", (DL_FUNC)&pansari, 4, pansari_t},
                                        {"pkendall", (DL_FUNC)&pkendall, 3, pkendall_t},
                                        {"pkstwo", (DL_FUNC)&pkstwo, 3, pkstwo_t},
                                        {"prho", (DL_FUNC)&prho, 5, prho_t},
                                        {"psmirnov2x", (DL_FUNC)&psmirnov2x, 3, psmirnov2x_t},
                                        {"qansari", (DL_FUNC)&qansari, 4, qansari_t},
                                        {"swilk", (DL_FUNC)&swilk, 9, swilk_t},
                                        {"BDRksmooth", (DL_FUNC)&BDRksmooth, 8},
                                        {"loess_raw", (DL_FUNC)&loess_raw, 24},
                                        {"loess_dfit", (DL_FUNC)&loess_dfit, 13},
                                        {"loess_dfitse", (DL_FUNC)&loess_dfitse, 16},
                                        {"loess_ifit", (DL_FUNC)&loess_ifit, 8},
                                        {"loess_ise", (DL_FUNC)&loess_ise, 15},
                                        {"Srunmed", (DL_FUNC)&Srunmed, 6, Srunmed_t},
                                        {"Trunmed", (DL_FUNC)&Trunmed, 9, Trunmed_t},
                                        {NULL, NULL, 0}};

static R_CallMethodDef CallEntries[] = {{"R_isoreg", (DL_FUNC)&R_isoreg, 1}, {NULL, NULL, 0}};

static R_FortranMethodDef FortEntries[] = {{"lowesw", (DL_FUNC)&F77_SUB(lowesw), 4},
                                           {"lowesp", (DL_FUNC)&F77_SUB(lowesp), 7},
                                           {"setppr", (DL_FUNC)&F77_SUB(setppr), 6},
                                           {"smart", (DL_FUNC)&F77_SUB(smart), 16},
                                           {"pppred", (DL_FUNC)&F77_SUB(pppred), 5},
                                           {"qsbart", (DL_FUNC)&F77_SUB(qsbart), 21},
                                           {"bvalus", (DL_FUNC)&F77_SUB(bvalus), 7},
                                           {"supsmu", (DL_FUNC)&F77_SUB(supsmu), 10},
                                           {NULL, NULL, 0}};

void R_init_stats(DllInfo *dll)
{
    R_useDynamicSymbols(dll, FALSE);
    R_registerRoutines(dll, cMethods, CallEntries, FortEntries, NULL);
}
