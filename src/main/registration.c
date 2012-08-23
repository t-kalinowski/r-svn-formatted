/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2002-2012	The R Core Team.
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
 *
 */
/*
 *  This file replaces the previously used ROUTINES file and is used to
 *  explicitly register native routines that are located in the R
 *  executable (e.g. R.bin, Rgui.exe) but which are intended to be
 *  accessible to S code via .C(), .Fortran(), .Call(), .External().
 *  The approach we use here is the regular registration mechanism that
 *  packages can use to explicitly list the symbols to be exported.
 *  For .C() and .Call() routines, we give the number of arguments
 *  expected.
 *  For .C() routines, we also specify the types of the arguments.
 *  For .Fortran() and .External() routines, we specify only the name
 *  and symbol.

 *  To add an entry, first determine by which interface the routine will
 *  be accessed:
 *   .C, .Call, .External or .Fortran
 *  Then add an entry to
 *    cMethods, callMethods, externalMethods, or fortranMethods
 *  respectively
 *
 *  DTL 14-Dec-2002
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Defn.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>
#include <R_ext/Linpack.h>

/*  These get the declarations of some routines refernced here but
    not explicitly declared.    This is necessary when we link with
    a C++ compiler because the linkage changes as the declarations
    are (currently) within extern "C" blocks.
*/
#include <R_ext/Callbacks.h>
#include <Rdynpriv.h>

#include "basedecl.h"

/* FIXME: bincode is no longer used in R, but is still used by
   packages mixOmics spam .  Remove after R 2.15.2.
*/
void bincode(double *x, int *n, double *breaks, int *nb, int *code, int *right, int *include_border, int *naok);

static R_NativePrimitiveArgType bincode_t[] = {REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, LGLSXP, LGLSXP, LGLSXP};

static R_NativePrimitiveArgType R_pretty_t[] = {REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP};

static R_NativePrimitiveArgType Rsockconnect_t[] = {INTSXP, STRSXP};
static R_NativePrimitiveArgType Rsockopen_t[] = {INTSXP};
static R_NativePrimitiveArgType Rsocklisten_t[] = {INTSXP, STRSXP, INTSXP};
static R_NativePrimitiveArgType Rsockclose_t[] = {INTSXP};
static R_NativePrimitiveArgType Rsockread_t[] = {INTSXP, STRSXP, INTSXP};
static R_NativePrimitiveArgType Rsockwrite_t[] = {INTSXP, STRSXP, INTSXP, INTSXP, INTSXP};

#define CDEF(name)                                                                                                     \
    {                                                                                                                  \
#name, (DL_FUNC)&name, sizeof(name##_t) / sizeof(name##_t[0]), name##_t                                        \
    }

static R_CMethodDef cMethods[] = {CDEF(bincode), // remove after R 2.15.2
                                  CDEF(R_pretty),
                                  {"str_signif", (DL_FUNC)&str_signif, 8, NULL}, // mutable first arg

                                  /* Sockets */
                                  CDEF(Rsockconnect),
                                  CDEF(Rsockopen),
                                  CDEF(Rsocklisten),
                                  CDEF(Rsockclose),
                                  CDEF(Rsockread),
                                  CDEF(Rsockwrite),

                                  {NULL, NULL, 0}};

#define CALLDEF(name, n)                                                                                               \
    {                                                                                                                  \
#name, (DL_FUNC)&name, n                                                                                       \
    }

static R_CallMethodDef callMethods[] = {
    /* lapack */
    CALLDEF(La_svd, 7),
    CALLDEF(La_rs, 2),
    CALLDEF(La_rg, 2),
    CALLDEF(La_dlange, 2),
    CALLDEF(La_dgecon, 2),
    CALLDEF(La_dtrcon, 2),
    CALLDEF(La_zgecon, 2),
    CALLDEF(La_ztrcon, 2),
    CALLDEF(La_zgesv, 2),
    CALLDEF(La_zgeqp3, 1),
    CALLDEF(qr_coef_cmplx, 2),
    CALLDEF(qr_qy_cmplx, 3),
    CALLDEF(La_svd_cmplx, 6),
    CALLDEF(La_rs_cmplx, 2),
    CALLDEF(La_rg_cmplx, 2),
    CALLDEF(La_chol2inv, 2),
    CALLDEF(La_chol, 1),
    CALLDEF(La_dgesv, 3),
    CALLDEF(La_dgeqp3, 1),
    CALLDEF(qr_coef_real, 2),
    CALLDEF(qr_qy_real, 3),
    CALLDEF(det_ge_real, 2),

    /* In ../main/unique.c to use hashing. */
    CALLDEF(Rrowsum_matrix, 5),
    CALLDEF(Rrowsum_df, 5),

    /* Top-level task callbacks */
    CALLDEF(R_getTaskCallbackNames, 0),
    CALLDEF(R_removeTaskCallback, 1),
    CALLDEF(R_addTaskCallback, 4),

    /* Methods related routines. */
    CALLDEF(R_isMethodsDispatchOn, 1),
    CALLDEF(R_traceOnOff, 1),
    CALLDEF(R_isS4Object, 1),
    CALLDEF(R_setS4Object, 3),
    CALLDEF(R_do_new_object, 1),

    /* compression and serialization routines */
    CALLDEF(R_compress1, 1),
    CALLDEF(R_decompress1, 1),
    CALLDEF(R_serializeb, 5),
    CALLDEF(R_serialize, 5),
    CALLDEF(R_unserialize, 2),

    /* lazy loading support */
    CALLDEF(R_getVarsFromFrame, 3),
    CALLDEF(R_lazyLoadDBinsertValue, 5),
    CALLDEF(R_lazyLoadDBflush, 1),

    CALLDEF(R_getbcprofcounts, 0),
    CALLDEF(R_startbcprof, 0),
    CALLDEF(R_stopbcprof, 0),

    CALLDEF(bitwiseNot, 1),
    CALLDEF(bitwiseAnd, 2),
    CALLDEF(bitwiseOr, 2),
    CALLDEF(bitwiseXor, 2),

    CALLDEF(crc64ToString, 1),
    CALLDEF(BinCode, 4),
    CALLDEF(R_Tabulate, 2),
    CALLDEF(FindIntervVec, 4),

    {NULL, NULL, 0}};

#define FDEF(name)                                                                                                     \
    {                                                                                                                  \
#name, (DL_FUNC)&F77_SYMBOL(name), -1, NULL                                                                    \
    }
static R_FortranMethodDef fortranMethods[] = {FDEF(ch2inv),
                                              FDEF(chol),
                                              FDEF(cg),
                                              FDEF(ch),
                                              FDEF(rg),
                                              FDEF(rs),
                                              /* Linpack */
                                              FDEF(dchdc),
                                              //    FDEF(dpbfa),
                                              //    FDEF(dpbsl),
                                              //    FDEF(dpoco),
                                              //    FDEF(dpodi),
                                              //    FDEF(dpofa),
                                              //    FDEF(dposl),
                                              FDEF(dqrcf),
                                              //    FDEF(dqrdc),
                                              FDEF(dqrdc2),
                                              FDEF(dqrls), // for historical reasons
                                              FDEF(dqrqty),
                                              FDEF(dqrqy),
                                              FDEF(dqrrsd),
                                              FDEF(dqrsl),
                                              FDEF(dqrxb),
                                              FDEF(dsvdc),
                                              //    FDEF(dtrsl),
                                              FDEF(dtrco),
                                              {NULL, NULL, 0}};

void attribute_hidden R_init_base(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, callMethods, fortranMethods, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
