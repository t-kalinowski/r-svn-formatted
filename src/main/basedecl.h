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
 *  This header contains the declarations of code to be used by
 *  .C, .Fortran, .Call or .External within the base package.
 *  These routines are `registered' in registration.c.
 */

void Rsockconnect(int *, char **);
void Rsockopen(int *);
void Rsocklisten(int *, char **, int *);
void Rsockclose(int *);
void Rsockread(int *, char **, int *);
void Rsockwrite(int *, char **, int *, int *, int *);
SEXP La_svd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP La_rs(SEXP, SEXP);
SEXP La_rg(SEXP, SEXP);
SEXP La_dlange(SEXP, SEXP);
SEXP La_dgecon(SEXP, SEXP);
SEXP La_dtrcon(SEXP, SEXP);
SEXP La_zgecon(SEXP, SEXP);
SEXP La_ztrcon(SEXP, SEXP);
SEXP La_zgesv(SEXP, SEXP);
SEXP La_zgeqp3(SEXP);
SEXP qr_coef_cmplx(SEXP, SEXP);
SEXP qr_qy_cmplx(SEXP, SEXP, SEXP);
SEXP La_svd_cmplx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP La_rs_cmplx(SEXP, SEXP);
SEXP La_rg_cmplx(SEXP, SEXP);
SEXP La_chol(SEXP);
SEXP La_chol2inv(SEXP, SEXP);
SEXP La_dgesv(SEXP, SEXP, SEXP);
SEXP La_dgeqp3(SEXP);
SEXP qr_coef_real(SEXP, SEXP);
SEXP qr_qy_real(SEXP, SEXP, SEXP);
SEXP det_ge_real(SEXP, SEXP);

SEXP R_getTaskCallbackNames(void);
SEXP R_removeTaskCallback(SEXP);
SEXP R_addTaskCallback(SEXP, SEXP, SEXP, SEXP);

SEXP R_traceOnOff(SEXP);
SEXP R_setS4Object(SEXP, SEXP);

SEXP Rrowsum_matrix(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rrowsum_df(SEXP, SEXP, SEXP, SEXP, SEXP);

void F77_SYMBOL(dchdc)(double *, int *, int *, double *, int *, int *, int *);
void F77_SYMBOL(dpbfa)(double *, int *, int *, int *, int *);
void F77_SYMBOL(dpbsl)(double *, int *, int *, int *, double *);

SEXP R_compress1(SEXP);
SEXP R_decompress1(SEXP);

SEXP R_serializeb(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_serialize(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_unserialize(SEXP, SEXP);

//SEXP R_getVarsFromFrame(SEXP, SEXP, SEXP);
//SEXP R_lazyLoadDBinsertValue(SEXP, SEXP, SEXP, SEXP, SEXP);
//SEXP R_lazyLoadDBflush(SEXP);

SEXP R_getbcprofcounts(void);
SEXP R_startbcprof(void);
SEXP R_stopbcprof(void);

SEXP bitwiseNot(SEXP);
SEXP bitwiseAnd(SEXP, SEXP);
SEXP bitwiseOr(SEXP, SEXP);
SEXP bitwiseXor(SEXP, SEXP);
SEXP crc64ToString(SEXP);

SEXP BinCode(SEXP x, SEXP breaks, SEXP right, SEXP lowest);
SEXP R_Tabulate(SEXP in, SEXP nbin);
SEXP FindIntervVec(SEXP xt, SEXP x, SEXP right, SEXP inside);

