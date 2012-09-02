/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2012   The R Core Team.
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

SEXP getListElement(SEXP list, char *str);

/* Declarations for .Call entry points */

SEXP Cdqrls(SEXP x, SEXP y, SEXP tol);
SEXP Cdist(SEXP x, SEXP method, SEXP attrs, SEXP p);
SEXP r2dtable(SEXP n, SEXP r, SEXP c);
SEXP cor(SEXP x, SEXP y, SEXP na_method, SEXP method);
SEXP cov(SEXP x, SEXP y, SEXP na_method, SEXP method);
SEXP updateform(SEXP old, SEXP new);
SEXP fft(SEXP z, SEXP inverse);
SEXP mvfft(SEXP z, SEXP inverse);
SEXP nextn(SEXP n, SEXP factors);

SEXP cfilter(SEXP sx, SEXP sfilter, SEXP ssides, SEXP scircular);
SEXP rfilter(SEXP x, SEXP filter, SEXP out);
SEXP lowess(SEXP x, SEXP y, SEXP sf, SEXP siter, SEXP sdelta);
SEXP DoubleCentre(SEXP A);
SEXP BinDist(SEXP x, SEXP weights, SEXP slo, SEXP sup, SEXP sn);

/* Declarations for .External[2] entry points */

SEXP compcases(SEXP args);
SEXP doD(SEXP args);
SEXP deriv(SEXP args);
SEXP modelframe(SEXP call, SEXP op, SEXP args, SEXP rho);
SEXP modelmatrix(SEXP call, SEXP op, SEXP args, SEXP rho);
SEXP termsform(SEXP args);
SEXP do_fmin(SEXP call, SEXP op, SEXP args, SEXP rho);
SEXP nlm(SEXP call, SEXP op, SEXP args, SEXP rho);
SEXP zeroin2(SEXP call, SEXP op, SEXP args, SEXP rho);
SEXP optim(SEXP call, SEXP op, SEXP args, SEXP rho);
SEXP optimhess(SEXP call, SEXP op, SEXP args, SEXP rho);
SEXP Rmultinom(SEXP args);
SEXP call_dqagi(SEXP);
SEXP call_dqags(SEXP);
SEXP Random1(SEXP args);
SEXP Random2(SEXP args);
SEXP Random3(SEXP args);
SEXP distn2(SEXP args);
SEXP distn3(SEXP args);
SEXP distn4(SEXP args);

SEXP Rsm(SEXP x, SEXP stype, SEXP send);
SEXP tukeyline(SEXP x, SEXP y, SEXP call);
SEXP runmed(SEXP x, SEXP stype, SEXP sk, SEXP end, SEXP print_level);
SEXP influence(SEXP mqr, SEXP do_coef, SEXP e, SEXP stol);
