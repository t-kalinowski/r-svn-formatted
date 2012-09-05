/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001-12   The R Core Team.
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

#ifndef R_CTEST_H
#define R_CTEST_H

#include <R.h>

void chisqsim(int *nrow, int *ncol, int *nrowt, int *ncolt, int *n,
	      int *b, double *expected, int *observed, double *fact,
	      int *jwork, double *results);
void fisher_sim(int *nrow, int *ncol, int *nrowt, int *ncolt, int *n,
		int *b, int *observed, double *fact,
		int *jwork, double *results);
void d2x2xk(int *k, double *m, double *n, double *t, double *d);
void fexact(int *nrow, int *ncol, int *table, int *ldtabl,
	    double *expect, double *percnt, double *emin, double *prt,
	    double *pre, int *workspace, int *mult);
void swilk(double *x, int *n, int *n1, double *w, double *pw, int *ifault);

void rcont2(int *nrow, int *ncol, int *nrowt, int *ncolt, int *ntotal,
	    double *fact, int *jwork, int *matrix);
#endif
