/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2003    Robert Gentleman, Ross Ihaka
 *                             and the R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 * Generally useful  UTILITIES  *NOT* relying on R internals (from Defn.h)
 */

#ifndef R_EXT_UTILS_H_
#define R_EXT_UTILS_H_

#include <R_ext/Boolean.h>
#include <R_ext/Complex.h>

#define revsort       Rf_revsort
#define iPsort        Rf_iPsort
#define rPsort        Rf_rPsort
#define cPsort        Rf_cPsort
#define IndexWidth    Rf_IndexWidth
#define setIVector    Rf_setIVector
#define setRVector    Rf_setRVector
#define StringFalse   Rf_StringFalse
#define StringTrue    Rf_StringTrue
#define isBlankString Rf_isBlankString
#define hsv2rgb       Rf_hsv2rgb

#ifdef  __cplusplus
extern "C" {
#endif

/* ../../main/sort.c : */
void	R_isort(int*, int);
void	R_rsort(double*, int);
void	R_csort(Rcomplex*, int);
void    rsort_with_index(double *, int *, int);
void	revsort(double*, int*, int);/* reverse; sort i[] alongside */
void	iPsort(int*,    int, int);
void	rPsort(double*, int, int);
void	cPsort(Rcomplex*, int, int);

/* ../../main/qsort.c : */
void R_qsort    (double *v,         int i, int j);
void R_qsort_I  (double *v, int *I, int i, int j);
void R_qsort_int  (int *iv,         int i, int j);
void R_qsort_int_I(int *iv, int *I, int i, int j);


/* ../../main/printutils.c : */
int	IndexWidth(int);
int	Rstrlen(char*, int);
/* ../../main/util.c  and others : */
char*	R_ExpandFileName(char*);
void	setIVector(int*, int, int);
void	setRVector(double*, int, double);
Rboolean StringFalse(char*);
Rboolean StringTrue(char*);
Rboolean isBlankString(char *);

double 	R_strtod(const char *c, char **end);
char *  R_tmpnam(const char *prefix, const char *tempdir);

void	hsv2rgb(double *h, double *s, double *v,/* in */
		double *r, double *g, double *b);/* out */

void R_CheckUserInterrupt(void);

#ifdef  __cplusplus
}
#endif

#endif /* R_EXT_UTILS_H_ */
