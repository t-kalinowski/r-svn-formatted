/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997-2003   The R Development Core Team.
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

#ifndef PRINT_H_
#define PRINT_H_

#include <R_ext/PrtUtil.h>
#include "Defn.h"

#define formatFactor        Rf_formatFactor
#define formatRaw           Rf_formatRaw
#define formatString        Rf_formatString
#define EncodeFactor        Rf_EncodeFactor
#define EncodeElement       Rf_EncodeElement
#define printArray          Rf_printArray
#define printMatrix         Rf_printMatrix
#define printNamedVector    Rf_printNamedVector
#define printVector         Rf_printVector

typedef struct {
    int width;
    int na_width;
    int na_width_noquote;
    int digits;
    int scipen;
    int gap;
    int quote;
    int right;
    SEXP na_string;
    SEXP na_string_noquote;
} R_print_par_t;
extern R_print_par_t R_print;

/* Computation of printing formats */
void formatFactor(int *, int, int *, SEXP, int);
void formatRaw(Rbyte *, int, int *);
void formatString(SEXP*, int, int*, int);

/* Formating of values */
char *EncodeFactor(int, int, int, SEXP);
char *EncodeElement(SEXP, int, int);

/* Printing */
void MatrixColumnLabel(SEXP, int, int);
void RightMatrixColumnLabel(SEXP, int, int);
void LeftMatrixColumnLabel(SEXP, int, int);
void MatrixRowLabel(SEXP, int, int, int);

/* In Rinternals.h (and MUST be there):
   CustomPrintValue,  PrintValue, PrintValueRec */
void printArray(SEXP, SEXP, int, int, SEXP);
void printMatrix(SEXP, int, SEXP, int, int, SEXP, SEXP, char*, char*);
void printNamedVector(SEXP, SEXP, int, char*);
void printVector(SEXP, int, int);

/* Utilities for S compatibility and debuggging */
int F77_SYMBOL(dblepr0)(char *, int *, double *, int *);
int F77_SYMBOL(intpr0) (char *, int *, int *, int *);
int F77_SYMBOL(realpr0)(char *, int *, float *, int *);
void R_PV(SEXP s);

/* Offset for rowlabels if there are named dimnames */
#define R_MIN_LBLOFF 2

#define R_MIN_WIDTH_OPT		10
#define R_MAX_WIDTH_OPT		10000
#define R_MIN_DIGITS_OPT	1
#define R_MAX_DIGITS_OPT	22

#endif
