/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998 ff  Robert Gentleman, Ross Ihaka and the R core team
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
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *
 * Generally useful  UTILITIES  *NOT* relying on R internals (from Defn.h)
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "Complex.h"

/* ../main/sort.c : */
void	isort(int*,     int);
void	rsort(double*, int);
void	csort(complex*, int);
void	revsort(double*, int*, int);/* reverse; sort i[] alongside */
void	iPsort(int*,    int, int);
void	rPsort(double*, int, int);
void	cPsort(complex*, int, int);

int	IndexWidth(int);
int	Rstrlen(char*);
char*	R_ExpandFileName(char*);
void	setIVector(int*, int, int);
void	setRVector(double*, int, double);
int	StringFalse(char*);
int	StringTrue(char*);

void	hsv2rgb(double *h, double *s, double *v,/* in */
		double *r, double *g, double *b);/* out */

#endif
