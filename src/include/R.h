/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2000 The R Development Core Team.
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

#ifndef R_R_H
#define R_R_H

#ifndef USING_R
#define USING_R
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#ifdef Macintosh
#include <fp.h>
#else
#include <math.h>
#endif

#include "Rconfig.h"
#include "Rversion.h"         /* R_VERSION */
#include "R_ext/Arith.h"      /* R_FINITE, ISNAN, ... */
#include "R_ext/Boolean.h"    /* Rboolean type */
#include "R_ext/Complex.h"    /* Rcomplex type */
#include "R_ext/Constants.h"  /* PI, DOUBLE_EPS, etc */
#include "R_ext/Error.h"      /* error and warning */
#include "R_ext/Memory.h"     /* R_alloc and S_alloc */
#include "R_ext/Random.h"     /* RNG interface */
#include "R_ext/Utils.h"      /* sort routines */
#include "R_ext/RS.h"
/* for PROBLEM ... Calloc, Realloc, Free, Memcpy, F77_xxxx */


typedef double Sfloat;
typedef int Sint;
#define SINT_MAX INT_MAX
#define SINT_MIN INT_MIN

#ifdef __cplusplus
}
#endif

#endif /* !R_R_H */
