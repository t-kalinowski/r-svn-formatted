/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999, The R Development Core Team.
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
 */

#ifndef R_DEFINES_H
#define R_DEFINES_H

/*
 *  Much is from John Chambers' "Programming With Data".
 *  Some of this is from Doug Bates.
 *
 *  It is presented here to support a joint programming style which
 *  will work in both R and S.  In particular it helps with:
 *
 *    1. Duncan Temple Lang's CORBA code.
 *    2. John Chambers' Java Code.
 *
 *  And to hide some internal nastiness.
 */

#define NULL_USER_OBJECT	R_NilValue

#define AS_LOGICAL(x)		coerceVector(x,LGLSXP)
#define AS_INTEGER(x)		coerceVector(x,INTSXP)
#define AS_NUMERIC(x)		coerceVector(x,REALSXP)
#define AS_CHARACTER(x)		coerceVector(x,STRSXP)
#define AS_COMPLEX(x)		coerceVector(x,CPLXSXP)
#define AS_VECTOR(x)		coerceVector(x,VECSXP)

#define NEW_LIST(n)		allocVector(VECSXP,n)
#define NEW_LOGICAL(n)		allocVector(LGLSXP,n)
#define NEW_INTEGER(n)		allocVector(INTSXP,n)
#define NEW_NUMERIC(n)		allocVector(REALSXP,n)
#define NEW_COMPLEX(n)		allocVector(CPLXSXP,n)
#define NEW_CHARACTER(n)	allocVector(STRSXP,n)
#define NEW_STRING(n)		NEW_CHARACTER(n)

#define GET_LENGTH(x)		length(x)

#define LOGICAL_DATA(x)		LOGICAL(x)
#define INTEGER_DATA(x)		INTEGER(x)
#define NUMERIC_DATA(x)		REAL(x)
#define COMPLEX_DATA(x)		COMPLEX(x)
#define STRING_DATA(x)		STRING(x)

#define R_PROBLEM_BUFSIZE	4096
#define PROBLEM			{char R_problem_buf[R_PROBLEM_BUFSIZE];sprintf(R_problem_buf,
#define ERROR			),error(R_problem_buf);}
#define RECOVER(x)		),error(R_problem_buf);}
#define WARNING(x)		),warning(R_problem_buf);}
#define LOCAL_EVALUATOR		/**/
#define NULL_ENTRY		/**/
#define WARN			WARNING(NULL)

#ifdef NEW_GC
#define COPY_TO_USER_STRING(x)	mkStringElement(x)
#define CREATE_STRING_VECTOR(x)	mkStringElement(x)
#else
#define COPY_TO_USER_STRING(x)	mkChar(x)
#define CREATE_STRING_VECTOR(x)	mkChar(x)
#endif

#define RECURSIVE_DATA(x)	VECTOR(x)
#define CHARACTER_DATA(x)	STRING(x)

#define CREATE_FUNCTION_CALL(name, argList) createFunctionCall(name, argList)

#define EVAL(x)			eval(x,R_GlobalEnv)

/* S Like Memory Management */

#define Calloc(n, t)   (t *) R_chk_calloc( (size_t) (n), sizeof(t) )
#define Realloc(p,n,t) (t *) R_chk_realloc( (void *)(p), (size_t)((n) * sizeof(t)) )
#define Free(p)        R_chk_free( (void *)(p) )
#define Memcpy(p,q,n)  memcpy( p, q, (size_t)( (n) * sizeof(*p) ) )

/* S Like Fortran Interface */

#define F77_CALL(x)    F77_SYMBOL(x)
#define F77_NAME(x)    F77_SYMBOL(x)

#endif
