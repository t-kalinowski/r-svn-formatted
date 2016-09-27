/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2016  The R Core Team
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
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *
 * Memory Allocation (garbage collected) --- INCLUDING S compatibility ---
 */

/* Included by R.h: API */

#ifndef R_EXT_MEMORY_H_
#define R_EXT_MEMORY_H_

#if defined(__cplusplus) && !defined(DO_NOT_USE_CXX_HEADERS)
# include <cstddef>
# define R_SIZE_T std::size_t
#else
# include <stddef.h> /* for size_t */
# define R_SIZE_T size_t
#endif

#ifdef  __cplusplus
extern "C" {
#endif

void*	vmaxget(void);
void	vmaxset(const void *);

void	R_gc(void);
int	R_gc_running();

char*	R_alloc(R_SIZE_T, int);
long double *R_allocLD(R_SIZE_T nelem);
char*	S_alloc(long, int);
char*	S_realloc(char *, long, long, int);

#ifdef  __cplusplus
}
#endif

#endif /* R_EXT_MEMORY_H_ */
