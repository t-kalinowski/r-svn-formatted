/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2003   The R Development Core Team.
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

#include <R.h>
#include "tools.h"
#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
    {"delim_match", (DL_FUNC)&delim_match, 2}, {"Rmd5", (DL_FUNC)&Rmd5, 1}, {NULL, NULL, 0}};

void R_init_tools(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
}
