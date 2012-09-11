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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <R.h>
#include <Rinternals.h>

#include "utils.h"
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)                                                                                               \
    {                                                                                                                  \
#name, (DL_FUNC)&name, n                                                                                       \
    }

static const R_CallMethodDef CallEntries[] = {CALLDEF(crc64, 1),
                                              CALLDEF(menu, 1),
                                              CALLDEF(nsl, 1),
                                              CALLDEF(objectSize, 1),

                                              /* Sockets */
                                              CALLDEF(sockconnect, 2),
                                              CALLDEF(sockread, 2),
                                              CALLDEF(sockclose, 1),
                                              CALLDEF(sockopen, 1),
                                              CALLDEF(socklisten, 1),
                                              CALLDEF(sockwrite, 2),

                                              {NULL, NULL, 0}};

#define EXTDEF(name, n)                                                                                                \
    {                                                                                                                  \
#name, (DL_FUNC)&name, n                                                                                       \
    }

static const R_ExternalMethodDef ExtEntries[] = {
    EXTDEF(unzip, 7),       EXTDEF(Rprof, 4),         EXTDEF(Rprofmem, 3),

    EXTDEF(countfields, 6), EXTDEF(readtablehead, 6), EXTDEF(typeconvert, 4), EXTDEF(writetable, 11),

    EXTDEF(addhistory, 1),  EXTDEF(loadhistory, 1),   EXTDEF(savehistory, 1),

    {NULL, NULL, 0}};

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
    __attribute__((visibility("default")))
#endif
    R_init_utils(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, ExtEntries);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
