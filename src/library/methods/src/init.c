/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2005   The R Development Core Team.
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
#include <Rinternals.h>

#include "methods.h"
#include <R_ext/Rdynload.h>

void methods_init(char **path)
{
#ifdef ENABLE_NLS
    char localedir[PATH_MAX];

    strcpy(localedir, path[0]);
    strcat(localedir, "/po");
    bindtextdomain("methods", localedir);
#endif
}

static const R_CMethodDef CEntries[] = {{"methods_init", (DL_FUNC)&methods_init, 1}, {NULL, NULL, 0}};

#define CALLDEF(name, n)                                                                                               \
    {                                                                                                                  \
#name, (DL_FUNC)&name, n                                                                                       \
    }

static R_CallMethodDef CallEntries[] = {CALLDEF(R_M_setPrimitiveMethods, 5),
                                        CALLDEF(R_clear_method_selection, 0),
                                        CALLDEF(R_dummy_extern_place, 0),
                                        CALLDEF(R_el_named, 2),
                                        CALLDEF(R_externalptr_prototype_object, 0),
                                        CALLDEF(R_getGeneric, 3),
                                        CALLDEF(R_get_slot, 2),
                                        CALLDEF(R_identC, 2),
                                        CALLDEF(R_initMethodDispatch, 1),
                                        CALLDEF(R_methodsPackageMetaName, 2),
                                        CALLDEF(R_methods_test_MAKE_CLASS, 1),
                                        CALLDEF(R_methods_test_NEW, 1),
                                        CALLDEF(R_missingArg, 2),
                                        CALLDEF(R_nextMethodCall, 2),
                                        CALLDEF(R_quick_method_check, 3),
                                        CALLDEF(R_selectMethod, 4),
                                        CALLDEF(R_set_el_named, 3),
                                        CALLDEF(R_set_slot, 3),
                                        CALLDEF(R_standardGeneric, 3),
                                        CALLDEF(do_substitute_direct, 2),
                                        {NULL, NULL, 0}};

void R_init_methods(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
