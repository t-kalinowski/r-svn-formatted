/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995-1996 Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997-2002 The R Development Core Team
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

/*
  This is an effort to merge the 3 different dynload.c files in the
  distribution from the unix/, macintosh/dll/ and gnuwin32/ directories.
  The aim is to consolidate these different implementations into
      i) a generic or platform-independent common core
     ii) platform-dependent routines that are registered
         as function pointers.
  The reason for using function pointers rather than explicit
  linking of symbols is
     a) to avoid confusion in the linking
     b) to allow for easily overriding these in embedded applications
        in which a host application needs to control how R finds
        symbols. This may be necessary for security reasons.
 */

/*  Dynamic Loading Support
 *
 *  This module provides support for run-time loading of shared libraries
 *  access to symbols within such libraries via .C and .Fortran.  This is
 *  done under Unix with dlopen, dlclose and dlsym (the exception is
 *  hpux, where we use compatibility code provided by Luke Tierney.
 *  There are two cases:
 *
 *
 *  1. The dlopen interface is available.
 *
 *  In this case all symbol location in packages is done using the dlopen routines.
 *  We maintain a list of currently loaded shared libraries in an array
 *  called "LoadedDLL" with the number of currenly loaded libraries
 *  being "CountDLL".  To locate a symbol, we probe the loaded libraries
 *  in order until the symbol is located.  If we do not find a symbol
 *  in the loaded libraries, we search the executable itself.  This
 *  search is not very efficient, but this probably pales into
 *  insignificance when compared with the inefficiencies in the R
 *  interpreter.
 *
 *  Loading and unloading of shared libraries is done via the routines
 *  AddDLL and DeleteDLL.  These routines maintain the list of currently
 *  loaded libraries.  When a library is added, any existing reference
 *  to that library are deleted and then the library is inserted at the
 *  start of the search list.  This way, symbols in more recently loaded
 *  libraries are found first.
 *
 *
 *  Accessing native routines in base (the R executable).
 *
 *  In this case, we use the registration mechanism and the DllInfo array
 *  in ../main/Rdynload.c to locate functions in the executable. We do this
 *  by straight linear search through the table.
 *  Note that the base routines registered are listed in
 *               ../main/registration.c
 *  and are registered during the initialization of the R engine.
 *  (This replaces the previous mechanism that built a table from ../appl/ROUTINES
 *  using Perl/sed).
 *
 *
 *  If speed is ever an issue in the lookup of registered symbols, we can
 *  store the registered routines in a hashtable or binary tree as they
 *  are being registered.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* Need "" rather than <> for the Macintosh, apparently. */
#include "Defn.h"
#include <Rmath.h>

#include <Rdynpriv.h>

#ifdef Unix
/* HP-UX 11.0 has dlfcn.h, but according to libtool as of Dec 2001
   this support is broken. So we force use of shlib even when dlfcn.h
   is available */
#ifdef __hpux
#ifdef HAVE_DL_H
#define HAVE_DYNAMIC_LOADING
#endif
#else
#ifdef HAVE_DLFCN_H
#define HAVE_DYNAMIC_LOADING
#endif
#endif /* __hpux */
#ifndef HAVE_NO_SYMBOL_UNDERSCORE
#ifdef HAVE_ELF_H
#define HAVE_NO_SYMBOL_UNDERSCORE
#endif /* HAVE_ELF_H */
#endif /* HAVE_NO_SYMBOL_UNDERSCORE */
#endif

#ifdef Macintosh
#define HAVE_NO_SYMBOL_UNDERSCORE
#define HAVE_DYNAMIC_LOADING
#endif

#ifdef Win32
#define HAVE_DYNAMIC_LOADING
#endif

/* The following code loads in a compatibility module written by Luke
   Tierney to support S version 4 on Hewlett-Packard machines.	The
   relevant defines are set up by autoconf. */

#ifdef HAVE_DYNAMIC_LOADING

#ifdef CACHE_DLL_SYM
/* keep a record of symbols that have been found */
R_CPFun CPFun[100];
int nCPFun = 0;
#endif

#define MAX_NUM_DLLS 100

static int CountDLL = 0;

#include <R_ext/Rdynload.h>

static DllInfo LoadedDLL[MAX_NUM_DLLS];

int addDLL(char *dpath, char *name, HINSTANCE handle);

OSDynSymbol Rf_osDynSymbol;
OSDynSymbol *R_osDynSymbol = &Rf_osDynSymbol;

void R_init_base(DllInfo *); /* In Registration.c */

void InitDynload()
{
    DllInfo *dll;
    int which = addDLL(strdup("base"), "base", NULL);
    dll = &LoadedDLL[which];
    R_init_base(dll);
    InitFunctionHashing();
}

DllInfo *getBaseDllInfo()
{
    return (&LoadedDLL[0]);
}

Rboolean R_useDynamicSymbols(DllInfo *info, Rboolean value)
{
    Rboolean old;
    old = info->useDynamicLookup;
    info->useDynamicLookup = value;

    return (old);
}

void R_addCRoutine(DllInfo *info, const R_CMethodDef *const croutine, Rf_DotCSymbol *sym);
void R_addCallRoutine(DllInfo *info, const R_CallMethodDef *const croutine, Rf_DotCallSymbol *sym);
void R_addFortranRoutine(DllInfo *info, const R_FortranMethodDef *const croutine, Rf_DotFortranSymbol *sym);
void R_addExternalRoutine(DllInfo *info, const R_ExternalMethodDef *const croutine, Rf_DotExternalSymbol *sym);

/*
 Returns a reference to the DllInfo object associated with the dynamic library
 with the path name `path'. This ensures uniqueness rather than having the
 undesirable situation of two libraries with the same name but in different
 directories.
 This is available so that it can be called from arbitrary C routines
 that need to call R_registerRoutines(). The initialization routine
 R_init_<library name> is passed the DllInfo reference as an argument.
 Other routines must explicitly request it using this routine.
 */
DllInfo *R_getDllInfo(const char *path)
{
    int i;
    for (i = 0; i < CountDLL; i++)
    {
        if (strcmp(LoadedDLL[i].path, path) == 0)
            return (&LoadedDLL[i]);
    }
    return ((DllInfo *)NULL);
}

/*
  Explicitly register the native routines for use in .Call(), .C() and .Fortran()
  functions. These registered values are used to resolve symbols in a library
  that makes a call to this routine, rather than the usual dynamic resolution
  done by dlsym() or the equivalent on the different platforms.
 */
int R_registerRoutines(DllInfo *info, const R_CMethodDef *const croutines, const R_CallMethodDef *const callRoutines,
                       const R_FortranMethodDef *const fortranRoutines,
                       const R_ExternalMethodDef *const externalRoutines)
{
    int i, num;

    if (info == NULL)
        error("R_RegisterRoutines called with invalid DllInfo object.");

    info->useDynamicLookup = TRUE; /* Default is to look in registered and then dynamic.
                                      Potentially change in the future to be only registered
                                      if there are any registered values.
                                    */

    if (croutines)
    {
        for (num = 0; croutines[num].name != NULL; num++)
        {
            ;
        }
        info->CSymbols = (Rf_DotCSymbol *)calloc(num, sizeof(Rf_DotCSymbol));
        info->numCSymbols = num;
        for (i = 0; i < num; i++)
        {
            R_addCRoutine(info, croutines + i, info->CSymbols + i);
        }
    }

    if (fortranRoutines)
    {
        for (num = 0; fortranRoutines[num].name != NULL; num++)
        {
            ;
        }
        info->FortranSymbols = (Rf_DotFortranSymbol *)calloc(num, sizeof(Rf_DotFortranSymbol));
        info->numFortranSymbols = num;

        for (i = 0; i < num; i++)
        {
            R_addFortranRoutine(info, fortranRoutines + i, info->FortranSymbols + i);
        }
    }

    if (callRoutines)
    {
        for (num = 0; callRoutines[num].name != NULL; num++)
        {
            ;
        }
        info->CallSymbols = (Rf_DotCallSymbol *)calloc(num, sizeof(Rf_DotCallSymbol));
        info->numCallSymbols = num;
        for (i = 0; i < num; i++)
        {
            R_addCallRoutine(info, callRoutines + i, info->CallSymbols + i);
        }
    }

    if (externalRoutines)
    {
        for (num = 0; externalRoutines[num].name != NULL; num++)
        {
            ;
        }
        info->ExternalSymbols = (Rf_DotExternalSymbol *)calloc(num, sizeof(Rf_DotExternalSymbol));
        info->numExternalSymbols = num;

        for (i = 0; i < num; i++)
        {
            R_addExternalRoutine(info, externalRoutines + i, info->ExternalSymbols + i);
        }
    }

    return (1);
}

static void R_setPrimitiveArgTypes(const R_FortranMethodDef *const croutine, Rf_DotFortranSymbol *sym)
{
    sym->types = (R_NativePrimitiveArgType *)malloc(sizeof(R_NativePrimitiveArgType) * croutine->numArgs);
    if (sym->types)
        memcpy(sym->types, croutine->types, sizeof(R_NativePrimitiveArgType) * croutine->numArgs);
}

static void R_setArgStyles(const R_FortranMethodDef *const croutine, Rf_DotFortranSymbol *sym)
{
    sym->styles = (R_NativeArgStyle *)malloc(sizeof(R_NativeArgStyle) * croutine->numArgs);
    if (sym->styles)
        memcpy(sym->styles, croutine->styles, sizeof(R_NativeArgStyle) * croutine->numArgs);
}

void R_addFortranRoutine(DllInfo *info, const R_FortranMethodDef *const croutine, Rf_DotFortranSymbol *sym)
{
    sym->name = strdup(croutine->name);
    sym->fun = croutine->fun;
    sym->numArgs = croutine->numArgs > -1 ? croutine->numArgs : -1;
    if (croutine->types)
        R_setPrimitiveArgTypes(croutine, sym);
    if (croutine->styles)
        R_setArgStyles(croutine, sym);
}

void R_addExternalRoutine(DllInfo *info, const R_ExternalMethodDef *const croutine, Rf_DotExternalSymbol *sym)
{
    sym->name = strdup(croutine->name);
    sym->fun = croutine->fun;
    sym->numArgs = croutine->numArgs > -1 ? croutine->numArgs : -1;
}

void R_addCRoutine(DllInfo *info, const R_CMethodDef *const croutine, Rf_DotCSymbol *sym)
{
    sym->name = strdup(croutine->name);
    sym->fun = croutine->fun;
    sym->numArgs = croutine->numArgs > -1 ? croutine->numArgs : -1;
    if (croutine->types)
        R_setPrimitiveArgTypes(croutine, sym);
    if (croutine->styles)
        R_setArgStyles(croutine, sym);
}

void R_addCallRoutine(DllInfo *info, const R_CallMethodDef *const croutine, Rf_DotCallSymbol *sym)
{
    sym->name = strdup(croutine->name);
    sym->fun = croutine->fun;
    sym->numArgs = croutine->numArgs > -1 ? croutine->numArgs : -1;
}

void Rf_freeCSymbol(Rf_DotCSymbol *sym)
{
    free(sym->name);
}

void Rf_freeCallSymbol(Rf_DotCallSymbol *sym)
{
    free(sym->name);
}

void Rf_freeFortranSymbol(Rf_DotFortranSymbol *sym)
{
    free(sym->name);
}

void Rf_freeDllInfo(DllInfo *info)
{
    int i;
    free(info->name);
    free(info->path);
    if (info->CSymbols)
    {
        for (i = 0; i < info->numCSymbols; i++)
            Rf_freeCSymbol(info->CSymbols + i);
        free(info->CSymbols);
    }
    if (info->CallSymbols)
    {
        for (i = 0; i < info->numCallSymbols; i++)
            Rf_freeCallSymbol(info->CallSymbols + i);
        free(info->CallSymbols);
    }
    if (info->FortranSymbols)
    {
        for (i = 0; i < info->numFortranSymbols; i++)
            Rf_freeFortranSymbol(info->FortranSymbols + i);
        free(info->FortranSymbols);
    }
}

/* Remove the specified DLL from the current DLL list */
/* Returns 1 if the DLL was found and removed from */
/* the list and returns 0 otherwise. */

static int DeleteDLL(char *path)
{
    int i, loc;

    for (i = 0; i < CountDLL; i++)
    {
        if (!strcmp(path, LoadedDLL[i].path))
        {
            loc = i;
            goto found;
        }
    }
    return 0;
found:
#ifdef CACHE_DLL_SYM
    if (R_osDynSymbol->deleteCachedSymbols)
        R_osDynSymbol->deleteCachedSymbols(&LoadedDLL[loc]);
#endif
    R_osDynSymbol->closeLibrary(LoadedDLL[loc].handle);
    Rf_freeDllInfo(LoadedDLL + loc);
    for (i = loc + 1; i < CountDLL; i++)
    {
        LoadedDLL[i - 1].path = LoadedDLL[i].path;
        LoadedDLL[i - 1].name = LoadedDLL[i].name;
        LoadedDLL[i - 1].handle = LoadedDLL[i].handle;
        LoadedDLL[i - 1].numCSymbols = LoadedDLL[i].numCSymbols;
        LoadedDLL[i - 1].numCallSymbols = LoadedDLL[i].numCallSymbols;
        LoadedDLL[i - 1].numFortranSymbols = LoadedDLL[i].numFortranSymbols;
        LoadedDLL[i - 1].CSymbols = LoadedDLL[i].CSymbols;
        LoadedDLL[i - 1].CallSymbols = LoadedDLL[i].CallSymbols;
        LoadedDLL[i - 1].FortranSymbols = LoadedDLL[i].FortranSymbols;
    }
    CountDLL--;
    return 1;
}

DL_FUNC Rf_lookupCachedSymbol(const char *name, const char *pkg, int all)
{
#ifdef CACHE_DLL_SYM
    int i;
#ifdef Macintosh
    all = 0;
#endif
    for (i = 0; i < nCPFun; i++)
        if (!strcmp(name, CPFun[i].name) && (all || !strcmp(pkg, CPFun[i].pkg)))
            return CPFun[i].func;
#endif

    return ((DL_FUNC)NULL);
}

#ifdef Win32
#define DLLerrBUFSIZE 4000
#else /* Not Windows */
#define DLLerrBUFSIZE 1000
#endif

static char DLLerror[DLLerrBUFSIZE] = "";

/* the error message; length taken from ERRBUFSIZE in ./hpdlfcn.c  */

/* Inserts the specified DLL at the head of the DLL list */
/* Returns 1 if the library was successfully added */
/* and returns 0 if the library table is full or */
/* or if dlopen fails for some reason. */

static DL_FUNC R_dlsym(DllInfo *dll, char const *name, R_RegisteredNativeSymbol *symbol);

static int AddDLL(char *path, int asLocal, int now)
{
    HINSTANCE handle;
    DllInfo *info;

    DeleteDLL(path);
    if (CountDLL == MAX_NUM_DLLS)
    {
        strcpy(DLLerror, "Maximal number of DLLs reached...");
        return 0;
    }

    handle = R_osDynSymbol->loadLibrary(path, asLocal, now);

    if (handle == NULL)
    {
        R_osDynSymbol->getError(DLLerror, DLLerrBUFSIZE);
        return 0;
    }

    info = R_RegisterDLL(handle, path);

    /* Now look for an initializing routine named R_init_<library name>.
       If it is present, we invoke it. It should take a reference to the
       DllInfo object currently being initialized.
    */
    if (info)
    {
        char *tmp;
        DL_FUNC f;
#ifdef HAVE_NO_SYMBOL_UNDERSCORE
        tmp = (char *)malloc(sizeof(char) * (strlen("R_init_") + strlen(info->name) + 1));
        sprintf(tmp, "%s%s", "R_init_", info->name);
#else
        tmp = (char *)malloc(sizeof(char) * (strlen("R_init_") + strlen(info->name) + 2));
        sprintf(tmp, "_%s%s", "R_init_", info->name);
#endif
        f = (DL_FUNC)R_osDynSymbol->dlsym(info, tmp);
        free(tmp);
        if (f)
            f(info);
    }

    return 1;
}

DllInfo *R_RegisterDLL(HINSTANCE handle, const char *path)
{
    char *dpath, DLLname[PATH_MAX], *p;
    DllInfo *info;

    info = &LoadedDLL[CountDLL];
    /* default is to use old-style dynamic lookup. Library's
       initialization routine can limit access by setting this to FALSE.
    */
    info->useDynamicLookup = TRUE;

    dpath = malloc(strlen(path) + 1);
    if (dpath == NULL)
    {
        strcpy(DLLerror, "Couldn't allocate space for 'path'");
        R_osDynSymbol->closeLibrary(handle);
        return 0;
    }
    strcpy(dpath, path);

    if (R_osDynSymbol->fixPath)
        R_osDynSymbol->fixPath(dpath);

    /* keep only basename from path */
    p = strrchr(dpath, FILESEP[0]);
    if (!p)
        p = dpath;
    else
        p++;
    if (strlen(p) < PATH_MAX)
        strcpy(DLLname, p);
    else
        error("DLLname %s is too long", p);

    /* FIXME: didn't work on Mac, unsafe
    p = strchr(DLLname, '.');
    if(p) *p = '\0'; */

    /* remove SHLIB_EXT if present */
    p = DLLname + strlen(DLLname) - strlen(SHLIB_EXT);
#ifdef Win32 /* case-insensitive file system */
    if (p > DLLname && stricmp(p, SHLIB_EXT) == 0)
        *p = '\0';
#else
    if (p > DLLname && strcmp(p, SHLIB_EXT) == 0)
        *p = '\0';
#endif

    addDLL(dpath, DLLname, handle);

    return (info);
}

int addDLL(char *dpath, char *DLLname, HINSTANCE handle)
{
    int ans = CountDLL;
    char *name = malloc(strlen(DLLname) + 1);
    if (name == NULL)
    {
        strcpy(DLLerror, "Couldn't allocate space for 'name'");
        if (handle)
            R_osDynSymbol->closeLibrary(handle);
        free(dpath);
        return 0;
    }

    strcpy(name, DLLname);
    LoadedDLL[CountDLL].path = dpath;
    LoadedDLL[CountDLL].name = name;
    LoadedDLL[CountDLL].handle = handle;

    LoadedDLL[CountDLL].numCSymbols = 0;
    LoadedDLL[CountDLL].numCallSymbols = 0;
    LoadedDLL[CountDLL].numFortranSymbols = 0;
    LoadedDLL[CountDLL].CSymbols = NULL;
    LoadedDLL[CountDLL].CallSymbols = NULL;
    LoadedDLL[CountDLL].FortranSymbols = NULL;
    CountDLL++;

    return (ans);
}

Rf_DotCSymbol *Rf_lookupRegisteredCSymbol(DllInfo *info, const char *name)
{
    int i;
    for (i = 0; i < info->numCSymbols; i++)
    {
        if (strcmp(name, info->CSymbols[i].name) == 0)
            return (&(info->CSymbols[i]));
    }

    return (NULL);
}

Rf_DotFortranSymbol *Rf_lookupRegisteredFortranSymbol(DllInfo *info, const char *name)
{
    int i;
    for (i = 0; i < info->numFortranSymbols; i++)
    {
        if (strcmp(name, info->FortranSymbols[i].name) == 0)
            return (&(info->FortranSymbols[i]));
    }

    return ((Rf_DotFortranSymbol *)NULL);
}

Rf_DotCallSymbol *Rf_lookupRegisteredCallSymbol(DllInfo *info, const char *name)
{
    int i;

    for (i = 0; i < info->numCallSymbols; i++)
    {
        if (strcmp(name, info->CallSymbols[i].name) == 0)
            return (&(info->CallSymbols[i]));
    }
    return ((Rf_DotCallSymbol *)NULL);
}

Rf_DotExternalSymbol *Rf_lookupRegisteredExternalSymbol(DllInfo *info, const char *name)
{
    int i;

    for (i = 0; i < info->numExternalSymbols; i++)
    {
        if (strcmp(name, info->ExternalSymbols[i].name) == 0)
            return (&(info->ExternalSymbols[i]));
    }
    return ((Rf_DotExternalSymbol *)NULL);
}

DL_FUNC R_getDLLRegisteredSymbol(DllInfo *info, const char *name, R_RegisteredNativeSymbol *symbol)
{
    int fail = 0;
    NativeSymbolType purpose = R_ANY_SYM;

    if (symbol)
    {
        purpose = symbol->type;
    }
    if ((purpose == R_ANY_SYM || purpose == R_C_SYM) && info->numCSymbols > 0)
    {
        Rf_DotCSymbol *sym;
        sym = Rf_lookupRegisteredCSymbol(info, name);
        if (sym)
        {
            if (symbol)
            {
                symbol->type = R_C_SYM;
                symbol->symbol.c = sym;
                symbol->dll = info;
            }

            return ((DL_FUNC)sym->fun);
        }
        fail = 1;
    }

    if ((purpose == R_ANY_SYM || purpose == R_CALL_SYM) && info->numCallSymbols > 0)
    {
        Rf_DotCallSymbol *sym;
        sym = Rf_lookupRegisteredCallSymbol(info, name);
        if (sym)
        {
            if (symbol)
            {
                symbol->type = R_CALL_SYM;
                symbol->symbol.call = sym;
                symbol->dll = info;
            }
            return ((DL_FUNC)sym->fun);
        }
        fail = 1;
    }

    if ((purpose == R_ANY_SYM || purpose == R_FORTRAN_SYM) && info->numCallSymbols > 0)
    {
        Rf_DotFortranSymbol *sym;
        sym = Rf_lookupRegisteredFortranSymbol(info, name);
        if (sym)
        {
            if (symbol)
            {
                symbol->type = R_FORTRAN_SYM;
                symbol->symbol.fortran = sym;
                symbol->dll = info;
            }
            return ((DL_FUNC)sym->fun);
        }
        fail = 1;
    }

    if ((purpose == R_ANY_SYM || purpose == R_EXTERNAL_SYM) && info->numExternalSymbols > 0)
    {
        Rf_DotExternalSymbol *sym;
        sym = Rf_lookupRegisteredExternalSymbol(info, name);
        if (sym)
        {
            if (symbol)
            {
                symbol->type = R_EXTERNAL_SYM;
                symbol->symbol.external = sym;
                symbol->dll = info;
            }
            return ((DL_FUNC)sym->fun);
        }
        fail = 1;
    }

    return ((DL_FUNC)NULL);
}

static DL_FUNC R_dlsym(DllInfo *info, char const *name, R_RegisteredNativeSymbol *symbol)
{
    char buf[MAXIDSIZE + 1];
    DL_FUNC f;

    f = R_getDLLRegisteredSymbol(info, name, symbol);
    if (f)
        return (f);

    if (info->useDynamicLookup == FALSE)
        return (NULL);

#ifdef HAVE_NO_SYMBOL_UNDERSCORE
    sprintf(buf, "%s", name);
#else
    sprintf(buf, "_%s", name);
#endif

#ifdef HAVE_F77_UNDERSCORE
    if (symbol && symbol->type == R_FORTRAN_SYM)
    {
        buf[strlen(buf) + 1] = '\0';
        buf[strlen(buf)] = '_';
    }
#endif

    return (DL_FUNC)R_osDynSymbol->dlsym(info, buf);
}

/* R_FindSymbol checks whether one of the libraries */
/* that have been loaded contains the symbol name and */
/* returns a pointer to that symbol upon success. */

/*
  In the future, this will receive an additional argument
  which will specify the nature of the symbol expected by the
  caller, specifically whether it is for a .C(), .Call(),
  .Fortran(), .External(), generic, etc. invocation. This will
  reduce the pool of possible symbols in the case of a library
  that registers its routines.
 */

DL_FUNC R_FindSymbol(char const *name, char const *pkg, R_RegisteredNativeSymbol *symbol)
{
    DL_FUNC fcnptr = (DL_FUNC)NULL;
#ifndef Macintosh
    int i, all = (strlen(pkg) == 0), doit;
#else /* cannot load locally */
    int i, all = (strlen("") == 0), doit;
#endif

    if (R_osDynSymbol->lookupCachedSymbol)
        fcnptr = R_osDynSymbol->lookupCachedSymbol(name, pkg, all);

    if (fcnptr)
        return (fcnptr);

    /* The following is not legal ANSI C. */
    /* It is only meant to be used in systems supporting */
    /* the dlopen() interface, in which systems data and  */
    /* function pointers _are_ the same size and _can_   */
    /* be cast without loss of information.		     */

    for (i = CountDLL - 1; i >= 0; i--)
    {
        doit = all;
        if (!doit && !strcmp(pkg, LoadedDLL[i].name))
            doit = 2;
        if (doit)
        {
            fcnptr = R_dlsym(&LoadedDLL[i], name, symbol); /* R_osDynSymbol->dlsym */
            if (fcnptr != (DL_FUNC)NULL)
            {
                if (symbol)
                    symbol->dll = LoadedDLL + i;
#ifdef CACHE_DLL_SYM
                if (strlen(pkg) <= 20 && strlen(name) <= 20 && nCPFun < 100)
                {
                    strcpy(CPFun[nCPFun].pkg, LoadedDLL[i].name);
                    strcpy(CPFun[nCPFun].name, name);
                    CPFun[nCPFun++].func = fcnptr;
                }
#endif
                return fcnptr;
            }
        }
        if (doit > 1)
            return (DL_FUNC)NULL; /* Only look in the first-matching DLL */
    }

    return (DL_FUNC)NULL;
}

static void GetFullDLLPath(SEXP call, char *buf, char *path)
{
    R_osDynSymbol->getFullDLLPath(call, buf, path);
}

/* do_dynload implements the R-Interface for the */
/* loading of shared libraries */

/*
  Extended to support 2 additional arguments (3 in total).
  First argument is the name of the library.
  Second argument is a logical indicating whether we
  want the symbols to be kept in their own local symbol table
  or added to the global symbol table of the application.
  Third argument is a logical indicating whether the
  dynamic loading should relocate all routine symbols
  now and signal any errors immediately or lazily relocate
  the symbols as they are invoked. This is useful for
  developers so that they can ensure that all the symbols
  are available before they release, and allows users to
  call routines from "incomplete" libraries.
 */

SEXP do_dynload(SEXP call, SEXP op, SEXP args, SEXP env)
{
    char buf[2 * PATH_MAX];

    checkArity(op, args);
    if (!isString(CAR(args)) || length(CAR(args)) < 1)
        errorcall(call, "character argument expected");
    GetFullDLLPath(call, buf, CHAR(STRING_ELT(CAR(args), 0)));
    /* AddDLL does this DeleteDLL(buf); */
    if (!AddDLL(buf, LOGICAL(CADR(args))[0], LOGICAL(CADDR(args))[0]))
        errorcall(call, "unable to load shared library \"%s\":\n  %s", buf, DLLerror);
    return R_NilValue;
}

SEXP do_dynunload(SEXP call, SEXP op, SEXP args, SEXP env)
{
    char buf[2 * PATH_MAX];

    checkArity(op, args);
    if (!isString(CAR(args)) || length(CAR(args)) < 1)
        errorcall(call, "character argument expected");
    GetFullDLLPath(call, buf, CHAR(STRING_ELT(CAR(args), 0)));
    if (!DeleteDLL(buf))
        errorcall(call, "dynamic/shared library \"%s\" was not loaded", buf);
    return R_NilValue;
}

int moduleCdynload(char *module, int local, int now)
{
#ifndef Macintosh
    char dllpath[PATH_MAX], *p = getenv("R_HOME");
#else
    char dllpath[PATH_MAX], *p = R_Home;
#endif
    if (!p)
        return 0;
    sprintf(dllpath, "%s%smodules%s%s%s", p, FILESEP, FILESEP, module, SHLIB_EXT);
    return AddDLL(dllpath, local, now);
}

/**
  Creates an R object representing the value of the
  function pointer given by `f'. This object has class
  NativeSymbol and can be used to relay symbols from
  one library to another.
 */
SEXP Rf_MakeNativeSymbolRef(DL_FUNC f)
{
    SEXP ref, klass;

    PROTECT(ref = R_MakeExternalPtr((void *)f, Rf_install("native symbol"), NULL));

    PROTECT(klass = allocVector(STRSXP, 1));
    SET_STRING_ELT(klass, 0, mkChar("NativeSymbol"));
    setAttrib(ref, R_ClassSymbol, klass);

    UNPROTECT(2);
    return (ref);
}

/**
 Creates an R object representing the public DLL information stored in
 info. Currently this is only the short and the long, fully qualified
 name of the DLL and whether we only look for symbols that have been
 registered in this DLL or do we also use dynamic lookup.
 */
SEXP Rf_MakeDLLInfo(DllInfo *info)
{
    SEXP ref, elNames, tmp;
    int i, n;
    const char *const names[] = {"name", "path", "dynamicLookup"};

    n = sizeof(names) / sizeof(names[0]);

    PROTECT(ref = allocVector(VECSXP, n));
    SET_VECTOR_ELT(ref, 0, tmp = allocVector(STRSXP, 1));
    if (info->name)
        SET_STRING_ELT(tmp, 0, mkChar(info->name));
    SET_VECTOR_ELT(ref, 1, tmp = allocVector(STRSXP, 1));
    if (info->path)
        SET_STRING_ELT(tmp, 0, mkChar(info->path));
    SET_VECTOR_ELT(ref, 2, ScalarLogical(info->useDynamicLookup));

    PROTECT(elNames = allocVector(STRSXP, n));
    for (i = 0; i < n; i++)
        SET_STRING_ELT(elNames, i, mkChar(names[i]));
    setAttrib(ref, R_NamesSymbol, elNames);
    UNPROTECT(2);

    return (ref);
}

/**
  This is the routine associated with the getNativeSymbolInfo()
  function and it takes the name of a symbol and optionally a
  library identifier (package usually) in which to restrict the search
  for this symbol. It resolves the symbol and returns it to the caller
  giving the symbol address, the package information (i.e. name and
  fully qualified shared library name). If the symbol was explicitly
  registered (rather than dynamically resolved by R), then we pass
  back that information also, giving the number of arguments it
  expects and the interface by which it should be called.
  The returned object has class NativeSymbol. If the symbol was
  registered, we add a class identifying the interface type
  for which it is intended (i.e. .C(), .Call(), etc.)
 */
SEXP R_getSymbolInfo(SEXP sname, SEXP spackage)
{
    char *package, *name;
    R_RegisteredNativeSymbol symbol = {R_ANY_SYM, {NULL}, NULL};
    SEXP sym = R_NilValue;
    DL_FUNC f;

    name = CHAR(STRING_ELT(sname, 0));
    if (length(spackage))
        package = CHAR(STRING_ELT(spackage, 0));
    else
        package = "";
    f = R_FindSymbol(name, package, &symbol);
    if (f)
    {
        SEXP tmp, klass;
        int n = (symbol.type != R_ANY_SYM) ? 4 : 3;
        PROTECT(sym = allocVector(VECSXP, n));

        SET_VECTOR_ELT(sym, 0, sname);
        SET_VECTOR_ELT(sym, 1, Rf_MakeNativeSymbolRef(f));
        if (symbol.dll)
            SET_VECTOR_ELT(sym, 2, Rf_MakeDLLInfo(symbol.dll));

        PROTECT(klass = allocVector(STRSXP, (symbol.type != R_ANY_SYM ? 2 : 1)));
        SET_STRING_ELT(klass, length(klass) - 1, mkChar("NativeSymbolInfo"));

        if (n > 3)
        {
            /* Add the registration information: the number of arguments and the classname. */
            int nargs = -1;
            char *className = "";
            switch (symbol.type)
            {
            case R_C_SYM:
                nargs = symbol.symbol.c->numArgs;
                className = "CRoutine";
                break;
            case R_CALL_SYM:
                nargs = symbol.symbol.call->numArgs;
                className = "CallRoutine";
                break;
            case R_FORTRAN_SYM:
                nargs = symbol.symbol.fortran->numArgs;
                className = "FortranRoutine";
                break;
            case R_EXTERNAL_SYM:
                nargs = symbol.symbol.external->numArgs;
                className = "ExternalRoutine";
                break;
            default:
                /* Something unintended has happened if we get here. */
                break;
            }
            SET_VECTOR_ELT(sym, 3, tmp = ScalarInteger(nargs));
            SET_STRING_ELT(klass, 0, mkChar(className));
        }
        setAttrib(sym, R_ClassSymbol, klass);
        UNPROTECT(2);
    }

    return (sym);
}

#else /* no dyn.load support */

void InitFunctionHashing()
{
#ifdef OLD
    NaokSymbol = install("NAOK");
    DupSymbol = install("DUP");
#endif
}

DL_FUNC R_FindSymbol(char const *name, char const *pkg, R_RegisteredNativeSymbol *symbol)
{
    int i;
    for (i = 0; CFunTab[i].name; i++)
        if (!strcmp(name, CFunTab[i].name))
            return CFunTab[i].func;
    return (DL_FUNC)0;
}

SEXP do_dynload(SEXP call, SEXP op, SEXP args, SEXP env)
{
    error("no dyn.load support in this R version");
    return (R_NilValue);
}

SEXP do_dynunload(SEXP call, SEXP op, SEXP args, SEXP env)
{
    error("no dyn.load support in this R version");
    return (R_NilValue);
}

SEXP R_getSymbolInfo(SEXP sname, SEXP spackage)
{
    error("no dyn.load support in this R version");
}
#endif
