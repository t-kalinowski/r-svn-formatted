/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995-1996 Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997-2001 The R Development Core Team
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

/*  Dynamic Loading Support: See ../main/Rdynload.c and ../include/Rdynpriv.h
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <stdlib.h>
#include <Defn.h>
#include <Rmath.h>
#include <direct.h>
#include <windows.h>

#include <R_ext/Rdynload.h>
#include <Rdynpriv.h>

/* Inserts the specified DLL at the head of the DLL list */
/* Returns 1 if the library was successfully added */
/* and returns 0 if the library table is full or */
/* or if LoadLibrary fails for some reason. */

static void fixPath(char *path)
{
    char *p;
    for (p = path; *p != '\0'; p++)
        if (*p == '\\')
            *p = '/';
}

static HINSTANCE R_loadLibrary(const char *path, int asLocal, int now);
static DL_FUNC getRoutine(DllInfo *info, char const *name);
static void R_deleteCachedSymbols(DllInfo *dll);

static void R_getDLLError(char *buf, int len);
static void GetFullDLLPath(SEXP call, char *buf, char *path);

static void closeLibrary(HINSTANCE handle)
{
    FreeLibrary(handle);
}

void InitFunctionHashing()
{
    R_osDynSymbol->loadLibrary = R_loadLibrary;
    R_osDynSymbol->dlsym = getRoutine;
    R_osDynSymbol->closeLibrary = closeLibrary;
    R_osDynSymbol->getError = R_getDLLError;

    R_osDynSymbol->deleteCachedSymbols = R_deleteCachedSymbols;
    R_osDynSymbol->lookupCachedSymbol = Rf_lookupCachedSymbol;

    R_osDynSymbol->fixPath = fixPath;
    R_osDynSymbol->getFullDLLPath = GetFullDLLPath;
}

static void R_deleteCachedSymbols(DllInfo *dll)
{
    int i;
    for (i = nCPFun - 1; i >= 0; i--)
        if (!strcmp(CPFun[i].pkg, dll->name))
        {
            if (i < nCPFun - 1)
            {
                strcpy(CPFun[i].name, CPFun[--nCPFun].name);
                strcpy(CPFun[i].pkg, CPFun[nCPFun].pkg);
                CPFun[i].func = CPFun[nCPFun].func;
            }
            else
                nCPFun--;
        }
}

HINSTANCE R_loadLibrary(const char *path, int asLocal, int now)
{
    HINSTANCE tdlh;
    unsigned int dllcw, rcw;

    rcw = _controlfp(0, 0);
    _clearfp();
    tdlh = LoadLibrary(path);
    dllcw = _controlfp(0, 0);

    if (dllcw != rcw)
    {
        warning("DLL attempted to change FPU control word from %x to %x", rcw, dllcw);
        _controlfp(rcw, _MCW_EM | _MCW_IC | _MCW_RC | _MCW_PC);
    }
    return (tdlh);
}

static DL_FUNC getRoutine(DllInfo *info, char const *name)
{
    DL_FUNC f;
    f = (DL_FUNC)GetProcAddress(info->handle, name);
    return (f);
}

static void R_getDLLError(char *buf, int len)
{
    LPVOID lpMsgBuf;
    FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, NULL,
                  GetLastError(), MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPTSTR)&lpMsgBuf, 0, NULL);
    strcpy(buf, "LoadLibrary failure:  ");
    strcat(buf, lpMsgBuf);
    LocalFree(lpMsgBuf);
}

static void GetFullDLLPath(SEXP call, char *buf, char *path)
{
    char *p;

    if ((path[0] != '/') && (path[0] != '\\') && (path[1] != ':'))
    {
        if (!getcwd(buf, MAX_PATH))
            errorcall(call, "can't get working directory!");
        strcat(buf, "\\");
        strcat(buf, path);
    }
    else
        strcpy(buf, path);
    /* fix slashes to allow inconsistent usage later */
    for (p = buf; *p; p++)
        if (*p == '\\')
            *p = '/';
}
