/*
 *  R : A Computer Language for Statistical Data Analysis
 *  file extra.c
 *  Copyright (C) 1998--2000  Guido Masarotto and Brian Ripley
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

/* extra commands for R */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "Defn.h"
#include "Fileio.h"
#include <io.h>
#include <direct.h>
#include <time.h>
#include <windows.h>
#include "graphapp/ga.h"
#include "rui.h"

char *Rwin32_tmpnam(char *prefix)
{
    char *tmp, tm[MAX_PATH], tmp1[MAX_PATH], *p, *res;
    int hasspace = 0;
    unsigned int n, done = 0;
    WIN32_FIND_DATA fd;
    HANDLE h;

    tmp = getenv("TMP");
    if (!tmp)
        tmp = getenv("TEMP");
    if (!tmp)
        tmp = getenv("R_USER"); /* this one will succeed */
    /* make sure no spaces in path */
    for (p = tmp; *p; p++)
        if (isspace(*p))
        {
            hasspace = 1;
            break;
        }
    if (hasspace)
        GetShortPathName(tmp, tmp1, MAX_PATH);
    else
        strcpy(tmp1, tmp);
    for (n = 0; n < 100; n++)
    {
        /* try a random number at the end */
        sprintf(tm, "%s\\%s%d", tmp1, prefix, rand());
        if ((h = FindFirstFile(tm, &fd)) == INVALID_HANDLE_VALUE)
        {
            done = 1;
            break;
        }
        FindClose(h);
        tm[0] = '\0';
    }
    if (!done)
        error("cannot find unused tempfile name");
    res = (char *)malloc((strlen(tm) + 1) * sizeof(char));
    strcpy(res, tm);
    return res;
}

SEXP do_tempfile(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans;
    char *tn, *tm;
    int i, slen = 0 /* -Wall */;

    checkArity(op, args);
    if (!isString(CAR(args)) || (slen = LENGTH(CAR(args))) < 1)
        errorcall(call, "invalid file name argument");
    PROTECT(ans = allocVector(STRSXP, slen));
    for (i = 0; i < slen; i++)
    {
        tn = CHAR(STRING_ELT(CAR(args), i));
        /* try to get a new file name */
        tm = Rwin32_tmpnam(tn);
        SET_STRING_ELT(ans, i, mkChar(tm));
        free(tm);
    }
    UNPROTECT(1);
    return (ans);
}

SEXP do_dircreate(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP path, ans;
    char *p, dir[MAX_PATH];
    int res;

    checkArity(op, args);
    path = CAR(args);
    if (!isString(path) || length(path) != 1)
        errorcall(call, "invalid path argument");
    strcpy(dir, CHAR(STRING_ELT(path, 0)));
    for (p = dir; *p != '\0'; p++)
        if (*p == '/')
            *p = '\\';
    res = mkdir(dir);
    PROTECT(ans = allocVector(LGLSXP, 1));
    LOGICAL(ans)[0] = (res == 0);
    UNPROTECT(1);
    return (ans);
}

#include <sys/types.h>
#include <sys/stat.h>

static int R_unlink(char *names, int recursive);

static int R_unlink_one(char *dir, char *name, int recursive)
{
    char tmp[MAX_PATH];

    if (strcmp(name, ".") == 0)
        return 0;
    if (strcmp(name, "..") == 0)
        return 0;
    strcpy(tmp, dir);
    strcat(tmp, "\\");
    strcat(tmp, name);
    return (recursive ? R_unlink(tmp, 1) : unlink(tmp)) != 0;
}

static int R_unlink(char *names, int recursive)
{
    int failures = 0;
    char *p, tmp[MAX_PATH], dir[MAX_PATH];
    WIN32_FIND_DATA find_data;
    HANDLE fh;
    struct stat sb;

    strcpy(tmp, names);
    for (p = tmp; *p != '\0'; p++)
        if (*p == '/')
            *p = '\\';
    if (stat(tmp, &sb) == 0)
    {
        /* Is this a directory? */
        if (sb.st_mode & _S_IFDIR)
        {
            if (recursive)
            {
                strcpy(dir, tmp);
                strcat(tmp, "\\*");
                fh = FindFirstFile(tmp, &find_data);
                if (fh != INVALID_HANDLE_VALUE)
                {
                    failures += R_unlink_one(dir, find_data.cFileName, 1);
                    while (FindNextFile(fh, &find_data))
                        failures += R_unlink_one(dir, find_data.cFileName, 1);
                    FindClose(fh);
                }
                if (rmdir(dir))
                    failures++;
            }
            else
                failures++; /* don't try to delete dirs */
        }
        else
        { /* Regular file (or several) */
            strcpy(dir, tmp);
            if ((p = strrchr(dir, '\\')))
                *(++p) = '\0';
            else
                *dir = '\0';
            /* check for wildcard matches */
            fh = FindFirstFile(tmp, &find_data);
            if (fh != INVALID_HANDLE_VALUE)
            {
                failures += R_unlink_one(dir, find_data.cFileName, 0);
                while (FindNextFile(fh, &find_data))
                    failures += R_unlink_one(dir, find_data.cFileName, 0);
                FindClose(fh);
            }
        }
    }
    return failures;
}

SEXP do_unlink(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP fn, ans;
    int i, nfiles, failures = 0, recursive;

    checkArity(op, args);
    fn = CAR(args);
    nfiles = length(fn);
    if (!isString(fn) || nfiles < 1)
        errorcall(call, "invalid file name argument");
    recursive = asLogical(CADR(args));
    if (recursive == NA_LOGICAL)
        errorcall(call, "invalid recursive argument");
    for (i = 0; i < nfiles; i++)
        failures += R_unlink(CHAR(STRING_ELT(fn, i)), recursive);
    PROTECT(ans = allocVector(INTSXP, 1));
    if (!failures)
        INTEGER(ans)[0] = 0;
    else
        INTEGER(ans)[0] = 1;
    UNPROTECT(1);
    return (ans);
}

SEXP do_helpstart(SEXP call, SEXP op, SEXP args, SEXP env)
{
    char *home, buf[MAX_PATH];
    FILE *ff;

    checkArity(op, args);
    home = getenv("R_HOME");
    if (home == NULL)
        error("R_HOME not set");
    sprintf(buf, "%s\\doc\\html\\rwin.html", home);
    ff = fopen(buf, "r");
    if (!ff)
    {
        sprintf(buf, "%s\\doc\\html\\rwin.htm", home);
        ff = fopen(buf, "r");
        if (!ff)
        {
            sprintf(buf, "%s\\doc\\html\\rwin.htm[l] not found", home);
            error(buf);
        }
    }
    fclose(ff);
    ShellExecute(NULL, "open", buf, NULL, home, SW_SHOW);
    return R_NilValue;
}

static int nhfiles = 0;
static char *hfiles[50];

SEXP do_helpitem(SEXP call, SEXP op, SEXP args, SEXP env)
{
    /*
     * type = 1: launch html file.
     *        2: "topic", 2, Windows help file.
     *        3: notify are finished with the help file.
     */

    char *item, *hfile;
    char *home, buf[MAX_PATH];
    FILE *ff;
    int type;

    checkArity(op, args);
    if (!isString(CAR(args)))
        errorcall(call, "invalid topic argument");
    item = CHAR(STRING_ELT(CAR(args), 0));
    type = asInteger(CADR(args));
    if (type == 1)
    {
        ff = fopen(item, "r");
        if (!ff)
        {
            sprintf(buf, "%s not found", item);
            error(buf);
        }
        fclose(ff);
        home = getenv("R_HOME");
        if (home == NULL)
            error("R_HOME not set");
        ShellExecute(NULL, "open", item, NULL, home, SW_SHOW);
    }
    else if (type == 2)
    {
        if (!isString(CADDR(args)))
            errorcall(call, "invalid hlpfile argument");
        hfile = CHAR(STRING_ELT(CADDR(args), 0));
        if (!WinHelp((HWND)0, hfile, HELP_KEY, (DWORD)item))
            warning("WinHelp call failed");
        else
        {
            if (nhfiles >= 50)
                error("too many .hlp files opened");
            hfiles[nhfiles] = malloc(strlen(hfile) * sizeof(char));
            strcpy(hfiles[nhfiles++], hfile);
        }
    }
    else if (type == 3)
    {
        if (!isString(CADDR(args)))
            warningcall(call, "invalid hlpfile argument");
        hfile = CHAR(STRING_ELT(CADDR(args), 0));
        if (!WinHelp((HWND)0, hfile, HELP_QUIT, (DWORD)0))
            error("WinHelp call failed");
    }
    else
        warning("type not yet implemented");
    return R_NilValue;
}

void closeAllHlpFiles()
{
    int i;

    for (i = nhfiles - 1; i >= 0; i--)
        WinHelp((HWND)0, hfiles[i], HELP_QUIT, (DWORD)0);
}

SEXP do_flushconsole(SEXP call, SEXP op, SEXP args, SEXP env)
{
    R_FlushConsole();
    return R_NilValue;
}

#include <winbase.h>
/* typedef struct _OSVERSIONINFO{
    DWORD dwOSVersionInfoSize;
    DWORD dwMajorVersion;
    DWORD dwMinorVersion;
    DWORD dwBuildNumber;
    DWORD dwPlatformId;
    TCHAR szCSDVersion[ 128 ];
    } OSVERSIONINFO; */

typedef struct _OSVERSIONINFOEXA
{
    DWORD dwOSVersionInfoSize;
    DWORD dwMajorVersion;
    DWORD dwMinorVersion;
    DWORD dwBuildNumber;
    DWORD dwPlatformId;
    CHAR szCSDVersion[128];
    WORD wServicePackMajor;
    WORD wServicePackMinor;
    WORD wReserved[2];
} OSVERSIONINFOEXA;

SEXP do_winver(SEXP call, SEXP op, SEXP args, SEXP env)
{
    char isNT[8] = "??", ver[256];
    SEXP ans;
    OSVERSIONINFO verinfo;
    OSVERSIONINFOEXA verinfoex;

    checkArity(op, args);
    verinfo.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
    GetVersionEx(&verinfo);
    switch (verinfo.dwPlatformId)
    {
    case VER_PLATFORM_WIN32_NT:
        strcpy(isNT, "NT");
        break;
    case VER_PLATFORM_WIN32_WINDOWS:
        strcpy(isNT, "9x");
        break;
    case VER_PLATFORM_WIN32s:
        strcpy(isNT, "win32s");
        break;
    default:
        sprintf(isNT, "ID=%d", (int)verinfo.dwPlatformId);
        break;
    }

    if ((int)verinfo.dwMajorVersion >= 5)
    {
        verinfo.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEXA);
        GetVersionEx(&verinfoex);
        sprintf(ver, "Windows %d.%d (build %d) Service Pack %d.%d (%s)", (int)verinfoex.dwMajorVersion,
                (int)verinfoex.dwMinorVersion, LOWORD(verinfoex.dwBuildNumber), (int)verinfoex.wServicePackMajor,
                (int)verinfoex.wServicePackMinor, verinfoex.szCSDVersion);
    }
    else
    {
        sprintf(ver, "Windows %s %d.%d (build %d) %s", isNT, (int)verinfo.dwMajorVersion, (int)verinfo.dwMinorVersion,
                LOWORD(verinfo.dwBuildNumber), verinfo.szCSDVersion);
    }

    PROTECT(ans = allocVector(STRSXP, 1));
    SET_STRING_ELT(ans, 0, mkChar(ver));
    UNPROTECT(1);
    return (ans);
}

void internal_shellexec(char *file)
{
    char *home, buf[MAX_PATH];

    home = getenv("R_HOME");
    if (home == NULL)
        error("R_HOME not set");
    strcpy(buf, file);
    ShellExecute(NULL, "open", buf, NULL, home, SW_SHOW);
}

SEXP do_shellexec(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP file;

    checkArity(op, args);
    file = CAR(args);
    if (!isString(file) || length(file) != 1)
        errorcall(call, "invalid file argument");
    internal_shellexec(CHAR(STRING_ELT(file, 0)));
    return R_NilValue;
}

int check_doc_file(char *file)
{
    char *home, path[MAX_PATH];
    struct stat sb;

    home = getenv("R_HOME");
    if (home == NULL)
        error("R_HOME not set");
    strcpy(path, home);
    strcat(path, "/");
    strcat(path, file);
    return stat(path, &sb) == 0;
}

SEXP do_windialog(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP message, ans;
    char *type;
    int res = YES;

    checkArity(op, args);
    type = CHAR(STRING_ELT(CAR(args), 0));
    message = CADR(args);
    if (!isString(message) || length(message) != 1)
        error("invalid `message' argument");
    if (strcmp(type, "ok") == 0)
    {
        askok(CHAR(STRING_ELT(message, 0)));
        res = 10;
    }
    else if (strcmp(type, "okcancel") == 0)
    {
        res = askokcancel(CHAR(STRING_ELT(message, 0)));
        if (res == YES)
            res = 2;
    }
    else if (strcmp(type, "yesno") == 0)
    {
        res = askyesno(CHAR(STRING_ELT(message, 0)));
    }
    else if (strcmp(type, "yesnocancel") == 0)
    {
        res = askyesnocancel(CHAR(STRING_ELT(message, 0)));
    }
    else
        errorcall(call, "unknown type");
    ans = allocVector(INTSXP, 1);
    INTEGER(ans)[0] = res;
    return (ans);
}

SEXP do_windialogstring(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP message, def, ans;
    char *string;

    checkArity(op, args);
    message = CAR(args);
    if (!isString(message) || length(message) != 1)
        error("invalid `message' argument");
    def = CADR(args);
    if (!isString(def) || length(def) != 1)
        error("invalid `default' argument");
    string = askstring(CHAR(STRING_ELT(message, 0)), CHAR(STRING_ELT(def, 0)));
    if (string)
    {
        ans = allocVector(STRSXP, 1);
        SET_STRING_ELT(ans, 0, mkChar(string));
        return (ans);
    }
    else
        return (R_NilValue);
}

#include "Startup.h"
extern UImode CharacterMode;
static char msgbuf[256];

SEXP do_winmenuadd(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP smenu, sitem;
    int res;
    char errmsg[50];

    checkArity(op, args);
    if (CharacterMode != RGui)
        errorcall(call, "Menu functions can only be used in the GUI");
    smenu = CAR(args);
    if (!isString(smenu) || length(smenu) != 1)
        error("invalid `menuname' argument");
    sitem = CADR(args);
    if (isNull(sitem))
    { /* add a menu */
        res = winaddmenu(CHAR(STRING_ELT(smenu, 0)), errmsg);
        if (res > 0)
        {
            sprintf(msgbuf, "unable to add menu (%s)", errmsg);
            errorcall(call, msgbuf);
        }
    }
    else
    { /* add an item */
        if (!isString(sitem) || length(sitem) != 1)
            error("invalid `itemname' argument");
        res = winaddmenuitem(CHAR(STRING_ELT(sitem, 0)), CHAR(STRING_ELT(smenu, 0)), CHAR(STRING_ELT(CADDR(args), 0)),
                             errmsg);
        if (res > 0)
        {
            sprintf(msgbuf, "unable to add menu item (%s)", errmsg);
            errorcall(call, msgbuf);
        }
    }
    return (R_NilValue);
}

SEXP do_winmenudel(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP smenu, sitem;
    int res;
    char errmsg[50];

    checkArity(op, args);
    if (CharacterMode != RGui)
        errorcall(call, "Menu functions can only be used in the GUI");
    smenu = CAR(args);
    if (!isString(smenu) || length(smenu) != 1)
        error("invalid `menuname' argument");
    sitem = CADR(args);
    if (isNull(sitem))
    { /* delete a menu */
        res = windelmenu(CHAR(STRING_ELT(smenu, 0)), errmsg);
        if (res > 0)
            errorcall(call, "menu does not exist");
    }
    else
    { /* delete an item */
        if (!isString(sitem) || length(sitem) != 1)
            error("invalid `itemname' argument");
        res = windelmenuitem(CHAR(STRING_ELT(sitem, 0)), CHAR(STRING_ELT(smenu, 0)), errmsg);
        if (res > 0)
        {
            sprintf(msgbuf, "unable to delete menu item (%s)", errmsg);
            errorcall(call, msgbuf);
        }
    }
    return (R_NilValue);
}

void Rwin_fpset()
{
    _fpreset();
    _controlfp(_MCW_EM, _MCW_EM);
}

#include "getline/getline.h" /* for gl_load/savehistory */
SEXP do_savehistory(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP sfile;

    checkArity(op, args);
    sfile = CAR(args);
    if (!isString(sfile) || LENGTH(sfile) < 1)
        errorcall(call, "invalid file argument");
    if (CharacterMode == RGui || (R_Interactive && CharacterMode == RTerm))
        gl_savehistory(CHAR(STRING_ELT(sfile, 0)));
    else
        errorcall(call, "savehistory can only be used in Rgui and Rterm");
    return R_NilValue;
}

SEXP do_loadhistory(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP sfile;

    checkArity(op, args);
    sfile = CAR(args);
    if (!isString(sfile) || LENGTH(sfile) < 1)
        errorcall(call, "invalid file argument");
    if (CharacterMode == RGui || (R_Interactive && CharacterMode == RTerm))
        gl_loadhistory(CHAR(STRING_ELT(sfile, 0)));
    else
        errorcall(call, "savehistory can only be used in Rgui and Rterm");
    return R_NilValue;
}

#include <lmcons.h>

SEXP do_sysinfo(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP ans, ansnames;
    OSVERSIONINFO verinfo;
    char isNT[8] = "??", ver[256], name[MAX_COMPUTERNAME_LENGTH + 1], user[UNLEN + 1];
    DWORD namelen = MAX_COMPUTERNAME_LENGTH + 1, userlen = UNLEN + 1;

    checkArity(op, args);
    PROTECT(ans = allocVector(STRSXP, 7));
    verinfo.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
    GetVersionEx(&verinfo);
    switch (verinfo.dwPlatformId)
    {
    case VER_PLATFORM_WIN32_NT:
        strcpy(isNT, "NT");
        break;
    case VER_PLATFORM_WIN32_WINDOWS:
        strcpy(isNT, "9x");
        break;
    case VER_PLATFORM_WIN32s:
        strcpy(isNT, "win32s");
        break;
    default:
        sprintf(isNT, "ID=%d", (int)verinfo.dwPlatformId);
        break;
    }

    SET_STRING_ELT(ans, 0, mkChar("Windows"));
    sprintf(ver, "%s %d.%d", isNT, (int)verinfo.dwMajorVersion, (int)verinfo.dwMinorVersion);
    SET_STRING_ELT(ans, 1, mkChar(ver));
    sprintf(ver, "(build %d) %s", LOWORD(verinfo.dwBuildNumber), verinfo.szCSDVersion);
    SET_STRING_ELT(ans, 2, mkChar(ver));
    GetComputerName(name, &namelen);
    SET_STRING_ELT(ans, 3, mkChar(name));
    SET_STRING_ELT(ans, 4, mkChar("x86"));
    GetUserName(user, &userlen);
    SET_STRING_ELT(ans, 5, mkChar(user));
    SET_STRING_ELT(ans, 6, STRING_ELT(ans, 5));
    PROTECT(ansnames = allocVector(STRSXP, 7));
    SET_STRING_ELT(ansnames, 0, mkChar("sysname"));
    SET_STRING_ELT(ansnames, 1, mkChar("release"));
    SET_STRING_ELT(ansnames, 2, mkChar("version"));
    SET_STRING_ELT(ansnames, 3, mkChar("nodename"));
    SET_STRING_ELT(ansnames, 4, mkChar("machine"));
    SET_STRING_ELT(ansnames, 5, mkChar("login"));
    SET_STRING_ELT(ansnames, 6, mkChar("user"));
    setAttrib(ans, R_NamesSymbol, ansnames);
    UNPROTECT(2);
    return ans;
}

void R_ProcessEvents(void); /* from system.c */

SEXP do_syssleep(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    DWORD mtime;
    int ntime;
    double time;

    checkArity(op, args);
    time = asReal(CAR(args));
    if (ISNAN(time) || time < 0)
        errorcall(call, "invalid time value");
    ntime = 1000 * (time) + 0.5;
    while (ntime > 0)
    {
        mtime = min(500, ntime);
        ntime -= mtime;
        Sleep(mtime);
        R_ProcessEvents();
    }
    return R_NilValue;
}
