/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997-1999   Robert Gentleman, Ross Ihaka
 *                            and the R Development Core Team
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

/* See ../unix/system.txt for a description of functions */

/* Windows analogue of unix/sys-unix.c: often rather similar */

#ifdef HAVE_CONFIG_H
#include <Rconfig.h>
#endif

#include "Defn.h"
#include "Fileio.h"
#include "Startup.h"

extern int LoadInitFile;
extern UImode CharacterMode;

/*
 *  4) INITIALIZATION AND TERMINATION ACTIONS
 */

FILE *R_OpenInitFile(void)
{
    char buf[256];
    FILE *fp;

    fp = NULL;
    if (LoadInitFile)
    {
        if ((fp = R_fopen(".Rprofile", "r")))
            return fp;
        sprintf(buf, "%s/.Rprofile", getenv("R_USER"));
        if ((fp = R_fopen(buf, "r")))
            return fp;
    }
    return fp;
}
/*
 *  5) FILESYSTEM INTERACTION
 */

static int HaveHOME = -1;
static char UserHOME[PATH_MAX];
static char newFileName[PATH_MAX];
char *R_ExpandFileName(char *s)
{
    char *p;

    if (s[0] != '~')
        return s;
    if (HaveHOME < 0)
    {
        HaveHOME = 0;
        p = getenv("HOME");
        if (p && strlen(p))
        {
            strcpy(UserHOME, p);
            HaveHOME = 1;
        }
        else
        {
            p = getenv("HOMEDIR");
            if (p)
            {
                strcpy(UserHOME, p);
                p = getenv("HOMEPATH");
                if (p)
                {
                    strcat(UserHOME, p);
                    HaveHOME = 1;
                }
            }
        }
    }
    if (HaveHOME > 0)
    {
        strcpy(newFileName, UserHOME);
        strcat(newFileName, s + 1);
        return newFileName;
    }
    else
        return s;
}

/*
 *  7) PLATFORM DEPENDENT FUNCTIONS
 */

SEXP do_getenv(SEXP call, SEXP op, SEXP args, SEXP env)
{
    int i, j;
    char *s;
    char **e;
    SEXP ans;

    char *_env[1];

    _env[0] = NULL;
    checkArity(op, args);

    if (!isString(CAR(args)))
        errorcall(call, "wrong type for argument\n");

    i = LENGTH(CAR(args));
    if (i == 0)
    {
        for (i = 0, e = _env; *e != NULL; i++, e++)
            ;
        PROTECT(ans = allocVector(STRSXP, i));
        for (i = 0, e = _env; *e != NULL; i++, e++)
            STRING(ans)[i] = mkChar(*e);
    }
    else
    {
        PROTECT(ans = allocVector(STRSXP, i));
        for (j = 0; j < i; j++)
        {
            s = getenv(CHAR(STRING(CAR(args))[j]));
            if (s == NULL)
                STRING(ans)[j] = mkChar("");
            else
                STRING(ans)[j] = mkChar(s);
        }
    }
    UNPROTECT(1);
    return (ans);
}

SEXP do_interactive(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP rval;

    rval = allocVector(LGLSXP, 1);
    LOGICAL(rval)[0] = R_Interactive;
    return rval;
}

SEXP do_machine(SEXP call, SEXP op, SEXP args, SEXP env)
{
    return mkString("Win32");
}

#ifdef HAVE_TIMES

#include <windows.h> /* for DWORD, GetTickCount */
static DWORD StartTime;

void setStartTime(void)
{
    StartTime = GetTickCount();
}

SEXP do_proctime(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans;
    long elapsed;

    elapsed = (GetTickCount() - StartTime) / 10;
    ans = allocVector(REALSXP, 5);
    REAL(ans)[0] = R_NaReal;
    REAL(ans)[1] = R_NaReal;
    REAL(ans)[2] = (double)elapsed / 100.0;
    REAL(ans)[3] = R_NaReal;
    REAL(ans)[4] = R_NaReal;
    return ans;
}
#endif /* HAVE_TIMES */

/*
 * flag =0 don't wait/ignore stdout
 * flag =1 wait/ignore stdout
 * flag =2 wait/copy stdout to the console
 * flag =3 wait/return stdout
 * Add 10 to minimize application
 * Add 20 to make application "invisible"
 */

#include "run.h"

SEXP do_system(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    rpipe *fp;
    char buf[120];
    int vis = 0, flag = 2, i = 0, j, ll;
    SEXP tlist = R_NilValue, tchar, rval;

    checkArity(op, args);
    if (!isString(CAR(args)))
        errorcall(call, "character string expected as first argument\n");
    if (isInteger(CADR(args)))
        flag = INTEGER(CADR(args))[0];
    if (flag > 20)
    {
        vis = -1;
        flag -= 20;
    }
    else if (flag > 10)
    {
        vis = 0;
        flag -= 10;
    }
    else
        vis = 1;
    if (!isString(CADDR(args)))
        errorcall(call, "character string expected as third argument\n");
    if ((CharacterMode != RGui) && (flag == 2))
        flag = 1;
    if (CharacterMode == RGui)
    {
        SetStdHandle(STD_INPUT_HANDLE, INVALID_HANDLE_VALUE);
        SetStdHandle(STD_OUTPUT_HANDLE, INVALID_HANDLE_VALUE);
        SetStdHandle(STD_ERROR_HANDLE, INVALID_HANDLE_VALUE);
    }
    if (flag < 2)
    {
        ll = runcmd(CHAR(STRING(CAR(args))[0]), flag, vis, CHAR(STRING(CADDR(args))[0]));
        if (ll == NOLAUNCH)
            warning(runerror());
    }
    else
    {
        fp = rpipeOpen(CHAR(STRING(CAR(args))[0]), vis, CHAR(STRING(CADDR(args))[0]));
        if (!fp)
        {
            /* If we are returning standard output generate an error */
            if (flag == 3)
                error(runerror());
            warning(runerror());
            ll = NOLAUNCH;
        }
        else
        {
            if (flag == 3)
                PROTECT(tlist);
            for (i = 0; rpipeGets(fp, buf, 120); i++)
            {
                if (flag == 3)
                {
                    ll = strlen(buf) - 1;
                    if ((ll >= 0) && (buf[ll] == '\n'))
                        buf[ll] = '\0';
                    tchar = mkChar(buf);
                    UNPROTECT(1);
                    PROTECT(tlist = CONS(tchar, tlist));
                }
                else
                    R_WriteConsole(buf, strlen(buf));
            }
            ll = rpipeClose(fp);
        }
    }
    if (flag == 3)
    {
        rval = allocVector(STRSXP, i);
        ;
        for (j = (i - 1); j >= 0; j--)
        {
            STRING(rval)[j] = CAR(tlist);
            tlist = CDR(tlist);
        }
        UNPROTECT(1);
        return (rval);
    }
    else
    {
        tlist = allocVector(INTSXP, 1);
        INTEGER(tlist)[0] = ll;
        R_Visible = 0;
        return tlist;
    }
}
