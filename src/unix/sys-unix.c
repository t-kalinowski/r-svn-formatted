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
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* See system.txt for a description of functions */

#ifdef HAVE_CONFIG_H
#include <Rconfig.h>
#endif

#include "Defn.h"
#include "Fileio.h"

#ifdef HAVE_LIBREADLINE
#ifdef HAVE_READLINE_READLINE_H
#include <readline/readline.h>
#endif
#ifdef HAVE_READLINE_HISTORY_H
#include <readline/history.h>
#endif
#endif

extern int LoadInitFile;

/*
 *  4) INITIALIZATION AND TERMINATION ACTIONS
 */

FILE *R_OpenInitFile(void)
{
    char buf[256], *home;
    FILE *fp;

    fp = NULL;
    if (LoadInitFile)
    {
        if ((fp = R_fopen(".Rprofile", "r")))
            return fp;
        if ((home = getenv("HOME")) == NULL)
            return NULL;
        sprintf(buf, "%s/.Rprofile", home);
        if ((fp = R_fopen(buf, "r")))
            return fp;
    }
    return fp;
}
/*
 *   R_CleanUp is interface-specific
 */

/*
 *  5) FILESYSTEM INTERACTION
 */

/*
 *   R_ShowFiles is interface-specific
 */

/*
 *   R_ChooseFile is interface-specific
 */

#ifdef HAVE_LIBREADLINE
char *tilde_expand(char *);

char *R_ExpandFileName(char *s)
{
    return (tilde_expand(s));
}
#else
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
        p = getenv("HOME");
        if (p && strlen(p))
        {
            strcpy(UserHOME, p);
            HaveHOME = 1;
        }
        else
            HaveHOME = 0;
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
#endif

/*
 *  7) PLATFORM DEPENDENT FUNCTIONS
 */

extern char **environ;

SEXP do_getenv(SEXP call, SEXP op, SEXP args, SEXP env)
{
    int i, j;
    char *s;
    char **e;
    SEXP ans;

    checkArity(op, args);

    if (!isString(CAR(args)))
        errorcall(call, "wrong type for argument\n");

    i = LENGTH(CAR(args));
    if (i == 0)
    {
        for (i = 0, e = environ; *e != NULL; i++, e++)
            ;
        PROTECT(ans = allocVector(STRSXP, i));
        for (i = 0, e = environ; *e != NULL; i++, e++)
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
    if (R_Interactive)
        LOGICAL(rval)[0] = 1;
    else
        LOGICAL(rval)[0] = 0;
    return rval;
}

SEXP do_machine(SEXP call, SEXP op, SEXP args, SEXP env)
{
    return mkString("Unix");
}

#ifdef HAVE_TIMES
#include <sys/times.h>
#ifndef CLK_TCK
/* this is in ticks/second, generally 60 on BSD style Unix, 100? on SysV */
#ifdef HZ
#define CLK_TCK HZ
#else
#define CLK_TCK 60
#endif
#endif /* CLK_TCK */

static clock_t StartTime;
static struct tms timeinfo;

void setStartTime()
{
    StartTime = times(&timeinfo);
}

SEXP do_proctime(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans;
    clock_t elapsed;
    elapsed = (times(&timeinfo) - StartTime) / (double)CLK_TCK;
    ans = allocVector(REALSXP, 5);
    REAL(ans)[0] = timeinfo.tms_utime / (double)CLK_TCK;
    REAL(ans)[1] = timeinfo.tms_stime / (double)CLK_TCK;
    REAL(ans)[2] = elapsed;
    REAL(ans)[3] = timeinfo.tms_cutime / (double)CLK_TCK;
    REAL(ans)[4] = timeinfo.tms_cstime / (double)CLK_TCK;
    return ans;
}
#endif /* HAVE_TIMES */

SEXP do_system(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    FILE *fp;
    char *x = "r", buf[120];
    int read = 0, i, j;
    SEXP tlist = R_NilValue, tchar, rval;

    checkArity(op, args);
    if (!isString(CAR(args)))
        errorcall(call, "character argument expected\n");
    if (isLogical(CADR(args)))
        read = INTEGER(CADR(args))[0];
    if (read)
    {
        PROTECT(tlist);
        fp = popen(CHAR(STRING(CAR(args))[0]), x);
        for (i = 0; fgets(buf, 120, fp); i++)
        {
            read = strlen(buf);
            buf[read - 1] = '\0';
            tchar = mkChar(buf);
            UNPROTECT(1);
            PROTECT(tlist = CONS(tchar, tlist));
        }
        pclose(fp);
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
        fflush(stdout);
        INTEGER(tlist)[0] = system(CHAR(STRING(CAR(args))[0]));
        R_Visible = 0;
        return tlist;
    }
}

/*
 *  helpers for start-up code
 */

#ifdef __FreeBSD__
#ifdef HAVE_FLOATINGPOINT_H
#include <floatingpoint.h>
#endif
#endif

#ifdef linux
#ifdef HAVE_FPU_CONTROL_H
#include <fpu_control.h>
#endif
#endif

void fpu_setup(int start)
{
    if (start)
    {
#ifdef __FreeBSD__
        fpsetmask(0);
#endif

#ifdef NEED___SETFPUCW
        __setfpucw(_FPU_IEEE);
#endif
    }
    else
    {
#ifdef __FreeBSD__
        fpsetmask(~0);
#endif

#ifdef NEED___SETFPUCW
        __setfpucw(_FPU_DEFAULT);
#endif
    }
}
