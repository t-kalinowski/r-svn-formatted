/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997-2000   Robert Gentleman, Ross Ihaka
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
#include <config.h>
#endif

#include "Defn.h"
#include "Fileio.h"
#include "Runix.h"

/* HP-UX headers need this before CLK_TCK */
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_LIBREADLINE
#ifdef HAVE_READLINE_READLINE_H
#include <readline/readline.h>
#endif
#ifdef HAVE_READLINE_HISTORY_H
#include <readline/history.h>
#endif
#endif

extern Rboolean LoadInitFile;

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

SEXP do_machine(SEXP call, SEXP op, SEXP args, SEXP env)
{
    return mkString("Unix");
}

#ifdef HAVE_TIMES
#include <time.h>
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

void R_setStartTime(void)
{
    StartTime = times(&timeinfo);
}

void R_getProcTime(double *data)
{
    double elapsed;
    elapsed = (times(&timeinfo) - StartTime) / (double)CLK_TCK;
    data[0] = timeinfo.tms_utime / (double)CLK_TCK;
    data[1] = timeinfo.tms_stime / (double)CLK_TCK;
    data[2] = elapsed;
    data[3] = timeinfo.tms_cutime / (double)CLK_TCK;
    data[4] = timeinfo.tms_cstime / (double)CLK_TCK;
}

double R_getClockIncrement(void)
{
    return 1.0 / (double)CLK_TCK;
}

SEXP do_proctime(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans = allocVector(REALSXP, 5);
    R_getProcTime(REAL(ans));
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
    if (!isValidStringF(CAR(args)))
        errorcall(call, "non-empty character argument expected");
    if (isLogical(CADR(args)))
        read = INTEGER(CADR(args))[0];
    if (read)
    {
#ifdef HAVE_POPEN
        PROTECT(tlist);
        fp = popen(CHAR(STRING_ELT(CAR(args), 0)), x);
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
            SET_STRING_ELT(rval, j, CAR(tlist));
            tlist = CDR(tlist);
        }
        UNPROTECT(1);
        return (rval);
#else
        errorcall(call, "intern=TRUE is not implemented on this platform");
        return R_NilValue;
#endif
    }
    else
    {
        tlist = allocVector(INTSXP, 1);
        fflush(stdout);
        INTEGER(tlist)[0] = system(CHAR(STRING_ELT(CAR(args), 0)));
        R_Visible = 0;
        return tlist;
    }
}

static char *Runix_tmpnam(char *prefix)
{
    char *tmp, tm[PATH_MAX], tmp1[PATH_MAX], *res;
    unsigned int n, done = 0, pid;

    tmp = getenv("TMP");
    if (!tmp)
        tmp = getenv("TEMP");
    if (tmp)
        strcpy(tmp1, tmp);
    else
        strcpy(tmp1, "/tmp");
    pid = (unsigned int)getpid();
    for (n = 0; n < 100; n++)
    {
        /* try a random number at the end */
        sprintf(tm, "%s/%sR%xS%x", tmp1, prefix, pid, rand());
        if (!R_FileExists(tm))
        {
            done = 1;
            break;
        }
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
        tm = Runix_tmpnam(tn);
        SET_STRING_ELT(ans, i, mkChar(tm));
        free(tm);
    }
    UNPROTECT(1);
    return (ans);
}

#ifdef HAVE_SYS_UTSNAME_H
#include <sys/utsname.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_PWD_H
#include <pwd.h>
#endif

SEXP do_sysinfo(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP ans, ansnames;
    struct utsname name;
    char *login;

    checkArity(op, args);
    PROTECT(ans = allocVector(STRSXP, 7));
    if (uname(&name) == -1)
    {
        UNPROTECT(1);
        return R_NilValue;
    }
    SET_STRING_ELT(ans, 0, mkChar(name.sysname));
    SET_STRING_ELT(ans, 1, mkChar(name.release));
    SET_STRING_ELT(ans, 2, mkChar(name.version));
    SET_STRING_ELT(ans, 3, mkChar(name.nodename));
    SET_STRING_ELT(ans, 4, mkChar(name.machine));
    login = getlogin();
    SET_STRING_ELT(ans, 5, login ? mkChar(login) : mkChar("unknown"));
#if defined(HAVE_PWD_H) && defined(HAVE_GETPWUID) && defined(HAVE_GETUID)
    {
        struct passwd *stpwd;
        stpwd = getpwuid(getuid());
        SET_STRING_ELT(ans, 6, stpwd ? mkChar(stpwd->pw_name) : mkChar("unknown"));
    }
#else
    SET_STRING_ELT(ans, 6, mkChar("unknown"));
#endif
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
#else
SEXP do_sysinfo(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    warning("Sys,info is not implemented on this system");
    return R_NilValue; /* -Wall */
}
#endif

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

void fpu_setup(Rboolean start)
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
