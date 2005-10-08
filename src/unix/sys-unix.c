/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997--2004  Robert Gentleman, Ross Ihaka
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

/* <UTF8>
   char here is mainly handled as a whole string.
   Does handle file names.
   Chopping final \n is OK in UTF-8.
 */

/* See system.txt for a description of functions */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Defn.h"
#include "Fileio.h"
#include "Runix.h"
#include <sys/stat.h> /* for mkdir */

/* HP-UX headers need this before CLK_TCK */
#ifdef HAVE_UNISTD_H
#include <unistd.h>
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

char *R_ExpandFileName_readline(char *s, char *buff); /* sys-std.c */

static char newFileName[PATH_MAX];
static int HaveHOME = -1;
static char UserHOME[PATH_MAX];

/* Only interpret inputs of the form ~ and ~/... */
static char *R_ExpandFileName_unix(char *s, char *buff)
{
    char *p;

    if (s[0] != '~')
        return s;
    if (strlen(s) > 1 && s[1] != '/')
        return s;
    if (HaveHOME < 0)
    {
        p = getenv("HOME");
        if (p && strlen(p) && (strlen(p) < PATH_MAX))
        {
            strcpy(UserHOME, p);
            HaveHOME = 1;
        }
        else
            HaveHOME = 0;
    }
    if (HaveHOME > 0 && (strlen(UserHOME) + strlen(s + 1) < PATH_MAX))
    {
        strcpy(buff, UserHOME);
        strcat(buff, s + 1);
        return buff;
    }
    else
        return s;
}

/* tilde_expand (in libreadline) mallocs storage for its return value.
   The R entry point does not require that storage to be freed, so we
   copy the value to a static buffer, to void a memory leak in R<=1.6.0.

   This is not thread-safe, but as R_ExpandFileName is a public entry
   point (in R-exts.texi) it will need to deprecated and replaced by a
   version which takes a buffer as an argument.

   BDR 10/2002
*/

extern Rboolean UsingReadline;

char *R_ExpandFileName(char *s)
{
#ifdef HAVE_LIBREADLINE
    if (UsingReadline)
    {
        char *c = R_ExpandFileName_readline(s, newFileName);
        /* we can return the result only if tilde_expand is not broken */
        if (!c || c[0] != '~' || (c[1] != '\0' && c[1] != '/'))
            return c;
    }
#endif
    return R_ExpandFileName_unix(s, newFileName);
}

/*
 *  7) PLATFORM DEPENDENT FUNCTIONS
 */

SEXP do_machine(SEXP call, SEXP op, SEXP args, SEXP env)
{
    return mkString("Unix");
}

#ifdef _R_HAVE_TIMING_
#include <time.h>
#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif
#ifndef CLK_TCK
/* this is in ticks/second, generally 60 on BSD style Unix, 100? on SysV
 */
#ifdef HZ
#define CLK_TCK HZ
#else
#define CLK_TCK 60
#endif
#endif /* not CLK_TCK */

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
#else  /* not _R_HAVE_TIMING_ */
void R_setStartTime(void)
{
}

SEXP do_proctime(SEXP call, SEXP op, SEXP args, SEXP env)
{
    error(_("proc.time() is not implemented on this system"));
    return R_NilValue; /* -Wall */
}
#endif /* not _R_HAVE_TIMING_ */

#define INTERN_BUFSIZE 8096
SEXP do_system(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    FILE *fp;
    char *x = "r", buf[INTERN_BUFSIZE];
    int read = 0, i, j;
    SEXP tlist = R_NilValue, tchar, rval;

    checkArity(op, args);
    if (!isValidStringF(CAR(args)))
        errorcall(call, _("non-empty character argument expected"));
    if (isLogical(CADR(args)))
        read = INTEGER(CADR(args))[0];
    if (read)
    {
#ifdef HAVE_POPEN
        PROTECT(tlist);
        fp = R_popen(CHAR(STRING_ELT(CAR(args), 0)), x);
        for (i = 0; fgets(buf, INTERN_BUFSIZE, fp); i++)
        {
            read = strlen(buf);
            if (read > 0 && buf[read - 1] == '\n')
                buf[read - 1] = '\0'; /* chop final CR */
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
#else  /* not HAVE_POPEN */
        errorcall(call, _("intern=TRUE is not implemented on this platform"));
        return R_NilValue;
#endif /* not HAVE_POPEN */
    }
    else
    {
#ifdef HAVE_AQUA
        R_Busy(1);
#endif
        tlist = allocVector(INTSXP, 1);
        fflush(stdout);
        INTEGER(tlist)[0] = R_system(CHAR(STRING_ELT(CAR(args), 0)));
#ifdef HAVE_AQUA
        R_Busy(0);
#endif
        R_Visible = 0;
        return tlist;
    }
}

void InitTempDir()
{
    char *tmp, *tm, tmp1[PATH_MAX + 10], *p;
    int len;
#ifndef HAVE_MKDTEMP
    int res;
#endif

    tmp = getenv("R_SESSION_TMPDIR");
    if (!tmp)
    {
        /* This looks like it will only be called in the embedded case
           since this is done in the script. Also should test if directory
           exists rather than just attempting to remove it.
    */
        char *buf;
        tm = getenv("TMPDIR");
        if (!tm)
            tm = getenv("TMP");
        if (!tm)
            tm = getenv("TEMP");
        if (!tm)
            tm = "/tmp";
#ifdef HAVE_MKDTEMP
        sprintf(tmp1, "%s/RtmpXXXXXX", tm);
        tmp = mkdtemp(tmp1);
        if (!tmp)
            R_Suicide(_("cannot mkdir R_TempDir"));
#else
        sprintf(tmp1, "rm -rf %s/Rtmp%u", tm, (unsigned int)getpid());
        R_system(tmp1);
        sprintf(tmp1, "%s/Rtmp%u", tm, (unsigned int)getpid());
        res = mkdir(tmp1, 0755);
        if (res)
        {
            /* Try one more time, in case a dir left around from that
               process number from another user */
            sprintf(tmp1, "rm -rf %s/Rtmp%u-%d", tm, (unsigned int)getpid(), rand() % 1000);
            R_system(tmp1);
            sprintf(tmp1, "%s/Rtmp%u-%d", tm, (unsigned int)getpid(), rand() % 1000);
            res = mkdir(tmp1, 0755);
        }
        if (res)
            R_Suicide(_("cannot mkdir R_TempDir"));
#endif
        tmp = tmp1;
        buf = (char *)malloc((strlen(tmp) + 20) * sizeof(char));
        if (buf)
        {
            sprintf(buf, "R_SESSION_TMPDIR=%s", tmp);
            putenv(buf);
            /* no free here: storage remains in use */
        }
    }

    len = strlen(tmp) + 1;
    p = (char *)malloc(len);
    if (!p)
        R_Suicide(_("cannot allocate R_TempDir"));
    else
    {
        R_TempDir = p;
        strcpy(R_TempDir, tmp);
    }
}

char *R_tmpnam(const char *prefix, const char *tempdir)
{
    char tm[PATH_MAX], tmp1[PATH_MAX], *res;
    unsigned int n, done = 0;

    if (!prefix)
        prefix = ""; /* NULL */
    strcpy(tmp1, tempdir);
    for (n = 0; n < 100; n++)
    {
        /* try a random number at the end */
        sprintf(tm, "%s/%s%x", tmp1, prefix, rand());
        if (!R_FileExists(tm))
        {
            done = 1;
            break;
        }
    }
    if (!done)
        error(_("cannot find unused tempfile name"));
    res = (char *)malloc((strlen(tm) + 1) * sizeof(char));
    strcpy(res, tm);
    return res;
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
#else  /* not HAVE_SYS_UTSNAME_H */
SEXP do_sysinfo(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    warning(_("Sys.info() is not implemented on this system"));
    return R_NilValue; /* -Wall */
}
#endif /* not HAVE_SYS_UTSNAME_H */

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

/* patch from Ei-ji Nakama for Intel compilers on ix86.
   From http://www.nakama.ne.jp/memo/ia32_linux/R-2.1.1.iccftzdaz.patch.txt.
   Since updated to include x86_64.
 */
#if (defined(__i386) || defined(__x86_64)) && defined(__INTEL_COMPILER) && __INTEL_COMPILER > 800
#include <xmmintrin.h>
#include <pmmintrin.h>
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
#if (defined(__i386) || defined(__x86_64)) && defined(__INTEL_COMPILER) && __INTEL_COMPILER > 800
        _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_OFF);
        _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_OFF);
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
