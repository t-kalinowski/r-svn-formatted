/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998, 2000, 2003  The R Development Core Team
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Defn.h"
#include "Print.h"
#include "Fileio.h"
#include "IOStuff.h"
#include "Parse.h"

#include "Runix.h"

#include <stdio.h>
#ifdef Win32
#include "run.h"
int Rgui_Edit(char *filename, char *title, int modal);
#endif

#ifdef HAVE_AQUA
extern DL_FUNC ptr_Raqua_Edit;

int Raqua_Edit(char *filename)
{
    ptr_Raqua_Edit(filename);
}
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h> /* for unlink() */
#endif

/*
 * ed, vi etc have 3 parameters. the data, a file and an editor
 *
 * If `file' is specified then the given file is used (and not removed on
 * exit). If `file' is not specified then a temporary file is used; since
 * only one temporary file is used for an entire session previous
 * editing is lost. That file is removed at the end of the R session.
 *
 * If `data' is specified then it is passed out to be edited; if `data' is not
 * specified then either `file' (if specified) or the temporary file is used
 * (thus errors can be re-edited by calling edit a second time with no
 * arguments).
 *
 * If the editor is specified then the specified editor is invoked if
 * possible and an error message reported otherwise
 */

static char *DefaultFileName;
static int EdFileUsed = 0;

void InitEd()
{
#ifdef Win32
    DefaultFileName = R_tmpnam("Redit", R_TempDir);
#else
    DefaultFileName = R_tmpnam(NULL, R_TempDir);
#endif
}

void CleanEd()
{
    if (EdFileUsed)
        unlink(DefaultFileName);
}

SEXP do_edit(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    int i, rc;
    ParseStatus status;
    SEXP x, fn, envir, ti, ed, t;
    char *filename, *title, *editcmd, *vmaxsave, *cmd;
    FILE *fp;

    checkArity(op, args);

    vmaxsave = vmaxget();

    x = CAR(args);
    args = CDR(args);
    if (TYPEOF(x) == CLOSXP)
        envir = CLOENV(x);
    else
        envir = R_NilValue;
    PROTECT(envir);

    fn = CAR(args);
    args = CDR(args);
    if (!isString(fn))
        error("invalid argument to edit()");

    if (LENGTH(STRING_ELT(fn, 0)) > 0)
    {
        filename = R_alloc(strlen(CHAR(STRING_ELT(fn, 0))), sizeof(char));
        strcpy(filename, CHAR(STRING_ELT(fn, 0)));
    }
    else
        filename = DefaultFileName;

    if (x != R_NilValue)
    {

        if ((fp = R_fopen(R_ExpandFileName(filename), "w")) == NULL)
            errorcall(call, "unable to open file");
        if (LENGTH(STRING_ELT(fn, 0)) == 0)
            EdFileUsed++;
        if (TYPEOF(x) != CLOSXP || isNull(t = getAttrib(x, R_SourceSymbol)))
            t = deparse1(x, 0);
        for (i = 0; i < LENGTH(t); i++)
            fprintf(fp, "%s\n", CHAR(STRING_ELT(t, i)));
        fclose(fp);
    }
    ti = CAR(args);
    args = CDR(args);
    ed = CAR(args);
    if (!isString(ed))
        errorcall(call, "argument `editor' type not valid");
    cmd = CHAR(STRING_ELT(ed, 0));
    if (strlen(cmd) == 0)
        errorcall(call, "argument `editor' is not set");
    editcmd = R_alloc(strlen(cmd) + strlen(filename) + 6, sizeof(char));
#ifdef Win32
    if (!strcmp(cmd, "internal"))
    {
        if (!isString(ti))
            error("title must be a string");
        if (LENGTH(STRING_ELT(ti, 0)) > 0)
        {
            title = R_alloc(strlen(CHAR(STRING_ELT(ti, 0))) + 1, sizeof(char));
            strcpy(title, CHAR(STRING_ELT(ti, 0)));
        }
        else
        {
            title = R_alloc(strlen(filename) + 1, sizeof(char));
            strcpy(title, filename);
        }
        Rgui_Edit(filename, title, 1);
    }
    else
    {
        /* Quote path if necessary */
        if (cmd[0] != '"' && strchr(cmd, ' '))
            sprintf(editcmd, "\"%s\" \"%s\"", cmd, filename);
        else
            sprintf(editcmd, "%s \"%s\"", cmd, filename);
        rc = runcmd(editcmd, 1, 1, "");
        if (rc == NOLAUNCH)
            errorcall(call, "unable to run editor %s", cmd);
        if (rc != 0)
            warningcall(call, "editor ran but returned error status");
    }
#else
#if defined(HAVE_AQUA)
    if (!strcmp(R_GUIType, "AQUA"))
        rc = Raqua_Edit(filename);
    else
    {
        sprintf(editcmd, "%s %s", cmd, filename);
        rc = R_system(editcmd);
    }
#else
    sprintf(editcmd, "%s %s", cmd, filename);
    rc = R_system(editcmd);
#endif
    if (rc != 0)
        errorcall(call, "problem with running editor %s", cmd);
#endif

    if ((fp = R_fopen(R_ExpandFileName(filename), "r")) == NULL)
        errorcall(call, "unable to open file to read");
    R_ParseCnt = 0;
    x = PROTECT(R_ParseFile(fp, -1, &status));
    fclose(fp);
    if (status != PARSE_OK)
        errorcall(call, "An error occurred on line %d\n use a command like\n x <- edit()\n to recover", R_ParseError);
    R_ResetConsole();
    { /* can't just eval(x) here */
        int j, n;
        SEXP tmp = R_NilValue;

        n = LENGTH(x);
        for (j = 0; j < n; j++)
            tmp = eval(VECTOR_ELT(x, j), R_GlobalEnv);
        x = tmp;
    }
    if (TYPEOF(x) == CLOSXP && envir != R_NilValue)
        SET_CLOENV(x, envir);
    UNPROTECT(2);
    vmaxset(vmaxsave);
    return (x);
}
