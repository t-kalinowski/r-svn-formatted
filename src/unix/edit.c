/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998, 1999  The R Development Core Team
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
#include <stdio.h>
#ifdef Win32
#include "run.h"
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
    DefaultFileName = tmpnam(NULL);
}

void CleanEd()
{
    if (EdFileUsed)
        unlink(DefaultFileName);
}

SEXP do_edit(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    int i, rc, status;
    SEXP x, fn, envir, ed, t;
    char *filename, *editcmd, *vmaxsave;
    FILE *fp;

    checkArity(op, args);

    vmaxsave = vmaxget();

    x = CAR(args);
    if (TYPEOF(x) == CLOSXP)
        envir = CLOENV(x);
    else
        envir = R_NilValue;
    PROTECT(envir);

    fn = CADR(args);
    if (!isString(fn))
        error("invalid argument to edit()");

    if (LENGTH(STRING(fn)[0]) > 0)
    {
        filename = R_alloc(strlen(CHAR(STRING(fn)[0])), sizeof(char));
        strcpy(filename, CHAR(STRING(fn)[0]));
    }
    else
        filename = DefaultFileName;

    if (x != R_NilValue)
    {

        if ((fp = R_fopen(R_ExpandFileName(filename), "w")) == NULL)
            errorcall(call, "unable to open file");
        if (LENGTH(STRING(fn)[0]) == 0)
            EdFileUsed++;
        if (TYPEOF(x) != CLOSXP || isNull(t = getAttrib(x, R_SourceSymbol)))
            t = deparse1(x, 0);
        for (i = 0; i < LENGTH(t); i++)
            fprintf(fp, "%s\n", CHAR(STRING(t)[i]));
        fclose(fp);
    }

    ed = CAR(CDDR(args));
    if (!isString(ed))
        error("editor type not valid");
    editcmd = R_alloc(strlen(CHAR(STRING(ed)[0])) + strlen(filename) + 2, sizeof(char));
#ifdef Win32
    sprintf(editcmd, "%s \"%s\"", CHAR(STRING(ed)[0]), filename);
    rc = runcmd(editcmd, 1, 1, "");
    if (rc == NOLAUNCH)
        errorcall(call, "unable to run editor");
    if (rc != 0)
        warningcall(call, "editor ran but returned error status");
#else
    sprintf(editcmd, "%s %s", CHAR(STRING(ed)[0]), filename);
    rc = system(editcmd);
#endif

    if ((fp = R_fopen(R_ExpandFileName(filename), "r")) == NULL)
        errorcall(call, "unable to open file to read");
    R_ParseCnt = 0;
    x = PROTECT(R_ParseFile(fp, -1, &status));
    if (status != PARSE_OK)
        errorcall(call, "An error occurred on line %d\n use a command like\n x <- vi()\n to recover", R_ParseError);
    else
        fclose(fp);
    R_ResetConsole();
    { /* can't just eval(x) here */
        int i, n;
        SEXP tmp = R_NilValue;

        n = LENGTH(x);
        for (i = 0; i < n; i++)
            tmp = eval(VECTOR(x)[i], R_GlobalEnv);
        x = tmp;
    }
    if (TYPEOF(x) == CLOSXP && envir != R_NilValue)
        CLOENV(x) = envir;
    UNPROTECT(2);
    vmaxset(vmaxsave);
    return (x);
}
