/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998, 1999  Guido Masarotto
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

#include "Defn.h"
#include "Print.h"
#include "Fileio.h"
#include "IOStuff.h"
#include "Parse.h"
#include <stdio.h>
#include "run.h"

/*
 * ed, vi etc have 3 parameters. the data, a file and an editor
 *
 * if file is specified then the given file is used (and not removed on exit) if
 * file is not specified then a temporary file is used; since only one
 * temporary file is used for an entire session previous editing is lost
 *
 * if data is specified then it is passed out to be edited; if data is not
 * specified then either file (if specified) or the temporary file is used
 * (thus errors can be re-editied by calling edit a second time with no
 * arguments).
 *
 * if the editor is specified then the specified editor is invoked if possible
 * and an error message reported otherwise
 */

static char DefaultFileName[MAX_PATH];

void InitEd()
{
    char *tmp;

    tmp = getenv("TMP");
    if (!tmp)
        tmp = getenv("TEMP");
    if (!tmp)
        getenv("R_USER");
    sprintf(DefaultFileName, "%s/XXXXXX", tmp);
    mktemp(DefaultFileName);
}

SEXP do_edit(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    int i, rc, status;
    SEXP x, fn, envir, ed;
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
        error("invalid argument to edit()\n");

    if (LENGTH(STRING(fn)[0]) > 0)
    {
        filename = R_alloc(strlen(CHAR(STRING(fn)[0])), sizeof(char));
        strcpy(filename, CHAR(STRING(fn)[0]));
    }
    else
        filename = DefaultFileName;

    if (x != R_NilValue)
    {
        if ((fp = R_fopen(filename, "w")) == NULL)
            errorcall(call, "unable to open file\n");
        x = deparse1(x, 0);
        for (i = 0; i < LENGTH(x); i++)
            fprintf(fp, "%s\n", CHAR(STRING(x)[i]));
        fclose(fp);
    }
    ed = CAR(CDDR(args));
    if (!isString(ed))
        error("editor type not valid\n");
    editcmd = R_alloc(strlen(CHAR(STRING(ed)[0])) + strlen(filename) + 2, sizeof(char));
    sprintf(editcmd, "%s %s", CHAR(STRING(ed)[0]), filename);
    rc = runcmd(editcmd, 1, 1, "");
    if (rc == NOLAUNCH)
        errorcall(call, "unable to run editor\n");
    if (rc != 0)
        warningcall(call, "editor ran but returned error status\n");

    if ((fp = R_fopen(filename, "r")) == NULL)
        errorcall(call, "unable to open file to read\n");
    R_ParseCnt = 0;
    x = PROTECT(R_ParseFile(fp, -1, &status));
    if (status != PARSE_OK)
        errorcall(call, "An error occurred on line %d\n use a command like\n x <- vi()\n to recover\n", R_ParseError);
    else
        fclose(fp);
    R_ResetConsole();
    { /* can't just eval(x) here */
        int i, n;
        SEXP tmp;

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
