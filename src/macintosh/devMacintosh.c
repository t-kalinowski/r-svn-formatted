/*
 *  R : A Computer Language for Statistical Data Analysis
 *  file devMacintosh.c
 *  Copyright (C) 1998-1999  Ross Ihaka
 *                2000-2001  Stefano M. Iacus and the R core team
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

#include <RCarbon.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <Defn.h>
#include <Graphics.h>
#include <Rdevices.h>

static SEXP gcall;
static char *SaveString(SEXP sxp, int offset)
{
    char *s;
    if (!isString(sxp) || length(sxp) <= offset)
        errorcall(gcall, "invalid string argument");
    s = R_alloc(strlen(CHAR(STRING_ELT(sxp, offset))) + 1, sizeof(char));
    strcpy(s, CHAR(STRING_ELT(sxp, offset)));
    return s;
}

void InitEd()
{
}

/*  Macintosh Device Driver Parameters:
 *  -----------------		cf with ../unix/X11/devX11.c
 *  display	= display
 *  width	= width in inches
 *  height	= height in inches
 *  ps		= pointsize
 */
SEXP do_Macintosh(SEXP call, SEXP op, SEXP args, SEXP env)
{
    DevDesc *dd;
    char *display, *vmax;
    double height, width, ps;
    gcall = call;
    vmax = vmaxget();
    display = SaveString(CAR(args), 0);
    args = CDR(args);
    width = asReal(CAR(args));
    args = CDR(args);
    height = asReal(CAR(args));
    args = CDR(args);
    if (width <= 0 || height <= 0)
        errorcall(call, "invalid width or height");
    ps = asReal(CAR(args));

    R_CheckDeviceAvailable();
    /* Allocate and initialize the device driver data */
    BEGIN_SUSPEND_INTERRUPTS
    {
        if (!(dd = (DevDesc *)malloc(sizeof(DevDesc))))
            return 0;
        /* Do this for early redraw attempts */
        dd->displayList = R_NilValue;
        GInit(&dd->dp);
        if (!MacDeviceDriver(dd, width, height, ps))
        {
            free(dd);
            errorcall(call, "unable to start device Macintosh\n");
        }
        gsetVar(install(".Device"), mkString("Macintosh"), R_NilValue);
        addDevice(dd);
        initDisplayList(dd);
    }
    END_SUSPEND_INTERRUPTS;
    vmaxset(vmax);
    return R_NilValue;
}
