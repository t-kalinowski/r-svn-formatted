/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998	Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998-2000   The R Development Core Team.
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
#include "Mathlib.h"
#include "Graphics.h"

/* Return a non-relocatable copy of a string */

static SEXP gcall;
static char *SaveString(SEXP sxp, int offset)
{
    char *s;
    if (!isString(sxp) || length(sxp) <= offset)
        errorcall(gcall, "invalid string argument");
    s = R_alloc(strlen(CHAR(STRING(sxp)[offset])) + 1, sizeof(char));
    strcpy(s, CHAR(STRING(sxp)[offset]));
    return s;
}

#include "../unix/devX11.h"

/*  X11 Device Driver Parameters:
 *  -----------------		--> ../unix/devX11.c
 *  display	= display
 *  width	= width in inches
 *  height	= height in inches
 *  ps		= pointsize
 *  gamma       = gamma correction
 *  colormodel  = color model
 */

SEXP do_X11(SEXP call, SEXP op, SEXP args, SEXP env)
{
    DevDesc *dd;
    char *display, *vmax, *cname;
    double height, width, ps, gamma;
    int colormodel, maxcubesize;
    gcall = call;
    vmax = vmaxget();

    /* Decode the arguments */
    display = SaveString(CAR(args), 0);
    args = CDR(args);
    width = asReal(CAR(args));
    args = CDR(args);
    height = asReal(CAR(args));
    args = CDR(args);
    if (width <= 0 || height <= 0)
        errorcall(call, "invalid width or height");
    ps = asReal(CAR(args));
    args = CDR(args);
    gamma = asReal(CAR(args));
    args = CDR(args);
    if (gamma < 0 || gamma > 100)
        errorcall(call, "invalid gamma value");

    if (!isValidString(CAR(args)))
        error("invalid colortype passed to X11 driver");
    cname = CHAR(STRING(CAR(args))[0]);
    if (strcmp(cname, "mono") == 0)
        colormodel = 0;
    else if (strcmp(cname, "gray") == 0 || strcmp(cname, "grey") == 0)
        colormodel = 1;
    else if (strcmp(cname, "pseudo.cube") == 0)
        colormodel = 2;
    else if (strcmp(cname, "pseudo") == 0)
        colormodel = 3;
    else if (strcmp(cname, "true") == 0)
        colormodel = 4;
    else
    {
        warningcall(call, "unknown X11 color/colour model -- using monochrome");
        colormodel = 0;
    }
    args = CDR(args);
    maxcubesize = asInteger(CAR(args));
    if (maxcubesize < 1 || maxcubesize > 256)
        maxcubesize = 256;

    /* Allocate and initialize the device driver data */
    if (!(dd = (DevDesc *)malloc(sizeof(DevDesc))))
        return 0;
    /* Do this for early redraw attempts */
    dd->displayList = R_NilValue;
    GInit(&dd->dp);
    if (!X11DeviceDriver(dd, display, width, height, ps, gamma, colormodel, maxcubesize))
    {
        free(dd);
        errorcall(call, "unable to start device X11");
    }
    gsetVar(install(".Device"), mkString("X11"), R_NilValue);
    addDevice(dd);
    initDisplayList(dd);
    vmaxset(vmax);
    return R_NilValue;
}

SEXP do_Gnome(SEXP call, SEXP op, SEXP args, SEXP env)
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
    /* Allocate and initialize the device driver data */
    if (!(dd = (DevDesc *)malloc(sizeof(DevDesc))))
        return 0;
    /* Do this for early redraw attempts */
    dd->displayList = R_NilValue;
    GInit(&dd->dp);
    if (!GnomeDeviceDriver(dd, display, width, height, ps))
    {
        free(dd);
        errorcall(call, "unable to start device Gnome");
    }
    gsetVar(install(".Device"), mkString("Gnome"), R_NilValue);
    addDevice(dd);
    initDisplayList(dd);
    vmaxset(vmax);
    return R_NilValue;
}
