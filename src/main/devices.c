/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998	Robert Gentleman and Ross Ihaka
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
#include "Mathlib.h"
#include "Graphics.h"

static SEXP gcall;

/* Return a non-relocatable copy of a string */

static char *SaveString(SEXP sxp, int offset)
{
    char *s;
    if (!isString(sxp) || length(sxp) <= offset)
        errorcall(gcall, "invalid string argument\n");
    s = R_alloc(strlen(CHAR(STRING(sxp)[offset])) + 1, sizeof(char));
    strcpy(s, CHAR(STRING(sxp)[offset]));
    return s;
}

static void DeviceUnavailable(char *dev)
{
    error("The %s device driver is unavailable.\n", dev);
}

#ifdef Macintosh
/*  Macintosh Device Driver Parameters:
 *  -----------------		--> ../unix/devX11.c
 *  display	= display
 *  width	= width in inches
 *  height	= height in inches
 *  ps		= pointsize
 */
int MacDeviceDriver();

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
    /* Allocate and initialize the device driver data */
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
    vmaxset(vmax);
    return R_NilValue;
}
#else
SEXP do_Macintosh(SEXP call, SEXP op, SEXP args, SEXP env)
{
    gcall = call;
    DeviceUnavailable("Macintosh");
}
#endif

/*  PostScript Device Driver Parameters:
 *  ------------------------		--> ../unix/devPS.c
 *  file	= output filename
 *  paper	= paper type
 *  face	= typeface = "family"
 *  bg		= background color
 *  fg		= foreground color
 *  width	= width in inches
 *  height	= height in inches
 *  horizontal	= {TRUE: landscape; FALSE: portrait}
 *  ps		= pointsize
 */

SEXP do_PS(SEXP call, SEXP op, SEXP args, SEXP env)
{
    DevDesc *dd;
    char *vmax;
    char *file, *paper, *face, *bg, *fg;
    int horizontal;
    double height, width, ps;
    gcall = call;
    vmax = vmaxget();
    file = SaveString(CAR(args), 0);
    args = CDR(args);
    paper = SaveString(CAR(args), 0);
    args = CDR(args);
    face = SaveString(CAR(args), 0);
    args = CDR(args);
    bg = SaveString(CAR(args), 0);
    args = CDR(args);
    fg = SaveString(CAR(args), 0);
    args = CDR(args);
    width = asReal(CAR(args));
    args = CDR(args);
    height = asReal(CAR(args));
    args = CDR(args);
    horizontal = asLogical(CAR(args));
    args = CDR(args);
    if (horizontal == NA_LOGICAL)
        horizontal = 1;
    ps = asReal(CAR(args));

    if (!(dd = (DevDesc *)malloc(sizeof(DevDesc))))
        return 0;
    /* Do this for early redraw attempts */
    dd->displayList = R_NilValue;
    GInit(&dd->dp);
    if (!PSDeviceDriver(dd, file, paper, face, bg, fg, width, height, (double)horizontal, ps))
    {
        free(dd);
        errorcall(call, "unable to start device PostScript\n");
    }
    gsetVar(install(".Device"), mkString("postscript"), R_NilValue);
    addDevice(dd);
    initDisplayList(dd);
    vmaxset(vmax);
    return R_NilValue;
}

/*  PicTeX Device Driver Parameters
 *  --------------------		--> ../unix/devPicTeX.c
 *  file    = output filename
 *  bg	    = background color
 *  fg	    = foreground color
 *  width   = width in inches
 *  height  = height in inches
 *  debug   = int; if non-0, write TeX-Comments into output.
 */

SEXP do_PicTeX(SEXP call, SEXP op, SEXP args, SEXP env)
{
    DevDesc *dd;
    char *vmax;
    char *file, *bg, *fg;
    double height, width;
    int debug;
    gcall = call;
    vmax = vmaxget();
    file = SaveString(CAR(args), 0);
    args = CDR(args);
    bg = SaveString(CAR(args), 0);
    args = CDR(args);
    fg = SaveString(CAR(args), 0);
    args = CDR(args);
    width = asReal(CAR(args));
    args = CDR(args);
    height = asReal(CAR(args));
    args = CDR(args);
    debug = asInteger(CAR(args));
    args = CDR(args);
    if (!(dd = (DevDesc *)malloc(sizeof(DevDesc))))
        return 0;
    /* Do this for early redraw attempts */
    dd->displayList = R_NilValue;
    GInit(&dd->dp);
    if (!PicTeXDeviceDriver(dd, file, bg, fg, width, height, debug))
    {
        free(dd);
        errorcall(call, "unable to start device PicTeX\n");
    }
    gsetVar(install(".Device"), mkString("pictex"), R_NilValue);
    addDevice(dd);
    initDisplayList(dd);
    vmaxset(vmax);
    return R_NilValue;
}
#ifdef Unix
#include "../unix/devX11.h"
/*  X11 Device Driver Parameters:
 *  -----------------		--> ../unix/devX11.c
 *  display	= display
 *  width	= width in inches
 *  height	= height in inches
 *  ps		= pointsize
 */
SEXP do_X11(SEXP call, SEXP op, SEXP args, SEXP env)
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
    if (!X11DeviceDriver(dd, display, width, height, ps))
    {
        free(dd);
        errorcall(call, "unable to start device X11\n");
    }
    gsetVar(install(".Device"), mkString("X11"), R_NilValue);
    addDevice(dd);
    initDisplayList(dd);
    vmaxset(vmax);
    return R_NilValue;
}
#else
SEXP do_x11(SEXP call, SEXP op, SEXP args, SEXP env)
{
    gcall = call;
    DeviceUnavailable("X11");
}
#endif
