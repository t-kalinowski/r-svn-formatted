/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2003  Guido Masarotto and Brian Ripley
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

#include <Defn.h>
#include <Graphics.h>
#include <Rdevices.h>
#include "devga.h"

/* Return a non-relocatable copy of a string */

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

/* This is Guido's devga device. */

SEXP do_devga(SEXP call, SEXP op, SEXP args, SEXP env)
{
    NewDevDesc *dev;
    GEDevDesc *dd;
    char *display, *vmax;
    double height, width, ps, xpinch, ypinch, gamma;
    int recording = 0, resize = 1, canvas, xpos, ypos;
    SEXP sc;

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
    args = CDR(args);
    recording = asLogical(CAR(args));
    if (recording == NA_LOGICAL)
        errorcall(call, "invalid value of `recording'");
    args = CDR(args);
    resize = asInteger(CAR(args));
    if (resize == NA_INTEGER)
        errorcall(call, "invalid value of `resize'");
    args = CDR(args);
    xpinch = asReal(CAR(args));
    args = CDR(args);
    ypinch = asReal(CAR(args));
    args = CDR(args);
    sc = CAR(args);
    if (!isString(sc) && !isInteger(sc) && !isLogical(sc) && !isReal(sc))
        errorcall(call, "invalid value of `canvas'");
    canvas = RGBpar(sc, 0);
    args = CDR(args);
    gamma = asReal(CAR(args));
    args = CDR(args);
    xpos = asInteger(CAR(args));
    args = CDR(args);
    ypos = asInteger(CAR(args));

    R_CheckDeviceAvailable();
    BEGIN_SUSPEND_INTERRUPTS
    {
        /* Allocate and initialize the device driver data */
        if (!(dev = (NewDevDesc *)calloc(1, sizeof(NewDevDesc))))
            return 0;
        /* Do this for early redraw attempts */
        dev->displayList = R_NilValue;
        /* Make sure that this is initialised before a GC can occur.
         * This (and displayList) get protected during GC
         */
        dev->savedSnapshot = R_NilValue;
        GAsetunits(xpinch, ypinch);
        if (!GADeviceDriver(dev, display, width, height, ps, (Rboolean)recording, resize, canvas, gamma, xpos, ypos))
        {
            free(dev);
            errorcall(call, "unable to start device devga");
        }
        gsetVar(install(".Device"), mkString(display[0] ? display : "windows"), R_NilValue);
        dd = GEcreateDevDesc(dev);
        addDevice((DevDesc *)dd);
        GEinitDisplayList(dd);
    }
    END_SUSPEND_INTERRUPTS;
    vmaxset(vmax);
    return R_NilValue;
}
