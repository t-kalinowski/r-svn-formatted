/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2000   the R Development Core Team
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

/* this is just a stub to call the X11 version */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Defn.h"
#include "Graphics.h"

SEXP RX11_dataentry(SEXP call, SEXP op, SEXP args, SEXP rho);

SEXP do_dataentry(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    return RX11_dataentry(call, op, args, rho);
}

int stub_X11DeviceDriver(DevDesc *dd, char *display, double width, double height, double pointsize, double gamma,
                         int colormodel, int maxcube)
{
    error("no x11 device in this version of R");
    return 0;
}
