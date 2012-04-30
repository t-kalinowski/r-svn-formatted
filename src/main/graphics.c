/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997--2012  The R Core Team
 *  Copyright (C) 2002--2011  The R Foundation
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
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/


 *  This is an extensive reworking by Paul Murrell of an original
 *  quick hack by Ross Ihaka designed to give a superset of the
 *  functionality in the AT&T Bell Laboratories GRZ library.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Defn.h>
#include <float.h> /* for DBL_EPSILON etc */
#include <Graphics.h>
// --> R_ext/GraphicsEngine.h + Rgraphics.h
#include <GraphicsBase.h> /* setBaseDevice */
#include <Rmath.h>        /* eg. fmax2() */

/*--->> Documentation now in  ../include/Rgraphics.h  "API" ----- */

/*-------------------------------------------------------------------
 *
 *  TRANSFORMATIONS
 *
 *    There are five major regions on a device, for any
 *    particular figure:  the outer margins, which "stick"
 *    to the edges of the device;  the inner region, which
 *    is defined as the total device less the outer margins;
 *    the figure region, which defaults from the current
 *    layout (mfrow, mfcol, layout) unless the user specifies
 *    it directly (fig, fin);  the figure margins, which
 *    "stick" to the edges of the plot region;	and thed
 *    plot region, which is the figure region less the figure
 *    margins by default unless the user specifies it directly
 *    (plt, pin)
 *
 *  COORDINATE SYSTEMS
 *
 *    DEVICE  = devices natural coordinate system
 *		(e.g., pixels, 1/72", ...)
 *    NDC     = normalised device coordinates (0..1 on device)
 *    INCHES  = inches
 *    OMA1..4 = outer margin coordinates
 *    NIC     = normalised inner region coordinates
 *		(0..1 on inner region)
 *    NFC     = normalised figure coordinates
 *		(0..1 on figure region)
 *    MAR1..4 = figure margin coordinates
 *    NPC     = normalised plot coordinates
 *		(0..1 on plot region)
 *    USER    = world or data coordinates
 *
 *
 *  UNITS
 *
 *    All of the above, except OMA1..4 and MAR1..4, plus ...
 *
 *	LINES = line coordinates (lines of margin;  based on mex)
 *	CHARS = char coordinates (lines of text;  based on cex)
 *
 *    The function Convert(value, from, to) is provided
 *    to transform between any pair of coordinate systems
 *    (for transforming locations)
 *
 *    The functions ConvertXUnits(value, from, to) and
 *    ConvertYUnits(value, from, to) are provided to transform
 *    between any pair of units (for transforming dimensions)
 *
 *    IMPORTANT:  if user coordinates are logged, then the
 *    conversion to/from USER units will not work.  in this
 *    case it is necessary to use convert(x1) - convert(x2)
 *    rather than convert(x1 - x2)
 *
 */

double attribute_hidden R_Log10(double x)
{
    return (R_FINITE(x) && x > 0.0) ? log10(x) : NA_REAL;
}

/* In interpreted R, units are as follows:
 *	1 = "user"
 *	2 = "figure"
 *	3 = "inches"
 * the function GMapUnits provides a mapping
 * between interpreted units and internal units.
 */
GUnit GMapUnits(int Runits)
{
    switch (Runits)
    {
    case 1:
        return USER;
    case 2:
        return NFC;
    case 3:
        return INCHES;
    default:
        return 0;
    }
}

/* Conversions Between Units*/

/* Used to be global (non-static) -- but are nowhere declared.
 * The public interface is through G[XY]ConvertUnits() */

static double xNDCtoDevUnits(double x, pGEDevDesc dd)
{
    return x * fabs(gpptr(dd)->ndc2dev.bx);
}

static double yNDCtoDevUnits(double y, pGEDevDesc dd)
{
    return y * fabs(gpptr(dd)->ndc2dev.by);
}

static double xNICtoDevUnits(double x, pGEDevDesc dd)
{
    return x * fabs(gpptr(dd)->inner2dev.bx);
}

static double yNICtoDevUnits(double y, pGEDevDesc dd)
{
    return y * fabs(gpptr(dd)->inner2dev.by);
}

static double xNFCtoDevUnits(double x, pGEDevDesc dd)
{
    return x * fabs(gpptr(dd)->fig2dev.bx);
}

static double yNFCtoDevUnits(double y, pGEDevDesc dd)
{
    return y * fabs(gpptr(dd)->fig2dev.by);
}

static double xNPCtoDevUnits(double x, pGEDevDesc dd)
{
    return xNFCtoDevUnits(x * (gpptr(dd)->plt[1] - gpptr(dd)->plt[0]), dd);
}

static double yNPCtoDevUnits(double y, pGEDevDesc dd)
{
    return yNFCtoDevUnits(y * (gpptr(dd)->plt[3] - gpptr(dd)->plt[2]), dd);
}

static double xUsrtoDevUnits(double x, pGEDevDesc dd)
{
    return xNFCtoDevUnits(x * gpptr(dd)->win2fig.bx, dd);
}

static double yUsrtoDevUnits(double y, pGEDevDesc dd)
{
    return yNFCtoDevUnits(y * gpptr(dd)->win2fig.by, dd);
}

static double xInchtoDevUnits(double x, pGEDevDesc dd)
{
    return xNDCtoDevUnits(x * gpptr(dd)->xNDCPerInch, dd);
}

static double yInchtoDevUnits(double y, pGEDevDesc dd)
{
    return yNDCtoDevUnits(y * gpptr(dd)->yNDCPerInch, dd);
}

static double xLinetoDevUnits(double x, pGEDevDesc dd)
{
    return xNDCtoDevUnits(x * gpptr(dd)->xNDCPerLine, dd);
}

static double yLinetoDevUnits(double y, pGEDevDesc dd)
{
    return yNDCtoDevUnits(y * gpptr(dd)->yNDCPerLine, dd);
}

static double xChartoDevUnits(double x, pGEDevDesc dd)
{
    return xNDCtoDevUnits(x * gpptr(dd)->cex * gpptr(dd)->xNDCPerChar, dd);
}

static double yChartoDevUnits(double y, pGEDevDesc dd)
{
    return yNDCtoDevUnits(y * gpptr(dd)->cex * gpptr(dd)->yNDCPerChar, dd);
}

static double xDevtoNDCUnits(double x, pGEDevDesc dd)
{
    return x / fabs(gpptr(dd)->ndc2dev.bx);
}

static double yDevtoNDCUnits(double y, pGEDevDesc dd)
{
    return y / fabs(gpptr(dd)->ndc2dev.by);
}

static double xDevtoNICUnits(double x, pGEDevDesc dd)
{
    return x / fabs(gpptr(dd)->inner2dev.bx);
}

static double yDevtoNICUnits(double y, pGEDevDesc dd)
{
    return y / fabs(gpptr(dd)->inner2dev.by);
}

static double xDevtoNFCUnits(double x, pGEDevDesc dd)
{
    return x / fabs(gpptr(dd)->fig2dev.bx);
}

static double yDevtoNFCUnits(double y, pGEDevDesc dd)
{
    return y / fabs(gpptr(dd)->fig2dev.by);
}

static double xDevtoNPCUnits(double x, pGEDevDesc dd)
{
    return xDevtoNFCUnits(x, dd) / (gpptr(dd)->plt[1] - gpptr(dd)->plt[0]);
}

static double yDevtoNPCUnits(double y, pGEDevDesc dd)
{
    return yDevtoNFCUnits(y, dd) / (gpptr(dd)->plt[3] - gpptr(dd)->plt[2]);
}

static double xDevtoUsrUnits(double x, pGEDevDesc dd)
{
    return xDevtoNFCUnits(x, dd) / gpptr(dd)->win2fig.bx;
}

static double yDevtoUsrUnits(double y, pGEDevDesc dd)
{
    return yDevtoNFCUnits(y, dd) / gpptr(dd)->win2fig.by;
}

static double xDevtoInchUnits(double x, pGEDevDesc dd)
{
    return xDevtoNDCUnits(x, dd) / gpptr(dd)->xNDCPerInch;
}

static double yDevtoInchUnits(double y, pGEDevDesc dd)
{
    return yDevtoNDCUnits(y, dd) / gpptr(dd)->yNDCPerInch;
}

static double xDevtoLineUnits(double x, pGEDevDesc dd)
{
    return xDevtoNDCUnits(x, dd) / gpptr(dd)->xNDCPerLine;
}

static double yDevtoLineUnits(double y, pGEDevDesc dd)
{
    return yDevtoNDCUnits(y, dd) / gpptr(dd)->yNDCPerLine;
}

/* NOTE that use the _current_ gpptr(dd)->cex here */
/* the conversion for lines doesn't have to worry about */
/* this because gpptr(dd)->mex can only be set once per plot */

static double xDevtoCharUnits(double x, pGEDevDesc dd)
{
    return xDevtoNDCUnits(x, dd) / (gpptr(dd)->cex * gpptr(dd)->xNDCPerChar);
}

static double yDevtoCharUnits(double y, pGEDevDesc dd)
{
    return yDevtoNDCUnits(y, dd) / (gpptr(dd)->cex * gpptr(dd)->yNDCPerChar);
}

static void BadUnitsError(const char *where)
{
    error(_("bad units specified in '%s'"), where);
}

/* GConvertXUnits() and GConvertYUnits() convert
   a single value fromUnits toUnits : */

double GConvertXUnits(double x, GUnit fromUnits, GUnit toUnits, pGEDevDesc dd)
{
    double dev, final;
    switch (fromUnits)
    {
    case DEVICE:
        dev = x;
        break;
    case NDC:
        dev = xNDCtoDevUnits(x, dd);
        break;
    case NIC:
        dev = xNICtoDevUnits(x, dd);
        break;
    case NFC:
        dev = xNFCtoDevUnits(x, dd);
        break;
    case NPC:
        dev = xNPCtoDevUnits(x, dd);
        break;
    case USER:
        dev = xUsrtoDevUnits(x, dd);
        break;
    case INCHES:
        dev = xInchtoDevUnits(x, dd);
        break;
    case LINES:
        dev = xLinetoDevUnits(x, dd);
        break;
    case CHARS:
        dev = xChartoDevUnits(x, dd);
        break;
    default:
        dev = 0;
        BadUnitsError("GConvertXUnits");
    }
    switch (toUnits)
    {
    case DEVICE:
        final = dev;
        break;
    case NDC:
        final = xDevtoNDCUnits(dev, dd);
        break;
    case NIC:
        final = xDevtoNICUnits(dev, dd);
        break;
    case NFC:
        final = xDevtoNFCUnits(dev, dd);
        break;
    case NPC:
        final = xDevtoNPCUnits(dev, dd);
        break;
    case USER:
        final = xDevtoUsrUnits(dev, dd);
        break;
    case INCHES:
        final = xDevtoInchUnits(dev, dd);
        break;
    case LINES:
        final = xDevtoLineUnits(dev, dd);
        break;
    case CHARS:
        final = xDevtoCharUnits(dev, dd);
        break;
    default:
        final = 0;
        BadUnitsError("GConvertXUnits");
    }
    return final;
}

double GConvertYUnits(double y, GUnit fromUnits, GUnit toUnits, pGEDevDesc dd)
{
    double dev, final;
    switch (fromUnits)
    {
    case DEVICE:
        dev = y;
        break;
    case NDC:
        dev = yNDCtoDevUnits(y, dd);
        break;
    case NIC:
        dev = yNICtoDevUnits(y, dd);
        break;
    case NFC:
        dev = yNFCtoDevUnits(y, dd);
        break;
    case NPC:
        dev = yNPCtoDevUnits(y, dd);
        break;
    case USER:
        dev = yUsrtoDevUnits(y, dd);
        break;
    case INCHES:
        dev = yInchtoDevUnits(y, dd);
        break;
    case LINES:
        dev = yLinetoDevUnits(y, dd);
        break;
    case CHARS:
        dev = yChartoDevUnits(y, dd);
        break;
    default:
        dev = 0;
        BadUnitsError("GConvertYUnits");
    }
    switch (toUnits)
    {
    case DEVICE:
        final = dev;
        break;
    case NDC:
        final = yDevtoNDCUnits(dev, dd);
        break;
    case NIC:
        final = yDevtoNICUnits(dev, dd);
        break;
    case NFC:
        final = yDevtoNFCUnits(dev, dd);
        break;
    case NPC:
        final = yDevtoNPCUnits(dev, dd);
        break;
    case USER:
        final = yDevtoUsrUnits(dev, dd);
        break;
    case INCHES:
        final = yDevtoInchUnits(dev, dd);
        break;
    case LINES:
        final = yDevtoLineUnits(dev, dd);
        break;
    case CHARS:
        final = yDevtoCharUnits(dev, dd);
        break;
    default:
        final = 0;
        BadUnitsError("GConvertYUnits");
    }
    return final;
}

/* Functions to convert locations from one coordinate system to another */

/* OTHER coordinate systems to DEVICE */

/* Used to be global (non-static) -- but are nowhere declared.
 * The public interface is  GConvert(), GConvertX(), GConvertY() */
static double xNDCtoDev(double x, pGEDevDesc dd)
{
    return gpptr(dd)->ndc2dev.ax + x * gpptr(dd)->ndc2dev.bx;
}

static double yNDCtoDev(double y, pGEDevDesc dd)
{
    return gpptr(dd)->ndc2dev.ay + y * gpptr(dd)->ndc2dev.by;
}

static double xInchtoDev(double x, pGEDevDesc dd)
{
    return xNDCtoDev(x * gpptr(dd)->xNDCPerInch, dd);
}

static double yInchtoDev(double y, pGEDevDesc dd)
{
    return yNDCtoDev(y * gpptr(dd)->yNDCPerInch, dd);
}

static double xLinetoDev(double x, pGEDevDesc dd)
{
    return xNDCtoDev(x * gpptr(dd)->xNDCPerLine, dd);
}

static double yLinetoDev(double y, pGEDevDesc dd)
{
    return yNDCtoDev(y * gpptr(dd)->yNDCPerLine, dd);
}

static double xNICtoDev(double x, pGEDevDesc dd)
{
    return gpptr(dd)->inner2dev.ax + x * gpptr(dd)->inner2dev.bx;
}

static double yNICtoDev(double y, pGEDevDesc dd)
{
    return gpptr(dd)->inner2dev.ay + y * gpptr(dd)->inner2dev.by;
}
/* NOTE that an x-coordinate in OMA2 or OMA4 converts to a */
/* y-coordinate in Dev and a y-coordinate in OMA2 or OMA4 */
/* converts to an x-coordinate in Dev */

static double xOMA1toDev(double x, pGEDevDesc dd)
{
    return xNICtoDev(x, dd);
}

static double yOMA1toDev(double y, pGEDevDesc dd)
{
    return yLinetoDev((gpptr(dd)->oma[0] - y), dd);
}

static double xOMA2toyDev(double x, pGEDevDesc dd)
{
    return yNICtoDev(x, dd);
}

static double yOMA2toxDev(double y, pGEDevDesc dd)
{
    return xLinetoDev((gpptr(dd)->oma[1] - y), dd);
}

static double xOMA3toDev(double x, pGEDevDesc dd)
{
    return xNICtoDev(x, dd);
}

static double yOMA3toDev(double y, pGEDevDesc dd)
{
    double ndc = 1.0 - yDevtoNDC(yLinetoDev((gpptr(dd)->oma[2] - y), dd), dd);
    return yNDCtoDev(ndc, dd);
}

static double xOMA4toyDev(double x, pGEDevDesc dd)
{
    return yNICtoDev(x, dd);
}

static double yOMA4toxDev(double y, pGEDevDesc dd)
{
    double ndc = 1.0 - xDevtoNDC(xLinetoDev(gpptr(dd)->oma[3] - y, dd), dd);
    return xNDCtoDev(ndc, dd);
}

static double xNFCtoDev(double x, pGEDevDesc dd)
{
    return gpptr(dd)->fig2dev.ax + x * gpptr(dd)->fig2dev.bx;
}

static double yNFCtoDev(double y, pGEDevDesc dd)
{
    return gpptr(dd)->fig2dev.ay + y * gpptr(dd)->fig2dev.by;
}

static double xNPCtoDev(double x, pGEDevDesc dd)
{
    return xNFCtoDev(gpptr(dd)->plt[0] + x * (gpptr(dd)->plt[1] - gpptr(dd)->plt[0]), dd);
}

static double yNPCtoDev(double y, pGEDevDesc dd)
{
    return yNFCtoDev(gpptr(dd)->plt[2] + y * (gpptr(dd)->plt[3] - gpptr(dd)->plt[2]), dd);
}

static double xUsrtoDev(double x, pGEDevDesc dd)
{
    if (gpptr(dd)->xlog)
        x = R_Log10(x);
    return xNFCtoDev(gpptr(dd)->win2fig.ax + x * gpptr(dd)->win2fig.bx, dd);
}

static double yUsrtoDev(double y, pGEDevDesc dd)
{
    if (gpptr(dd)->ylog)
        y = R_Log10(y);
    return yNFCtoDev(gpptr(dd)->win2fig.ay + y * gpptr(dd)->win2fig.by, dd);
}

/* NOTE that an x-coordinate in MAR2 or MAR4 converts to a */
/* y-coordinate in Dev and a y-coordinate in MAR2 or MAR4 */
/* converts to an x-coordinate in Dev */

static double xMAR1toDev(double x, pGEDevDesc dd)
{
    return xUsrtoDev(x, dd);
}

static double yMAR1toDev(double y, pGEDevDesc dd)
{
    double nfc = GConvertYUnits(y, LINES, NFC, dd);
    return yNFCtoDev(gpptr(dd)->plt[2] - nfc, dd);
}

static double xMAR2toyDev(double x, pGEDevDesc dd)
{
    return yUsrtoDev(x, dd);
}

static double yMAR2toxDev(double y, pGEDevDesc dd)
{
    double nfc = GConvertXUnits(y, LINES, NFC, dd);
    return xNFCtoDev(gpptr(dd)->plt[0] - nfc, dd);
}

static double xMAR3toDev(double x, pGEDevDesc dd)
{
    return xUsrtoDev(x, dd);
}

static double yMAR3toDev(double y, pGEDevDesc dd)
{
    double nfc = GConvertYUnits(y, LINES, NFC, dd);
    return yNFCtoDev(gpptr(dd)->plt[3] + nfc, dd);
}

static double xMAR4toyDev(double x, pGEDevDesc dd)
{
    return yUsrtoDev(x, dd);
}

static double yMAR4toxDev(double y, pGEDevDesc dd)
{
    double nfc = GConvertXUnits(y, LINES, NFC, dd);
    return xNFCtoDev(gpptr(dd)->plt[1] + nfc, dd);
}

/* DEVICE coordinates to OTHER */

double xDevtoNDC(double x, pGEDevDesc dd)
{
    return (x - gpptr(dd)->ndc2dev.ax) / gpptr(dd)->ndc2dev.bx;
}

double yDevtoNDC(double y, pGEDevDesc dd)
{
    return (y - gpptr(dd)->ndc2dev.ay) / gpptr(dd)->ndc2dev.by;
}

static double xDevtoInch(double x, pGEDevDesc dd)
{
    return xDevtoNDC(x, dd) / gpptr(dd)->xNDCPerInch;
}

static double yDevtoInch(double y, pGEDevDesc dd)
{
    return yDevtoNDC(y, dd) / gpptr(dd)->yNDCPerInch;
}

static double xDevtoLine(double x, pGEDevDesc dd)
{
    return xDevtoNDC(x, dd) / gpptr(dd)->xNDCPerLine;
}

static double yDevtoLine(double y, pGEDevDesc dd)
{
    return yDevtoNDC(y, dd) / gpptr(dd)->yNDCPerLine;
}

static double xDevtoNIC(double x, pGEDevDesc dd)
{
    return (x - gpptr(dd)->inner2dev.ax) / gpptr(dd)->inner2dev.bx;
}

static double yDevtoNIC(double y, pGEDevDesc dd)
{
    return (y - gpptr(dd)->inner2dev.ay) / gpptr(dd)->inner2dev.by;
}

static double xDevtoOMA1(double x, pGEDevDesc dd)
{
    return xDevtoNIC(x, dd);
}

static double yDevtoOMA1(double y, pGEDevDesc dd)
{
    return gpptr(dd)->oma[0] - yDevtoLine(y, dd);
}

static double xDevtoyOMA2(double x, pGEDevDesc dd)
{
    return gpptr(dd)->oma[1] - xDevtoLine(x, dd);
}

static double yDevtoxOMA2(double y, pGEDevDesc dd)
{
    return yDevtoNIC(y, dd);
}

static double xDevtoOMA3(double x, pGEDevDesc dd)
{
    return xDevtoNIC(x, dd);
}

static double yDevtoOMA3(double y, pGEDevDesc dd)
{
    double line = (1.0 - yDevtoNDC(y, dd)) / gpptr(dd)->yNDCPerLine;
    return gpptr(dd)->oma[2] - line;
}

static double xDevtoyOMA4(double x, pGEDevDesc dd)
{
    double line = (1.0 - xDevtoNDC(x, dd)) / gpptr(dd)->xNDCPerLine;
    return gpptr(dd)->oma[3] - line;
}

static double yDevtoxOMA4(double y, pGEDevDesc dd)
{
    return yDevtoNIC(y, dd);
}

double xDevtoNFC(double x, pGEDevDesc dd)
{
    return (x - gpptr(dd)->fig2dev.ax) / gpptr(dd)->fig2dev.bx;
}

double yDevtoNFC(double y, pGEDevDesc dd)
{
    return (y - gpptr(dd)->fig2dev.ay) / gpptr(dd)->fig2dev.by;
}

double xDevtoNPC(double x, pGEDevDesc dd)
{
    return (xDevtoNFC(x, dd) - gpptr(dd)->plt[0]) / (gpptr(dd)->plt[1] - gpptr(dd)->plt[0]);
}

double yDevtoNPC(double y, pGEDevDesc dd)
{
    return (yDevtoNFC(y, dd) - gpptr(dd)->plt[2]) / (gpptr(dd)->plt[3] - gpptr(dd)->plt[2]);
}

/* a special case (NPC = normalised plot region coordinates) */

double xNPCtoUsr(double x, pGEDevDesc dd)
{
    if (gpptr(dd)->xlog)
        return pow(10., gpptr(dd)->logusr[0] + x * (gpptr(dd)->logusr[1] - gpptr(dd)->logusr[0]));
    else
        return gpptr(dd)->usr[0] + x * (gpptr(dd)->usr[1] - gpptr(dd)->usr[0]);
}

double yNPCtoUsr(double y, pGEDevDesc dd)
{
    if (gpptr(dd)->ylog)
        return pow(10., gpptr(dd)->logusr[2] + y * (gpptr(dd)->logusr[3] - gpptr(dd)->logusr[2]));
    else
        return gpptr(dd)->usr[2] + y * (gpptr(dd)->usr[3] - gpptr(dd)->usr[2]);
}

double xDevtoUsr(double x, pGEDevDesc dd)
{
    double nfc = xDevtoNFC(x, dd);
    if (gpptr(dd)->xlog)
        return pow(10., (nfc - gpptr(dd)->win2fig.ax) / gpptr(dd)->win2fig.bx);
    else
        return (nfc - gpptr(dd)->win2fig.ax) / gpptr(dd)->win2fig.bx;
}

double yDevtoUsr(double y, pGEDevDesc dd)
{
    double nfc = yDevtoNFC(y, dd);
    if (gpptr(dd)->ylog)
        return pow(10., (nfc - gpptr(dd)->win2fig.ay) / gpptr(dd)->win2fig.by);
    else
        return (nfc - gpptr(dd)->win2fig.ay) / gpptr(dd)->win2fig.by;
}

static double xDevtoMAR1(double x, pGEDevDesc dd)
{
    return xDevtoUsr(x, dd);
}

static double yDevtoMAR1(double y, pGEDevDesc dd)
{
    return gpptr(dd)->oma[0] + gpptr(dd)->mar[0] - yDevtoLine(y, dd);
}

static double xDevtoyMAR2(double x, pGEDevDesc dd)
{
    return gpptr(dd)->oma[1] + gpptr(dd)->mar[1] - xDevtoLine(x, dd);
}

static double yDevtoxMAR2(double y, pGEDevDesc dd)
{
    return yDevtoUsr(y, dd);
}

static double xDevtoMAR3(double x, pGEDevDesc dd)
{
    return xDevtoUsr(x, dd);
}

static double yDevtoMAR3(double y, pGEDevDesc dd)
{
    double line = GConvertYUnits(1.0 - yDevtoNFC(y, dd), NFC, LINES, dd);
    return gpptr(dd)->mar[2] - line;
}

static double xDevtoyMAR4(double x, pGEDevDesc dd)
{
    double line = GConvertXUnits(1.0 - xDevtoNFC(x, dd), NFC, LINES, dd);
    return gpptr(dd)->mar[3] - line;
}

static double yDevtoxMAR4(double y, pGEDevDesc dd)
{
    return yDevtoUsr(y, dd);
}

/* the Convert function converts a LOCATION in the FROM coordinate */
/* system to a LOCATION in the TO coordinate system */

void GConvert(double *x, double *y, GUnit from, GUnit to, pGEDevDesc dd)
{
    double devx, devy;

    switch (from)
    {
    case DEVICE:
        devx = *x;
        devy = *y;
        break;
    case NDC:
        devx = xNDCtoDev(*x, dd);
        devy = yNDCtoDev(*y, dd);
        break;
    case INCHES:
        devx = xInchtoDev(*x, dd);
        devy = yInchtoDev(*y, dd);
        break;
    case OMA1:
        devx = xOMA1toDev(*x, dd);
        devy = yOMA1toDev(*y, dd);
        break;
    case OMA2:
        devx = yOMA2toxDev(*y, dd);
        devy = xOMA2toyDev(*x, dd);
        break;
    case OMA3:
        devx = xOMA3toDev(*x, dd);
        devy = yOMA3toDev(*y, dd);
        break;
    case OMA4:
        devx = yOMA4toxDev(*y, dd);
        devy = xOMA4toyDev(*x, dd);
        break;
    case NIC:
        devx = xNICtoDev(*x, dd);
        devy = yNICtoDev(*y, dd);
        break;
    case NFC:
        devx = xNFCtoDev(*x, dd);
        devy = yNFCtoDev(*y, dd);
        break;
    case MAR1:
        devx = xMAR1toDev(*x, dd);
        devy = yMAR1toDev(*y, dd);
        break;
    case MAR2:
        devx = yMAR2toxDev(*y, dd);
        devy = xMAR2toyDev(*x, dd);
        break;
    case MAR3:
        devx = xMAR3toDev(*x, dd);
        devy = yMAR3toDev(*y, dd);
        break;
    case MAR4:
        devx = yMAR4toxDev(*y, dd);
        devy = xMAR4toyDev(*x, dd);
        break;
    case NPC:
        devx = xNPCtoDev(*x, dd);
        devy = yNPCtoDev(*y, dd);
        break;
    case USER:
        devx = xUsrtoDev(*x, dd);
        devy = yUsrtoDev(*y, dd);
        break;
    default:
        devx = 0; /* for -Wall */
        devy = 0;
        BadUnitsError("GConvert");
    }

    switch (to)
    {
    case DEVICE:
        *x = devx;
        *y = devy;
        break;
    case NDC:
        *x = xDevtoNDC(devx, dd);
        *y = yDevtoNDC(devy, dd);
        break;
    case INCHES:
        *x = xDevtoInch(devx, dd);
        *y = yDevtoInch(devy, dd);
        break;
    case LINES:
        *x = xDevtoLine(devx, dd);
        *y = yDevtoLine(devy, dd);
        break;
    case NIC:
        *x = xDevtoNIC(devx, dd);
        *y = yDevtoNIC(devy, dd);
        break;
    case OMA1:
        *x = xDevtoOMA1(devx, dd);
        *y = yDevtoOMA1(devy, dd);
        break;
    case OMA2:
        *x = yDevtoxOMA2(devy, dd);
        *y = xDevtoyOMA2(devx, dd);
        break;
    case OMA3:
        *x = xDevtoOMA3(devx, dd);
        *y = yDevtoOMA3(devy, dd);
        break;
    case OMA4:
        *x = yDevtoxOMA4(devy, dd);
        *y = xDevtoyOMA4(devx, dd);
        break;
    case NFC:
        *x = xDevtoNFC(devx, dd);
        *y = yDevtoNFC(devy, dd);
        break;
    case NPC:
        *x = xDevtoNPC(devx, dd);
        *y = yDevtoNPC(devy, dd);
        break;
    case USER:
        *x = xDevtoUsr(devx, dd);
        *y = yDevtoUsr(devy, dd);
        break;
    case MAR1:
        *x = xDevtoMAR1(devx, dd);
        *y = yDevtoMAR1(devy, dd);
        break;
    case MAR2:
        *x = yDevtoxMAR2(devy, dd);
        *y = xDevtoyMAR2(devx, dd);
        break;
    case MAR3:
        *x = xDevtoMAR3(devx, dd);
        *y = yDevtoMAR3(devy, dd);
        break;
    case MAR4:
        *x = yDevtoxMAR4(devy, dd);
        *y = xDevtoyMAR4(devx, dd);
        break;
    default:
        BadUnitsError("GConvert");
    }
}

double GConvertX(double x, GUnit from, GUnit to, pGEDevDesc dd)
{
    double devx;
    switch (from)
    {
    case DEVICE:
        devx = x;
        break;
    case NDC:
        devx = xNDCtoDev(x, dd);
        break;
    case INCHES:
        devx = xInchtoDev(x, dd);
        break;
    case LINES:
        devx = xLinetoDev(x, dd);
        break;
    case OMA1:
        devx = xOMA1toDev(x, dd);
        break;
    /*case OMA2:	x <--> y */
    case OMA3:
        devx = xOMA3toDev(x, dd);
        break;
    /*case OMA4:	x <--> y */
    case NIC:
        devx = xNICtoDev(x, dd);
        break;
    case NFC:
        devx = xNFCtoDev(x, dd);
        break;
    case MAR1:
        devx = xMAR1toDev(x, dd);
        break;
    /*case MAR2:	x <--> y */
    case MAR3:
        devx = xMAR3toDev(x, dd);
        break;
    /*case MAR4:	x <--> y */
    case NPC:
        devx = xNPCtoDev(x, dd);
        break;
    case USER:
        devx = xUsrtoDev(x, dd);
        break;
    default:
        devx = 0; /* for -Wall */
        BadUnitsError("GConvertX");
    }

    switch (to)
    {
    case DEVICE:
        x = devx;
        break;
    case NDC:
        x = xDevtoNDC(devx, dd);
        break;
    case INCHES:
        x = xDevtoInch(devx, dd);
        break;
    case LINES:
        x = xDevtoLine(devx, dd);
        break;
    case NIC:
        x = xDevtoNIC(devx, dd);
        break;
    case OMA1:
        x = xDevtoOMA1(devx, dd);
        break;
    /*case OMA2:	x <--> y */
    case OMA3:
        x = xDevtoOMA3(devx, dd);
        break;
    /*case OMA4:	x <--> y */
    case NFC:
        x = xDevtoNFC(devx, dd);
        break;
    case USER:
        x = xDevtoUsr(devx, dd);
        break;
    case MAR1:
        x = xDevtoMAR1(devx, dd);
        break;
    /*case MAR2:	x <--> y */
    case MAR3:
        x = xDevtoMAR3(devx, dd);
        break;
    /*case MAR4:	x <--> y */
    case NPC:
        x = xDevtoNPC(devx, dd);
        break;
    default:
        BadUnitsError("GConvertX");
    }
    return x;
}

double GConvertY(double y, GUnit from, GUnit to, pGEDevDesc dd)
{
    double devy;
    switch (from)
    {
    case DEVICE:
        devy = y;
        break;
    case NDC:
        devy = yNDCtoDev(y, dd);
        break;
    case INCHES:
        devy = yInchtoDev(y, dd);
        break;
    case LINES:
        devy = yLinetoDev(y, dd);
        break;
    case OMA1:
        devy = yOMA1toDev(y, dd);
        break;
    /*case OMA2:	x <--> y */
    case OMA3:
        devy = yOMA3toDev(y, dd);
        break;
    /*case OMA4:	x <--> y */
    case NIC:
        devy = yNICtoDev(y, dd);
        break;
    case NFC:
        devy = yNFCtoDev(y, dd);
        break;
    case MAR1:
        devy = yMAR1toDev(y, dd);
        break;
    /*case MAR2:	x <--> y */
    case MAR3:
        devy = yMAR3toDev(y, dd);
        break;
    /*case MAR4:	x <--> y */
    case NPC:
        devy = yNPCtoDev(y, dd);
        break;
    case USER:
        devy = yUsrtoDev(y, dd);
        break;
    default:
        devy = 0; /* for -Wall */
        BadUnitsError("GConvertY");
    }

    switch (to)
    {
    case DEVICE:
        y = devy;
        break;
    case NDC:
        y = yDevtoNDC(devy, dd);
        break;
    case INCHES:
        y = yDevtoInch(devy, dd);
        break;
    case LINES:
        y = yDevtoLine(devy, dd);
        break;
    case NIC:
        y = yDevtoNIC(devy, dd);
        break;
    case OMA1:
        y = yDevtoOMA1(devy, dd);
        break;
    /*case OMA2:	x <--> y */
    case OMA3:
        y = yDevtoOMA3(devy, dd);
        break;
    /*case OMA4:	x <--> y */
    case NFC:
        y = yDevtoNFC(devy, dd);
        break;
    case USER:
        y = yDevtoUsr(devy, dd);
        break;
    case MAR1:
        y = yDevtoMAR1(devy, dd);
        break;
    /*case MAR2:	x <--> y */
    case MAR3:
        y = yDevtoMAR3(devy, dd);
        break;
    /*case MAR4:	x <--> y */
    case NPC:
        y = yDevtoNPC(devy, dd);
        break;
    default:
        BadUnitsError("GConvertY");
    }
    return y;
}

/* Code for layouts */

static double sum(double values[], int n, int cmValues[], int cmSum)
{
    int i;
    double s = 0;
    for (i = 0; i < n; i++)
        if ((cmSum && cmValues[i]) || (!cmSum && !cmValues[i]))
            s = s + values[i];
    return s;
}

static double sumWidths(pGEDevDesc dd)
{
    return sum(gpptr(dd)->widths, gpptr(dd)->numcols, gpptr(dd)->cmWidths, 0);
}

static double sumCmWidths(pGEDevDesc dd)
{
    return sum(gpptr(dd)->widths, gpptr(dd)->numcols, gpptr(dd)->cmWidths, 1);
}

static double sumHeights(pGEDevDesc dd)
{
    return sum(gpptr(dd)->heights, gpptr(dd)->numrows, gpptr(dd)->cmHeights, 0);
}

static double sumCmHeights(pGEDevDesc dd)
{
    return sum(gpptr(dd)->heights, gpptr(dd)->numrows, gpptr(dd)->cmHeights, 1);
}

static int tallLayout(double cmWidth, double cmHeight, pGEDevDesc dd)
{
    return (cmHeight / sumHeights(dd)) > (cmWidth / sumWidths(dd));
}

static void figureExtent(int *minCol, int *maxCol, int *minRow, int *maxRow, int figureNum, pGEDevDesc dd)
{
    int minc = -1;
    int maxc = -1;
    int minr = -1;
    int maxr = -1;
    int i, j;
    int nr = gpptr(dd)->numrows;
    for (i = 0; i < nr; i++)
        for (j = 0; j < gpptr(dd)->numcols; j++)
            if (gpptr(dd)->order[i + j * nr] == figureNum)
            {
                if ((minc == -1) || (j < minc))
                    minc = j;
                if ((maxc == -1) || (j > maxc))
                    maxc = j;
                if ((minr == -1) || (i < minr))
                    minr = i;
                if ((maxr == -1) || (i > maxr))
                    maxr = i;
            }
    *minCol = minc;
    *maxCol = maxc;
    *minRow = minr;
    *maxRow = maxr;
}

static double sumRegions(double regions[], int from, int to)
{
    int i;
    double s = 0;
    for (i = from; i < to + 1; i++)
        s = s + regions[i];
    return s;
}

static void largestRegion(double *width, double *height, double layoutAspectRatio, double innerAspectRatio)
{
    if (layoutAspectRatio < innerAspectRatio)
    {
        *width = 1.0;
        *height = layoutAspectRatio / innerAspectRatio;
    }
    else
    {
        *width = innerAspectRatio / layoutAspectRatio;
        *height = 1.0;
    }
}

static void layoutRegion(double *width, double *height, double widths[], double heights[], double cmWidth,
                         double cmHeight, pGEDevDesc dd)
{
    largestRegion(width, height,
                  sum(heights, gpptr(dd)->numrows, gpptr(dd)->cmHeights, 0) /
                      sum(widths, gpptr(dd)->numcols, gpptr(dd)->cmWidths, 0),
                  cmHeight / cmWidth);
}

/* allocate one dimension (width or height) for either */
/* relative or cm units */

static void allocDimension(double dimensions[], double sumDimensions, int n, int cmDimensions[], int cmDimension)
{
    int i;
    for (i = 0; i < n; i++)
        if ((cmDimension && cmDimensions[i]) || (!cmDimension && !cmDimensions[i]))
            dimensions[i] = dimensions[i] / sumDimensions;
}

static void allCmRegions(double widths[], double heights[], double cmWidth, double cmHeight, pGEDevDesc dd)
{
    allocDimension(widths, cmWidth, gpptr(dd)->numcols, gpptr(dd)->cmWidths, 1);
    allocDimension(heights, cmHeight, gpptr(dd)->numrows, gpptr(dd)->cmHeights, 1);
}

static void modifyDimension(double dimension[], double multiplier, double n, int cmDimensions[])
{
    int i;
    for (i = 0; i < n; i++)
        if (!cmDimensions[i])
            dimension[i] = dimension[i] * multiplier;
}

static void modifyRegions(double widths[], double heights[], double colMultiplier, double rowMultiplier, pGEDevDesc dd)
{
    modifyDimension(widths, colMultiplier, gpptr(dd)->numcols, gpptr(dd)->cmWidths);
    modifyDimension(heights, rowMultiplier, gpptr(dd)->numrows, gpptr(dd)->cmHeights);
}

static void regionsWithoutRespect(double widths[], double heights[], pGEDevDesc dd)
{
    allocDimension(widths, sum(widths, gpptr(dd)->numcols, gpptr(dd)->cmWidths, 0), gpptr(dd)->numcols,
                   gpptr(dd)->cmWidths, 0);
    allocDimension(heights, sum(heights, gpptr(dd)->numrows, gpptr(dd)->cmHeights, 0), gpptr(dd)->numrows,
                   gpptr(dd)->cmHeights, 0);
}

static void regionsWithRespect(double widths[], double heights[], double cmWidth, double cmHeight, pGEDevDesc dd)
{
    double cm, rm;
    layoutRegion(&cm, &rm, widths, heights, cmWidth, cmHeight, dd);
    regionsWithoutRespect(widths, heights, dd);
    modifyRegions(widths, heights, cm, rm, dd);
}

static void widthsRespectingHeights(double widths[], double cmWidth, double cmHeight, pGEDevDesc dd)
{
    int i, j;
    int respectedCols[MAX_LAYOUT_COLS];
    double widthLeft;
    double disrespectedWidth = 0;
    int nr = gpptr(dd)->numrows;
    for (j = 0; j < gpptr(dd)->numcols; j++)
    {
        respectedCols[j] = 0;
        widths[j] = gpptr(dd)->widths[j];
    }
    for (i = 0; i < nr; i++)
        for (j = 0; j < gpptr(dd)->numcols; j++)
            if (gpptr(dd)->respect[i + j * nr] && !gpptr(dd)->cmWidths[j])
                respectedCols[j] = 1;
    for (j = 0; j < gpptr(dd)->numcols; j++)
        if (!respectedCols[j])
            disrespectedWidth += gpptr(dd)->widths[j];
    widthLeft = sumHeights(dd) * cmWidth / cmHeight - sumWidths(dd) + disrespectedWidth;
    for (j = 0; j < gpptr(dd)->numcols; j++)
        if (!respectedCols[j])
            widths[j] = widthLeft * widths[j] / disrespectedWidth;
}

static void regionsRespectingHeight(double widths[], double heights[], double cmWidth, double cmHeight, pGEDevDesc dd)
{
    widthsRespectingHeights(widths, cmWidth, cmHeight, dd);
    regionsWithRespect(widths, heights, cmWidth, cmHeight, dd);
}

static void heightsRespectingWidths(double heights[], double cmWidth, double cmHeight, pGEDevDesc dd)
{
    int i, j;
    int respectedRows[MAX_LAYOUT_ROWS];
    double heightLeft;
    double disrespectedHeight = 0;
    int nr = gpptr(dd)->numrows;
    for (i = 0; i < nr; i++)
    {
        respectedRows[i] = 0;
        heights[i] = gpptr(dd)->heights[i];
    }
    for (i = 0; i < nr; i++)
        for (j = 0; j < gpptr(dd)->numcols; j++)
            if (gpptr(dd)->respect[i + j * nr] && !gpptr(dd)->cmHeights[i])
                respectedRows[i] = 1;
    for (i = 0; i < gpptr(dd)->numrows; i++)
        if (!respectedRows[i])
            disrespectedHeight += gpptr(dd)->heights[i];
    heightLeft = sumWidths(dd) * cmHeight / cmWidth - sumHeights(dd) + disrespectedHeight;
    for (i = 0; i < gpptr(dd)->numrows; i++)
        if (!respectedRows[i])
            heights[i] = heightLeft * heights[i] / disrespectedHeight;
}

static void regionsRespectingWidth(double widths[], double heights[], double cmWidth, double cmHeight, pGEDevDesc dd)
{
    heightsRespectingWidths(heights, cmWidth, cmHeight, dd);
    regionsWithRespect(widths, heights, cmWidth, cmHeight, dd);
}

static void noCmRegions(double widths[], double heights[], double cmWidth, double cmHeight, pGEDevDesc dd)
{
    switch (gpptr(dd)->rspct)
    {
    case 0:
        regionsWithoutRespect(widths, heights, dd);
        break;
    case 1:
        regionsWithRespect(widths, heights, cmWidth, cmHeight, dd);
        break;
    case 2:
        if (tallLayout(cmWidth, cmHeight, dd))
            regionsRespectingWidth(widths, heights, cmWidth, cmHeight, dd);
        else
            regionsRespectingHeight(widths, heights, cmWidth, cmHeight, dd);
    }
}

static void notAllCmRegions(double widths[], double heights[], double cmWidth, double cmHeight, pGEDevDesc dd)
{
    double newCmWidth, newCmHeight;
    newCmWidth = cmWidth - sumCmWidths(dd);
    newCmHeight = cmHeight - sumCmHeights(dd);
    noCmRegions(widths, heights, newCmWidth, newCmHeight, dd);
    allocDimension(widths, cmWidth, gpptr(dd)->numcols, gpptr(dd)->cmWidths, 1);
    allocDimension(heights, cmHeight, gpptr(dd)->numrows, gpptr(dd)->cmHeights, 1);
    modifyDimension(widths, newCmWidth / cmWidth, gpptr(dd)->numcols, gpptr(dd)->cmWidths);
    modifyDimension(heights, newCmHeight / cmHeight, gpptr(dd)->numrows, gpptr(dd)->cmHeights);
}

static void widthCmRegions(double widths[], double heights[], double cmWidth, double cmHeight, pGEDevDesc dd)
{
    allocDimension(widths, cmWidth, gpptr(dd)->numcols, gpptr(dd)->cmWidths, 1);
    allocDimension(heights, sumHeights(dd), gpptr(dd)->numrows, gpptr(dd)->cmHeights, 0);
    modifyDimension(heights, (cmHeight - sumCmHeights(dd)) / cmHeight, gpptr(dd)->numrows, gpptr(dd)->cmHeights);
    allocDimension(heights, cmHeight, gpptr(dd)->numrows, gpptr(dd)->cmHeights, 1);
}

static void heightCmRegions(double widths[], double heights[], double cmWidth, double cmHeight, pGEDevDesc dd)
{
    allocDimension(heights, cmHeight, gpptr(dd)->numrows, gpptr(dd)->cmHeights, 1);
    allocDimension(widths, sumWidths(dd), gpptr(dd)->numcols, gpptr(dd)->cmWidths, 0);
    modifyDimension(widths, (cmWidth - sumCmWidths(dd)) / cmWidth, gpptr(dd)->numcols, gpptr(dd)->cmWidths);
    allocDimension(widths, cmWidth, gpptr(dd)->numcols, gpptr(dd)->cmWidths, 1);
}

static Rboolean allCmWidths(pGEDevDesc dd)
{
    int j;
    for (j = 0; j < gpptr(dd)->numcols; j++)
        if (!gpptr(dd)->cmWidths[j])
            return FALSE;
    return TRUE;
}

static Rboolean allCmHeights(pGEDevDesc dd)
{
    int i;
    for (i = 0; i < gpptr(dd)->numrows; i++)
        if (!gpptr(dd)->cmHeights[i])
            return FALSE;
    return TRUE;
}

static Rboolean noCmWidths(pGEDevDesc dd)
{
    int j;
    for (j = 0; j < gpptr(dd)->numcols; j++)
        if (gpptr(dd)->cmWidths[j])
            return FALSE;
    return TRUE;
}

static Rboolean noCmHeights(pGEDevDesc dd)
{
    int i;
    for (i = 0; i < gpptr(dd)->numrows; i++)
        if (gpptr(dd)->cmHeights[i])
            return FALSE;
    return TRUE;
}

static void someCmRegions(double widths[], double heights[], double cmWidth, double cmHeight, pGEDevDesc dd)
{
    if (allCmWidths(dd))
        widthCmRegions(widths, heights, cmWidth, cmHeight, dd);
    else if (allCmHeights(dd))
        heightCmRegions(widths, heights, cmWidth, cmHeight, dd);
    else
        notAllCmRegions(widths, heights, cmWidth, cmHeight, dd);
}

static Rboolean allCm(pGEDevDesc dd)
{
    return allCmWidths(dd) && allCmHeights(dd);
}

static Rboolean noCm(pGEDevDesc dd)
{
    return noCmWidths(dd) && noCmHeights(dd);
}

static void layoutRegions(double widths[], double heights[], double cmWidth, double cmHeight, pGEDevDesc dd)
{
    int i, j;
    for (j = 0; j < gpptr(dd)->numcols; j++)
        widths[j] = gpptr(dd)->widths[j];
    for (i = 0; i < gpptr(dd)->numrows; i++)
        heights[i] = gpptr(dd)->heights[i];

    if (allCm(dd))
        allCmRegions(widths, heights, cmWidth, cmHeight, dd);
    else if (noCm(dd))
        noCmRegions(widths, heights, cmWidth, cmHeight, dd);
    else
        someCmRegions(widths, heights, cmWidth, cmHeight, dd);
}

static void subRegion(double *left, double *right, double *bottom, double *top, int mincol, int maxcol, int minrow,
                      int maxrow, double widths[], double heights[], pGEDevDesc dd)
{
    double totalWidth = sumRegions(widths, 0, gpptr(dd)->numcols - 1);
    double totalHeight = sumRegions(heights, 0, gpptr(dd)->numrows - 1);
    *left = (0.5 - totalWidth / 2) + sumRegions(widths, 0, mincol - 1);
    *right = (0.5 - totalWidth / 2) + sumRegions(widths, 0, maxcol);
    *bottom = (0.5 - totalHeight / 2) + totalHeight - sumRegions(heights, 0, maxrow);
    *top = (0.5 - totalHeight / 2) + totalHeight - sumRegions(heights, 0, minrow - 1);
}

/* a fudge for backwards compatibility (of sorts) with par(mfg) */
/* return the top-left-most row/col that the current figure */
/* occupies in the current layout */

void currentFigureLocation(int *row, int *col, pGEDevDesc dd)
{
    int maxcol, maxrow;
    if (gpptr(dd)->layout)
        figureExtent(col, &maxcol, row, &maxrow, gpptr(dd)->currentFigure, dd);
    else if (gpptr(dd)->mfind)
    { /* mfcol */
        *row = (gpptr(dd)->currentFigure - 1) % gpptr(dd)->numrows;
        *col = (gpptr(dd)->currentFigure - 1) / gpptr(dd)->numrows;
    }
    else
    { /* mfrow */
        *row = (gpptr(dd)->currentFigure - 1) / gpptr(dd)->numcols;
        *col = (gpptr(dd)->currentFigure - 1) % gpptr(dd)->numcols;
    }
}

/*  mapNDC2Dev -- transformation from NDC to Dev  */
/*  Use this coordinate system for outer margin coordinates  */
/*  This must be called if the device is resized */

static void mapNDC2Dev(pGEDevDesc dd)
{
    /* For new devices, have to check the device's idea of its size
     * in case there has been a resize.
     */
    double asp = dd->dev->ipr[1] / dd->dev->ipr[0];

    gpptr(dd)->ndc2dev.bx = dpptr(dd)->ndc2dev.bx = dd->dev->right - dd->dev->left;
    gpptr(dd)->ndc2dev.ax = dpptr(dd)->ndc2dev.ax = dd->dev->left;
    gpptr(dd)->ndc2dev.by = dpptr(dd)->ndc2dev.by = dd->dev->top - dd->dev->bottom;
    gpptr(dd)->ndc2dev.ay = dpptr(dd)->ndc2dev.ay = dd->dev->bottom;
    /* Units Conversion */

    gpptr(dd)->xNDCPerInch = dpptr(dd)->xNDCPerInch = 1.0 / fabs(gpptr(dd)->ndc2dev.bx * dd->dev->ipr[0]);
    gpptr(dd)->yNDCPerInch = dpptr(dd)->yNDCPerInch = 1.0 / fabs(gpptr(dd)->ndc2dev.by * dd->dev->ipr[1]);
    gpptr(dd)->xNDCPerChar = dpptr(dd)->xNDCPerChar =
        fabs(gpptr(dd)->cexbase * gpptr(dd)->scale * dd->dev->cra[1] * asp / gpptr(dd)->ndc2dev.bx);
    gpptr(dd)->yNDCPerChar = dpptr(dd)->yNDCPerChar =
        fabs(gpptr(dd)->cexbase * gpptr(dd)->scale * dd->dev->cra[1] / gpptr(dd)->ndc2dev.by);
    gpptr(dd)->xNDCPerLine = dpptr(dd)->xNDCPerLine =
        fabs(gpptr(dd)->mex * gpptr(dd)->cexbase * gpptr(dd)->scale * dd->dev->cra[1] * asp / gpptr(dd)->ndc2dev.bx);
    gpptr(dd)->yNDCPerLine = dpptr(dd)->yNDCPerLine =
        fabs(gpptr(dd)->mex * gpptr(dd)->cexbase * gpptr(dd)->scale * dd->dev->cra[1] / gpptr(dd)->ndc2dev.by);
}

static void updateOuterMargins(pGEDevDesc dd)
{
    switch (gpptr(dd)->oUnits)
    {
    case LINES:
        gpptr(dd)->omi[0] = dpptr(dd)->omi[0] = GConvertYUnits(gpptr(dd)->oma[0], LINES, INCHES, dd);
        gpptr(dd)->omi[1] = dpptr(dd)->omi[1] = GConvertXUnits(gpptr(dd)->oma[1], LINES, INCHES, dd);
        gpptr(dd)->omi[2] = dpptr(dd)->omi[2] = GConvertYUnits(gpptr(dd)->oma[2], LINES, INCHES, dd);
        gpptr(dd)->omi[3] = dpptr(dd)->omi[3] = GConvertXUnits(gpptr(dd)->oma[3], LINES, INCHES, dd);
        gpptr(dd)->omd[0] = dpptr(dd)->omd[0] = GConvertXUnits(gpptr(dd)->oma[1], LINES, NDC, dd);
        gpptr(dd)->omd[1] = dpptr(dd)->omd[1] = 1 - GConvertXUnits(gpptr(dd)->oma[3], LINES, NDC, dd);
        gpptr(dd)->omd[2] = dpptr(dd)->omd[2] = GConvertYUnits(gpptr(dd)->oma[0], LINES, NDC, dd);
        gpptr(dd)->omd[3] = dpptr(dd)->omd[3] = 1 - GConvertYUnits(gpptr(dd)->oma[2], LINES, NDC, dd);
        break;
    case INCHES:
        gpptr(dd)->oma[0] = dpptr(dd)->oma[0] = GConvertYUnits(gpptr(dd)->omi[0], INCHES, LINES, dd);
        gpptr(dd)->oma[1] = dpptr(dd)->oma[1] = GConvertXUnits(gpptr(dd)->omi[1], INCHES, LINES, dd);
        gpptr(dd)->oma[2] = dpptr(dd)->oma[2] = GConvertYUnits(gpptr(dd)->omi[2], INCHES, LINES, dd);
        gpptr(dd)->oma[3] = dpptr(dd)->oma[3] = GConvertXUnits(gpptr(dd)->omi[3], INCHES, LINES, dd);
        gpptr(dd)->omd[0] = dpptr(dd)->omd[0] = GConvertXUnits(gpptr(dd)->omi[1], INCHES, NDC, dd);
        gpptr(dd)->omd[1] = dpptr(dd)->omd[1] = 1 - GConvertXUnits(gpptr(dd)->omi[3], INCHES, NDC, dd);
        gpptr(dd)->omd[2] = dpptr(dd)->omd[2] = GConvertYUnits(gpptr(dd)->omi[0], INCHES, NDC, dd);
        gpptr(dd)->omd[3] = dpptr(dd)->omd[3] = 1 - GConvertYUnits(gpptr(dd)->omi[2], INCHES, NDC, dd);
        break;
    case NDC:
        gpptr(dd)->oma[0] = dpptr(dd)->oma[0] = GConvertYUnits(gpptr(dd)->omd[2], NDC, LINES, dd);
        gpptr(dd)->oma[1] = dpptr(dd)->oma[1] = GConvertXUnits(gpptr(dd)->omd[0], NDC, LINES, dd);
        gpptr(dd)->oma[2] = dpptr(dd)->oma[2] = GConvertYUnits(1 - gpptr(dd)->omd[3], NDC, LINES, dd);
        gpptr(dd)->oma[3] = dpptr(dd)->oma[3] = GConvertXUnits(1 - gpptr(dd)->omd[1], NDC, LINES, dd);
        gpptr(dd)->omi[0] = dpptr(dd)->omi[0] = GConvertYUnits(gpptr(dd)->omd[2], NDC, INCHES, dd);
        gpptr(dd)->omi[1] = dpptr(dd)->omi[1] = GConvertXUnits(gpptr(dd)->omd[0], NDC, INCHES, dd);
        gpptr(dd)->omi[2] = dpptr(dd)->omi[2] = GConvertYUnits(1 - gpptr(dd)->omd[3], NDC, INCHES, dd);
        gpptr(dd)->omi[3] = dpptr(dd)->omi[3] = GConvertXUnits(1 - gpptr(dd)->omd[1], NDC, INCHES, dd);
        break;
    default:
        break; /*nothing (-Wall) */
    }
}

/*  mapInner2Dev -- transformation from NIC to Dev  */
/*  Use this coordinate system for setting up multiple figures	*/
/*  This is also used when specifying the figure region directly  */
/*  Note that this is incompatible with S which uses then entire  */
/*  device surface for such a plot  */
/*  This must be called per DevNewPlot, if the NDCtoDev transformation */
/*  changes, and if oma changes */

static void mapInner2Dev(pGEDevDesc dd)
{
    double x0, x1, y0, y1;
    x0 = xLinetoDev(gpptr(dd)->oma[1], dd);
    y0 = yLinetoDev(gpptr(dd)->oma[0], dd);
    x1 = GConvertXUnits(gpptr(dd)->oma[3], LINES, NDC, dd);
    x1 = xNDCtoDev(1.0 - x1, dd);
    y1 = GConvertYUnits(gpptr(dd)->oma[2], LINES, NDC, dd);
    y1 = yNDCtoDev(1.0 - y1, dd);
    gpptr(dd)->inner2dev.bx = dpptr(dd)->inner2dev.bx = x1 - x0;
    gpptr(dd)->inner2dev.ax = dpptr(dd)->inner2dev.ax = x0;
    gpptr(dd)->inner2dev.by = dpptr(dd)->inner2dev.by = y1 - y0;
    gpptr(dd)->inner2dev.ay = dpptr(dd)->inner2dev.ay = y0;
}

/* mapFigureRegion -- calculate figure region in NIC  */

static void mapFigureRegion(pGEDevDesc dd)
{
    int mincol, maxcol, minrow, maxrow;
    double x0, x1, y0, y1;
    double widths[MAX_LAYOUT_COLS], heights[MAX_LAYOUT_ROWS];
    if (gpptr(dd)->layout)
    {
        layoutRegions(widths, heights, GConvertXUnits(1.0, NIC, INCHES, dd) * 2.54,
                      GConvertYUnits(1.0, NIC, INCHES, dd) * 2.54, dd);
        figureExtent(&mincol, &maxcol, &minrow, &maxrow, gpptr(dd)->currentFigure, dd);
        subRegion(&x0, &x1, &y0, &y1, mincol, maxcol, minrow, maxrow, widths, heights, dd);
    }
    else
    {
        int row, col;
        if (gpptr(dd)->mfind)
        {
            col = (gpptr(dd)->currentFigure - 1) / gpptr(dd)->numrows + 1;
            row = gpptr(dd)->currentFigure - (col - 1) * gpptr(dd)->numrows;
        }
        else
        {
            row = (gpptr(dd)->currentFigure - 1) / gpptr(dd)->numcols + 1;
            col = gpptr(dd)->currentFigure - (row - 1) * gpptr(dd)->numcols;
        }
        x0 = (double)(col - 1) / gpptr(dd)->numcols;
        x1 = (double)col / gpptr(dd)->numcols;
        y0 = (double)(gpptr(dd)->numrows - row) / gpptr(dd)->numrows;
        y1 = (double)(gpptr(dd)->numrows - row + 1) / gpptr(dd)->numrows;
    }
    gpptr(dd)->fig[0] = dpptr(dd)->fig[0] = x0;
    gpptr(dd)->fig[1] = dpptr(dd)->fig[1] = x1;
    gpptr(dd)->fig[2] = dpptr(dd)->fig[2] = y0;
    gpptr(dd)->fig[3] = dpptr(dd)->fig[3] = y1;
    gpptr(dd)->fUnits = dpptr(dd)->fUnits = NIC;
}

static void updateFigureRegion(pGEDevDesc dd)
{
    double nicWidth, nicHeight;
    switch (gpptr(dd)->fUnits)
    {
    case NIC:
        gpptr(dd)->fin[0] = dpptr(dd)->fin[0] = GConvertXUnits(gpptr(dd)->fig[1] - gpptr(dd)->fig[0], NIC, INCHES, dd);
        gpptr(dd)->fin[1] = dpptr(dd)->fin[1] = GConvertYUnits(gpptr(dd)->fig[3] - gpptr(dd)->fig[2], NIC, INCHES, dd);
        break;
    case INCHES:
        nicWidth = GConvertXUnits(gpptr(dd)->fin[0], INCHES, NIC, dd);
        nicHeight = GConvertYUnits(gpptr(dd)->fin[1], INCHES, NIC, dd);
        gpptr(dd)->fig[0] = dpptr(dd)->fig[0] = 0.5 - nicWidth / 2;
        gpptr(dd)->fig[1] = dpptr(dd)->fig[1] = gpptr(dd)->fig[0] + nicWidth;
        gpptr(dd)->fig[2] = dpptr(dd)->fig[2] = 0.5 - nicHeight / 2;
        gpptr(dd)->fig[3] = dpptr(dd)->fig[3] = gpptr(dd)->fig[2] + nicHeight;
        break;
    default: /*nothing*/
        break;
    }
}

/*  mapFig2Dev -- Transformation from NFC to Dev  */
/* This must be called per plot.new and if the NICtoDev transformation */
/* changes */

static void mapFig2Dev(pGEDevDesc dd)
{
    double x0, x1, y0, y1;
    y0 = yNICtoDev(gpptr(dd)->fig[2], dd);
    y1 = yNICtoDev(gpptr(dd)->fig[3], dd);
    x0 = xNICtoDev(gpptr(dd)->fig[0], dd);
    x1 = xNICtoDev(gpptr(dd)->fig[1], dd);
    gpptr(dd)->fig2dev.bx = dpptr(dd)->fig2dev.bx = x1 - x0;
    gpptr(dd)->fig2dev.ax = dpptr(dd)->fig2dev.ax = x0;
    gpptr(dd)->fig2dev.by = dpptr(dd)->fig2dev.by = y1 - y0;
    gpptr(dd)->fig2dev.ay = dpptr(dd)->fig2dev.ay = y0;
}

static void updateFigureMargins(pGEDevDesc dd)
{
    switch (gpptr(dd)->mUnits)
    {
    case LINES:
        gpptr(dd)->mai[0] = dpptr(dd)->mai[0] = GConvertYUnits(gpptr(dd)->mar[0], LINES, INCHES, dd);
        gpptr(dd)->mai[1] = dpptr(dd)->mai[1] = GConvertXUnits(gpptr(dd)->mar[1], LINES, INCHES, dd);
        gpptr(dd)->mai[2] = dpptr(dd)->mai[2] = GConvertYUnits(gpptr(dd)->mar[2], LINES, INCHES, dd);
        gpptr(dd)->mai[3] = dpptr(dd)->mai[3] = GConvertXUnits(gpptr(dd)->mar[3], LINES, INCHES, dd);
        break;
    case INCHES:
        gpptr(dd)->mar[0] = dpptr(dd)->mar[0] = GConvertYUnits(gpptr(dd)->mai[0], INCHES, LINES, dd);
        gpptr(dd)->mar[1] = dpptr(dd)->mar[1] = GConvertXUnits(gpptr(dd)->mai[1], INCHES, LINES, dd);
        gpptr(dd)->mar[2] = dpptr(dd)->mar[2] = GConvertYUnits(gpptr(dd)->mai[2], INCHES, LINES, dd);
        gpptr(dd)->mar[3] = dpptr(dd)->mar[3] = GConvertXUnits(gpptr(dd)->mai[3], INCHES, LINES, dd);
        break;
    default: /*nothing*/
        break;
    }
}

/* mapPlotRegion -- plot region in NFC */

static void mapPlotRegion(pGEDevDesc dd)
{
    double x0, x1, y0, y1;
    x0 = GConvertXUnits(gpptr(dd)->mar[1], LINES, NFC, dd);
    y0 = GConvertYUnits(gpptr(dd)->mar[0], LINES, NFC, dd);
    x1 = 1.0 - GConvertXUnits(gpptr(dd)->mar[3], LINES, NFC, dd);
    y1 = 1.0 - GConvertYUnits(gpptr(dd)->mar[2], LINES, NFC, dd);
    if (gpptr(dd)->pty == 's')
    {
        /* maximal plot size in inches */
        double center, width, height;
        double inchWidth = GConvertXUnits(x1 - x0, NFC, INCHES, dd);
        double inchHeight = GConvertYUnits(y1 - y0, NFC, INCHES, dd);
        /* shrink the longer side */
        if (inchWidth > inchHeight)
        {
            width = 0.5 * GConvertXUnits(inchHeight, INCHES, NFC, dd);
            center = 0.5 * (x1 + x0);
            x0 = center - width;
            x1 = center + width;
        }
        else
        {
            height = 0.5 * GConvertYUnits(inchWidth, INCHES, NFC, dd);
            center = 0.5 * (y1 + y0);
            y0 = center - height;
            y1 = center + height;
        }
    }
    gpptr(dd)->plt[0] = dpptr(dd)->plt[0] = x0;
    gpptr(dd)->plt[1] = dpptr(dd)->plt[1] = x1;
    gpptr(dd)->plt[2] = dpptr(dd)->plt[2] = y0;
    gpptr(dd)->plt[3] = dpptr(dd)->plt[3] = y1;
    gpptr(dd)->pUnits = dpptr(dd)->pUnits = NFC;
}

static void updatePlotRegion(pGEDevDesc dd)
{
    double nfcWidth, nfcHeight;
    switch (gpptr(dd)->pUnits)
    {
    case NFC:
        gpptr(dd)->pin[0] = dpptr(dd)->pin[0] = GConvertXUnits(gpptr(dd)->plt[1] - gpptr(dd)->plt[0], NFC, INCHES, dd);
        gpptr(dd)->pin[1] = dpptr(dd)->pin[1] = GConvertYUnits(gpptr(dd)->plt[3] - gpptr(dd)->plt[2], NFC, INCHES, dd);
        break;
    case INCHES:
        nfcWidth = GConvertXUnits(gpptr(dd)->pin[0], INCHES, NFC, dd);
        nfcHeight = GConvertYUnits(gpptr(dd)->pin[1], INCHES, NFC, dd);
        gpptr(dd)->plt[0] = dpptr(dd)->plt[0] = 0.5 - nfcWidth / 2;
        gpptr(dd)->plt[1] = dpptr(dd)->plt[1] = gpptr(dd)->plt[0] + nfcWidth;
        gpptr(dd)->plt[2] = dpptr(dd)->plt[2] = 0.5 - nfcHeight / 2;
        gpptr(dd)->plt[3] = dpptr(dd)->plt[3] = gpptr(dd)->plt[2] + nfcHeight;
        break;
    default: /*nothing*/
        break;
    }
}

/*  GMapWin2Fig -- transformation from Usr to NFC */

void GMapWin2Fig(pGEDevDesc dd)
{
    if (gpptr(dd)->xlog)
    {
        gpptr(dd)->win2fig.bx = dpptr(dd)->win2fig.bx =
            (gpptr(dd)->plt[1] - gpptr(dd)->plt[0]) / (gpptr(dd)->logusr[1] - gpptr(dd)->logusr[0]);
        gpptr(dd)->win2fig.ax = dpptr(dd)->win2fig.ax =
            gpptr(dd)->plt[0] - gpptr(dd)->win2fig.bx * gpptr(dd)->logusr[0];
    }
    else
    {
        gpptr(dd)->win2fig.bx = dpptr(dd)->win2fig.bx =
            (gpptr(dd)->plt[1] - gpptr(dd)->plt[0]) / (gpptr(dd)->usr[1] - gpptr(dd)->usr[0]);
        gpptr(dd)->win2fig.ax = dpptr(dd)->win2fig.ax = gpptr(dd)->plt[0] - gpptr(dd)->win2fig.bx * gpptr(dd)->usr[0];
    }
    if (gpptr(dd)->ylog)
    {
        gpptr(dd)->win2fig.by = dpptr(dd)->win2fig.by =
            (gpptr(dd)->plt[3] - gpptr(dd)->plt[2]) / (gpptr(dd)->logusr[3] - gpptr(dd)->logusr[2]);
        gpptr(dd)->win2fig.ay = dpptr(dd)->win2fig.ay =
            gpptr(dd)->plt[2] - gpptr(dd)->win2fig.by * gpptr(dd)->logusr[2];
    }
    else
    {
        gpptr(dd)->win2fig.by = dpptr(dd)->win2fig.by =
            (gpptr(dd)->plt[3] - gpptr(dd)->plt[2]) / (gpptr(dd)->usr[3] - gpptr(dd)->usr[2]);
        gpptr(dd)->win2fig.ay = dpptr(dd)->win2fig.ay = gpptr(dd)->plt[2] - gpptr(dd)->win2fig.by * gpptr(dd)->usr[2];
    }
}

/*  mapping -- Set up mappings between coordinate systems  */
/*  This is the user's interface to the mapping routines above */

static void mapping(pGEDevDesc dd, int which)
{
    switch (which)
    {
    case 0:
        mapNDC2Dev(dd);
    case 1:
        updateOuterMargins(dd);
        mapInner2Dev(dd);
    case 2:
        if (gpptr(dd)->defaultFigure)
            mapFigureRegion(dd);
        updateFigureRegion(dd);
        mapFig2Dev(dd);
    case 3:
        updateFigureMargins(dd);
        if (gpptr(dd)->defaultPlot)
            mapPlotRegion(dd);
        updatePlotRegion(dd);
    }
}

/*  GReset -- Reset coordinate systems mappings and unit yardsticks */

void GReset(pGEDevDesc dd)
{
    /* Character extents are based on the raster size */
    gpptr(dd)->mkh = gpptr(dd)->scale * dd->dev->cra[0] * dd->dev->ipr[0];

    /* Recompute Mappings */
    mapping(dd, 0);
}

// used in GScale(), but also ../library/grDevices/src/axis_scales.c :
// (usr, log, n_inp) |--> (axp, n_out) :
void GAxisPars(double *min, double *max, int *n, Rboolean log, int axis)
{
#define EPS_FAC_2 100
    Rboolean swap = *min > *max;
    double t_, min_o, max_o;

    if (swap)
    { /* Feature: in R, something like  xlim = c(100,0)  just works */
        t_ = *min;
        *min = *max;
        *max = t_;
    }
    /* save only for the extreme case (EPS_FAC_2): */
    min_o = *min;
    max_o = *max;

    if (log)
    {
        /* Avoid infinities */
        if (*max > 308)
            *max = 308;
        if (*min < -307)
            *min = -307;
        *min = pow(10., *min);
        *max = pow(10., *max);
        GLPretty(min, max, n);
    }
    else
        GEPretty(min, max, n);

    double tmp2 = EPS_FAC_2 * DBL_EPSILON; /* << prevent overflow in product below */
    if (fabs(*max - *min) < (t_ = fmax2(fabs(*max), fabs(*min))) * tmp2)
    {
        /* Treat this case somewhat similar to the (min ~= max) case above */
        /* Too much accuracy here just shows machine differences */
        warning(_("relative range of values =%4.0f * EPS, is small (axis %d)")
                /*"to compute accurately"*/,
                fabs(*max - *min) / (t_ * DBL_EPSILON), axis);

        /* No pretty()ing anymore */
        *min = min_o;
        *max = max_o;
        double eps = .005 * fabs(*max - *min); /* .005: not to go to DBL_MIN/MAX */
        *min += eps;
        *max -= eps;
        if (log)
        {
            *min = pow(10., *min);
            *max = pow(10., *max);
        }
        *n = 1;
    }
    if (swap)
    {
        t_ = *min;
        *min = *max;
        *max = t_;
    }
}

#define LPR_SMALL 2
#define LPR_MEDIUM 3

void GLPretty(double *ul, double *uh, int *n)
{
    /* Generate pretty tick values --	LOGARITHMIC scale
     * __ ul < uh __
     * This only does a very simple setup.
     * The real work happens when the axis is drawn. */
    int p1, p2;
    double dl = *ul, dh = *uh;
    p1 = (int)ceil(log10(dl));
    p2 = (int)floor(log10(dh));
    if (p2 <= p1 && dh / dl > 10.0)
    {
        p1 = (int)ceil(log10(dl) - 0.5);
        p2 = (int)floor(log10(dh) + 0.5);
    }

    if (p2 <= p1)
    { /* floor(log10(uh)) <= ceil(log10(ul))
       * <==>	 log10(uh) - log10(ul) < 2
       * <==>		uh / ul	       < 100 */
        /* Very small range : Use tickmarks from a LINEAR scale
         *		      Splus uses n = 9 here, but that is dumb */
        GPretty(ul, uh, n);
        *n = -*n;
    }
    else
    { /* extra tickmarks --> CreateAtVector() in ./plot.c */
        /* round to nice "1e<N>" */
        *ul = pow(10., (double)p1);
        *uh = pow(10., (double)p2);
        if (p2 - p1 <= LPR_SMALL)
            *n = 3; /* Small range :	Use 1,2,5,10 times 10^k tickmarks */
        else if (p2 - p1 <= LPR_MEDIUM)
            *n = 2; /* Medium range :	Use 1,5 times 10^k tickmarks */
        else
            *n = 1; /* Large range :	Use 10^k tickmarks
                     *			But decimate, when there are too many*/
    }
}

void GPretty(double *lo, double *up, int *ndiv)
{
    GEPretty(lo, up, ndiv);
}

/*-------------------------------------------------------------------
 *
 *  GPAR FUNCTIONS
 *
 */

/* Set default graphics parameter values in a GPar.
 * This initialises the plot state, plus the graphical
 * parameters that are not the responsibility of the device initialisation.

 * Called from baseCallback.
 */

void attribute_hidden GInit(GPar *dp)
{
    dp->state = 0;
    dp->valid = FALSE;

    dp->ann = TRUE;
    dp->err = 0;
    dp->bty = 'o';

    dp->mkh = .001; /* dummy value > 0  --- set in GReset : unused in R */
    dp->cex = 1.0;
    dp->lheight = 1.0;
    dp->cexbase = 1.0;
    dp->cexmain = 1.2;
    dp->cexlab = 1.0;
    dp->cexsub = 1.0;
    dp->cexaxis = 1.0;

    dp->col = R_RGB(0, 0, 0);
    dp->colmain = R_RGB(0, 0, 0);
    dp->collab = R_RGB(0, 0, 0);
    dp->colsub = R_RGB(0, 0, 0);
    dp->colaxis = R_RGB(0, 0, 0);
    dp->gamma = 1;

    dp->scale = 1.0;
    strcpy(dp->family, "");
    dp->font = 1;
    dp->fontmain = 2;
    dp->fontlab = 1;
    dp->fontsub = 1;
    dp->fontaxis = 1;

    dp->pch = 1;
    dp->lty = LTY_SOLID;
    dp->lend = GE_ROUND_CAP;
    dp->ljoin = GE_ROUND_JOIN;
    dp->lmitre = 10.0;
    dp->smo = 1;

    /* String Adjustment and rotation */
    dp->adj = 0.5;
    dp->crt = 0.0;
    dp->srt = 0.0;

    /* Positioning of margin text */
    dp->mgp[0] = 3;
    dp->mgp[1] = 1;
    dp->mgp[2] = 0;

    /* Axis annotation parameters */
    dp->lab[0] = 5;
    dp->lab[1] = 5;
    dp->lab[2] = 7;
    dp->las = 0;
    dp->tck = NA_REAL;
    dp->tcl = -0.5;
    dp->xaxp[0] = 0.0;
    dp->xaxp[1] = 1.0;
    dp->xaxp[2] = 5.0;
    dp->xaxs = 'r';
    dp->xaxt = 's';
    dp->xlog = FALSE;
    dp->xpd = 0;
    dp->oldxpd = -99;
    dp->yaxp[0] = 0.0;
    dp->yaxp[1] = 1.0;
    dp->yaxp[2] = 5.0;
    dp->yaxs = 'r';
    dp->yaxt = 's';
    dp->ylog = FALSE;

    /* Outer Margins */
    dp->mex = 1.0;
    dp->oma[0] = 0.0;
    dp->oma[1] = 0.0;
    dp->oma[2] = 0.0;
    dp->oma[3] = 0.0;
    dp->oUnits = LINES;
    dp->fig[0] = 0.0;
    dp->fig[1] = 1.0;
    dp->fig[2] = 0.0;
    dp->fig[3] = 1.0;
    dp->fUnits = NIC;
    dp->defaultFigure = TRUE; /* the figure region is calculated from */
    /* the layout by default */
    dp->pUnits = NFC;
    dp->defaultPlot = TRUE; /* the plot region is calculated as */
    /* figure-margin by default */

    /* Inner Margins */
    dp->mar[0] = 5.1;
    dp->mar[1] = 4.1;
    dp->mar[2] = 4.1;
    dp->mar[3] = 2.1;
    dp->mUnits = LINES;

    /* Multi-figure parameters */
    dp->layout = FALSE;
    dp->mfind = 0;

    dp->numrows = 1;
    dp->numcols = 1;
    dp->currentFigure = 1;
    dp->lastFigure = 1;
    dp->heights[0] = 1;
    dp->widths[0] = 1;
    dp->cmHeights[0] = 0;
    dp->cmWidths[0] = 0;
    dp->order[0] = 1;
    dp->rspct = 0;
    dp->respect[0] = 0;

    /* Misc plotting parameters */
    dp->new = FALSE;
    dp->devmode = -99;
    dp->pty = 'm';
    dp->lwd = 1;

    /* Data window */
    dp->usr[0] = 0.0;
    dp->usr[1] = 1.0;
    dp->usr[2] = 0.0;
    dp->usr[3] = 1.0;
}

/* Copy a GPar structure from source to dest. */
void copyGPar(GPar *source, GPar *dest)
{
    memcpy(dest, source, sizeof(GPar));
}
