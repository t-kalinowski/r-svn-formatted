/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997 Robert Gentleman, Ross Ihaka and the R Core Team
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
 *
 *
 *  This is an extensive reworking by Paul Murrell of an original
 *  quick hack by Ross Ihaka designed to give a superset of the
 *  functionality in the AT&T Bell Laboratories GRZ library.
 *
 */

#include "Defn.h"
#include "Graphics.h"
#include "Errormsg.h"
#include "Arith.h"
#include "Platform.h"
#include "Mathlib.h" /* for floor(), fmax2(),.. in GPretty(.) */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/********************************************************/
/*							*/
/* 	TRANSFORMATIONS are concerned with converting 	*/
/*	locations between coordinate systems and 	*/
/*	dimensions between different units		*/
/*							*/
/********************************************************/

int GMapUnits(int);
/* Convert an R unit (e.g., "user") into */
/* an internal unit (e.g., USER) 	 */

void GConvert(double *, double *, int, int, DevDesc *);
/* Convert a location from one coordinate */
/* system to another */

double GConvertXUnits(double, int, int, DevDesc *);
double GConvertYUnits(double, int, int, DevDesc *);
/* Convert an x/y-dimension from one set of */
/* units to another */

void GReset(DevDesc *);
/* Set up the different regions on a device */
/* (i.e., inner region, figure region,      */
/* plot region) and transformations for     */
/* associated coordinate systems	    */
/* (called whenever anything that affects   */
/* the coordinate transformations changes)  */

void GMapWin2Fig(DevDesc *);
/* Set up the user coordinate transformations */

void GNewPlot(DevDesc *, int);
/* Set up the device for a new plot by */
/* Resetting graphics parameters and   */
/* Resetting the regions and coordinate */
/* Systems */

void GScale(double, double, int, DevDesc *);
/* Set up the user coordinates based on */
/* the axis limits */

void GSetupAxis(int, DevDesc *);
/* Set up the axis limits based on the */
/* user coordinates */

void currentFigureLocation(int *, int *, DevDesc *);
/* Return row and column of current figure */
/* in the layout matrix */

/********************************************************/
/*							*/
/*	GPAR FUNCTIONS are concerned with operations 	*/
/*	on the entire set of graphics parameters for 	*/
/* 	a device (e.g., initialisation, saving, and	*/
/* 	restoring)					*/
/*							*/
/********************************************************/

void GInit(GPar *);
/* Default the settings for general     */
/* graphical parameters (i.e., defaults */
/* that do not depend on the device     */
/* type */

void GRestore(DevDesc *);
/* Reset the current graphical parameters */
/* from the default graphical parameters */

void GSavePars(DevDesc *);
/* Make a temporary copy of the current */
/* parameters */

void GRestorePars(DevDesc *);
/* Restore the temporary copy saved by */
/* GSavePars */

/********************************************************/
/*							*/
/* 	DEVICE STATE FUNCTIONS are concerned with 	*/
/* 	getting and setting the current state of 	*/
/* 	the device;  is it ready to be drawn into ?	*/
/* 	(i.e., has plot.new been called yet ?)		*/
/*							*/
/********************************************************/

void GCheckState(DevDesc *);
/* has plot.new been called yet ? */

void GSetState(int, DevDesc *);
/* set to 1 when plot.new succeeds */
/* set to 0 when don't want drawing to */
/* go ahead */

/********************************************************/
/*							*/
/* 	GRAPHICAL PRIMITIVES are the generic front-end	*/
/* 	for the functions that every device driver must */
/* 	provide						*/
/* 	NOTE that locations supplied to these functions */
/*	may be in any of the valid coordinate systems	*/
/* 	(each function takes a "coords" parameter to 	*/
/*	indicate the coordinate system);  the device-	*/
/* 	specific version of the function is responsible */
/* 	for calling GConvert to get the location into 	*/
/* 	device coordinates				*/
/*							*/
/********************************************************/

void GCircle(double, double, int, double, int, int, DevDesc *);
/* draw a circle, centred on (x,y) with */
/* radius r (in inches) */

void GClip(DevDesc *);
/* set clipping region (based on current */
/* setting of dd->gp.xpd) */
/* only clip if new clipping region is */
/* different from current clipping region */

void GForceClip(DevDesc *);
/* always clips */

void GLine(double, double, double, double, int, DevDesc *);
/* draw a line from x1,y1 to x2,y2 */

int GLocator(double *, double *, int, DevDesc *);
/* return the location of the next mouse click */

void GMetricInfo(int, double *, double *, double *, int, DevDesc *);
/* return the height, depth, and width of */
/* the specified character in the specified */
/* units */

void GMode(DevDesc *, int);
/* set device "mode" (drawing or not drawing) */
/* here for windows and mac drivers */

void GPolygon(int, double *, double *, int, int, int, DevDesc *);
/* draw a polygon using the specified lists */
/* of x and y values */

void GPolyline(int, double *, double *, int, DevDesc *);
/* draw series of straight lines using the */
/* specified lists of x and y values */

void GRect(double, double, double, double, int, int, int, DevDesc *);
/* draw a rectangle given two opposite corners */

double GStrHeight(char *, int, DevDesc *);
/* return the height of the specified string */
/* in the specified units */

double GStrWidth(char *, int, DevDesc *);
/* return the width of the specified string */
/* in the specified units */

void GText(double, double, int, char *, double, double, double, DevDesc *);
/* draw the specified text at location (x,y) */
/* with the specified rotation and justification */

/********************************************************/
/*							*/
/* 	GRAPHICAL UTILITIES are functions that produce	*/
/* 	graphical output using the graphical primitives */
/* 	(i.e., they are generic - NOT device-specific)	*/
/*							*/
/********************************************************/

void GArrow(double, double, double, double, int, double, double, int, DevDesc *);
/* draw a line from (x1,y1) to (x2,y2)  */
/* with an arrow head at either or both */
/* ends */

void GBox(int, DevDesc *);
/* draw a box around specified region */
/* 1=plot region, 2=figure region,    */
/* 3=inner region, 4=device */

void GLPretty(double *, double *, int *);
/* return a "nice" min, max and number */
/* of intervals for a given range */
/* on a log scale */

void GMtext(char *, int, double, int, double, int, DevDesc *);
/* draw text in margins */

void GPretty(double *, double *, int *);
/* return a "nice" min, max and number */
/* of intervals for a given range on a linear scale */

void GSymbol(double, double, int, int, DevDesc *);
/* draw one of the predefined symbols */
/* (circle, square, diamond, ...) */

/********************************************************/
/*							*/
/* 	COLOUR CODE is concerned with the internals of 	*/
/* 	R colour representation 			*/
/*							*/
/********************************************************/

unsigned RGBpar(SEXP, int, DevDesc *);
/* convert an R colour specification (which */
/* might be a number or a string) */
/* into an internal colour specification */

/********************************************************/
/*							*/
/* 	LINE TEXTURE CODE is concerned with the		*/
/*	internals of R line texture representation	*/
/*							*/
/********************************************************/

unsigned int LTYpar(SEXP, int);

/********************************************************/
/*							*/
/* 	DEVICE FUNCTIONS are concerned with the		*/
/*	creation and destruction of devices		*/
/*							*/
/********************************************************/

DevDesc *CurrentDevice();
/* return a pointer to the current device */

DevDesc *GetDevice(int);
/* return a pointer to a device which is */
/* identified by number */

void InitGraphics(void);
/* initialise internal device structures */

void KillDevice(DevDesc *);
/* kill device which is identified by number */

void KillAllDevices();
/* kill all active devices (used at shutdown) */

int NoDevices();
/* is the null device the current device ? */

int NumDevices();
/* how many devices exist ? (>= 1) */

int deviceNumber(DevDesc *);
/* get the index of the specified device */

int StartDevice(SEXP, SEXP, int, SEXP, int);
/* create a new device */

/********************************************************/
/*							*/
/* 	DEVICE UTILITIES are concerned with providing  	*/
/* 	information for R interpreted functions		*/
/*							*/
/********************************************************/

int curDevice();
/* return the number of the current device */

int nextDevice(int);
/* return the number of the next device */

int prevDevice(int);
/* return the number of the previous device */

int selectDevice(int);
/* make the specified device (specified by number) */
/* the current device */

void killDevice(int);
/* kill device which is identified by number */

/***************************************************************
 *
 * TRANSFORMATIONS
 *
 *    There are five major regions on a device, for any
 *    particular figure:  the outer margins, which "stick"
 *    to the edges of the device;  the inner region, which
 *    is defined as the total device less the outer margins;
 *    the figure region, which defaults from the current
 *    layout (mfrow, mfcol, layout) unless the user specifies
 *    it directly (fig, fin);  the figure margins, which
 *    "stick" to the edges of the plot region;  and thed
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
 *      LINES = line coordinates (lines of margin;  based on mex)
 *      CHARS = char coordinates (lines of text;  based on cex)
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
 **************************************************************/

#define DEG2RAD 0.01745329251994329576

double Log10(double x)
{
    return (FINITE(x) && x > 0.0) ? log10(x) : NA_REAL;
}

int GMapUnits(int Runits)
{
    /* within interpreted R, units are as follows:	*/
    /* 	1 = "user"				*/
    /*	2 = "figure"				*/
    /*	3 = "inches"				*/
    /* the function mapUnits provides a mapping	*/
    /* between interpreted units and internal	*/
    /* units.					*/

    switch (Runits)
    {
    case 1:
        return USER;
    case 2:
        return NFC;
    case 3:
        return INCHES;
    }
    return 0; /* added to keep -Wall happy */
}

/* functions to convert from one set of units to another */

/* forward declarations */

double xDevtoNFC(double x, DevDesc *dd);
double yDevtoNFC(double y, DevDesc *dd);

double xNDCtoDevUnits(double x, DevDesc *dd)
{
    return x * fabs(dd->gp.ndc2dev.bx);
}
double yNDCtoDevUnits(double y, DevDesc *dd)
{
    return y * fabs(dd->gp.ndc2dev.by);
}
double xNICtoDevUnits(double x, DevDesc *dd)
{
    return x * fabs(dd->gp.inner2dev.bx);
}
double yNICtoDevUnits(double y, DevDesc *dd)
{
    return y * fabs(dd->gp.inner2dev.by);
}
double xNFCtoDevUnits(double x, DevDesc *dd)
{
    return x * fabs(dd->gp.fig2dev.bx);
}
double yNFCtoDevUnits(double y, DevDesc *dd)
{
    return y * fabs(dd->gp.fig2dev.by);
}
double xNPCtoDevUnits(double x, DevDesc *dd)
{
    return xNFCtoDevUnits(x * (dd->gp.plt[1] - dd->gp.plt[0]), dd);
}
double yNPCtoDevUnits(double y, DevDesc *dd)
{
    return yNFCtoDevUnits(y * (dd->gp.plt[3] - dd->gp.plt[2]), dd);
}
double xUsrtoDevUnits(double x, DevDesc *dd)
{
    return xNFCtoDevUnits(x * dd->gp.win2fig.bx, dd);
}
double yUsrtoDevUnits(double y, DevDesc *dd)
{
    return yNFCtoDevUnits(y * dd->gp.win2fig.by, dd);
}
double xInchtoDevUnits(double x, DevDesc *dd)
{
    return xNDCtoDevUnits(x * dd->gp.xNDCPerInch, dd);
}
double yInchtoDevUnits(double y, DevDesc *dd)
{
    return yNDCtoDevUnits(y * dd->gp.yNDCPerInch, dd);
}
double xLinetoDevUnits(double x, DevDesc *dd)
{
    return xNDCtoDevUnits(x * dd->gp.xNDCPerLine, dd);
}
double yLinetoDevUnits(double y, DevDesc *dd)
{
    return yNDCtoDevUnits(y * dd->gp.yNDCPerLine, dd);
}
double xChartoDevUnits(double x, DevDesc *dd)
{
    return xNDCtoDevUnits(x * dd->gp.cex * dd->gp.xNDCPerChar, dd);
}
double yChartoDevUnits(double y, DevDesc *dd)
{
    return yNDCtoDevUnits(y * dd->gp.cex * dd->gp.yNDCPerChar, dd);
}

double xDevtoNDCUnits(double x, DevDesc *dd)
{
    return x / fabs(dd->gp.ndc2dev.bx);
}
double yDevtoNDCUnits(double y, DevDesc *dd)
{
    return y / fabs(dd->gp.ndc2dev.by);
}
double xDevtoNICUnits(double x, DevDesc *dd)
{
    return x / fabs(dd->gp.inner2dev.bx);
}
double yDevtoNICUnits(double y, DevDesc *dd)
{
    return y / fabs(dd->gp.inner2dev.by);
}
double xDevtoNFCUnits(double x, DevDesc *dd)
{
    return x / fabs(dd->gp.fig2dev.bx);
}
double yDevtoNFCUnits(double y, DevDesc *dd)
{
    return y / fabs(dd->gp.fig2dev.by);
}
double xDevtoNPCUnits(double x, DevDesc *dd)
{
    return xDevtoNFCUnits(x, dd) / (dd->gp.plt[1] - dd->gp.plt[0]);
}
double yDevtoNPCUnits(double y, DevDesc *dd)
{
    return yDevtoNFCUnits(y, dd) / (dd->gp.plt[3] - dd->gp.plt[2]);
}
double xDevtoUsrUnits(double x, DevDesc *dd)
{
    return xDevtoNFCUnits(x, dd) / dd->gp.win2fig.bx;
}
double yDevtoUsrUnits(double y, DevDesc *dd)
{
    return yDevtoNFCUnits(y, dd) / dd->gp.win2fig.by;
}
double xDevtoInchUnits(double x, DevDesc *dd)
{
    return xDevtoNDCUnits(x, dd) / dd->gp.xNDCPerInch;
}
double yDevtoInchUnits(double y, DevDesc *dd)
{
    return yDevtoNDCUnits(y, dd) / dd->gp.yNDCPerInch;
}
double xDevtoLineUnits(double x, DevDesc *dd)
{
    return xDevtoNDCUnits(x, dd) / dd->gp.xNDCPerLine;
}
double yDevtoLineUnits(double y, DevDesc *dd)
{
    return yDevtoNDCUnits(y, dd) / dd->gp.yNDCPerLine;
}

/* NOTE that use the _current_ dd->gp.cex here */
/* the conversion for lines doesn't have to worry about */
/* this because dd->gp.mex can only be set once per plot */

double xDevtoCharUnits(double x, DevDesc *dd)
{
    return xDevtoNDCUnits(x, dd) / (dd->gp.cex * dd->gp.xNDCPerChar);
}
double yDevtoCharUnits(double y, DevDesc *dd)
{
    return yDevtoNDCUnits(y, dd) / (dd->gp.cex * dd->gp.yNDCPerChar);
}

static void BadUnitsError(char *where)
{
    error("Bad units specified in %s, please report!\n", where);
}

/* the functions GConvertXUnits and ConvertYUnits convert a single */
/* value fromUnits toUnits */

double GConvertXUnits(double x, int fromUnits, int toUnits, DevDesc *dd)
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

double GConvertYUnits(double y, int fromUnits, int toUnits, DevDesc *dd)
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

/* functions to convert locations from one coordinate system */
/* to another                                                */

/* forward declarations */

double xDevtoNDC(double x, DevDesc *dd);
double yDevtoNDC(double y, DevDesc *dd);
double xDevtoNFC(double x, DevDesc *dd);
double yDevtoNFC(double y, DevDesc *dd);

/* OTHER coordinate systems to DEVICE */

double xNDCtoDev(double x, DevDesc *dd)
{
    return dd->gp.ndc2dev.ax + x * dd->gp.ndc2dev.bx;
}
double yNDCtoDev(double y, DevDesc *dd)
{
    return dd->gp.ndc2dev.ay + y * dd->gp.ndc2dev.by;
}
double xInchtoDev(double x, DevDesc *dd)
{
    return xNDCtoDev(x * dd->gp.xNDCPerInch, dd);
}
double yInchtoDev(double y, DevDesc *dd)
{
    return yNDCtoDev(y * dd->gp.yNDCPerInch, dd);
}
double xLinetoDev(double x, DevDesc *dd)
{
    return xNDCtoDev(x * dd->gp.xNDCPerLine, dd);
}
double yLinetoDev(double y, DevDesc *dd)
{
    return yNDCtoDev(y * dd->gp.yNDCPerLine, dd);
}
double xNICtoDev(double x, DevDesc *dd)
{
    return dd->gp.inner2dev.ax + x * dd->gp.inner2dev.bx;
}
double yNICtoDev(double y, DevDesc *dd)
{
    return dd->gp.inner2dev.ay + y * dd->gp.inner2dev.by;
}
/* NOTE that an x-coordinate in OMA2 or OMA4 converts to a */
/* y-coordinate in Dev and a y-coordinate in OMA2 or OMA4 */
/* converts to an x-coordinate in Dev */
double xOMA1toDev(double x, DevDesc *dd)
{
    return xNICtoDev(x, dd);
}
double yOMA1toDev(double y, DevDesc *dd)
{
    return yLinetoDev((dd->gp.oma[0] - y - 1), dd);
}
double xOMA2toyDev(double x, DevDesc *dd)
{
    return yNICtoDev(x, dd);
}
double yOMA2toxDev(double y, DevDesc *dd)
{
    return xLinetoDev((dd->gp.oma[1] - y), dd);
}
double xOMA3toDev(double x, DevDesc *dd)
{
    return xNICtoDev(x, dd);
}
double yOMA3toDev(double y, DevDesc *dd)
{
    double ndc = 1.0 - yDevtoNDC(yLinetoDev((dd->gp.oma[2] - y), dd), dd);
    return yNDCtoDev(ndc, dd);
}
double xOMA4toyDev(double x, DevDesc *dd)
{
    return yNICtoDev(x, dd);
}
double yOMA4toxDev(double y, DevDesc *dd)
{
    double ndc = 1.0 - xDevtoNDC(xLinetoDev(dd->gp.oma[3] - y - 1, dd), dd);
    return xNDCtoDev(ndc, dd);
}
double xNFCtoDev(double x, DevDesc *dd)
{
    return dd->gp.fig2dev.ax + x * dd->gp.fig2dev.bx;
}
double yNFCtoDev(double y, DevDesc *dd)
{
    return dd->gp.fig2dev.ay + y * dd->gp.fig2dev.by;
}
double xNPCtoDev(double x, DevDesc *dd)
{
    return xNFCtoDev(dd->gp.plt[0] + x * (dd->gp.plt[1] - dd->gp.plt[0]), dd);
}
double yNPCtoDev(double y, DevDesc *dd)
{
    return yNFCtoDev(dd->gp.plt[2] + y * (dd->gp.plt[3] - dd->gp.plt[2]), dd);
}
double xUsrtoDev(double x, DevDesc *dd)
{
    if (dd->gp.xlog)
        x = Log10(x);
    return xNFCtoDev(dd->gp.win2fig.ax + x * dd->gp.win2fig.bx, dd);
}
double yUsrtoDev(double y, DevDesc *dd)
{
    if (dd->gp.ylog)
        y = Log10(y);
    return yNFCtoDev(dd->gp.win2fig.ay + y * dd->gp.win2fig.by, dd);
}
/* NOTE that an x-coordinate in MAR2 or MAR4 converts to a */
/* y-coordinate in Dev and a y-coordinate in MAR2 or MAR4 */
/* converts to an x-coordinate in Dev */
double xMAR1toDev(double x, DevDesc *dd)
{
    return xUsrtoDev(x, dd);
}
double yMAR1toDev(double y, DevDesc *dd)
{
    double nfc = GConvertYUnits(y, LINES, NFC, dd);
    return yNFCtoDev(dd->gp.plt[2] - nfc, dd);
}
double xMAR2toyDev(double x, DevDesc *dd)
{
    return yUsrtoDev(x, dd);
}
double yMAR2toxDev(double y, DevDesc *dd)
{
    double nfc = GConvertXUnits(y, LINES, NFC, dd);
    return xNFCtoDev(dd->gp.plt[0] - nfc, dd);
}
double xMAR3toDev(double x, DevDesc *dd)
{
    return xUsrtoDev(x, dd);
}
double yMAR3toDev(double y, DevDesc *dd)
{
    double nfc = GConvertYUnits(y, LINES, NFC, dd);
    return yNFCtoDev(dd->gp.plt[3] + nfc, dd);
}
double xMAR4toyDev(double x, DevDesc *dd)
{
    return yUsrtoDev(x, dd);
}
double yMAR4toxDev(double y, DevDesc *dd)
{
    double nfc = GConvertXUnits(y, LINES, NFC, dd);
    return xNFCtoDev(dd->gp.plt[1] + nfc, dd);
}

/* DEVICE coordinates to OTHER */

double xDevtoNDC(double x, DevDesc *dd)
{
    return (x - dd->gp.ndc2dev.ax) / dd->gp.ndc2dev.bx;
}
double yDevtoNDC(double y, DevDesc *dd)
{
    return (y - dd->gp.ndc2dev.ay) / dd->gp.ndc2dev.by;
}
double xDevtoInch(double x, DevDesc *dd)
{
    return xDevtoNDC(x, dd) / dd->gp.xNDCPerInch;
}
double yDevtoInch(double y, DevDesc *dd)
{
    return yDevtoNDC(y, dd) / dd->gp.yNDCPerInch;
}
double xDevtoLine(double x, DevDesc *dd)
{
    return xDevtoNDC(x, dd) / dd->gp.xNDCPerLine;
}
double yDevtoLine(double y, DevDesc *dd)
{
    return yDevtoNDC(y, dd) / dd->gp.yNDCPerLine;
}
double xDevtoNIC(double x, DevDesc *dd)
{
    return (x - dd->gp.inner2dev.ax) / dd->gp.inner2dev.bx;
}
double yDevtoNIC(double y, DevDesc *dd)
{
    return (y - dd->gp.inner2dev.ay) / dd->gp.inner2dev.by;
}
double xDevtoOMA1(double x, DevDesc *dd)
{
    return xDevtoNIC(x, dd);
}
double yDevtoOMA1(double y, DevDesc *dd)
{
    return dd->gp.oma[0] - yDevtoLine(y, dd);
}
double xDevtoyOMA2(double x, DevDesc *dd)
{
    return dd->gp.oma[1] - xDevtoLine(x, dd);
}
double yDevtoxOMA2(double y, DevDesc *dd)
{
    return yDevtoNIC(y, dd);
}
double xDevtoOMA3(double x, DevDesc *dd)
{
    return xDevtoNIC(x, dd);
}
double yDevtoOMA3(double y, DevDesc *dd)
{
    double line = (1.0 - yDevtoNDC(y, dd)) / dd->gp.yNDCPerLine;
    return dd->gp.oma[2] - line;
}
double xDevtoyOMA4(double x, DevDesc *dd)
{
    double line = (1.0 - xDevtoNDC(x, dd)) / dd->gp.xNDCPerLine;
    return dd->gp.oma[3] - line;
}
double yDevtoxOMA4(double y, DevDesc *dd)
{
    return yDevtoNIC(y, dd);
}
double xDevtoNFC(double x, DevDesc *dd)
{
    return (x - dd->gp.fig2dev.ax) / dd->gp.fig2dev.bx;
}
double yDevtoNFC(double y, DevDesc *dd)
{
    return (y - dd->gp.fig2dev.ay) / dd->gp.fig2dev.by;
}
double xDevtoNPC(double x, DevDesc *dd)
{
    return (xDevtoNFC(x, dd) - dd->gp.plt[0]) / (dd->gp.plt[1] - dd->gp.plt[0]);
}
double yDevtoNPC(double y, DevDesc *dd)
{
    return (yDevtoNFC(y, dd) - dd->gp.plt[2]) / (dd->gp.plt[3] - dd->gp.plt[2]);
}

/* a special case (NPC = normalised plot region coordinates) */

double xNPCtoUsr(double x, DevDesc *dd)
{
    if (dd->gp.xlog)
        return pow(10, dd->gp.logusr[0] + x * (dd->gp.logusr[1] - dd->gp.logusr[0]));
    else
        return dd->gp.usr[0] + x * (dd->gp.usr[1] - dd->gp.usr[0]);
}
double yNPCtoUsr(double y, DevDesc *dd)
{
    if (dd->gp.ylog)
        return pow(10, dd->gp.logusr[2] + y * (dd->gp.logusr[3] - dd->gp.logusr[2]));
    else
        return dd->gp.usr[2] + y * (dd->gp.usr[3] - dd->gp.usr[2]);
}

double xDevtoUsr(double x, DevDesc *dd)
{
    double nfc = xDevtoNFC(x, dd);
    return (nfc - dd->gp.win2fig.ax) / dd->gp.win2fig.bx;
}
double yDevtoUsr(double y, DevDesc *dd)
{
    double nfc = yDevtoNFC(y, dd);
    return (nfc - dd->gp.win2fig.ay) / dd->gp.win2fig.by;
}
double xDevtoMAR1(double x, DevDesc *dd)
{
    return xDevtoUsr(x, dd);
}
double yDevtoMAR1(double y, DevDesc *dd)
{
    return dd->gp.mar[0] - yDevtoLine(y, dd);
}
double xDevtoyMAR2(double x, DevDesc *dd)
{
    return dd->gp.mar[1] - xDevtoLine(x, dd);
}
double yDevtoxMAR2(double y, DevDesc *dd)
{
    return yDevtoUsr(y, dd);
}
double xDevtoMAR3(double x, DevDesc *dd)
{
    return xDevtoUsr(x, dd);
}
double yDevtoMAR3(double y, DevDesc *dd)
{
    double line = GConvertYUnits(1.0 - yDevtoNFC(y, dd), NFC, LINES, dd);
    return dd->gp.mar[2] - line;
}
double xDevtoyMAR4(double x, DevDesc *dd)
{
    double line = GConvertXUnits(1.0 - xDevtoNFC(x, dd), NFC, LINES, dd);
    return dd->gp.mar[3] - line;
}
double yDevtoxMAR4(double y, DevDesc *dd)
{
    return yDevtoUsr(y, dd);
}

/* the Convert function converts a LOCATION in the FROM coordinate */
/* system to a LOCATION in the TO coordinate system */

void GConvert(double *x, double *y, int from, int to, DevDesc *dd)
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

void NewFrameConfirm(); /* bring this into this file. */

/* code for layouts */

static double sum(double values[], int n, int cmValues[], int cmSum)
{
    int i;
    double s = 0;
    for (i = 0; i < n; i++)
        if ((cmSum && cmValues[i]) || (!cmSum && !cmValues[i]))
            s = s + values[i];
    return s;
}

static double sumWidths(DevDesc *dd)
{
    return sum(dd->gp.widths, dd->gp.numcols, dd->gp.cmWidths, 0);
}
static double sumCmWidths(DevDesc *dd)
{
    return sum(dd->gp.widths, dd->gp.numcols, dd->gp.cmWidths, 1);
}
static double sumHeights(DevDesc *dd)
{
    return sum(dd->gp.heights, dd->gp.numrows, dd->gp.cmHeights, 0);
}
static double sumCmHeights(DevDesc *dd)
{
    return sum(dd->gp.heights, dd->gp.numrows, dd->gp.cmHeights, 1);
}

static int tallLayout(double cmWidth, double cmHeight, DevDesc *dd)
{
    return (cmHeight / sumHeights(dd)) > (cmWidth / sumWidths(dd));
}

static void figureExtent(int *minCol, int *maxCol, int *minRow, int *maxRow, int figureNum, DevDesc *dd)
{
    int minc = -1;
    int maxc = -1;
    int minr = -1;
    int maxr = -1;
    int i, j;
    for (i = 0; i < dd->gp.numrows; i++)
        for (j = 0; j < dd->gp.numcols; j++)
            if (dd->gp.order[i][j] == figureNum)
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
                         double cmHeight, DevDesc *dd)
{
    largestRegion(width, height,
                  sum(heights, dd->gp.numrows, dd->gp.cmHeights, 0) / sum(widths, dd->gp.numcols, dd->gp.cmWidths, 0),
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

static void allCmRegions(double widths[], double heights[], double cmWidth, double cmHeight, DevDesc *dd)
{
    allocDimension(widths, cmWidth, dd->gp.numcols, dd->gp.cmWidths, 1);
    allocDimension(heights, cmHeight, dd->gp.numrows, dd->gp.cmHeights, 1);
}

static void modifyDimension(double dimension[], double multiplier, double n, int cmDimensions[])
{
    int i;
    for (i = 0; i < n; i++)
        if (!cmDimensions[i])
            dimension[i] = dimension[i] * multiplier;
}

static void modifyRegions(double widths[], double heights[], double colMultiplier, double rowMultiplier, DevDesc *dd)
{
    modifyDimension(widths, colMultiplier, dd->gp.numcols, dd->gp.cmWidths);
    modifyDimension(heights, rowMultiplier, dd->gp.numrows, dd->gp.cmHeights);
}

static void regionsWithoutRespect(double widths[], double heights[], DevDesc *dd)
{
    allocDimension(widths, sum(widths, dd->gp.numcols, dd->gp.cmWidths, 0), dd->gp.numcols, dd->gp.cmWidths, 0);
    allocDimension(heights, sum(heights, dd->gp.numrows, dd->gp.cmHeights, 0), dd->gp.numrows, dd->gp.cmHeights, 0);
}

static void regionsWithRespect(double widths[], double heights[], double cmWidth, double cmHeight, DevDesc *dd)
{
    double cm, rm;
    layoutRegion(&cm, &rm, widths, heights, cmWidth, cmHeight, dd);
    regionsWithoutRespect(widths, heights, dd);
    modifyRegions(widths, heights, cm, rm, dd);
}

static void widthsRespectingHeights(double widths[], double cmWidth, double cmHeight, DevDesc *dd)
{
    int i, j;
    int respectedCols[MAX_LAYOUT_COLS];
    double widthLeft;
    double disrespectedWidth = 0;
    for (j = 0; j < dd->gp.numcols; j++)
    {
        respectedCols[j] = 0;
        widths[j] = dd->gp.widths[j];
    }
    for (i = 0; i < dd->gp.numrows; i++)
        for (j = 0; j < dd->gp.numcols; j++)
            if (dd->gp.respect[i][j] && !dd->gp.cmWidths[j])
                respectedCols[j] = 1;
    for (j = 0; j < dd->gp.numcols; j++)
        if (!respectedCols[j])
            disrespectedWidth += dd->gp.widths[j];
    widthLeft = sumHeights(dd) * cmWidth / cmHeight - sumWidths(dd) + disrespectedWidth;
    for (j = 0; j < dd->gp.numcols; j++)
        if (!respectedCols[j])
            widths[j] = widthLeft * widths[j] / disrespectedWidth;
}

static void regionsRespectingHeight(double widths[], double heights[], double cmWidth, double cmHeight, DevDesc *dd)
{
    widthsRespectingHeights(widths, cmWidth, cmHeight, dd);
    regionsWithRespect(widths, heights, cmWidth, cmHeight, dd);
}

static void heightsRespectingWidths(double heights[], double cmWidth, double cmHeight, DevDesc *dd)
{
    int i, j;
    int respectedRows[MAX_LAYOUT_ROWS];
    double heightLeft;
    double disrespectedHeight = 0;
    for (i = 0; i < dd->gp.numrows; i++)
    {
        respectedRows[i] = 0;
        heights[i] = dd->gp.heights[i];
    }
    for (i = 0; i < dd->gp.numrows; i++)
        for (j = 0; j < dd->gp.numcols; j++)
            if (dd->gp.respect[i][j] && !dd->gp.cmHeights[i])
                respectedRows[i] = 1;
    for (i = 0; i < dd->gp.numrows; i++)
        if (!respectedRows[i])
            disrespectedHeight += dd->gp.heights[i];
    heightLeft = sumWidths(dd) * cmHeight / cmWidth - sumHeights(dd) + disrespectedHeight;
    for (i = 0; i < dd->gp.numrows; i++)
        if (!respectedRows[i])
            heights[i] = heightLeft * heights[i] / disrespectedHeight;
}

static void regionsRespectingWidth(double widths[], double heights[], double cmWidth, double cmHeight, DevDesc *dd)
{
    heightsRespectingWidths(heights, cmWidth, cmHeight, dd);
    regionsWithRespect(widths, heights, cmWidth, cmHeight, dd);
}

static void noCmRegions(double widths[], double heights[], double cmWidth, double cmHeight, DevDesc *dd)
{
    switch (dd->gp.rspct)
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

static void notAllCmRegions(double widths[], double heights[], double cmWidth, double cmHeight, DevDesc *dd)
{
    double newCmWidth, newCmHeight;
    newCmWidth = cmWidth - sumCmWidths(dd);
    newCmHeight = cmHeight - sumCmHeights(dd);
    noCmRegions(widths, heights, newCmWidth, newCmHeight, dd);
    allocDimension(widths, cmWidth, dd->gp.numcols, dd->gp.cmWidths, 1);
    allocDimension(heights, cmHeight, dd->gp.numrows, dd->gp.cmHeights, 1);
    modifyDimension(widths, newCmWidth / cmWidth, dd->gp.numcols, dd->gp.cmWidths);
    modifyDimension(heights, newCmHeight / cmHeight, dd->gp.numrows, dd->gp.cmHeights);
}

static void widthCmRegions(double widths[], double heights[], double cmWidth, double cmHeight, DevDesc *dd)
{
    allocDimension(widths, cmWidth, dd->gp.numcols, dd->gp.cmWidths, 1);
    allocDimension(heights, sumHeights(dd), dd->gp.numrows, dd->gp.cmHeights, 0);
    modifyDimension(heights, (cmHeight - sumCmHeights(dd)) / cmHeight, dd->gp.numrows, dd->gp.cmHeights);
    allocDimension(heights, cmHeight, dd->gp.numrows, dd->gp.cmHeights, 1);
}

static void heightCmRegions(double widths[], double heights[], double cmWidth, double cmHeight, DevDesc *dd)
{
    allocDimension(heights, cmHeight, dd->gp.numrows, dd->gp.cmHeights, 1);
    allocDimension(widths, sumWidths(dd), dd->gp.numcols, dd->gp.cmWidths, 0);
    modifyDimension(widths, (cmWidth - sumCmWidths(dd)) / cmWidth, dd->gp.numcols, dd->gp.cmWidths);
    allocDimension(widths, cmWidth, dd->gp.numcols, dd->gp.cmWidths, 1);
}

static int allCmWidths(DevDesc *dd)
{
    int j;
    int allCm = 1;
    for (j = 0; j < dd->gp.numcols; j++)
        if (!dd->gp.cmWidths[j])
            allCm = 0;
    return allCm;
}

static int allCmHeights(DevDesc *dd)
{
    int i;
    int allCm = 1;
    for (i = 0; i < dd->gp.numrows; i++)
        if (!dd->gp.cmHeights[i])
            allCm = 0;
    return allCm;
}

static int noCmWidths(DevDesc *dd)
{
    int j;
    int noCm = 1;
    for (j = 0; j < dd->gp.numcols; j++)
        if (dd->gp.cmWidths[j])
            noCm = 0;
    return noCm;
}

static int noCmHeights(DevDesc *dd)
{
    int i;
    int noCm = 1;
    for (i = 0; i < dd->gp.numrows; i++)
        if (dd->gp.cmHeights[i])
            noCm = 0;
    return noCm;
}

static void someCmRegions(double widths[], double heights[], double cmWidth, double cmHeight, DevDesc *dd)
{
    if (allCmWidths(dd))
        widthCmRegions(widths, heights, cmWidth, cmHeight, dd);
    else if (allCmHeights(dd))
        heightCmRegions(widths, heights, cmWidth, cmHeight, dd);
    else
        notAllCmRegions(widths, heights, cmWidth, cmHeight, dd);
}

static int allCm(DevDesc *dd)
{
    return allCmWidths(dd) && allCmHeights(dd);
}

static int noCm(DevDesc *dd)
{
    return noCmWidths(dd) && noCmHeights(dd);
}

static void layoutRegions(double widths[], double heights[], double cmWidth, double cmHeight, DevDesc *dd)
{
    int i, j;
    for (j = 0; j < dd->gp.numcols; j++)
        widths[j] = dd->gp.widths[j];
    for (i = 0; i < dd->gp.numrows; i++)
        heights[i] = dd->gp.heights[i];

    if (allCm(dd))
        allCmRegions(widths, heights, cmWidth, cmHeight, dd);
    else if (noCm(dd))
        noCmRegions(widths, heights, cmWidth, cmHeight, dd);
    else
        someCmRegions(widths, heights, cmWidth, cmHeight, dd);
}

static void subRegion(double *left, double *right, double *bottom, double *top, double mincol, double maxcol,
                      double minrow, double maxrow, double widths[], double heights[], DevDesc *dd)
{
    double totalWidth = sumRegions(widths, 0, dd->gp.numcols - 1);
    double totalHeight = sumRegions(heights, 0, dd->gp.numrows - 1);
    *left = (0.5 - totalWidth / 2) + sumRegions(widths, 0, mincol - 1);
    *right = (0.5 - totalWidth / 2) + sumRegions(widths, 0, maxcol);
    *bottom = (0.5 - totalHeight / 2) + totalHeight - sumRegions(heights, 0, maxrow);
    *top = (0.5 - totalHeight / 2) + totalHeight - sumRegions(heights, 0, minrow - 1);
}

/* a fudge for backwards compatibility (of sorts) with par(mfg) */
/* return the top-left-most row/col that the current figure */
/* occupies in the current layout */

void currentFigureLocation(int *row, int *col, DevDesc *dd)
{
    int maxcol, maxrow;
    if (dd->gp.layout)
        figureExtent(col, &maxcol, row, &maxrow, dd->gp.currentFigure, dd);
    else if (dd->gp.mfind)
    {
        *row = (dd->gp.currentFigure - 1) % dd->gp.numrows;
        *col = (dd->gp.currentFigure - 1) / dd->gp.numcols;
    }
    else
    {
        *row = (dd->gp.currentFigure - 1) / dd->gp.numrows;
        *col = (dd->gp.currentFigure - 1) % dd->gp.numcols;
    }
}

/*  mapNDC2Dev -- transformation from NDC to Dev  */
/*  Use this coordinate system for outer margin coordinates  */
/*  This must be called if the device is resized */

static void mapNDC2Dev(DevDesc *dd)
{
    dd->gp.ndc2dev.bx = dd->dp.ndc2dev.bx = (dd->gp.right - dd->gp.left);
    dd->gp.ndc2dev.ax = dd->dp.ndc2dev.ax = dd->gp.left;
    dd->gp.ndc2dev.by = dd->dp.ndc2dev.by = (dd->gp.top - dd->gp.bottom);
    dd->gp.ndc2dev.ay = dd->dp.ndc2dev.ay = dd->gp.bottom;

    /* Units Conversion */

    dd->gp.xNDCPerInch = dd->dp.xNDCPerInch = 1.0 / fabs(dd->gp.ndc2dev.bx * dd->gp.ipr[0]);
    dd->gp.yNDCPerInch = dd->dp.yNDCPerInch = 1.0 / fabs(dd->gp.ndc2dev.by * dd->gp.ipr[1]);
    dd->gp.xNDCPerChar = dd->dp.xNDCPerChar = fabs(dd->gp.cexbase * dd->gp.cra[1] * dd->gp.asp / dd->gp.ndc2dev.bx);
    dd->gp.yNDCPerChar = dd->dp.yNDCPerChar = fabs(dd->gp.cexbase * dd->gp.cra[1] / dd->gp.ndc2dev.by);
    dd->gp.xNDCPerLine = dd->dp.xNDCPerLine =
        fabs(dd->gp.mex * dd->gp.cexbase * dd->gp.cra[1] * dd->gp.asp / dd->gp.ndc2dev.bx);
    dd->gp.yNDCPerLine = dd->dp.yNDCPerLine = fabs(dd->gp.mex * dd->gp.cexbase * dd->gp.cra[1] / dd->gp.ndc2dev.by);
}

static void updateOuterMargins(DevDesc *dd)
{
    switch (dd->gp.oUnits)
    {
    case LINES:
        dd->gp.omi[0] = dd->dp.omi[0] = GConvertYUnits(dd->gp.oma[0], LINES, INCHES, dd);
        dd->gp.omi[1] = dd->dp.omi[1] = GConvertXUnits(dd->gp.oma[1], LINES, INCHES, dd);
        dd->gp.omi[2] = dd->dp.omi[2] = GConvertYUnits(dd->gp.oma[2], LINES, INCHES, dd);
        dd->gp.omi[3] = dd->dp.omi[3] = GConvertXUnits(dd->gp.oma[3], LINES, INCHES, dd);
        dd->gp.omd[0] = dd->dp.omd[0] = GConvertYUnits(dd->gp.oma[0], LINES, NDC, dd);
        dd->gp.omd[1] = dd->dp.omd[1] = GConvertXUnits(dd->gp.oma[1], LINES, NDC, dd);
        dd->gp.omd[2] = dd->dp.omd[2] = GConvertYUnits(dd->gp.oma[2], LINES, NDC, dd);
        dd->gp.omd[3] = dd->dp.omd[3] = GConvertXUnits(dd->gp.oma[3], LINES, NDC, dd);
        break;
    case INCHES:
        dd->gp.oma[0] = dd->dp.oma[0] = GConvertYUnits(dd->gp.omi[0], INCHES, LINES, dd);
        dd->gp.oma[1] = dd->dp.oma[1] = GConvertXUnits(dd->gp.omi[1], INCHES, LINES, dd);
        dd->gp.oma[2] = dd->dp.oma[2] = GConvertYUnits(dd->gp.omi[2], INCHES, LINES, dd);
        dd->gp.oma[3] = dd->dp.oma[3] = GConvertXUnits(dd->gp.omi[3], INCHES, LINES, dd);
        dd->gp.omd[0] = dd->dp.omd[0] = GConvertYUnits(dd->gp.omi[0], INCHES, NDC, dd);
        dd->gp.omd[1] = dd->dp.omd[1] = GConvertXUnits(dd->gp.omi[1], INCHES, NDC, dd);
        dd->gp.omd[2] = dd->dp.omd[2] = GConvertYUnits(dd->gp.omi[2], INCHES, NDC, dd);
        dd->gp.omd[3] = dd->dp.omd[3] = GConvertXUnits(dd->gp.omi[3], INCHES, NDC, dd);
        break;
    case NDC:
        dd->gp.oma[0] = dd->dp.oma[0] = GConvertYUnits(dd->gp.omd[0], NDC, LINES, dd);
        dd->gp.oma[1] = dd->dp.oma[1] = GConvertXUnits(dd->gp.omd[1], NDC, LINES, dd);
        dd->gp.oma[2] = dd->dp.oma[2] = GConvertYUnits(dd->gp.omd[2], NDC, LINES, dd);
        dd->gp.oma[3] = dd->dp.oma[3] = GConvertXUnits(dd->gp.omd[3], NDC, LINES, dd);
        dd->gp.omi[0] = dd->dp.omi[0] = GConvertYUnits(dd->gp.omd[0], NDC, INCHES, dd);
        dd->gp.omi[1] = dd->dp.omi[1] = GConvertXUnits(dd->gp.omd[1], NDC, INCHES, dd);
        dd->gp.omi[2] = dd->dp.omi[2] = GConvertYUnits(dd->gp.omd[2], NDC, INCHES, dd);
        dd->gp.omi[3] = dd->dp.omi[3] = GConvertXUnits(dd->gp.omd[3], NDC, INCHES, dd);
        break;
    }
}

/*  mapInner2Dev -- transformation from NIC to Dev  */
/*  Use this coordinate system for setting up multiple figures  */
/*  This is also used when specifying the figure region directly  */
/*  Note that this is incompatible with S which uses then entire  */
/*  device surface for such a plot  */
/*  This must be called per DevNewPlot, if the NDCtoDev transformation */
/*  changes, and if oma changes */

static void mapInner2Dev(DevDesc *dd)
{
    double x0, x1, y0, y1;
    x0 = xLinetoDev(dd->gp.oma[1], dd);
    y0 = yLinetoDev(dd->gp.oma[0], dd);
    x1 = GConvertXUnits(dd->gp.oma[3], LINES, NDC, dd);
    x1 = xNDCtoDev(1.0 - x1, dd);
    y1 = GConvertYUnits(dd->gp.oma[2], LINES, NDC, dd);
    y1 = yNDCtoDev(1.0 - y1, dd);
    dd->gp.inner2dev.bx = dd->dp.inner2dev.bx = x1 - x0;
    dd->gp.inner2dev.ax = dd->dp.inner2dev.ax = x0;
    dd->gp.inner2dev.by = dd->dp.inner2dev.by = y1 - y0;
    dd->gp.inner2dev.ay = dd->dp.inner2dev.ay = y0;
}

/* mapFigureRegion -- calculate figure region in NIC  */

static void mapFigureRegion(DevDesc *dd)
{
    int mincol, maxcol, minrow, maxrow;
    double x0, x1, y0, y1;
    double widths[MAX_LAYOUT_COLS], heights[MAX_LAYOUT_ROWS];
    if (dd->gp.layout)
    {
        layoutRegions(widths, heights, GConvertXUnits(1.0, NIC, INCHES, dd) * 2.54,
                      GConvertYUnits(1.0, NIC, INCHES, dd) * 2.54, dd);
        figureExtent(&mincol, &maxcol, &minrow, &maxrow, dd->gp.currentFigure, dd);
        subRegion(&x0, &x1, &y0, &y1, mincol, maxcol, minrow, maxrow, widths, heights, dd);
    }
    else
    {
        int row, col;
        if (dd->gp.mfind)
        {
            col = (dd->gp.currentFigure - 1) / dd->gp.numrows + 1;
            row = dd->gp.currentFigure - (col - 1) * dd->gp.numrows;
        }
        else
        {
            row = (dd->gp.currentFigure - 1) / dd->gp.numcols + 1;
            col = dd->gp.currentFigure - (row - 1) * dd->gp.numcols;
        }
        x0 = (double)(col - 1) / dd->gp.numcols;
        x1 = (double)col / dd->gp.numcols;
        y0 = (double)(dd->gp.numrows - row) / dd->gp.numrows;
        y1 = (double)(dd->gp.numrows - row + 1) / dd->gp.numrows;
    }
    dd->gp.fig[0] = dd->dp.fig[0] = x0;
    dd->gp.fig[1] = dd->dp.fig[1] = x1;
    dd->gp.fig[2] = dd->dp.fig[2] = y0;
    dd->gp.fig[3] = dd->dp.fig[3] = y1;
    dd->gp.fUnits = dd->dp.fUnits = NIC;
}

static void updateFigureRegion(DevDesc *dd)
{
    double nicWidth, nicHeight;
    switch (dd->gp.fUnits)
    {
    case NIC:
        dd->gp.fin[0] = dd->dp.fin[0] = GConvertXUnits(dd->gp.fig[1] - dd->gp.fig[0], NIC, INCHES, dd);
        dd->gp.fin[1] = dd->dp.fin[1] = GConvertYUnits(dd->gp.fig[3] - dd->gp.fig[2], NIC, INCHES, dd);
        break;
    case INCHES:
        nicWidth = GConvertXUnits(dd->gp.fin[0], INCHES, NIC, dd);
        nicHeight = GConvertYUnits(dd->gp.fin[1], INCHES, NIC, dd);
        dd->gp.fig[0] = dd->dp.fig[0] = 0.5 - nicWidth / 2;
        dd->gp.fig[1] = dd->dp.fig[1] = dd->gp.fig[0] + nicWidth;
        dd->gp.fig[2] = dd->dp.fig[2] = 0.5 - nicHeight / 2;
        dd->gp.fig[3] = dd->dp.fig[3] = dd->gp.fig[2] + nicHeight;
        break;
    }
}

/*  mapFig2Dev -- Transformation from NFC to Dev  */
/* This must be called per plot.new and if the NICtoDev transformation */
/* changes */

static void mapFig2Dev(DevDesc *dd)
{
    double x0, x1, y0, y1;
    y0 = yNICtoDev(dd->gp.fig[2], dd);
    y1 = yNICtoDev(dd->gp.fig[3], dd);
    x0 = xNICtoDev(dd->gp.fig[0], dd);
    x1 = xNICtoDev(dd->gp.fig[1], dd);
    dd->gp.fig2dev.bx = dd->dp.fig2dev.bx = x1 - x0;
    dd->gp.fig2dev.ax = dd->dp.fig2dev.ax = x0;
    dd->gp.fig2dev.by = dd->dp.fig2dev.by = y1 - y0;
    dd->gp.fig2dev.ay = dd->dp.fig2dev.ay = y0;
}

static void updateFigureMargins(DevDesc *dd)
{
    switch (dd->gp.mUnits)
    {
    case LINES:
        dd->gp.mai[0] = dd->dp.mai[0] = GConvertYUnits(dd->gp.mar[0], LINES, INCHES, dd);
        dd->gp.mai[1] = dd->dp.mai[1] = GConvertXUnits(dd->gp.mar[1], LINES, INCHES, dd);
        dd->gp.mai[2] = dd->dp.mai[2] = GConvertYUnits(dd->gp.mar[2], LINES, INCHES, dd);
        dd->gp.mai[3] = dd->dp.mai[3] = GConvertXUnits(dd->gp.mar[3], LINES, INCHES, dd);
        break;
    case INCHES:
        dd->gp.mar[0] = dd->dp.mar[0] = GConvertYUnits(dd->gp.mai[0], INCHES, LINES, dd);
        dd->gp.mar[1] = dd->dp.mar[1] = GConvertXUnits(dd->gp.mai[1], INCHES, LINES, dd);
        dd->gp.mar[2] = dd->dp.mar[2] = GConvertYUnits(dd->gp.mai[2], INCHES, LINES, dd);
        dd->gp.mar[3] = dd->dp.mar[3] = GConvertXUnits(dd->gp.mai[3], INCHES, LINES, dd);
        break;
    }
}

/* mapPlotRegion -- plot region in NFC */

static void mapPlotRegion(DevDesc *dd)
{
    double x0, x1, y0, y1;
    x0 = GConvertXUnits(dd->gp.mar[1], LINES, NFC, dd);
    y0 = GConvertYUnits(dd->gp.mar[0], LINES, NFC, dd);
    x1 = 1.0 - GConvertXUnits(dd->gp.mar[3], LINES, NFC, dd);
    y1 = 1.0 - GConvertYUnits(dd->gp.mar[2], LINES, NFC, dd);
    if (dd->gp.pty == 's')
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
    dd->gp.plt[0] = dd->dp.plt[0] = x0;
    dd->gp.plt[1] = dd->dp.plt[1] = x1;
    dd->gp.plt[2] = dd->dp.plt[2] = y0;
    dd->gp.plt[3] = dd->dp.plt[3] = y1;
    dd->gp.pUnits = dd->dp.pUnits = NFC;
}

static void updatePlotRegion(DevDesc *dd)
{
    double nfcWidth, nfcHeight;
    switch (dd->gp.pUnits)
    {
    case NFC:
        dd->gp.pin[0] = dd->dp.pin[0] = GConvertXUnits(dd->gp.plt[1] - dd->gp.plt[0], NFC, INCHES, dd);
        dd->gp.pin[1] = dd->dp.pin[1] = GConvertYUnits(dd->gp.plt[3] - dd->gp.plt[2], NFC, INCHES, dd);
        break;
    case INCHES:
        nfcWidth = GConvertXUnits(dd->gp.pin[0], INCHES, NFC, dd);
        nfcHeight = GConvertYUnits(dd->gp.pin[1], INCHES, NFC, dd);
        dd->gp.plt[0] = dd->dp.plt[0] = 0.5 - nfcWidth / 2;
        dd->gp.plt[1] = dd->dp.plt[1] = dd->gp.plt[0] + nfcWidth;
        dd->gp.plt[2] = dd->dp.plt[2] = 0.5 - nfcHeight / 2;
        dd->gp.plt[3] = dd->dp.plt[3] = dd->gp.plt[2] + nfcHeight;
        break;
    }
}

/*  GMapWin2Fig -- transformation from Usr to NFC */

void GMapWin2Fig(DevDesc *dd)
{
    if (dd->gp.xlog)
    {
        dd->gp.win2fig.bx = dd->dp.win2fig.bx = (dd->gp.plt[1] - dd->gp.plt[0]) / (dd->gp.logusr[1] - dd->gp.logusr[0]);
        dd->gp.win2fig.ax = dd->dp.win2fig.ax = dd->gp.plt[0] - dd->gp.win2fig.bx * dd->gp.logusr[0];
    }
    else
    {
        dd->gp.win2fig.bx = dd->dp.win2fig.bx = (dd->gp.plt[1] - dd->gp.plt[0]) / (dd->gp.usr[1] - dd->gp.usr[0]);
        dd->gp.win2fig.ax = dd->dp.win2fig.ax = dd->gp.plt[0] - dd->gp.win2fig.bx * dd->gp.usr[0];
    }
    if (dd->gp.ylog)
    {
        dd->gp.win2fig.by = dd->dp.win2fig.by = (dd->gp.plt[3] - dd->gp.plt[2]) / (dd->gp.logusr[3] - dd->gp.logusr[2]);
        dd->gp.win2fig.ay = dd->dp.win2fig.ay = dd->gp.plt[2] - dd->gp.win2fig.by * dd->gp.logusr[2];
    }
    else
    {
        dd->gp.win2fig.by = dd->dp.win2fig.by = (dd->gp.plt[3] - dd->gp.plt[2]) / (dd->gp.usr[3] - dd->gp.usr[2]);
        dd->gp.win2fig.ay = dd->dp.win2fig.ay = dd->gp.plt[2] - dd->gp.win2fig.by * dd->gp.usr[2];
    }
}

/*  mapping -- Set up mappings between coordinate systems  */
/*  This is the user's interface to the mapping routines above */

void mapping(DevDesc *dd, int which)
{
    switch (which)
    {
    case 0:
        mapNDC2Dev(dd);
    case 1:
        updateOuterMargins(dd);
        mapInner2Dev(dd);
    case 2:
        if (dd->gp.defaultFigure)
            mapFigureRegion(dd);
        updateFigureRegion(dd);
        mapFig2Dev(dd);
    case 3:
        updateFigureMargins(dd);
        if (dd->gp.defaultPlot)
            mapPlotRegion(dd);
        updatePlotRegion(dd);
    }
}

/*  GReset -- Reset coordinate systems mappings and unit yardsticks */

void GReset(DevDesc *dd)
{
    /* Character extents are based on the raster size */
    dd->gp.asp = dd->gp.ipr[1] / dd->gp.ipr[0];
    dd->gp.mkh = dd->gp.cra[0] * dd->gp.ipr[0];

    /* Recompute Mappings */
    mapping(dd, 0);
}

/*  Is the figure region too big ? */

int validFigureRegion(DevDesc *dd)
{
    return ((dd->gp.fig[0] > 0 - FLT_EPSILON) && (dd->gp.fig[1] < 1 + FLT_EPSILON) &&
            (dd->gp.fig[2] > 0 - FLT_EPSILON) && (dd->gp.fig[3] < 1 + FLT_EPSILON));
}

/*  Is the figure region too small ? */

int validOuterMargins(DevDesc *dd)
{
    return ((dd->gp.fig[0] < dd->gp.fig[1]) && (dd->gp.fig[2] < dd->gp.fig[3]));
}

/* Is the plot region too big ? */

int validPlotRegion(DevDesc *dd)
{
    return ((dd->gp.plt[0] > 0 - FLT_EPSILON) && (dd->gp.plt[1] < 1 + FLT_EPSILON) &&
            (dd->gp.plt[2] > 0 - FLT_EPSILON) && (dd->gp.plt[3] < 1 + FLT_EPSILON));
}

/* Is the plot region too small ? */

int validFigureMargins(DevDesc *dd)
{
    return ((dd->gp.plt[0] < dd->gp.plt[1]) && (dd->gp.plt[2] < dd->gp.plt[3]));
}

void initDisplayList();

static void invalidError(char *message, DevDesc *dd)
{
    dd->dp.currentFigure -= 1;
    if (dd->dp.currentFigure < 1)
        dd->dp.currentFigure = dd->dp.lastFigure;
    dd->gp.currentFigure = dd->dp.currentFigure;
    error(message);
}

/*  GNewPlot -- Begin a new plot (advance to new frame if needed)  */

void GNewPlot(DevDesc *dd, int recording)
{
    if (NoDevices())
        error("No graphics device is active\n");

    /* Restore Default Parameters */
    GRestore(dd);

    dd->dp.currentFigure += 1;
    dd->gp.currentFigure = dd->dp.currentFigure;
    if (dd->gp.currentFigure > dd->gp.lastFigure)
    {
        if (!dd->gp.new)
        {
            if (dd->gp.ask && recording)
                NewFrameConfirm();
            if (recording)
                initDisplayList(dd);
            dd->dp.newPage(dd);
        }
        dd->gp.new = dd->gp.new = 1;
        dd->dp.currentFigure = dd->gp.currentFigure = 1;
    }

    GReset(dd);
    GForceClip(dd);

    /* IF the division of the device into separate regions */
    /* has resulted in any invalid regions ... */
    /* IF this was a user command (i.e., we are recording) */
    /* send an error message to the command line */
    /* IF we are replaying then draw a message in the output */

    dd->dp.valid = dd->gp.valid = 0;
    if (!validOuterMargins(dd))
        if (recording)
            invalidError("Outer margins too large\n", dd);
        else
            GText(0.5, 0.5, NFC, "Outer margins too large", 0.5, 0.5, 0, dd);
    else if (!validFigureRegion(dd))
        if (recording)
            invalidError("Figure region too large\n", dd);
        else
            GText(0.5, 0.5, NFC, "Figure region too large", 0.5, 0.5, 0, dd);
    else if (!validFigureMargins(dd))
        if (recording)
            invalidError("Figure margins too large\n", dd);
        else
            GText(0.5, 0.5, NFC, "Figure margins too large", 0.5, 0.5, 0, dd);
    else if (!validPlotRegion(dd))
        if (recording)
            invalidError("Plot region too large\n", dd);
        else
            GText(0.5, 0.5, NFC, "Plot region too small / large", 0.5, 0.5, 0, dd);
    else
        dd->dp.valid = dd->gp.valid = 1;
}

/* used to default axis information */
/* i.e., if user hasn't specified par(usr=...) */

void GScale(double min, double max, int axis, DevDesc *dd)
{
    int log, n, style, swap;
    double temp;

    if (axis == 1 || axis == 3)
    {
        n = dd->gp.lab[0];
        style = dd->gp.xaxs;
        log = dd->gp.xlog;
    }
    else
    {
        n = dd->gp.lab[1];
        style = dd->gp.yaxs;
        log = dd->gp.ylog;
    }

    if (log)
    {
        min = log10(min);
        max = log10(max);
    }

    if (min == max)
    {
        if (min == 0)
        {
            min = -1;
            max = 1;
        }
        else
        {
            min = 0.6 * min;
            max = 1.4 * max;
        }
    }

    switch (style)
    {
    case 'r':
        temp = 0.04 * (max - min);
        min = min - temp;
        max = max + temp;
        break;
    case 'i':
        break;
    case 's':
    case 'e':
    default:
        error("axis style \"%c\" unimplemented\n", style);
    }

    if (axis == 1 || axis == 3)
    {
        if (log)
        {
            dd->gp.usr[0] = dd->dp.usr[0] = pow(10, min);
            dd->gp.usr[1] = dd->dp.usr[1] = pow(10, max);
            dd->gp.logusr[0] = dd->dp.logusr[0] = min;
            dd->gp.logusr[1] = dd->dp.logusr[1] = max;
        }
        else
        {
            dd->gp.usr[0] = dd->dp.usr[0] = min;
            dd->gp.usr[1] = dd->dp.usr[1] = max;
            dd->gp.logusr[0] = dd->dp.logusr[0] = log10(min);
            dd->gp.logusr[1] = dd->dp.logusr[1] = log10(max);
        }
    }
    else
    {
        if (log)
        {
            dd->gp.usr[2] = dd->dp.usr[2] = pow(10, min);
            dd->gp.usr[3] = dd->dp.usr[3] = pow(10, max);
            dd->gp.logusr[2] = dd->dp.logusr[2] = min;
            dd->gp.logusr[3] = dd->dp.logusr[3] = max;
        }
        else
        {
            dd->gp.usr[2] = dd->dp.usr[2] = min;
            dd->gp.usr[3] = dd->dp.usr[3] = max;
            dd->gp.logusr[2] = dd->dp.logusr[2] = log10(min);
            dd->gp.logusr[3] = dd->dp.logusr[3] = log10(max);
        }
    }

    if (min > max)
    {
        swap = 1;
        temp = min;
        min = max;
        max = temp;
    }
    else
        swap = 0;

    if (log)
    {
        min = pow(10.0, min);
        max = pow(10.0, max);
        GLPretty(&min, &max, &n);
    }
    else
        GPretty(&min, &max, &n);

    if (swap)
    {
        temp = min;
        min = max;
        max = temp;
    }

    if (axis == 1 || axis == 3)
    {
        dd->gp.xaxp[0] = dd->dp.xaxp[0] = min;
        dd->gp.xaxp[1] = dd->dp.xaxp[1] = max;
        dd->gp.xaxp[2] = dd->dp.xaxp[2] = n;
    }
    else
    {
        dd->gp.yaxp[0] = dd->dp.yaxp[0] = min;
        dd->gp.yaxp[1] = dd->dp.yaxp[1] = max;
        dd->gp.yaxp[2] = dd->dp.yaxp[2] = n;
    }
}

/*  GSetupAxis -- Set up the default axis information  */
/*  called when user specifies par(usr=...) */
/*  What should happen if logx = 1 or logy = 1 ? */

void GSetupAxis(int axis, DevDesc *dd)
{
    double min, max;
    int n;

    if (axis == 1 || axis == 3)
    {
        n = dd->gp.lab[0];
        min = dd->gp.usr[0];
        max = dd->gp.usr[1];
    }
    else
    {
        n = dd->gp.lab[1];
        min = dd->gp.usr[2];
        max = dd->gp.usr[3];
    }

    GPretty(&min, &max, &n);

    if (axis == 1 || axis == 3)
    {
        dd->gp.xaxp[0] = dd->dp.xaxp[0] = min;
        dd->gp.xaxp[1] = dd->dp.xaxp[1] = max;
        dd->gp.xaxp[2] = dd->dp.xaxp[2] = n;
    }
    else
    {
        dd->gp.yaxp[0] = dd->dp.yaxp[0] = min;
        dd->gp.yaxp[1] = dd->dp.yaxp[1] = max;
        dd->gp.yaxp[2] = dd->dp.yaxp[2] = n;
    }
}

/********************************************************/
/* GPAR FUNCTIONS					*/
/********************************************************/

void GInit(GPar *gp)
{
    /*  Set default graphics parameter values in  */
    /*  a GPar This initialises the plot state,   */
    /*  plus the other graphical parameters that  */
    /*  are not the responsibility of the device  */
    /*  initialisation                            */

    gp->state = 0;

    gp->ann = 1;
    gp->ask = 0;
    gp->err = 0;
    gp->bty = 'o';

    gp->mkh = .001; /* dummy value > 0  --- FIXME : */
    /* GREset has dd->gp.mkh = dd->gp.cra[0] * dd->gp.ipr[0]; */
    gp->cex = 1.0;
    gp->cexbase = 1.0;
    gp->cexmain = 1.2;
    gp->cexlab = 1.0;
    gp->cexsub = 1.0;
    gp->cexaxis = 1.0;

    gp->col = 0;
    gp->colmain = 0;
    gp->collab = 0;
    gp->colsub = 0;
    gp->colaxis = 0;
    gp->gamma = 1;

    /* gp->ps = 10; */ /* Device Specific */
    gp->metricInfo = 0;
    gp->font = 1;
    gp->fontmain = 2;
    gp->fontlab = 1;
    gp->fontsub = 1;
    gp->fontaxis = 1;

    gp->pch = 1;
    gp->lty = LTY_SOLID;
    gp->smo = 1;

    /* String Adjustment and rotation */
    gp->adj = 0.5;
    gp->crt = 0.0;
    gp->srt = 0.0;

    /* Positioning of margin text */
    gp->mgp[0] = 3;
    gp->mgp[1] = 1;
    gp->mgp[2] = 0;

    /* Axis annotation parameters */
    gp->lab[0] = 5;
    gp->lab[1] = 5;
    gp->lab[2] = 7;
    gp->las = 0;
    gp->tck = NA_REAL;
    gp->tcl = -0.5;
    gp->tmag = 1.2;
    gp->type = 'p';
    gp->xaxp[0] = 0.0;
    gp->xaxp[1] = 1.0;
    gp->xaxp[2] = 5.0;
    gp->xaxs = 'r';
    gp->xaxt = 's';
    gp->xlog = 0;
    gp->xpd = 0;
    gp->oldxpd = -99;
    gp->yaxp[0] = 0.0;
    gp->yaxp[1] = 1.0;
    gp->yaxp[2] = 5.0;
    gp->yaxs = 'r';
    gp->yaxt = 's';
    gp->ylog = 0;

    /* Outer Margins */
    gp->mex = 1.0;
    gp->oma[0] = 0.0;
    gp->oma[1] = 0.0;
    gp->oma[2] = 0.0;
    gp->oma[3] = 0.0;
    gp->oUnits = LINES;
    gp->fig[0] = 0.0;
    gp->fig[1] = 1.0;
    gp->fig[2] = 0.0;
    gp->fig[3] = 1.0;
    gp->fUnits = NIC;
    gp->defaultFigure = 1; /* the figure region is calculated from */
                           /* the layout by default */
    gp->pUnits = NFC;
    gp->defaultPlot = 1; /* the plot region is calculated as */
                         /* figure-margin by default */

    /* Inner Margins */
    gp->mar[0] = 5.1;
    gp->mar[1] = 4.1;
    gp->mar[2] = 4.1;
    gp->mar[3] = 2.1;
    gp->mUnits = LINES;

    /* Multi-figure parameters */
    gp->layout = 0;
    gp->mfind = 0;

    gp->numrows = 1;
    gp->numcols = 1;
    gp->currentFigure = 1;
    gp->lastFigure = 1;
    gp->heights[0] = 1;
    gp->widths[0] = 1;
    gp->cmHeights[0] = 0;
    gp->cmWidths[0] = 0;
    gp->order[0][0] = 1;
    gp->rspct = 0;
    gp->respect[0][0] = 0;

    /* Misc plotting parameters */
    gp->new = 0;
    gp->devmode = -99;
    gp->pty = 'm';
    gp->lwd = 1;

    /* Data window */
    gp->usr[0] = 0.0;
    gp->usr[1] = 1.0;
    gp->usr[2] = 0.0;
    gp->usr[3] = 1.0;
}

void copyGPar(GPar *source, GPar *dest)
{
    /* copy a GPar structure from source to dest */

    memcpy(dest, source, sizeof(GPar));
}

void GRestore(DevDesc *dd)
{
    /*  Restore the graphics parameters from the device copy  */

    if (NoDevices())
        error("No graphics device is active\n");
    copyGPar(&(dd->dp), &(dd->gp));
}

/*  Saving and restoring of "inline" graphical  */
/*  parameters.  These are the ones which can be  */
/*  specified as a arguments to high-level  */
/*  graphics functions.  */

static double adjsave;     /* adj */
static int annsave;        /* ann */
static int btysave;        /* bty */
static double cexsave;     /* cex */
static double cexbasesave; /* cexbase */
static double cexmainsave; /* cex.main */
static double cexlabsave;  /* cex.lab */
static double cexsubsave;  /* cex.sub */
static double cexaxissave; /* cex.axis */
static int colsave;        /* col */
static int fgsave;         /* fg */
static int bgsave;         /* bg */
static int colmainsave;    /* col.main */
static int collabsave;     /* col.lab */
static int colsubsave;     /* col.sub */
static int colaxissave;    /* col.axis */
static double crtsave;     /* character rotation */
static int fontsave;       /* font */
static int fontmainsave;   /* font.main */
static int fontlabsave;    /* font.lab */
static int fontsubsave;    /* font.sub */
static int fontaxissave;   /* font.axis */
#ifdef NO
static int csisave; /* line spacing in inches */
#endif
static int errsave;        /* error mode */
static int labsave[3];     /* axis labelling parameters */
static int lassave;        /* label style */
static int ltysave;        /* line type */
static double lwdsave;     /* line width */
static double mgpsave[3];  /* margin position for annotation */
static double mkhsave;     /* mark height */
static int pchsave;        /* plotting character */
static double srtsave;     /* string rotation */
static double tcksave;     /* tick mark length */
static double xaxpsave[3]; /* x axis parameters */
static int xaxssave;       /* x axis calculation style */
static int xaxtsave;       /* x axis type */
static int xpdsave;        /* clipping control */
static double yaxpsave[3]; /* y axis parameters */
static int yaxssave;       /* y axis calculation style */
static int yaxtsave;       /* y axis type */

void GSavePars(DevDesc *dd)
{
    /*  Make a temporary copy of the inline parameter values  */

    adjsave = dd->gp.adj;
    annsave = dd->gp.ann;
    btysave = dd->gp.bty;
    cexsave = dd->gp.cex;
    cexbasesave = dd->gp.cexbase;
    cexlabsave = dd->gp.cexlab;
    cexmainsave = dd->gp.cexmain;
    cexsubsave = dd->gp.cexsub;
    cexaxissave = dd->gp.cexaxis;
    colsave = dd->gp.col;
    fgsave = dd->gp.fg;
    bgsave = dd->gp.bg;
    collabsave = dd->gp.collab;
    colmainsave = dd->gp.colmain;
    colsubsave = dd->gp.colsub;
    colaxissave = dd->gp.colaxis;
    crtsave = dd->gp.crt;
    errsave = dd->gp.err;
    fontsave = dd->gp.font;
    fontmainsave = dd->gp.fontmain;
    fontlabsave = dd->gp.fontlab;
    fontsubsave = dd->gp.fontsub;
    fontaxissave = dd->gp.fontaxis;
    /* csisave = dd->gp.csi; */
    labsave[0] = dd->gp.lab[0];
    labsave[1] = dd->gp.lab[1];
    labsave[2] = dd->gp.lab[2];
    lassave = dd->gp.las;
    ltysave = dd->gp.lty;
    lwdsave = dd->gp.lwd;
    mgpsave[0] = dd->gp.mgp[0];
    mgpsave[1] = dd->gp.mgp[1];
    mgpsave[2] = dd->gp.mgp[2];
    mkhsave = dd->gp.mkh;
    pchsave = dd->gp.pch;
    srtsave = dd->gp.srt;
    tcksave = dd->gp.tck;
    xaxpsave[0] = dd->gp.xaxp[0];
    xaxpsave[1] = dd->gp.xaxp[1];
    xaxpsave[2] = dd->gp.xaxp[2];
    xaxssave = dd->gp.xaxs;
    xaxtsave = dd->gp.xaxt;
    xpdsave = dd->gp.xpd;
    yaxpsave[0] = dd->gp.yaxp[0];
    yaxpsave[1] = dd->gp.yaxp[1];
    yaxpsave[2] = dd->gp.yaxp[2];
    yaxssave = dd->gp.yaxs;
    yaxtsave = dd->gp.yaxt;
}

void GRestorePars(DevDesc *dd)
{
    /*  Restore temorarily saved inline parameter values  */

    dd->gp.adj = adjsave;
    dd->gp.ann = annsave;
    dd->gp.bty = btysave;
    dd->gp.cex = cexsave;
    dd->gp.cexbase = cexbasesave;
    dd->gp.cexlab = cexlabsave;
    dd->gp.cexmain = cexmainsave;
    dd->gp.cexsub = cexsubsave;
    dd->gp.cexaxis = cexaxissave;
    dd->gp.col = colsave;
    dd->gp.fg = fgsave;
    dd->gp.bg = bgsave;
    dd->gp.collab = collabsave;
    dd->gp.colmain = colmainsave;
    dd->gp.colsub = colsubsave;
    dd->gp.colaxis = colaxissave;
    dd->gp.crt = crtsave;
    dd->gp.err = errsave;
    dd->gp.font = fontsave;
    dd->gp.fontmain = fontmainsave;
    dd->gp.fontlab = fontlabsave;
    dd->gp.fontsub = fontsubsave;
    dd->gp.fontaxis = fontaxissave;
    /* dd->gp.csi = csisave; */
    dd->gp.lab[0] = labsave[0];
    dd->gp.lab[1] = labsave[1];
    dd->gp.lab[2] = labsave[2];
    dd->gp.las = lassave;
    dd->gp.lty = ltysave;
    dd->gp.lwd = lwdsave;
    dd->gp.mgp[0] = mgpsave[0];
    dd->gp.mgp[1] = mgpsave[1];
    dd->gp.mgp[2] = mgpsave[2];
    dd->gp.mkh = mkhsave;
    dd->gp.pch = pchsave;
    dd->gp.srt = srtsave;
    dd->gp.tck = tcksave;
    dd->gp.xaxp[0] = xaxpsave[0];
    dd->gp.xaxp[1] = xaxpsave[1];
    dd->gp.xaxp[2] = xaxpsave[2];
    dd->gp.xaxs = xaxssave;
    dd->gp.xaxt = xaxtsave;
    dd->gp.xpd = xpdsave;
    dd->gp.yaxp[0] = yaxpsave[0];
    dd->gp.yaxp[1] = yaxpsave[1];
    dd->gp.yaxp[2] = yaxpsave[2];
    dd->gp.yaxs = yaxssave;
    dd->gp.yaxt = yaxtsave;
}

/********************************************************/
/* DEVICE STATE FUNCTIONS				*/
/********************************************************/

void GSetState(int newstate, DevDesc *dd)
{
    /*  This records whether GNewPlot has been called  */

    dd->dp.state = dd->gp.state = newstate;
}

void GCheckState(DevDesc *dd)
{
    /*  Enquire whether GNewPlot has been called  */

    if (dd->gp.state == 0)
        error("plot.new has not been called yet\n");
    if (dd->gp.valid == 0)
        onintr();
}

/********************************************************/
/* GRAPHICAL PRIMITIVES 				*/
/********************************************************/

void GCircle(double x, double y, int coords, double radius, int col, int border, DevDesc *dd)
{
    /*  Draw a circle (radius is given in inches)  */

    double ir;
#ifdef POINTS
    ir = radius / (72.0 * dd->gp.ipr[0]);
#endif
    ir = radius / dd->gp.ipr[0];
    dd->dp.circle(x, y, coords, ir, col, border, dd);
}

void GClip(DevDesc *dd)
{
    /*  Update the device clipping region (depends on GP->xpd)  */

    if (dd->gp.xpd != dd->gp.oldxpd)
    {
        if (dd->gp.xpd)
            dd->dp.clip(dd->gp.left, dd->gp.right, dd->gp.bottom, dd->gp.top, dd);
        else
            dd->dp.clip(xNFCtoDev(dd->gp.plt[0], dd), xNFCtoDev(dd->gp.plt[1], dd), yNFCtoDev(dd->gp.plt[2], dd),
                        yNFCtoDev(dd->gp.plt[3], dd), dd);
        dd->gp.oldxpd = dd->gp.xpd;
    }
}

void GForceClip(DevDesc *dd)
{
    /*  Forced update of the device clipping region  */

    if (dd->gp.xpd)
        dd->dp.clip(dd->gp.left, dd->gp.right, dd->gp.bottom, dd->gp.top, dd);
    else
        dd->dp.clip(xNFCtoDev(dd->gp.plt[0], dd), xNFCtoDev(dd->gp.plt[1], dd), yNFCtoDev(dd->gp.plt[2], dd),
                    yNFCtoDev(dd->gp.plt[3], dd), dd);
}

/* Draw Line Segments, Clipping to the Viewport */
/* Cohen-Sutherland Algorithm */
/* Unneeded if the device can do the clipping */

#define CS_BOTTOM 001
#define CS_LEFT 002
#define CS_TOP 004
#define CS_RIGHT 010

static int clipcode(double x, double y, double left, double right, double bottom, double top, DevDesc *dd)
{
    int c = 0;
    if (x < left)
        c |= CS_LEFT;
    else if (x > right)
        c |= CS_RIGHT;
    if (y < bottom)
        c |= CS_BOTTOM;
    else if (y > top)
        c |= CS_TOP;
    return c;
}

static int CSclipline(double *x1, double *y1, double *x2, double *y2, int *clipped1, int *clipped2, int coords,
                      DevDesc *dd)
{
    int c, c1, c2;
    double x, y, temp;
    double Clipxl = 0;
    double Clipxr = 1;
    double Clipyb = 0;
    double Clipyt = 1;

    if (dd->gp.xpd)
    {
        GConvert(&Clipxl, &Clipyb, NDC, coords, dd);
        GConvert(&Clipxr, &Clipyt, NDC, coords, dd);
    }
    else
    {
        GConvert(&Clipxl, &Clipyb, NPC, coords, dd);
        GConvert(&Clipxr, &Clipyt, NPC, coords, dd);
    }
    if (Clipxr < Clipxl)
    {
        temp = Clipxl;
        Clipxl = Clipxr;
        Clipxr = temp;
    }
    if (Clipyt < Clipyb)
    {
        temp = Clipyb;
        Clipyb = Clipyt;
        Clipyt = temp;
    }
    *clipped1 = 0;
    *clipped2 = 0;

    c1 = clipcode(*x1, *y1, Clipxl, Clipxr, Clipyb, Clipyt, dd);
    c2 = clipcode(*x2, *y2, Clipxl, Clipxr, Clipyb, Clipyt, dd);
    x = Clipxl; /* keep -Wall happy */
    y = Clipyb; /* keep -Wall happy */
    while (c1 || c2)
    {
        if (c1 & c2)
            return 0;
        if (c1)
            c = c1;
        else
            c = c2;
        if (c & CS_LEFT)
        {
            y = *y1 + (*y2 - *y1) * (Clipxl - *x1) / (*x2 - *x1);
            x = Clipxl;
        }
        else if (c & CS_RIGHT)
        {
            y = *y1 + (*y2 - *y1) * (Clipxr - *x1) / (*x2 - *x1);
            x = Clipxr;
        }
        else if (c & CS_BOTTOM)
        {
            x = *x1 + (*x2 - *x1) * (Clipyb - *y1) / (*y2 - *y1);
            y = Clipyb;
        }
        else if (c & CS_TOP)
        {
            x = *x1 + (*x2 - *x1) * (Clipyt - *y1) / (*y2 - *y1);
            y = Clipyt;
        }

        if (c == c1)
        {
            *x1 = x;
            *y1 = y;
            *clipped1 = 1;
            c1 = clipcode(x, y, Clipxl, Clipxr, Clipyb, Clipyt, dd);
        }
        else
        {
            *x2 = x;
            *y2 = y;
            *clipped2 = 1;
            c2 = clipcode(x, y, Clipxl, Clipxr, Clipyb, Clipyt, dd);
        }
    }

    return 1;
}

static void CScliplines(int n, double *x, double *y, int coords, DevDesc *dd)
{
    int lineStart, lineEnd;
    int firstPoint = 1;
    int count = 0;
    int i = 0;
    double *xx, *yy;
    double x1, y1, x2, y2;

    yy = (double *)NULL; /* keep -Wall happy */
    if (((xx = (double *)C_alloc(n, sizeof(double))) == NULL) || ((yy = (double *)C_alloc(n, sizeof(double))) == NULL))
        error("out of memory while clipping polyline\n");

    x1 = x[i];
    y1 = y[i++];
    while (i < n)
    {
        x2 = x[i];
        y2 = y[i++];
        if (CSclipline(&x1, &y1, &x2, &y2, &lineStart, &lineEnd, coords, dd))
        {
            if (firstPoint)
            {
                xx[count] = x1;
                yy[count++] = y1;
                firstPoint = 0;
            }
            else if (lineStart)
            {
                xx[count] = x1;
                yy[count++] = y1;
            }
            xx[count] = x2;
            yy[count++] = y2;
            if (lineEnd || (i == n))
            {
                if (count > 1)
                    dd->dp.polyline(count, xx, yy, coords, dd);
                count = 0;
            }
        }
        /*
                x1 = x[i-1]; y1 = y[i-1];
        */
        x1 = x2;
        y1 = y2;
    }

    C_free((char *)xx);
    C_free((char *)yy);
}

void GLine(double x1, double y1, double x2, double y2, int coords, DevDesc *dd)
{
    /*  Draw a line */

    int lineStart, lineEnd;

    if (dd->dp.canClip)
    {
        GClip(dd);
        dd->dp.line(x1, y1, x2, y2, coords, dd);
    }
    else if (CSclipline(&x1, &y1, &x2, &y2, &lineStart, &lineEnd, coords, dd))
        dd->dp.line(x1, y1, x2, y2, coords, dd);
}

int GLocator(double *x, double *y, int coords, DevDesc *dd)
{
    /*  Read the current pen position */

    if (!dd->dp.locator)
        error("no locator capability in device driver\n");
    if (dd->dp.locator(x, y, dd))
    {
        GConvert(x, y, DEVICE, coords, dd);
        return 1;
    }
    else
        return 0;
}

void GMetricInfo(int c, double *ascent, double *descent, double *width, int units, DevDesc *dd)
{
    /* Access character font metric information  */

    if (dd->dp.metricInfo)
        dd->dp.metricInfo(c, ascent, descent, width, dd);
    else
        error("detailed character metric information unavailable\n");

    if (units != DEVICE)
    {
        *ascent = GConvertYUnits(*ascent, DEVICE, units, dd);
        *descent = GConvertYUnits(*descent, DEVICE, units, dd);
        *width = GConvertXUnits(*width, DEVICE, units, dd);
    }
}

void GMode(DevDesc *dd, int mode)
{
    /*  Check that everything is initialized  */
    /*  Interpretation  */
    /*  mode = 0, graphics off */
    /*  mode = 1, graphics on */
    /*  mode = 2, graphical input on */
    /*  (Ignored by most drivers)  */

    if (NoDevices())
        error("No graphics device is active\n");
    if (mode != dd->gp.devmode)
        dd->dp.mode(mode);
    dd->gp.new = dd->dp.new = 0;
    dd->gp.devmode = dd->dp.devmode = mode;
}

void GPolygon(int n, double *x, double *y, int coords, int bg, int fg, DevDesc *dd)
{
    /*  GPolygon -- Draw a polygon  */
    /*  Filled with color bg and outlined with color fg  */
    /*  These may both be NA_INTEGER  */
    /*  If device can't clip we should use something */
    /*  like Sutherland-Hodgman here */

    if (dd->dp.canClip)
        GClip(dd);
    dd->dp.polygon(n, x, y, coords, bg, fg, dd);
}

void GPolyline(int n, double *x, double *y, int coords, DevDesc *dd)
{
    /*  Draw a series of line segments */

    if (dd->dp.canClip)
    {
        GClip(dd);
        dd->dp.polyline(n, x, y, coords, dd);
    }
    else
        CScliplines(n, x, y, coords, dd);
}

void GRect(double x0, double y0, double x1, double y1, int coords, int bg, int fg, DevDesc *dd)
{
    /*  Draw a rectangle  */
    /*  Filled with color bg and outlined with color fg  */
    /*  These may both be NA_INTEGER  */
    /*  This should be adjusted so that rectangles are clipped */

    dd->dp.rect(x0, y0, x1, y1, coords, bg, fg, dd);
}

double GStrWidth(char *str, int units, DevDesc *dd)
{
    /*  Compute string width  */

    double w = dd->dp.strWidth(str, dd);
    if (units != DEVICE)
        w = GConvertXUnits(w, DEVICE, units, dd);
    return w;
}

double GStrHeight(char *str, int units, DevDesc *dd)
{
    /*  Compute string height  */

    double h = dd->gp.cex * dd->gp.cra[1];
    if (units != DEVICE)
        h = GConvertYUnits(h, DEVICE, units, dd);
    return h;
}

void GText(double x, double y, int coords, char *str, double xc, double yc, double rot, DevDesc *dd)
{
    /*  Draw text in a plot  */

    if (str && *str)
        if (dd->dp.canClip)
        {
            GClip(dd);
            dd->dp.text(x, y, coords, str, xc, yc, rot, dd);
        }
        else
        {
            if (!dd->dp.xpd)
            {
                double xtest = x;
                double ytest = y;
                GConvert(&xtest, &ytest, coords, NFC, dd);
                if (xtest < 0 || ytest < 0 || xtest > 1 || ytest > 1)
                    return;
            }
            dd->dp.text(x, y, coords, str, xc, yc, rot, dd);
        }
}

/********************************************************/
/* GRAPHICAL UTILITIES					*/
/********************************************************/

void GArrow(double xfrom, double yfrom, double xto, double yto, int coords, double length, double angle, int code,
            DevDesc *dd)
{
    /*  GArrow -- Draw an arrow  */
    /*  NOTE that the length parameter is in inches */

    double xfromInch = xfrom;
    double yfromInch = yfrom;
    double xtoInch = xto;
    double ytoInch = yto;
    double rot, xc, yc;
    double x[3], y[3];

    GLine(xfrom, yfrom, xto, yto, coords, dd);

    GConvert(&xfromInch, &yfromInch, coords, INCHES, dd);
    GConvert(&xtoInch, &ytoInch, coords, INCHES, dd);
    if (code & 1)
    {
        xc = xtoInch - xfromInch;
        yc = ytoInch - yfromInch;
        rot = atan2(yc, xc);
        x[0] = xfromInch + length * cos(rot + angle * DEG2RAD);
        y[0] = yfromInch + length * sin(rot + angle * DEG2RAD);
        x[1] = xfromInch;
        y[1] = yfromInch;
        x[2] = xfromInch + length * cos(rot - angle * DEG2RAD);
        y[2] = yfromInch + length * sin(rot - angle * DEG2RAD);
        GPolyline(3, x, y, INCHES, dd);
    }
    if (code & 2)
    {
        xc = xfromInch - xtoInch;
        yc = yfromInch - ytoInch;
        rot = atan2(yc, xc);
        x[0] = xtoInch + length * cos(rot + angle * DEG2RAD);
        y[0] = ytoInch + length * sin(rot + angle * DEG2RAD);
        x[1] = xtoInch;
        y[1] = ytoInch;
        x[2] = xtoInch + length * cos(rot - angle * DEG2RAD);
        y[2] = ytoInch + length * sin(rot - angle * DEG2RAD);
        GPolyline(3, x, y, INCHES, dd);
    }
}

void GBox(int which, DevDesc *dd)
{
    /*  Draw a box about one of several regions  */

    double x[5], y[5];
    GClip(dd);
    if (which == 1)
    {
        x[0] = dd->gp.plt[1];
        y[0] = dd->gp.plt[2];
        x[1] = dd->gp.plt[1];
        y[1] = dd->gp.plt[3];
        x[2] = dd->gp.plt[0];
        y[2] = dd->gp.plt[3];
        x[3] = dd->gp.plt[0];
        y[3] = dd->gp.plt[2];
        x[4] = dd->gp.plt[1];
        y[4] = dd->gp.plt[2];
    }
    else
    {
        x[0] = 0.0;
        y[0] = 0.0;
        x[1] = 1.0;
        y[1] = 0.0;
        x[2] = 1.0;
        y[2] = 1.0;
        x[3] = 0.0;
        y[3] = 1.0;
    }
    switch (which)
    {
    case 1: /* Plot */
        switch (dd->gp.bty)
        {
        case 'o':
        case 'O':
            GPolygon(4, x, y, NFC, NA_INTEGER, dd->gp.col, dd);
            break;
        case 'l':
        case 'L':
            GPolyline(3, x + 2, y + 2, NFC, dd);
            break;
        case '7':
            GPolyline(3, x, y, NFC, dd);
            break;
        case 'c':
        case 'C':
            GPolyline(4, x + 1, y + 1, NFC, dd);
            break;
        default:
            break;
        }
        break;
    case 2: /* Figure */
        GPolygon(4, x, y, NFC, NA_INTEGER, dd->gp.col, dd);
        break;
    case 3: /* Inner Region */
        GPolygon(4, x, y, NIC, NA_INTEGER, dd->gp.col, dd);
        break;
    case 4: /* Device border */
        GPolygon(4, x, y, NDC, NA_INTEGER, dd->gp.col, dd);
        break;
    default:
        error("invalid GBox argument\n");
    }
}

#ifdef OLD
void GLPretty(double *xmin, double *xmax, int *n)
{
    /*  Set scale and ticks for logarithmic scales  */
    /*  Note: 1 2 5 10 looks good on logarithmic scales  */

    double u1, u2, v1, v2, p1, p2;

    if ((*xmax) <= 7.8651 * (*xmin))
        GPretty(xmin, xmax, n);
    else
    {
        p1 = ceil(*xmin - 0.0001);
        p2 = floor(*xmax + 0.0001);
        u1 = pow(10.0, p1);
        v1 = pow(10.0, p2);
        *n = p2 - p1;

        if (*n <= 1)
        {
            *n = 2 * (*n);
            u2 = 0.5 * u1;
            if (u2 <= *xmin)
            {
                u2 = u1;
            }
            else
                (*n)++;

            v2 = 5.0 * v1;
            if (v2 >= *xmax)
            {
                v2 = v1;
            }
            else
                (*n)++;

            *xmin = u2;
            *xmax = v2;
        }
        else
        {
            *xmin = u1;
            *xmax = v1;
        }
    }
}
#else

#define LPR_SMALL 2
#define LPR_MEDIUM 3

void GLPretty(double *ul, double *uh, int *n)
{
    /* Generate pretty tick values - logarithmic scale */
    /* This only does a very simple setup.  The real */
    /* work happens when the axis is drawn. */

    int p1, p2;
    p1 = ceil(log10(*ul));
    p2 = floor(log10(*uh));

    if (p2 - p1 <= 0)
    {
        /* Very small range */
        /* Use tickmarks from a linear scale */
        /* Splus uses n = 9 here, but that is dumb */
        GPretty(ul, uh, n);
        *n = -*n;
    }
    else if (p2 - p1 <= LPR_SMALL)
    {
        /* Small range */
        /* Use 1,2,5,10 times 10^k tickmarks */
        *ul = pow(10.0, (double)p1);
        *uh = pow(10.0, (double)p2);
        *n = 3;
    }
    else if (p2 - p1 <= LPR_MEDIUM)
    {
        /* Medium range */
        /* Use 1,5 times 10^k tickmarks */
        *ul = pow(10.0, (double)p1);
        *uh = pow(10.0, (double)p2);
        *n = 2;
    }
    else
    {
        /* Large range */
        /* Use 10^k tickmarks */
        /* But decimate, when there are too many */
        *ul = pow(10.0, (double)p1);
        *uh = pow(10.0, (double)p2);
        *n = 1;
    }
}
#endif

void GPretty(double *lo, double *up, int *ndiv)
{
    /*  Set scale and ticks for linear scales.
     *  Pre:         x1 = lo < up = x2
     *  Post: x1 <= y1 := lo < up =: y2 <= x2;  ndiv >= 1
     */

    double dx, cell, unit, base, U;
    double x1, x2;
    int nd0;
    int ns, nu;
    short i_small;

    if (*lo == R_PosInf || *up == R_PosInf || *lo == R_NegInf || *up == R_NegInf || *ndiv == 0)
        error("infinite axis extents\n");

    x1 = *lo;
    x2 = *up;
    nd0 = *ndiv;
    dx = *up - *lo;
    /* cell := "scale"  here */
    if (dx == 0 && *up == 0)
    { /*  up == lo == 0  */
        cell = i_small = 1;
    }
    else
    {
        cell = fmax2(fabs(*lo), fabs(*up));
        i_small = dx < cell * 10 / (double)INT_MAX;
    }

    /*OLD: cell = FLT_EPSILON+ dx / *ndiv; FLT_EPSILON = 1.192e-07 */
    if (i_small)
        cell *= .25; /* shrink_sml in  pretty(.) */
    else
        cell = dx;
    cell /= *ndiv;

    base = pow(10, floor(log10(cell))); /* base <= cell < 10*base */

    /* unit :=  arg min _u { |u - cell| ;  u = c(1,2,5,10) * base } */
    unit = base;
    if (fabs((U = 2 * base) - cell) < fabs(unit - cell))
        unit = U;
    if (fabs((U = 5 * base) - cell) < fabs(unit - cell))
        unit = U;
    if (fabs((U = 10 * base) - cell) < fabs(unit - cell))
        unit = U;

    ns = floor(*lo / unit);
#ifdef OLD
    while (ns * unit > *lo * (1 - DBL_EPSILON))
        ns--;
#else
    while (ns * unit >= *lo * (1 - DBL_EPSILON))
        ns--;
#endif
    ns++;

    nu = ceil(*up / unit);
#ifdef OLD
    while (nu * unit < *up * (1 + DBL_EPSILON))
        nu++;
#else
    while (nu * unit <= *up * (1 + DBL_EPSILON))
        nu++;
#endif
    nu--;

    *lo = ns * unit;
    *up = nu * unit;
    *ndiv = nu - ns;

    if (*ndiv <= 0)
        printf("Gpretty(%g,%g,%d): cell=%g,  ndiv= %d <=0;\t\t(ns,nu)=(%d,%d);
               dx = % g,
               unit = % 3e,
               ismall = % 1d.\n ",
                        x1,
               x2, nd0, cell, *ndiv, ns, nu, dx, unit, (int)i_small);
}

#define SMALL 0.25
#ifdef OLD
#define RADIUS 0.425
#else
#define RADIUS 0.375
#endif
#define SQRC 0.88622692545275801364 /* sqrt(pi / 4) */
#define DMDC 1.25331413731550025119 /* sqrt(pi / 4) * sqrt(2) */
#define TRC0 1.55512030155621416073 /* sqrt(4 * pi/(3 * sqrt(3))) */
#define TRC1 1.34677368708859836060 /* TRC0 * sqrt(3) / 2 */
#define TRC2 0.77756015077810708036 /* TRC0 / 2 */
#define CMAG 1.1                    /* Circle magnifier */

void GSymbol(double x, double y, int coords, int pch, DevDesc *dd)
{
    /*  Draw one of the R special symbols  */

    double r, xc, yc;
    double xx[4], yy[4];
    char str[2];
    int ltysave;

    if (' ' <= pch && pch <= 255)
    {
        str[0] = pch;
        str[1] = '\0';
        GText(x, y, coords, str, dd->gp.xCharOffset, dd->gp.yCharOffset, 0.0, dd);
    }
    else
    {
        ltysave = dd->gp.lty;
        dd->gp.lty = LTY_SOLID;

        switch (pch)
        {

        case 0: /* S square */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            GRect(x - xc, y - xc, x + xc, y + xc, INCHES, NA_INTEGER, dd->gp.col, dd);
            break;

        case 1: /* S octahedron ( circle) */
            xc = CMAG * RADIUS * GStrWidth("0", INCHES, dd);
            GCircle(x, y, coords, xc, NA_INTEGER, dd->gp.col, dd);
            break;

        case 2: /* S triangle - point up */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            r = TRC0 * xc;
            yc = TRC2 * xc;
            xc = TRC1 * xc;
            xx[0] = x;
            yy[0] = y + r;
            xx[1] = x + xc;
            yy[1] = y - yc;
            xx[2] = x - xc;
            yy[2] = y - yc;
            GPolygon(3, xx, yy, INCHES, NA_INTEGER, dd->gp.col, dd);
            break;

        case 3: /* S plus */
            xc = sqrt(2.0) * RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            GLine(x - xc, y, x + xc, y, INCHES, dd);
            GLine(x, y - xc, x, y + xc, INCHES, dd);
            break;

        case 4: /* S times */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            GLine(x - xc, y - xc, x + xc, y + xc, INCHES, dd);
            GLine(x - xc, y + xc, x + xc, y - xc, INCHES, dd);
            break;

        case 5: /* S diamond */
            xc = sqrt(2.0) * RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            xx[0] = x - xc;
            yy[0] = y;
            xx[1] = x;
            yy[1] = y + xc;
            xx[2] = x + xc;
            yy[2] = y;
            xx[3] = x;
            yy[3] = y - xc;
            GPolygon(4, xx, yy, INCHES, NA_INTEGER, dd->gp.col, dd);
            break;

        case 6: /* S triangle - point down */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            r = TRC0 * xc;
            yc = TRC2 * xc;
            xc = TRC1 * xc;
            xx[0] = x;
            yy[0] = y - r;
            xx[1] = x + xc;
            yy[1] = y + yc;
            xx[2] = x - xc;
            yy[2] = y + yc;
            GPolygon(3, xx, yy, INCHES, NA_INTEGER, dd->gp.col, dd);
            break;

        case 7: /* S square and times superimposed */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            GLine(x - xc, y - xc, x + xc, y + xc, INCHES, dd);
            GLine(x - xc, y + xc, x + xc, y - xc, INCHES, dd);
            GRect(x - xc, y - xc, x + xc, y + xc, INCHES, NA_INTEGER, dd->gp.col, dd);
            break;

        case 8: /* S plus and times superimposed */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            GLine(x - xc, y - xc, x + xc, y + xc, INCHES, dd);
            GLine(x - xc, y + xc, x + xc, y - xc, INCHES, dd);
            xc = sqrt(2.0) * xc;
            GLine(x - xc, y, x + xc, y, INCHES, dd);
            GLine(x, y - xc, x, y + xc, INCHES, dd);
            break;

        case 9: /* S diamond and plus superimposed */
            xc = sqrt(2.0) * RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            xx[0] = x - xc;
            yy[0] = y;
            xx[1] = x;
            yy[1] = y + xc;
            xx[2] = x + xc;
            yy[2] = y;
            xx[3] = x;
            yy[3] = y - xc;
            GPolygon(4, xx, yy, INCHES, NA_INTEGER, dd->gp.col, dd);
            GLine(x - xc, y, x + xc, y, INCHES, dd);
            GLine(x, y - xc, x, y + xc, INCHES, dd);
            break;

        case 10: /* S hexagon (circle) and plus superimposed */
            xc = CMAG * RADIUS * GStrWidth("0", INCHES, dd);
            GCircle(x, y, coords, xc, NA_INTEGER, dd->gp.col, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            GLine(x - xc, y, x + xc, y, INCHES, dd);
            GLine(x, y - xc, x, y + xc, INCHES, dd);
            break;

        case 11: /* S superimposed triangles */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            r = TRC0 * xc;
            yc = TRC2 * xc;
            yc = 0.5 * (yc + r);
            xc = TRC1 * xc;
            xx[0] = x;
            yy[0] = y + r;
            xx[1] = x + xc;
            yy[1] = y - yc;
            xx[2] = x - xc;
            yy[2] = y - yc;
            GPolygon(3, xx, yy, INCHES, NA_INTEGER, dd->gp.col, dd);
            xx[0] = x;
            yy[0] = y - r;
            xx[1] = x + xc;
            yy[1] = y + yc;
            xx[2] = x - xc;
            yy[2] = y + yc;
            GPolygon(3, xx, yy, INCHES, NA_INTEGER, dd->gp.col, dd);
            break;

        case 12: /* S square and plus superimposed */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            GLine(x - xc, y, x + xc, y, INCHES, dd);
            GLine(x, y - xc, x, y + xc, INCHES, dd);
            GRect(x - xc, y - xc, x + xc, y + xc, INCHES, NA_INTEGER, dd->gp.col, dd);
            break;

        case 13: /* S octagon (circle) and times superimposed */
            xc = CMAG * RADIUS * GStrWidth("0", INCHES, dd);
            GCircle(x, y, coords, xc, NA_INTEGER, dd->gp.col, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            GLine(x - xc, y - xc, x + xc, y + xc, INCHES, dd);
            GLine(x - xc, y + xc, x + xc, y - xc, INCHES, dd);
            break;

        case 14: /* S square and point-up triangle superimposed */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            xx[0] = x;
            yy[0] = y + xc;
            xx[1] = x + xc;
            yy[1] = y - xc;
            xx[2] = x - xc;
            yy[2] = y - xc;
            GPolygon(3, xx, yy, INCHES, NA_INTEGER, dd->gp.col, dd);
            GRect(x - xc, y - xc, x + xc, y + xc, INCHES, NA_INTEGER, dd->gp.col, dd);
            break;

        case 15: /* S filled square */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            xx[0] = x - xc;
            yy[0] = y - xc;
            xx[1] = x + xc;
            yy[1] = y - xc;
            xx[2] = x + xc;
            yy[2] = y + xc;
            xx[3] = x - xc;
            yy[3] = y + xc;
            GPolygon(4, xx, yy, INCHES, dd->gp.col, NA_INTEGER, dd);
            break;

        case 16: /* S filled octagon (circle) */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GCircle(x, y, coords, xc, dd->gp.col, dd->gp.col, dd);
            break;

        case 17: /* S filled point-up triangle */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            r = TRC0 * xc;
            yc = TRC2 * xc;
            xc = TRC1 * xc;
            xx[0] = x;
            yy[0] = y + r;
            xx[1] = x + xc;
            yy[1] = y - yc;
            xx[2] = x - xc;
            yy[2] = y - yc;
            GPolygon(3, xx, yy, INCHES, dd->gp.col, NA_INTEGER, dd);
            break;

        case 18:
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            xx[0] = x;
            yy[0] = y - xc;
            xx[1] = x + xc;
            yy[1] = y;
            xx[2] = x;
            yy[2] = y + xc;
            xx[3] = x - xc;
            yy[3] = y;
            GPolygon(4, xx, yy, INCHES, dd->gp.col, NA_INTEGER, dd);
            break;

        case 19: /* R filled circle */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GCircle(x, y, coords, xc, dd->gp.col, dd->gp.col, dd);
            break;

        case 20: /* R Dot */
            xc = SMALL * GStrWidth("0", INCHES, dd);
            GCircle(x, y, coords, xc, dd->gp.col, dd->gp.col, dd);
            break;

        case 21: /* circles */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GCircle(x, y, coords, xc, dd->gp.bg, dd->gp.col, dd);
            break;

        case 22: /* squares */
            xc = RADIUS * SQRC * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            GRect(x - xc, y - xc, x + xc, y + xc, INCHES, dd->gp.bg, dd->gp.col, dd);
            break;

        case 23: /* diamonds */
            xc = RADIUS * DMDC * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            xx[0] = x;
            yy[0] = y - xc;
            xx[1] = x + xc;
            yy[1] = y;
            xx[2] = x;
            yy[2] = y + xc;
            xx[3] = x - xc;
            yy[3] = y;
            GPolygon(4, xx, yy, INCHES, dd->gp.bg, dd->gp.col, dd);
            break;

        case 24: /* triangle (point up) */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            r = TRC0 * xc;
            yc = TRC2 * xc;
            xc = TRC1 * xc;
            xx[0] = x;
            yy[0] = y + r;
            xx[1] = x + xc;
            yy[1] = y - yc;
            xx[2] = x - xc;
            yy[2] = y - yc;
            GPolygon(3, xx, yy, INCHES, dd->gp.bg, dd->gp.col, dd);
            break;

        case 25: /* triangle (point down) */
            xc = RADIUS * GStrWidth("0", INCHES, dd);
            GConvert(&x, &y, coords, INCHES, dd);
            r = TRC0 * xc;
            yc = TRC2 * xc;
            xc = TRC1 * xc;
            xx[0] = x;
            yy[0] = y - r;
            xx[1] = x + xc;
            yy[1] = y + yc;
            xx[2] = x - xc;
            yy[2] = y + yc;
            GPolygon(3, xx, yy, INCHES, dd->gp.bg, dd->gp.col, dd);
            break;
        }
        dd->gp.lty = ltysave;
    }
}

void GMtext(char *str, int side, double line, int outer, double at, int las, DevDesc *dd)
{
    /*  Draw text in plot margins  */

    double a, xadj, yadj;
    int coords;

    a = xadj = yadj = 0.0; /* to keep -Wall happy */
    coords = 0;            /* to keep -Wall happy */
    if (outer)
    {
        switch (side)
        {
        case 1:
            line = line + 1;
            a = 0.0;
            xadj = dd->gp.adj;
            yadj = 0.0;
            coords = OMA1;
            break;
        case 2:
            a = 90.0;
            xadj = dd->gp.adj;
            yadj = 0.0;
            coords = OMA2;
            break;
        case 3:
            a = 0.0;
            xadj = dd->gp.adj;
            yadj = 0.0;
            coords = OMA3;
            break;
        case 4:
            line = line + 1;
            a = 90.0;
            xadj = dd->gp.adj;
            yadj = 0.0;
            coords = OMA4;
            break;
        }
        GText(at, line, coords, str, xadj, yadj, a, dd);
    }
    else
    {
        switch (side)
        {
        case 1:
            if (las == 2)
            {
                at = at - GConvertXUnits(dd->gp.yLineBias, LINES, USER, dd);
                line = line + dd->gp.yLineBias;
                a = 90.0;
                xadj = 1.0;
                yadj = 0.5;
            }
            else
            {
                line = line + 1 - dd->gp.yLineBias;
                a = 0.0;
                xadj = dd->gp.adj;
                yadj = 0.0;
            }
            coords = MAR1;
            break;
        case 2:
            if (las == 1 || las == 2)
            {
                at = at + GConvertYUnits(dd->gp.yLineBias, LINES, USER, dd);
                line = line + dd->gp.yLineBias;
                a = 0.0;
                xadj = 1.0;
                yadj = 0.5;
            }
            else
            {
                line = line + dd->gp.yLineBias;
                a = 90.0;
                xadj = dd->gp.adj;
                yadj = 0.0;
            }
            coords = MAR2;
            break;
        case 3:
            if (las == 2)
            {
                at = at - GConvertXUnits(dd->gp.yLineBias, LINES, USER, dd);
                line = line + dd->gp.yLineBias;
                a = 90.0;
                xadj = 0.0;
                yadj = 0.5;
            }
            else
            {
                line = line + dd->gp.yLineBias;
                a = 0.0;
                xadj = dd->gp.adj;
                yadj = 0.0;
            }
            coords = MAR3;
            break;
        case 4:
            if (las == 1 || las == 2)
            {
                at = at + GConvertYUnits(dd->gp.yLineBias, LINES, USER, dd);
                line = line + dd->gp.yLineBias;
                a = 0.0;
                xadj = 0.0;
                yadj = 0.5;
            }
            else
            {
                line = line + 1 - dd->gp.yLineBias;
                a = 90.0;
                xadj = dd->gp.adj;
                yadj = 0.0;
            }
            coords = MAR4;
            break;
        }
        GText(at, line, coords, str, xadj, yadj, a, dd);
    }
}

/********************************************************/
/* COLOUR CODE						*/
/********************************************************/

void hsv2rgb(double h, double s, double v, double *r, double *g, double *b)
{
    /*  hsv2rgb -- HSV to RGB conversion  */
    /*  Based on HSV_TO_RGB from Foley and Van Dam First Ed. Page 616 */
    /*  See Alvy Ray Smith, Color Gamut Transform Pairs, SIGGRAPH '78 */

    double f, p, q, t;
    double modf();
    int i;

    h = 6 * modf(h, &f);
    i = floor(h);
    f = modf(h, &p);
    p = v * (1 - s);
    q = v * (1 - s * f);
    t = v * (1 - (s * (1 - f)));
    switch (i)
    {
    case 0:
        *r = v;
        *g = t;
        *b = p;
        break;
    case 1:
        *r = q;
        *g = v;
        *b = p;
        break;
    case 2:
        *r = p;
        *g = v;
        *b = t;
        break;
    case 3:
        *r = p;
        *g = q;
        *b = v;
        break;
    case 4:
        *r = t;
        *g = p;
        *b = v;
        break;
    case 5:
        *r = v;
        *g = p;
        *b = q;
        break;
    default:
        error("bad hsv to rgb color conversion\n");
    }
}

/*
 *	Color Specification
 *
 *	Colors are stored internally in integers.  Each integer is
 *	broken into four bytes.  The three least significant bytes
 *	are used to contain levels of red, green and blue.  These
 *	levels are integers in the range [0,255].
 *
 *      Externally, colors are specified either
 *	a) by name, using a large table of color names,
 *	b) by RGB values using a string of the form "#rrggbb"
 *	   where rr, gg and bb are hex integers giving the level
 *	   of red green and blue,
 *	c) as an index into a user setable palette of colors.
 */

/* Default Color Palette */

char *DefaultPalette[] = {"black", "red", "green3", "blue", "cyan", "magenta", "yellow", "white", NULL};

/* The Table of Known Color Names */
/* Adapted from the X11 RGB database */
/* Note: the color "white" is moved */
/* to the top of the database to avoid */
/* its being known as "gray100" */

static int ColorDataBaseSize;

ColorDataBaseEntry ColorDataBase[] = {{"white", "#FFFFFF", 0},
                                      {"aliceblue", "#F0F8FF", 0},
                                      {"antiquewhite", "#FAEBD7", 0},
                                      {"antiquewhite1", "#FFEFDB", 0},
                                      {"antiquewhite2", "#EEDFCC", 0},
                                      {"antiquewhite3", "#CDC0B0", 0},
                                      {"antiquewhite4", "#8B8378", 0},
                                      {"aquamarine", "#7FFFD4", 0},
                                      {"aquamarine1", "#7FFFD4", 0},
                                      {"aquamarine2", "#76EEC6", 0},
                                      {"aquamarine3", "#66CDAA", 0},
                                      {"aquamarine4", "#458B74", 0},
                                      {"azure", "#F0FFFF", 0},
                                      {"azure1", "#F0FFFF", 0},
                                      {"azure2", "#E0EEEE", 0},
                                      {"azure3", "#C1CDCD", 0},
                                      {"azure4", "#838B8B", 0},
                                      {"beige", "#F5F5DC", 0},
                                      {"bisque", "#FFE4C4", 0},
                                      {"bisque1", "#FFE4C4", 0},
                                      {"bisque2", "#EED5B7", 0},
                                      {"bisque3", "#CDB79E", 0},
                                      {"bisque4", "#8B7D6B", 0},
                                      {"black", "#000000", 0},
                                      {"blanchedalmond", "#FFEBCD", 0},
                                      {"blue", "#0000FF", 0},
                                      {"blue1", "#0000FF", 0},
                                      {"blue2", "#0000EE", 0},
                                      {"blue3", "#0000CD", 0},
                                      {"blue4", "#00008B", 0},
                                      {"blueviolet", "#8A2BE2", 0},
                                      {"brown", "#A52A2A", 0},
                                      {"brown1", "#FF4040", 0},
                                      {"brown2", "#EE3B3B", 0},
                                      {"brown3", "#CD3333", 0},
                                      {"brown4", "#8B2323", 0},
                                      {"burlywood", "#DEB887", 0},
                                      {"burlywood1", "#FFD39B", 0},
                                      {"burlywood2", "#EEC591", 0},
                                      {"burlywood3", "#CDAA7D", 0},
                                      {"burlywood4", "#8B7355", 0},
                                      {"cadetblue", "#5F9EA0", 0},
                                      {"cadetblue1", "#98F5FF", 0},
                                      {"cadetblue2", "#8EE5EE", 0},
                                      {"cadetblue3", "#7AC5CD", 0},
                                      {"cadetblue4", "#53868B", 0},
                                      {"chartreuse", "#7FFF00", 0},
                                      {"chartreuse1", "#7FFF00", 0},
                                      {"chartreuse2", "#76EE00", 0},
                                      {"chartreuse3", "#66CD00", 0},
                                      {"chartreuse4", "#458B00", 0},
                                      {"chocolate", "#D2691E", 0},
                                      {"chocolate1", "#FF7F24", 0},
                                      {"chocolate2", "#EE7621", 0},
                                      {"chocolate3", "#CD661D", 0},
                                      {"chocolate4", "#8B4513", 0},
                                      {"coral", "#FF7F50", 0},
                                      {"coral1", "#FF7256", 0},
                                      {"coral2", "#EE6A50", 0},
                                      {"coral3", "#CD5B45", 0},
                                      {"coral4", "#8B3E2F", 0},
                                      {"cornflowerblue", "#6495ED", 0},
                                      {"cornsilk", "#FFF8DC", 0},
                                      {"cornsilk1", "#FFF8DC", 0},
                                      {"cornsilk2", "#EEE8CD", 0},
                                      {"cornsilk3", "#CDC8B1", 0},
                                      {"cornsilk4", "#8B8878", 0},
                                      {"cyan", "#00FFFF", 0},
                                      {"cyan1", "#00FFFF", 0},
                                      {"cyan2", "#00EEEE", 0},
                                      {"cyan3", "#00CDCD", 0},
                                      {"cyan4", "#008B8B", 0},
                                      {"darkblue", "#00008B", 0},
                                      {"darkcyan", "#008B8B", 0},
                                      {"darkgoldenrod", "#B8860B", 0},
                                      {"darkgoldenrod1", "#FFB90F", 0},
                                      {"darkgoldenrod2", "#EEAD0E", 0},
                                      {"darkgoldenrod3", "#CD950C", 0},
                                      {"darkgoldenrod4", "#8B6508", 0},
                                      {"darkgray", "#A9A9A9", 0},
                                      {"darkgreen", "#006400", 0},
                                      {"darkgrey", "#A9A9A9", 0},
                                      {"darkkhaki", "#BDB76B", 0},
                                      {"darkmagenta", "#8B008B", 0},
                                      {"darkolivegreen", "#556B2F", 0},
                                      {"darkolivegreen1", "#CAFF70", 0},
                                      {"darkolivegreen2", "#BCEE68", 0},
                                      {"darkolivegreen3", "#A2CD5A", 0},
                                      {"darkolivegreen4", "#6E8B3D", 0},
                                      {"darkorange", "#FF8C00", 0},
                                      {"darkorange1", "#FF7F00", 0},
                                      {"darkorange2", "#EE7600", 0},
                                      {"darkorange3", "#CD6600", 0},
                                      {"darkorange4", "#8B4500", 0},
                                      {"darkorchid", "#9932CC", 0},
                                      {"darkorchid1", "#BF3EFF", 0},
                                      {"darkorchid2", "#B23AEE", 0},
                                      {"darkorchid3", "#9A32CD", 0},
                                      {"darkorchid4", "#68228B", 0},
                                      {"darkred", "#8B0000", 0},
                                      {"darksalmon", "#E9967A", 0},
                                      {"darkseagreen", "#8FBC8F", 0},
                                      {"darkseagreen1", "#C1FFC1", 0},
                                      {"darkseagreen2", "#B4EEB4", 0},
                                      {"darkseagreen3", "#9BCD9B", 0},
                                      {"darkseagreen4", "#698B69", 0},
                                      {"darkslateblue", "#483D8B", 0},
                                      {"darkslategray", "#2F4F4F", 0},
                                      {"darkslategray1", "#97FFFF", 0},
                                      {"darkslategray2", "#8DEEEE", 0},
                                      {"darkslategray3", "#79CDCD", 0},
                                      {"darkslategray4", "#528B8B", 0},
                                      {"darkslategrey", "#2F4F4F", 0},
                                      {"darkturquoise", "#00CED1", 0},
                                      {"darkviolet", "#9400D3", 0},
                                      {"deeppink", "#FF1493", 0},
                                      {"deeppink1", "#FF1493", 0},
                                      {"deeppink2", "#EE1289", 0},
                                      {"deeppink3", "#CD1076", 0},
                                      {"deeppink4", "#8B0A50", 0},
                                      {"deepskyblue", "#00BFFF", 0},
                                      {"deepskyblue1", "#00BFFF", 0},
                                      {"deepskyblue2", "#00B2EE", 0},
                                      {"deepskyblue3", "#009ACD", 0},
                                      {"deepskyblue4", "#00688B", 0},
                                      {"dimgray", "#696969", 0},
                                      {"dimgrey", "#696969", 0},
                                      {"dodgerblue", "#1E90FF", 0},
                                      {"dodgerblue1", "#1E90FF", 0},
                                      {"dodgerblue2", "#1C86EE", 0},
                                      {"dodgerblue3", "#1874CD", 0},
                                      {"dodgerblue4", "#104E8B", 0},
                                      {"firebrick", "#B22222", 0},
                                      {"firebrick1", "#FF3030", 0},
                                      {"firebrick2", "#EE2C2C", 0},
                                      {"firebrick3", "#CD2626", 0},
                                      {"firebrick4", "#8B1A1A", 0},
                                      {"floralwhite", "#FFFAF0", 0},
                                      {"forestgreen", "#228B22", 0},
                                      {"gainsboro", "#DCDCDC", 0},
                                      {"ghostwhite", "#F8F8FF", 0},
                                      {"gold", "#FFD700", 0},
                                      {"gold1", "#FFD700", 0},
                                      {"gold2", "#EEC900", 0},
                                      {"gold3", "#CDAD00", 0},
                                      {"gold4", "#8B7500", 0},
                                      {"goldenrod", "#DAA520", 0},
                                      {"goldenrod1", "#FFC125", 0},
                                      {"goldenrod2", "#EEB422", 0},
                                      {"goldenrod3", "#CD9B1D", 0},
                                      {"goldenrod4", "#8B6914", 0},
                                      {"gray", "#BEBEBE", 0},
                                      {"gray0", "#000000", 0},
                                      {"gray1", "#030303", 0},
                                      {"gray2", "#050505", 0},
                                      {"gray3", "#080808", 0},
                                      {"gray4", "#0A0A0A", 0},
                                      {"gray5", "#0D0D0D", 0},
                                      {"gray6", "#0F0F0F", 0},
                                      {"gray7", "#121212", 0},
                                      {"gray8", "#141414", 0},
                                      {"gray9", "#171717", 0},
                                      {"gray10", "#1A1A1A", 0},
                                      {"gray11", "#1C1C1C", 0},
                                      {"gray12", "#1F1F1F", 0},
                                      {"gray13", "#212121", 0},
                                      {"gray14", "#242424", 0},
                                      {"gray15", "#262626", 0},
                                      {"gray16", "#292929", 0},
                                      {"gray17", "#2B2B2B", 0},
                                      {"gray18", "#2E2E2E", 0},
                                      {"gray19", "#303030", 0},
                                      {"gray20", "#333333", 0},
                                      {"gray21", "#363636", 0},
                                      {"gray22", "#383838", 0},
                                      {"gray23", "#3B3B3B", 0},
                                      {"gray24", "#3D3D3D", 0},
                                      {"gray25", "#404040", 0},
                                      {"gray26", "#424242", 0},
                                      {"gray27", "#454545", 0},
                                      {"gray28", "#474747", 0},
                                      {"gray29", "#4A4A4A", 0},
                                      {"gray30", "#4D4D4D", 0},
                                      {"gray31", "#4F4F4F", 0},
                                      {"gray32", "#525252", 0},
                                      {"gray33", "#545454", 0},
                                      {"gray34", "#575757", 0},
                                      {"gray35", "#595959", 0},
                                      {"gray36", "#5C5C5C", 0},
                                      {"gray37", "#5E5E5E", 0},
                                      {"gray38", "#616161", 0},
                                      {"gray39", "#636363", 0},
                                      {"gray40", "#666666", 0},
                                      {"gray41", "#696969", 0},
                                      {"gray42", "#6B6B6B", 0},
                                      {"gray43", "#6E6E6E", 0},
                                      {"gray44", "#707070", 0},
                                      {"gray45", "#737373", 0},
                                      {"gray46", "#757575", 0},
                                      {"gray47", "#787878", 0},
                                      {"gray48", "#7A7A7A", 0},
                                      {"gray49", "#7D7D7D", 0},
                                      {"gray50", "#7F7F7F", 0},
                                      {"gray51", "#828282", 0},
                                      {"gray52", "#858585", 0},
                                      {"gray53", "#878787", 0},
                                      {"gray54", "#8A8A8A", 0},
                                      {"gray55", "#8C8C8C", 0},
                                      {"gray56", "#8F8F8F", 0},
                                      {"gray57", "#919191", 0},
                                      {"gray58", "#949494", 0},
                                      {"gray59", "#969696", 0},
                                      {"gray60", "#999999", 0},
                                      {"gray61", "#9C9C9C", 0},
                                      {"gray62", "#9E9E9E", 0},
                                      {"gray63", "#A1A1A1", 0},
                                      {"gray64", "#A3A3A3", 0},
                                      {"gray65", "#A6A6A6", 0},
                                      {"gray66", "#A8A8A8", 0},
                                      {"gray67", "#ABABAB", 0},
                                      {"gray68", "#ADADAD", 0},
                                      {"gray69", "#B0B0B0", 0},
                                      {"gray70", "#B3B3B3", 0},
                                      {"gray71", "#B5B5B5", 0},
                                      {"gray72", "#B8B8B8", 0},
                                      {"gray73", "#BABABA", 0},
                                      {"gray74", "#BDBDBD", 0},
                                      {"gray75", "#BFBFBF", 0},
                                      {"gray76", "#C2C2C2", 0},
                                      {"gray77", "#C4C4C4", 0},
                                      {"gray78", "#C7C7C7", 0},
                                      {"gray79", "#C9C9C9", 0},
                                      {"gray80", "#CCCCCC", 0},
                                      {"gray81", "#CFCFCF", 0},
                                      {"gray82", "#D1D1D1", 0},
                                      {"gray83", "#D4D4D4", 0},
                                      {"gray84", "#D6D6D6", 0},
                                      {"gray85", "#D9D9D9", 0},
                                      {"gray86", "#DBDBDB", 0},
                                      {"gray87", "#DEDEDE", 0},
                                      {"gray88", "#E0E0E0", 0},
                                      {"gray89", "#E3E3E3", 0},
                                      {"gray90", "#E5E5E5", 0},
                                      {"gray91", "#E8E8E8", 0},
                                      {"gray92", "#EBEBEB", 0},
                                      {"gray93", "#EDEDED", 0},
                                      {"gray94", "#F0F0F0", 0},
                                      {"gray95", "#F2F2F2", 0},
                                      {"gray96", "#F5F5F5", 0},
                                      {"gray97", "#F7F7F7", 0},
                                      {"gray98", "#FAFAFA", 0},
                                      {"gray99", "#FCFCFC", 0},
                                      {"gray100", "#FFFFFF", 0},
                                      {"green", "#00FF00", 0},
                                      {"green1", "#00FF00", 0},
                                      {"green2", "#00EE00", 0},
                                      {"green3", "#00CD00", 0},
                                      {"green4", "#008B00", 0},
                                      {"greenyellow", "#ADFF2F", 0},
                                      {"grey", "#BEBEBE", 0},
                                      {"grey0", "#000000", 0},
                                      {"grey1", "#030303", 0},
                                      {"grey2", "#050505", 0},
                                      {"grey3", "#080808", 0},
                                      {"grey4", "#0A0A0A", 0},
                                      {"grey5", "#0D0D0D", 0},
                                      {"grey6", "#0F0F0F", 0},
                                      {"grey7", "#121212", 0},
                                      {"grey8", "#141414", 0},
                                      {"grey9", "#171717", 0},
                                      {"grey10", "#1A1A1A", 0},
                                      {"grey11", "#1C1C1C", 0},
                                      {"grey12", "#1F1F1F", 0},
                                      {"grey13", "#212121", 0},
                                      {"grey14", "#242424", 0},
                                      {"grey15", "#262626", 0},
                                      {"grey16", "#292929", 0},
                                      {"grey17", "#2B2B2B", 0},
                                      {"grey18", "#2E2E2E", 0},
                                      {"grey19", "#303030", 0},
                                      {"grey20", "#333333", 0},
                                      {"grey21", "#363636", 0},
                                      {"grey22", "#383838", 0},
                                      {"grey23", "#3B3B3B", 0},
                                      {"grey24", "#3D3D3D", 0},
                                      {"grey25", "#404040", 0},
                                      {"grey26", "#424242", 0},
                                      {"grey27", "#454545", 0},
                                      {"grey28", "#474747", 0},
                                      {"grey29", "#4A4A4A", 0},
                                      {"grey30", "#4D4D4D", 0},
                                      {"grey31", "#4F4F4F", 0},
                                      {"grey32", "#525252", 0},
                                      {"grey33", "#545454", 0},
                                      {"grey34", "#575757", 0},
                                      {"grey35", "#595959", 0},
                                      {"grey36", "#5C5C5C", 0},
                                      {"grey37", "#5E5E5E", 0},
                                      {"grey38", "#616161", 0},
                                      {"grey39", "#636363", 0},
                                      {"grey40", "#666666", 0},
                                      {"grey41", "#696969", 0},
                                      {"grey42", "#6B6B6B", 0},
                                      {"grey43", "#6E6E6E", 0},
                                      {"grey44", "#707070", 0},
                                      {"grey45", "#737373", 0},
                                      {"grey46", "#757575", 0},
                                      {"grey47", "#787878", 0},
                                      {"grey48", "#7A7A7A", 0},
                                      {"grey49", "#7D7D7D", 0},
                                      {"grey50", "#7F7F7F", 0},
                                      {"grey51", "#828282", 0},
                                      {"grey52", "#858585", 0},
                                      {"grey53", "#878787", 0},
                                      {"grey54", "#8A8A8A", 0},
                                      {"grey55", "#8C8C8C", 0},
                                      {"grey56", "#8F8F8F", 0},
                                      {"grey57", "#919191", 0},
                                      {"grey58", "#949494", 0},
                                      {"grey59", "#969696", 0},
                                      {"grey60", "#999999", 0},
                                      {"grey61", "#9C9C9C", 0},
                                      {"grey62", "#9E9E9E", 0},
                                      {"grey63", "#A1A1A1", 0},
                                      {"grey64", "#A3A3A3", 0},
                                      {"grey65", "#A6A6A6", 0},
                                      {"grey66", "#A8A8A8", 0},
                                      {"grey67", "#ABABAB", 0},
                                      {"grey68", "#ADADAD", 0},
                                      {"grey69", "#B0B0B0", 0},
                                      {"grey70", "#B3B3B3", 0},
                                      {"grey71", "#B5B5B5", 0},
                                      {"grey72", "#B8B8B8", 0},
                                      {"grey73", "#BABABA", 0},
                                      {"grey74", "#BDBDBD", 0},
                                      {"grey75", "#BFBFBF", 0},
                                      {"grey76", "#C2C2C2", 0},
                                      {"grey77", "#C4C4C4", 0},
                                      {"grey78", "#C7C7C7", 0},
                                      {"grey79", "#C9C9C9", 0},
                                      {"grey80", "#CCCCCC", 0},
                                      {"grey81", "#CFCFCF", 0},
                                      {"grey82", "#D1D1D1", 0},
                                      {"grey83", "#D4D4D4", 0},
                                      {"grey84", "#D6D6D6", 0},
                                      {"grey85", "#D9D9D9", 0},
                                      {"grey86", "#DBDBDB", 0},
                                      {"grey87", "#DEDEDE", 0},
                                      {"grey88", "#E0E0E0", 0},
                                      {"grey89", "#E3E3E3", 0},
                                      {"grey90", "#E5E5E5", 0},
                                      {"grey91", "#E8E8E8", 0},
                                      {"grey92", "#EBEBEB", 0},
                                      {"grey93", "#EDEDED", 0},
                                      {"grey94", "#F0F0F0", 0},
                                      {"grey95", "#F2F2F2", 0},
                                      {"grey96", "#F5F5F5", 0},
                                      {"grey97", "#F7F7F7", 0},
                                      {"grey98", "#FAFAFA", 0},
                                      {"grey99", "#FCFCFC", 0},
                                      {"grey100", "#FFFFFF", 0},
                                      {"honeydew", "#F0FFF0", 0},
                                      {"honeydew1", "#F0FFF0", 0},
                                      {"honeydew2", "#E0EEE0", 0},
                                      {"honeydew3", "#C1CDC1", 0},
                                      {"honeydew4", "#838B83", 0},
                                      {"hotpink", "#FF69B4", 0},
                                      {"hotpink1", "#FF6EB4", 0},
                                      {"hotpink2", "#EE6AA7", 0},
                                      {"hotpink3", "#CD6090", 0},
                                      {"hotpink4", "#8B3A62", 0},
                                      {"indianred", "#CD5C5C", 0},
                                      {"indianred1", "#FF6A6A", 0},
                                      {"indianred2", "#EE6363", 0},
                                      {"indianred3", "#CD5555", 0},
                                      {"indianred4", "#8B3A3A", 0},
                                      {"ivory", "#FFFFF0", 0},
                                      {"ivory1", "#FFFFF0", 0},
                                      {"ivory2", "#EEEEE0", 0},
                                      {"ivory3", "#CDCDC1", 0},
                                      {"ivory4", "#8B8B83", 0},
                                      {"khaki", "#F0E68C", 0},
                                      {"khaki1", "#FFF68F", 0},
                                      {"khaki2", "#EEE685", 0},
                                      {"khaki3", "#CDC673", 0},
                                      {"khaki4", "#8B864E", 0},
                                      {"lavender", "#E6E6FA", 0},
                                      {"lavenderblush", "#FFF0F5", 0},
                                      {"lavenderblush1", "#FFF0F5", 0},
                                      {"lavenderblush2", "#EEE0E5", 0},
                                      {"lavenderblush3", "#CDC1C5", 0},
                                      {"lavenderblush4", "#8B8386", 0},
                                      {"lawngreen", "#7CFC00", 0},
                                      {"lemonchiffon", "#FFFACD", 0},
                                      {"lemonchiffon1", "#FFFACD", 0},
                                      {"lemonchiffon2", "#EEE9BF", 0},
                                      {"lemonchiffon3", "#CDC9A5", 0},
                                      {"lemonchiffon4", "#8B8970", 0},
                                      {"lightblue", "#ADD8E6", 0},
                                      {"lightblue1", "#BFEFFF", 0},
                                      {"lightblue2", "#B2DFEE", 0},
                                      {"lightblue3", "#9AC0CD", 0},
                                      {"lightblue4", "#68838B", 0},
                                      {"lightcoral", "#F08080", 0},
                                      {"lightcyan", "#E0FFFF", 0},
                                      {"lightcyan1", "#E0FFFF", 0},
                                      {"lightcyan2", "#D1EEEE", 0},
                                      {"lightcyan3", "#B4CDCD", 0},
                                      {"lightcyan4", "#7A8B8B", 0},
                                      {"lightgoldenrod", "#EEDD82", 0},
                                      {"lightgoldenrod1", "#FFEC8B", 0},
                                      {"lightgoldenrod2", "#EEDC82", 0},
                                      {"lightgoldenrod3", "#CDBE70", 0},
                                      {"lightgoldenrod4", "#8B814C", 0},
                                      {"lightgoldenrodyellow", "#FAFAD2", 0},
                                      {"lightgray", "#D3D3D3", 0},
                                      {"lightgreen", "#90EE90", 0},
                                      {"lightgrey", "#D3D3D3", 0},
                                      {"lightpink", "#FFB6C1", 0},
                                      {"lightpink1", "#FFAEB9", 0},
                                      {"lightpink2", "#EEA2AD", 0},
                                      {"lightpink3", "#CD8C95", 0},
                                      {"lightpink4", "#8B5F65", 0},
                                      {"lightsalmon", "#FFA07A", 0},
                                      {"lightsalmon1", "#FFA07A", 0},
                                      {"lightsalmon2", "#EE9572", 0},
                                      {"lightsalmon3", "#CD8162", 0},
                                      {"lightsalmon4", "#8B5742", 0},
                                      {"lightseagreen", "#20B2AA", 0},
                                      {"lightskyblue", "#87CEFA", 0},
                                      {"lightskyblue1", "#B0E2FF", 0},
                                      {"lightskyblue2", "#A4D3EE", 0},
                                      {"lightskyblue3", "#8DB6CD", 0},
                                      {"lightskyblue4", "#607B8B", 0},
                                      {"lightslateblue", "#8470FF", 0},
                                      {"lightslategray", "#778899", 0},
                                      {"lightslategrey", "#778899", 0},
                                      {"lightsteelblue", "#B0C4DE", 0},
                                      {"lightsteelblue1", "#CAE1FF", 0},
                                      {"lightsteelblue2", "#BCD2EE", 0},
                                      {"lightsteelblue3", "#A2B5CD", 0},
                                      {"lightsteelblue4", "#6E7B8B", 0},
                                      {"lightyellow", "#FFFFE0", 0},
                                      {"lightyellow1", "#FFFFE0", 0},
                                      {"lightyellow2", "#EEEED1", 0},
                                      {"lightyellow3", "#CDCDB4", 0},
                                      {"lightyellow4", "#8B8B7A", 0},
                                      {"limegreen", "#32CD32", 0},
                                      {"linen", "#FAF0E6", 0},
                                      {"magenta", "#FF00FF", 0},
                                      {"magenta1", "#FF00FF", 0},
                                      {"magenta2", "#EE00EE", 0},
                                      {"magenta3", "#CD00CD", 0},
                                      {"magenta4", "#8B008B", 0},
                                      {"maroon", "#B03060", 0},
                                      {"maroon1", "#FF34B3", 0},
                                      {"maroon2", "#EE30A7", 0},
                                      {"maroon3", "#CD2990", 0},
                                      {"maroon4", "#8B1C62", 0},
                                      {"mediumaquamarine", "#66CDAA", 0},
                                      {"mediumblue", "#0000CD", 0},
                                      {"mediumorchid", "#BA55D3", 0},
                                      {"mediumorchid1", "#E066FF", 0},
                                      {"mediumorchid2", "#D15FEE", 0},
                                      {"mediumorchid3", "#B452CD", 0},
                                      {"mediumorchid4", "#7A378B", 0},
                                      {"mediumpurple", "#9370DB", 0},
                                      {"mediumpurple1", "#AB82FF", 0},
                                      {"mediumpurple2", "#9F79EE", 0},
                                      {"mediumpurple3", "#8968CD", 0},
                                      {"mediumpurple4", "#5D478B", 0},
                                      {"mediumseagreen", "#3CB371", 0},
                                      {"mediumslateblue", "#7B68EE", 0},
                                      {"mediumspringgreen", "#00FA9A", 0},
                                      {"mediumturquoise", "#48D1CC", 0},
                                      {"mediumvioletred", "#C71585", 0},
                                      {"midnightblue", "#191970", 0},
                                      {"mintcream", "#F5FFFA", 0},
                                      {"mistyrose", "#FFE4E1", 0},
                                      {"mistyrose1", "#FFE4E1", 0},
                                      {"mistyrose2", "#EED5D2", 0},
                                      {"mistyrose3", "#CDB7B5", 0},
                                      {"mistyrose4", "#8B7D7B", 0},
                                      {"moccasin", "#FFE4B5", 0},
                                      {"navajowhite", "#FFDEAD", 0},
                                      {"navajowhite1", "#FFDEAD", 0},
                                      {"navajowhite2", "#EECFA1", 0},
                                      {"navajowhite3", "#CDB38B", 0},
                                      {"navajowhite4", "#8B795E", 0},
                                      {"navy", "#000080", 0},
                                      {"navyblue", "#000080", 0},
                                      {"oldlace", "#FDF5E6", 0},
                                      {"olivedrab", "#6B8E23", 0},
                                      {"olivedrab1", "#C0FF3E", 0},
                                      {"olivedrab2", "#B3EE3A", 0},
                                      {"olivedrab3", "#9ACD32", 0},
                                      {"olivedrab4", "#698B22", 0},
                                      {"orange", "#FFA500", 0},
                                      {"orange1", "#FFA500", 0},
                                      {"orange2", "#EE9A00", 0},
                                      {"orange3", "#CD8500", 0},
                                      {"orange4", "#8B5A00", 0},
                                      {"orangered", "#FF4500", 0},
                                      {"orangered1", "#FF4500", 0},
                                      {"orangered2", "#EE4000", 0},
                                      {"orangered3", "#CD3700", 0},
                                      {"orangered4", "#8B2500", 0},
                                      {"orchid", "#DA70D6", 0},
                                      {"orchid1", "#FF83FA", 0},
                                      {"orchid2", "#EE7AE9", 0},
                                      {"orchid3", "#CD69C9", 0},
                                      {"orchid4", "#8B4789", 0},
                                      {"palegoldenrod", "#EEE8AA", 0},
                                      {"palegreen", "#98FB98", 0},
                                      {"palegreen1", "#9AFF9A", 0},
                                      {"palegreen2", "#90EE90", 0},
                                      {"palegreen3", "#7CCD7C", 0},
                                      {"palegreen4", "#548B54", 0},
                                      {"paleturquoise", "#AFEEEE", 0},
                                      {"paleturquoise1", "#BBFFFF", 0},
                                      {"paleturquoise2", "#AEEEEE", 0},
                                      {"paleturquoise3", "#96CDCD", 0},
                                      {"paleturquoise4", "#668B8B", 0},
                                      {"palevioletred", "#DB7093", 0},
                                      {"palevioletred1", "#FF82AB", 0},
                                      {"palevioletred2", "#EE799F", 0},
                                      {"palevioletred3", "#CD6889", 0},
                                      {"palevioletred4", "#8B475D", 0},
                                      {"papayawhip", "#FFEFD5", 0},
                                      {"peachpuff", "#FFDAB9", 0},
                                      {"peachpuff1", "#FFDAB9", 0},
                                      {"peachpuff2", "#EECBAD", 0},
                                      {"peachpuff3", "#CDAF95", 0},
                                      {"peachpuff4", "#8B7765", 0},
                                      {"peru", "#CD853F", 0},
                                      {"pink", "#FFC0CB", 0},
                                      {"pink1", "#FFB5C5", 0},
                                      {"pink2", "#EEA9B8", 0},
                                      {"pink3", "#CD919E", 0},
                                      {"pink4", "#8B636C", 0},
                                      {"plum", "#DDA0DD", 0},
                                      {"plum1", "#FFBBFF", 0},
                                      {"plum2", "#EEAEEE", 0},
                                      {"plum3", "#CD96CD", 0},
                                      {"plum4", "#8B668B", 0},
                                      {"powderblue", "#B0E0E6", 0},
                                      {"purple", "#A020F0", 0},
                                      {"purple1", "#9B30FF", 0},
                                      {"purple2", "#912CEE", 0},
                                      {"purple3", "#7D26CD", 0},
                                      {"purple4", "#551A8B", 0},
                                      {"red", "#FF0000", 0},
                                      {"red1", "#FF0000", 0},
                                      {"red2", "#EE0000", 0},
                                      {"red3", "#CD0000", 0},
                                      {"red4", "#8B0000", 0},
                                      {"rosybrown", "#BC8F8F", 0},
                                      {"rosybrown1", "#FFC1C1", 0},
                                      {"rosybrown2", "#EEB4B4", 0},
                                      {"rosybrown3", "#CD9B9B", 0},
                                      {"rosybrown4", "#8B6969", 0},
                                      {"royalblue", "#4169E1", 0},
                                      {"royalblue1", "#4876FF", 0},
                                      {"royalblue2", "#436EEE", 0},
                                      {"royalblue3", "#3A5FCD", 0},
                                      {"royalblue4", "#27408B", 0},
                                      {"saddlebrown", "#8B4513", 0},
                                      {"salmon", "#FA8072", 0},
                                      {"salmon1", "#FF8C69", 0},
                                      {"salmon2", "#EE8262", 0},
                                      {"salmon3", "#CD7054", 0},
                                      {"salmon4", "#8B4C39", 0},
                                      {"sandybrown", "#F4A460", 0},
                                      {"seagreen", "#2E8B57", 0},
                                      {"seagreen1", "#54FF9F", 0},
                                      {"seagreen2", "#4EEE94", 0},
                                      {"seagreen3", "#43CD80", 0},
                                      {"seagreen4", "#2E8B57", 0},
                                      {"seashell", "#FFF5EE", 0},
                                      {"seashell1", "#FFF5EE", 0},
                                      {"seashell2", "#EEE5DE", 0},
                                      {"seashell3", "#CDC5BF", 0},
                                      {"seashell4", "#8B8682", 0},
                                      {"sienna", "#A0522D", 0},
                                      {"sienna1", "#FF8247", 0},
                                      {"sienna2", "#EE7942", 0},
                                      {"sienna3", "#CD6839", 0},
                                      {"sienna4", "#8B4726", 0},
                                      {"skyblue", "#87CEEB", 0},
                                      {"skyblue1", "#87CEFF", 0},
                                      {"skyblue2", "#7EC0EE", 0},
                                      {"skyblue3", "#6CA6CD", 0},
                                      {"skyblue4", "#4A708B", 0},
                                      {"slateblue", "#6A5ACD", 0},
                                      {"slateblue1", "#836FFF", 0},
                                      {"slateblue2", "#7A67EE", 0},
                                      {"slateblue3", "#6959CD", 0},
                                      {"slateblue4", "#473C8B", 0},
                                      {"slategray", "#708090", 0},
                                      {"slategray1", "#C6E2FF", 0},
                                      {"slategray2", "#B9D3EE", 0},
                                      {"slategray3", "#9FB6CD", 0},
                                      {"slategray4", "#6C7B8B", 0},
                                      {"slategrey", "#708090", 0},
                                      {"snow", "#FFFAFA", 0},
                                      {"snow1", "#FFFAFA", 0},
                                      {"snow2", "#EEE9E9", 0},
                                      {"snow3", "#CDC9C9", 0},
                                      {"snow4", "#8B8989", 0},
                                      {"springgreen", "#00FF7F", 0},
                                      {"springgreen1", "#00FF7F", 0},
                                      {"springgreen2", "#00EE76", 0},
                                      {"springgreen3", "#00CD66", 0},
                                      {"springgreen4", "#008B45", 0},
                                      {"steelblue", "#4682B4", 0},
                                      {"steelblue1", "#63B8FF", 0},
                                      {"steelblue2", "#5CACEE", 0},
                                      {"steelblue3", "#4F94CD", 0},
                                      {"steelblue4", "#36648B", 0},
                                      {"tan", "#D2B48C", 0},
                                      {"tan1", "#FFA54F", 0},
                                      {"tan2", "#EE9A49", 0},
                                      {"tan3", "#CD853F", 0},
                                      {"tan4", "#8B5A2B", 0},
                                      {"thistle", "#D8BFD8", 0},
                                      {"thistle1", "#FFE1FF", 0},
                                      {"thistle2", "#EED2EE", 0},
                                      {"thistle3", "#CDB5CD", 0},
                                      {"thistle4", "#8B7B8B", 0},
                                      {"tomato", "#FF6347", 0},
                                      {"tomato1", "#FF6347", 0},
                                      {"tomato2", "#EE5C42", 0},
                                      {"tomato3", "#CD4F39", 0},
                                      {"tomato4", "#8B3626", 0},
                                      {"turquoise", "#40E0D0", 0},
                                      {"turquoise1", "#00F5FF", 0},
                                      {"turquoise2", "#00E5EE", 0},
                                      {"turquoise3", "#00C5CD", 0},
                                      {"turquoise4", "#00868B", 0},
                                      {"violet", "#EE82EE", 0},
                                      {"violetred", "#D02090", 0},
                                      {"violetred1", "#FF3E96", 0},
                                      {"violetred2", "#EE3A8C", 0},
                                      {"violetred3", "#CD3278", 0},
                                      {"violetred4", "#8B2252", 0},
                                      {"wheat", "#F5DEB3", 0},
                                      {"wheat1", "#FFE7BA", 0},
                                      {"wheat2", "#EED8AE", 0},
                                      {"wheat3", "#CDBA96", 0},
                                      {"wheat4", "#8B7E66", 0},
                                      {"whitesmoke", "#F5F5F5", 0},
                                      {"yellow", "#FFFF00", 0},
                                      {"yellow1", "#FFFF00", 0},
                                      {"yellow2", "#EEEE00", 0},
                                      {"yellow3", "#CDCD00", 0},
                                      {"yellow4", "#8B8B00", 0},
                                      {"yellowgreen", "#9ACD32", 0},
                                      {NULL, NULL, 0}};

int ColorTableSize;
unsigned int ColorTable[COLOR_TABLE_SIZE];

static unsigned hexdigit(int digit)
{
    /* Hex Digit to Integer Conversion */

    if ('0' <= digit && digit <= '9')
        return digit - '0';
    else if ('A' <= digit && digit <= 'F')
        return 10 + digit - 'A';
    else if ('a' <= digit && digit <= 'f')
        return 10 + digit - 'a';
    else
        error("invalid hex digit in color\n");
    return digit - '0'; /* never occurs but avoid compiler warnings */
}

static char HexDigits[] = "0123456789ABCDEF";

#ifdef UNUSED
static unsigned digithex(int digit)
{
    /* Integer to Hex Digit */

    return HexDigits[abs(digit) % 16];
}
#endif

static int StrMatch(char *s, char *t)
{
    /* String Comparison Ignoring Case and Squeezing Out Blanks */

    for (;;)
    {
        if (*s == '\0' && *t == '\0')
        {
            return 1;
        }
        if (*s == ' ')
        {
            s++;
            continue;
        }
        if (*t == ' ')
        {
            t++;
            continue;
        }
        if (tolower(*s++) != tolower(*t++))
            return 0;
    }
}

unsigned int rgb2col(char *rgb)
{
    /* #RRGGBB String to Internal Color Code */

    unsigned int r, g, b;
    if (rgb[0] != '#' || strlen(rgb) != 7)
        error("invalid RGB specification\n");
    r = 16 * hexdigit(rgb[1]) + hexdigit(rgb[2]);
    g = 16 * hexdigit(rgb[3]) + hexdigit(rgb[4]);
    b = 16 * hexdigit(rgb[5]) + hexdigit(rgb[6]);
    return R_RGB(r, g, b);
}

unsigned int name2col(char *nm)
{
    /* External Color Name to Internal Color Code */

    int i;
    for (i = 0; ColorDataBase[i].name; i++)
    {
        if (StrMatch(ColorDataBase[i].name, nm))
            return ColorDataBase[i].code;
    }
    error("invalid color name\n");
    return 0; /* never occurs but avoid compiler warnings */
}

unsigned int number2col(char *nm, DevDesc *dd)
{
    /* Index (as string) to Internal Color Code */

    int index;
    char *ptr;
    index = strtod(nm, &ptr);
    if (*ptr)
        error("invalid color specification\n");
    if (index == 0)
        return dd->dp.bg;
    else
        return ColorTable[(index - 1) % ColorTableSize];
}

static char ColBuf[8];

char *RGB2rgb(unsigned int r, unsigned int g, unsigned int b)
{
    ColBuf[0] = '#';
    ColBuf[1] = HexDigits[(r >> 4) & 15];
    ColBuf[2] = HexDigits[r & 15];
    ColBuf[3] = HexDigits[(g >> 4) & 15];
    ColBuf[4] = HexDigits[g & 15];
    ColBuf[5] = HexDigits[(b >> 4) & 15];
    ColBuf[6] = HexDigits[b & 15];
    ColBuf[7] = '\0';
    return &ColBuf[0];
}

char *col2name(unsigned int col)
{
    /* Internal to External Color Representation */
    /* Search the color name database first */
    /* If this fails, create an #RRGGBB string */

    int i;
    for (i = 0; ColorDataBase[i].name; i++)
    {
        if (col == ColorDataBase[i].code)
            return ColorDataBase[i].name;
    }
    ColBuf[0] = '#';
    ColBuf[1] = HexDigits[(col >> 4) & 15];
    ColBuf[2] = HexDigits[(col)&15];
    ColBuf[3] = HexDigits[(col >> 12) & 15];
    ColBuf[4] = HexDigits[(col >> 8) & 15];
    ColBuf[5] = HexDigits[(col >> 20) & 15];
    ColBuf[6] = HexDigits[(col >> 16) & 15];
    ColBuf[7] = '\0';
    return &ColBuf[0];
}

unsigned str2col(char *s, DevDesc *dd)
{
    /* NOTE that this is called with dd == NULL by */
    /* the initialisation code in which case, str2col */
    /* assumes that `s' is a name */

    if (dd)
    {
        if (s[0] == '#')
            return rgb2col(s);
        else if (isdigit(s[0]))
            return number2col(s, dd);
        else
            return name2col(s);
    }
    else
        return name2col(s);
}

unsigned RGBpar(SEXP x, int i, DevDesc *dd)
{
    /* Convert a sexp element to an R  color desc */
    /* We Assume that Checks Have Been Done */

    int index;
    if (isString(x))
    {
        return str2col(CHAR(STRING(x)[i]), dd);
    }
    else if (isInteger(x) || isLogical(x))
    {
        if (INTEGER(x)[i] == NA_INTEGER)
            return NA_INTEGER;
        index = INTEGER(x)[i] - 1;
        if (index < 0)
            return dd->dp.bg;
        else
            return ColorTable[abs(index) % ColorTableSize];
    }
    else if (isReal(x))
    {
        if (!FINITE(REAL(x)[i]))
            return NA_INTEGER;
        index = REAL(x)[i] - 1;
        if (index < 0)
            return dd->dp.bg;
        else
            return ColorTable[abs(index) % ColorTableSize];
    }
    return 0; /* should not occur */
}

void InitColors()
{
    /* Initialize the Color Databases */

    int i;

    /* Initialize the Color Database */
    for (i = 0; ColorDataBase[i].name; i++)
        ColorDataBase[i].code = rgb2col(ColorDataBase[i].rgb);
    ColorDataBaseSize = i;

    /* Install Default Palette */
    for (i = 0; DefaultPalette[i]; i++)
        ColorTable[i] = str2col(DefaultPalette[i], NULL);
    ColorTableSize = i;
}

/********************************************************/
/* LINE TEXTURE CODE					*/
/********************************************************/

/*
 *	LINE TEXTURE SPECIFICATION
 *
 *	Linetypes are stored internally in integers.  An integer
 *	is interpreted as containing a sequence of 8 4-bit integers
 *	which give the lengths of up to 8 on-off line segments.
 *	The lengths are typically interpreted as pixels on a screen
 *	and as "points" in postscript.
 */

typedef struct
{
    char *name;
    unsigned int pattern;
} LineTYPE;

static LineTYPE linetype[] = {
    {"solid", LTY_SOLID}, {"dashed", LTY_DASHED}, {"dotted", LTY_DOTTED}, {"dotdash", LTY_DOTDASH}, {NULL, 0},
};

static int nlinetype = (sizeof(linetype) / sizeof(LineTYPE) - 1);

unsigned int LTYpar(SEXP value, int index)
{
    char *p;
    int i, code, shift, digit;

    if (isString(value))
    {
        for (i = 0; linetype[i].name; i++)
        {
            if (!strcmp(CHAR(STRING(value)[index]), linetype[i].name))
                return linetype[i].pattern;
        }
        code = 0;
        shift = 0;
        for (p = CHAR(STRING(value)[index]); *p; p++)
        {
            digit = hexdigit(*p);
            code = code | (digit << shift);
            shift = shift + 4;
        }
        return code;
    }
    else if (isInteger(value))
    {
        code = INTEGER(value)[index];
        if (code == NA_INTEGER || code <= 0)
            return NA_INTEGER;
        code = (code - 1) % nlinetype;
        return linetype[code].pattern;
    }
    else if (isReal(value))
    {
        code = REAL(value)[index];
        if (!FINITE(code) || code <= 0)
            return NA_INTEGER;
        code = (code - 1) % nlinetype;
        return linetype[code].pattern;
    }
    else
        error("invalid line type\n");
    /*NOTREACHED*/
    return 0; /* never occurs but avoid compiler warnings */
}

SEXP LTYget(unsigned int lty)
{
    SEXP ans;
    int i, ndash;
    char dash[8];
    unsigned int l;

    for (i = 0; linetype[i].name; i++)
    {
        if (linetype[i].pattern == lty)
            return mkString(linetype[i].name);
    }

    l = lty;
    ndash = 0;
    for (i = 0; i < 8 && l & 15; i++)
    {
        dash[ndash++] = l & 15;
        l = l >> 4;
    }
    PROTECT(ans = allocVector(STRSXP, 1));
    STRING(ans)[0] = allocString(ndash);
    for (i = 0; i < ndash; i++)
    {
        CHAR(STRING(ans)[0])[i] = dash[i] + '0';
    }
    UNPROTECT(1);
    return ans;
}

/***********************************************************************
 *
 *  DEVICE FUNCTIONS
 *
 *    R allows there to be (up to 64) multiple devices in
 *    existence at the same time.  Only one device is the
 *    active device and all drawing occurs in this device
 *
 *    Each device has its own set of graphics parameters
 *    so that switching between devices, switches between
 *    their graphical contexts (e.g., if you set the line
 *    width on one device then switch to another device,
 *    don't expect to be using the line width you just set !)
 *
 *    Each device has additional device-specific graphics
 *    parameters which the device driver (i.e., NOT this
 *    generic graphics code) is wholly responsible for
 *    maintaining (including creating and destroying special
 *    resources such as X11 windows).
 *
 *    Each device has a display list which records every
 *    graphical operation since the last dd->dp.newPage;
 *    this is used to redraw the output on the device
 *    when it is resized and to copy output from one device
 *    to another (this can be disabled, which is the default
 *    for postscript).
 *
 *    NOTE: that graphical operations should only be
 *    recorded in the displayList if they are "guaranteed"
 *    to succeed (to avoid heaps of error messages on a
 *    redraw) which means that the recording should be the
 *    last thing done in a graphical operation (see do_*
 *    in plot.c).
 *
 ***********************************************************************/

void DevNull()
{
}

static int R_CurrentDevice = 0;
int R_NumDevices = 1;
DevDesc *R_Devices[R_MaxDevices];

/* a dummy description to point to when there are no active devices */

DevDesc nullDevice;

#ifdef Unix
int PSDeviceDriver(DevDesc *, SEXP, int, SEXP, int);
int X11DeviceDriver(DevDesc *, SEXP, int, SEXP, int);
#ifdef SOON
int XFigDeviceDriver(char **, int, double *, int);
#endif
#endif

#ifdef Win32
int PSDeviceDriver(DevDesc *, SEXP, int, SEXP, int);
int WinDeviceDriver(char **, int, double *, int);
#endif

#ifdef Macintosh
int MacDeviceDriver(char **, int, double *, int);
#endif

void devError()
{
    error("No graphics device is active\n");
}

int NoDevices()
{
    return (R_NumDevices == 1 || R_CurrentDevice == 0);
}

int NumDevices()
{
    return R_NumDevices;
}

DevDesc *CurrentDevice()
{
    return R_Devices[R_CurrentDevice];
}

DevDesc *GetDevice(int i)
{
    return R_Devices[i];
}

void InitGraphics(void)
{
    int i;
    SEXP s, t;

    /* init R_Devices */
    R_Devices[0] = &nullDevice;
    for (i = 1; i < R_MaxDevices; i++)
        R_Devices[i] = NULL;

    /* init .Device and .Devices */
    PROTECT(s = mkString("null device"));
    gsetVar(install(".Device"), s, R_NilValue);
    PROTECT(t = mkString("null device"));
    gsetVar(install(".Devices"), CONS(t, R_NilValue), R_NilValue);
    UNPROTECT(2);
}

SEXP elt(SEXP list, int i)
{
    /* Get the i-th element of a list */

    int j;
    SEXP result = list;

    if ((i < 0) || (i > length(list)))
        return R_NilValue;
    else
        for (j = 0; j < i; j++)
            result = CDR(result);

    return CAR(result);
}

static SEXP getSymbolValue(char *symbolName)
{
    SEXP s, t;
    t = findVar(s = install(symbolName), R_NilValue);
    return t;
}

int curDevice()
{
    return R_CurrentDevice;
}

int nextDevice(int from)
{
    if (R_NumDevices == 1)
        return 0;
    else
    {
        int i = from;
        int nextDev = 0;
        while ((i < (R_MaxDevices - 1)) && (nextDev == 0))
            if (R_Devices[++i] != NULL)
                nextDev = i;
        if (nextDev == 0)
        {
            /* start again from 1 */
            i = 0;
            while (nextDev == 0)
                if (R_Devices[++i] != NULL)
                    nextDev = i;
        }
        return nextDev;
    }
}

int prevDevice(int from)
{
    if (R_NumDevices == 1)
        return 0;
    else
    {
        int i = from;
        int prevDev = 0;
        while ((i > 1) && (prevDev == 0))
            if (R_Devices[--i] != NULL)
                prevDev = i;
        if (prevDev == 0)
        {
            /* start again from R_MaxDevices */
            i = R_MaxDevices;
            while (prevDev == 0)
                if (R_Devices[--i] != NULL)
                    prevDev = i;
        }
        return prevDev;
    }
}

void addDevice(DevDesc *dd)
{
    int i, append;
    SEXP s, t;
    DevDesc *oldd;
    PROTECT(s = getSymbolValue(".Devices"));

    if (!NoDevices())
    {
        oldd = CurrentDevice();
        oldd->dp.deactivate(oldd);
    }

    /* find empty slot for new descriptor */
    i = 1;
    if (CDR(s) == R_NilValue)
        append = 1;
    else
    {
        s = CDR(s);
        append = 0;
    }
    while (R_Devices[i] != NULL)
    {
        i++;
        if (CDR(s) == R_NilValue)
            append = 1;
        else
            s = CDR(s);
    }
    R_CurrentDevice = i;
    R_NumDevices += 1;
    R_Devices[i] = dd;

    dd->dp.activate(dd);

    /* maintain .Devices (.Device has already been set) */
    PROTECT(t = mkString(CHAR(STRING(getSymbolValue(".Device"))[0])));
    if (append)
        CDR(s) = CONS(t, R_NilValue);
    else
        CAR(s) = t;

    UNPROTECT(2);

    copyGPar(&(dd->dp), &(dd->gp));
    GReset(dd);
}

int deviceNumber(DevDesc *dd)
{
    int i;
    for (i = 1; i < R_MaxDevices; i++)
        if (R_Devices[i] == dd)
            return i;
    return 0;
}

int selectDevice(int devNum)
{
    /* valid to select nullDevice */
    if ((devNum >= 0) && (devNum < R_MaxDevices) && (R_Devices[devNum] != NULL))
    {
        DevDesc *oldd, *dd;

        if (!NoDevices())
        {
            oldd = CurrentDevice();
            oldd->dp.deactivate(oldd);
        }

        R_CurrentDevice = devNum;

        /* maintain .Device */
        gsetVar(install(".Device"), elt(getSymbolValue(".Devices"), devNum), R_NilValue);

        dd = CurrentDevice();
        if (!NoDevices())
        {
            dd->dp.activate(dd);
        }
        copyGPar(&(dd->dp), &(dd->gp));
        GReset(dd);
        return devNum;
    }
    else
        return selectDevice(nextDevice(devNum));
}

void removeDevice(int devNum)
{
    if ((devNum > 0) && (devNum < R_MaxDevices) && (R_Devices[devNum] != NULL))
    {
        int i;
        SEXP s;

        free(R_Devices[devNum]);
        R_Devices[devNum] = NULL;
        R_NumDevices -= 1;

        /* maintain .Devices */
        PROTECT(s = getSymbolValue(".Devices"));
        for (i = 0; i < devNum; i++)
            s = CDR(s);
        CAR(s) = mkString("");
        UNPROTECT(1);

        /* determine new current device */
        if (devNum == R_CurrentDevice)
        {
            DevDesc *dd;

            R_CurrentDevice = nextDevice(R_CurrentDevice);

            /* maintain .Device */
            gsetVar(install(".Device"), elt(getSymbolValue(".Devices"), R_CurrentDevice), R_NilValue);

            dd = CurrentDevice();
            if (!NoDevices())
            {
                dd->dp.activate(dd);
            }
            copyGPar(&(dd->dp), &(dd->gp));
            GReset(dd);
        }
    }
}

void initDisplayList();

/*
int SetDevice(SEXP name, SEXP cpars, int ncpars, SEXP npars, int nnpars)
{
    DevDesc *dd;

    if (!(dd = (DevDesc *) malloc(sizeof(DevDesc))))
        return 0;
    GInit(&dd->dp);

#ifdef Macintosh
    if(!strcmp(CHAR(STRING(name)[0]), "Macintosh"))
        if(MacDeviceDriver(cpars, ncpars, npars, nnpars))
            goto have_device;
#endif

#ifndef Macintosh
    if(!strcmp(CHAR(STRING(name)[0]), "postscript")) {
        if (PSDeviceDriver(dd, cpars, ncpars, npars, nnpars))
            goto have_device;
    }
#endif

#ifdef Unix
    if(!strcmp(CHAR(STRING(name)[0]), "X11")) {
        if (X11DeviceDriver(dd, cpars, ncpars, npars, nnpars))
            goto have_device;
    }
#endif

#ifdef Win32
    if( !strcmp(CHAR(STRING(name)[0]), "Win32") )
        if(WinDeviceDriver(cpars, ncpars, npars, nnpars))
            goto have_device;
#endif

#ifdef Unix
#ifdef SOON
    if( !strcmp(name, "XFig") )
        if(XFigDeviceDriver(cpars, ncpars, npars, nnpars))
            goto have_device;
#endif
#endif

    return 0;

have_device:
    gsetVar(install(".Device"), mkString(CHAR(STRING(name)[0])), R_NilValue);
    addDevice(dd);
    initDisplayList(dd);

    return 1;
}
*/

void KillDevice(DevDesc *dd)
{
    dd->dp.close(dd);
    removeDevice(deviceNumber(dd));
}

void killDevice(int devNum)
{
    if (!NoDevices() && (devNum > 0) && (devNum < R_MaxDevices) && (R_Devices[devNum] != NULL))
    {
        R_Devices[devNum]->dp.close(R_Devices[devNum]);
        removeDevice(devNum);
    }
}

void KillAllDevices()
{
    /* don't try to close or remove the null device ! */
    while (R_NumDevices > 1)
        killDevice(R_CurrentDevice);
}

/* code for maintaining */
/* DISPLAY LISTS */

SEXP lastElt(SEXP list)
{
    /* Return the last element of a list */

    SEXP result = R_NilValue;
    while (list != R_NilValue)
    {
        result = list;
        list = CDR(list);
    }
    return result;
}

void initDisplayList(DevDesc *dd)
{
    /* init saveParams */
    copyGPar(&(dd->dp), &(dd->dpSaved));
    dd->displayList = R_NilValue;
}

void recordGraphicOperation(SEXP op, SEXP args, DevDesc *dd)
{
    SEXP lastOperation = lastElt(dd->displayList);
    if (dd->displayListOn)
    {
        SEXP newOperation = CONS(op, args);
        if (lastOperation == R_NilValue)
            dd->displayList = CONS(newOperation, R_NilValue);
        else
            CDR(lastOperation) = CONS(newOperation, R_NilValue);
    }
}

void restoredpSaved(DevDesc *dd)
{
    /* NOTE that not all params should be restored before playing */
    /* the display list (e.g., don't restore the device size) */

    int i, j;

    /* do NOT restore basic device driver properties;  they are */
    /* either meant to be different (e.g., left, right, bottom, top */
    /* changed because of window resize) or never change (e.g., ipr) */

    dd->dp.state = dd->dpSaved.state;
    dd->dp.adj = dd->dpSaved.adj;
    dd->dp.ann = dd->dpSaved.ann;
    dd->dp.bg = dd->dpSaved.bg;
    dd->dp.bty = dd->dpSaved.bty;
    dd->dp.cex = dd->dpSaved.cex;
    dd->dp.col = dd->dpSaved.col;
    dd->dp.crt = dd->dpSaved.crt;
    dd->dp.err = dd->dpSaved.err;
    dd->dp.fg = dd->dpSaved.fg;
    dd->dp.font = dd->dpSaved.font;
    dd->dp.gamma = dd->dpSaved.gamma;
    dd->dp.lab[0] = dd->dpSaved.lab[0];
    dd->dp.lab[1] = dd->dpSaved.lab[1];
    dd->dp.lab[2] = dd->dpSaved.lab[2];
    dd->dp.las = dd->dpSaved.las;
    dd->dp.lty = dd->dpSaved.lty;
    dd->dp.lwd = dd->dpSaved.lwd;
    dd->dp.mgp[0] = dd->dpSaved.mgp[0];
    dd->dp.mgp[1] = dd->dpSaved.mgp[1];
    dd->dp.mgp[2] = dd->dpSaved.mgp[2];
    dd->dp.mkh = dd->dpSaved.mkh;
    dd->dp.pch = dd->dpSaved.pch;
    dd->dp.ps = dd->dpSaved.ps;
    dd->dp.smo = dd->dpSaved.smo;
    dd->dp.srt = dd->dpSaved.srt;
    dd->dp.tck = dd->dpSaved.tck;
    dd->dp.tmag = dd->dpSaved.tmag;
    dd->dp.type = dd->dpSaved.type;
    dd->dp.xaxp[0] = dd->dpSaved.xaxp[0];
    dd->dp.xaxp[1] = dd->dpSaved.xaxp[1];
    dd->dp.xaxp[2] = dd->dpSaved.xaxp[2];
    dd->dp.xaxs = dd->dpSaved.xaxs;
    dd->dp.xaxt = dd->dpSaved.xaxt;
    dd->dp.xpd = dd->dpSaved.xpd;
    dd->dp.xlog = dd->dpSaved.xlog;
    dd->dp.yaxp[0] = dd->dpSaved.yaxp[0];
    dd->dp.yaxp[1] = dd->dpSaved.yaxp[1];
    dd->dp.yaxp[2] = dd->dpSaved.yaxp[2];
    dd->dp.yaxs = dd->dpSaved.yaxs;
    dd->dp.yaxt = dd->dpSaved.yaxt;
    dd->dp.ylog = dd->dpSaved.ylog;
    dd->dp.cexbase = dd->dpSaved.cexbase;
    dd->dp.cexmain = dd->dpSaved.cexmain;
    dd->dp.cexlab = dd->dpSaved.cexlab;
    dd->dp.cexsub = dd->dpSaved.cexsub;
    dd->dp.cexaxis = dd->dpSaved.cexaxis;
    dd->dp.fontmain = dd->dpSaved.fontmain;
    dd->dp.fontlab = dd->dpSaved.fontlab;
    dd->dp.fontsub = dd->dpSaved.fontsub;
    dd->dp.fontaxis = dd->dpSaved.fontaxis;
    dd->dp.colmain = dd->dpSaved.colmain;
    dd->dp.collab = dd->dpSaved.collab;
    dd->dp.colsub = dd->dpSaved.colsub;
    dd->dp.colaxis = dd->dpSaved.colaxis;

    /* must restore layout parameters;  the different graphics */
    /* regions and coordinate transformations will be recalculated */
    /* but they need all of the layout information restored for this */
    /* to happen correctly */

    dd->dp.devmode = dd->dpSaved.devmode;
    dd->dp.fig[0] = dd->dpSaved.fig[0];
    dd->dp.fig[1] = dd->dpSaved.fig[1];
    dd->dp.fig[2] = dd->dpSaved.fig[2];
    dd->dp.fig[3] = dd->dpSaved.fig[3];
    dd->dp.fin[0] = dd->dpSaved.fin[0];
    dd->dp.fin[1] = dd->dpSaved.fin[1];
    dd->dp.fin[2] = dd->dpSaved.fin[2];
    dd->dp.fin[3] = dd->dpSaved.fin[3];
    dd->dp.fUnits = dd->dpSaved.fUnits;
    dd->dp.defaultFigure = dd->dpSaved.defaultFigure;
    dd->dp.mar[0] = dd->dpSaved.mar[0];
    dd->dp.mar[1] = dd->dpSaved.mar[1];
    dd->dp.mar[2] = dd->dpSaved.mar[2];
    dd->dp.mar[3] = dd->dpSaved.mar[3];
    dd->dp.mai[0] = dd->dpSaved.mai[0];
    dd->dp.mai[1] = dd->dpSaved.mai[1];
    dd->dp.mai[2] = dd->dpSaved.mai[2];
    dd->dp.mai[3] = dd->dpSaved.mai[3];
    dd->dp.mUnits = dd->dpSaved.mUnits;
    dd->dp.mex = dd->dpSaved.mex;
    dd->dp.numrows = dd->dpSaved.numrows;
    dd->dp.numcols = dd->dpSaved.numcols;
    dd->dp.currentFigure = dd->dpSaved.currentFigure;
    dd->dp.lastFigure = dd->dpSaved.lastFigure;
    for (i = 0; i < dd->dpSaved.numrows; i++)
    {
        dd->dp.heights[i] = dd->dpSaved.heights[i];
        dd->dp.cmHeights[i] = dd->dpSaved.cmHeights[i];
    }
    for (j = 0; j < dd->dpSaved.numcols; j++)
    {
        dd->dp.widths[j] = dd->dpSaved.widths[j];
        dd->dp.cmWidths[j] = dd->dpSaved.cmWidths[j];
    }
    for (i = 0; i < dd->dpSaved.numrows; i++)
        for (j = 0; j < dd->dpSaved.numcols; j++)
        {
            dd->dp.order[i][j] = dd->dpSaved.order[i][j];
            dd->dp.respect[i][j] = dd->dpSaved.respect[i][j];
        }
    dd->dp.rspct = dd->dpSaved.rspct;
    dd->dp.layout = dd->dpSaved.layout;
    dd->dp.mfind = dd->dpSaved.mfind;
    dd->dp.new = dd->dpSaved.new;
    dd->dp.oma[0] = dd->dpSaved.oma[0];
    dd->dp.oma[1] = dd->dpSaved.oma[1];
    dd->dp.oma[2] = dd->dpSaved.oma[2];
    dd->dp.oma[3] = dd->dpSaved.oma[3];
    dd->dp.omi[0] = dd->dpSaved.omi[0];
    dd->dp.omi[1] = dd->dpSaved.omi[1];
    dd->dp.omi[2] = dd->dpSaved.omi[2];
    dd->dp.omi[3] = dd->dpSaved.omi[3];
    dd->dp.omd[0] = dd->dpSaved.omd[0];
    dd->dp.omd[1] = dd->dpSaved.omd[1];
    dd->dp.omd[2] = dd->dpSaved.omd[2];
    dd->dp.omd[3] = dd->dpSaved.omd[3];
    dd->dp.oUnits = dd->dpSaved.oUnits;
    dd->dp.plt[0] = dd->dpSaved.plt[0];
    dd->dp.plt[1] = dd->dpSaved.plt[1];
    dd->dp.plt[2] = dd->dpSaved.plt[2];
    dd->dp.plt[3] = dd->dpSaved.plt[3];
    dd->dp.pin[0] = dd->dpSaved.pin[0];
    dd->dp.pin[1] = dd->dpSaved.pin[1];
    dd->dp.pin[2] = dd->dpSaved.pin[2];
    dd->dp.pin[3] = dd->dpSaved.pin[3];
    dd->dp.pUnits = dd->dpSaved.pUnits;
    dd->dp.defaultPlot = dd->dpSaved.defaultPlot;
    dd->dp.pty = dd->dpSaved.pty;
    dd->dp.usr[0] = dd->dpSaved.usr[0];
    dd->dp.usr[1] = dd->dpSaved.usr[1];
    dd->dp.usr[2] = dd->dpSaved.usr[2];
    dd->dp.usr[3] = dd->dpSaved.usr[3];
    dd->dp.logusr[0] = dd->dpSaved.logusr[0];
    dd->dp.logusr[1] = dd->dpSaved.logusr[1];
    dd->dp.logusr[2] = dd->dpSaved.logusr[2];
    dd->dp.logusr[3] = dd->dpSaved.logusr[3];
}

void playDisplayList(DevDesc *dd)
{
    int ask, savedDevice;
    SEXP theList = dd->displayList;
    ask = dd->gp.ask;
    dd->gp.ask = 1;
    restoredpSaved(dd);
    copyGPar(&(dd->dp), &(dd->gp));
    GReset(dd);
    savedDevice = curDevice();
    selectDevice(deviceNumber(dd));
    while (theList != R_NilValue)
    {
        SEXP theOperation = CAR(theList);
        SEXP op = CAR(theOperation);
        SEXP args = CDR(theOperation);
        PRIMFUN(op)(R_NilValue, op, args, R_NilValue);
        theList = CDR(theList);
    }
    dd->gp.ask = ask;
    selectDevice(savedDevice);
}

void copyDisplayList(int fromDevice)
{
    DevDesc *dd = CurrentDevice();
    dd->displayList = R_Devices[fromDevice]->displayList;
    dd->dpSaved = R_Devices[fromDevice]->dpSaved;
    playDisplayList(dd);
    if (!dd->displayListOn)
        initDisplayList(dd);
}

void inhibitDisplayList(DevDesc *dd)
{
    initDisplayList(dd);
    dd->displayListOn = 0;
}
