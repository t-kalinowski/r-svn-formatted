/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998--2000  Robert Gentleman, Ross Ihaka and the
 *                            R Development Core Team
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

#include "dataentry.h"
#include <stdlib.h>
static Atom _XA_WM_PROTOCOLS, protocol;

static void clearwindow(void);
static int newcol;
static int xmaxused, ymaxused;
static int CellModified;

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif
#define BOXW(x) (x < 100 ? boxw[x] : box_w)

/*
   The spreadsheet function returns a list of vectors. The types of
   these vectors can be specified by the user as can their names. It
   the names are specified they are set during initialization. The
   user can change these via a menu interface, they can also change
   the type.

   The vectors are created too long and if they need to be increased
   this is done by using the next higher power of 2. They start 100
   long. To cut them to the correct length for return you need to know
   the largest row number that was assigned to. LEVELS (sxpinfo.gp) is
   used to keep track of this, separately for each vector. Vectors are
   initialized to NA when they are created so that NA is returned for
   any cell that was not set by the user.  So that coercion back and
   forth maintains values of ssNA_REAL and ssNA_STRING I have set
   ssNA_STRING to be coerceVector(ssNA_REAL), very weird but easy.

   In Macintosh we need to call the main event loop to get
   events. This ensures that the spreadsheet interacts well with the
   other windows. Under X windows we let the window manager handle
   those sorts of details.

 */

static char *menu_label[] = {
    "Real",
    "Character",
    "Change Name",
};

/*
   ssNewVector is just an interface to allocVector but it lets us
   set the fields to NA. We need to have a special NA for reals and
   strings so that we can differentiate between uninitialized elements
   in the vectors and user supplied NA's; hence ssNA_REAL and ssNA_STRING
 */

static SEXP ssNewVector(SEXPTYPE type, int vlen)
{
    SEXP tvec;
    int j;

    tvec = allocVector(type, vlen);
    for (j = 0; j < vlen; j++)
        if (type == REALSXP)
            REAL(tvec)[j] = ssNA_REAL;
        else if (type == STRSXP)
            STRING(tvec)[j] = STRING(ssNA_STRING)[0];
    LEVELS(tvec) = 0;
    return (tvec);
}

SEXP RX11_dataentry(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP tvec2, tvec, colmodes, indata;
    SEXPTYPE type;
    int i, j, len, nprotect, tmp;
    RCNTXT cntxt;

    nprotect = 0; /* count the PROTECT()s */

    PROTECT(indata = VectorToPairList(CAR(args)));
    nprotect++;
    PROTECT(colmodes = VectorToPairList(CADR(args)));
    nprotect++;

    if (!isList(indata) || !isList(colmodes))
        errorcall(call, "invalid argument");

    /* initialize the constants */

    bufp = buf;
    ne = 0;
    currentexp = 0;
    nneg = 0;
    ndecimal = 0;
    clength = 0;
    ccol = 1;
    crow = 1;
    colmin = 1;
    rowmin = 1;
    ssNA_REAL = -NA_REAL;
    tvec = allocVector(REALSXP, 1);
    REAL(tvec)[0] = ssNA_REAL;
    PROTECT(ssNA_STRING = coerceVector(tvec, STRSXP));
    nprotect++;
    bwidth = 5;
    hwidth = 30;

    /* setup inputlist  */

    if (indata != R_NilValue)
    {
        xmaxused = 0;
        ymaxused = 0;
        PROTECT(inputlist = duplicate(indata));
        nprotect++;
        for (tvec = inputlist, tvec2 = colmodes; tvec != R_NilValue; tvec = CDR(tvec), tvec2 = CDR(tvec2))
        {
            type = TYPEOF(CAR(tvec));
            xmaxused++;
            if (CAR(tvec2) != R_NilValue)
                type = str2type(CHAR(STRING(CAR(tvec2))[0]));
            if (type != STRSXP)
                type = REALSXP;
            if (CAR(tvec) == R_NilValue)
            {
                if (type == NILSXP)
                    type = REALSXP;
                CAR(tvec) = ssNewVector(type, 100);
                TAG(tvec) = install("var1");
                LEVELS(CAR(tvec)) = 0;
            }
            else if (!isVector(CAR(tvec)))
                errorcall(call, "invalid type for value");
            else
            {
                if (TYPEOF(CAR(tvec)) != type)
                    CAR(tvec) = coerceVector(CAR(tvec), type);
                tmp = LEVELS(CAR(tvec)) = LENGTH(CAR(tvec));
                ymaxused = max(tmp, ymaxused);
            }
        }
    }
    else if (colmodes == R_NilValue)
    {
        PROTECT(inputlist = allocList(1));
        nprotect++;
        CAR(inputlist) = ssNewVector(REALSXP, 100);
        TAG(inputlist) = install("var1");
        LEVELS(CAR(inputlist)) = 0;
    }
    else
    {
        errorcall(call, "invalid parameter(s) ");
    }

    /* start up the window, more initializing in here */
    if (initwin())
        errorcall(call, "invalid device");

    /* set up a context which will close the window if there is an error */
    begincontext(&cntxt, 8, R_NilValue, R_NilValue, R_NilValue, R_NilValue);
    cntxt.cend = &closewin;

    highlightrect();

    eventloop();

    endcontext(&cntxt);
    closewin();

    /* drop out unused columns */
    i = 0;
    for (tvec = inputlist; tvec != R_NilValue; tvec = CDR(tvec))
        if (CAR(tvec) == R_NilValue)
        {
            if (i == 0)
                inputlist = CDR(inputlist);
            else
            {
                tvec2 = nthcdr(inputlist, (i - 1));
                SETCDR(tvec2, CDR(tvec));
            }
        }
        else
            i++;

    for (tvec = inputlist; tvec != R_NilValue; tvec = CDR(tvec))
    {
        len = LEVELS(CAR(tvec));
        if (LENGTH(CAR(tvec)) != len)
        {
            tvec2 = ssNewVector(TYPEOF(CAR(tvec)), len);
            PROTECT(tvec);
            for (j = 0; j < len; j++)
                if (TYPEOF(CAR(tvec)) == REALSXP)
                {
                    if (REAL(CAR(tvec))[j] != ssNA_REAL)
                        REAL(tvec2)[j] = REAL(CAR(tvec))[j];
                    else
                        REAL(tvec2)[j] = NA_REAL;
                }
                else if (TYPEOF(CAR(tvec)) == STRSXP)
                {
                    if (!streql(CHAR(STRING(CAR(tvec))[j]), CHAR(STRING(ssNA_STRING)[0])))
                        STRING(tvec2)[j] = STRING(CAR(tvec))[j];
                    else
                        STRING(tvec2)[j] = NA_STRING;
                }
                else
                    error("dataentry: internal memory problem");
            CAR(tvec) = tvec2;
            UNPROTECT(1);
        }
    }

    UNPROTECT(nprotect);
    return PairToVectorList(inputlist);
}

/* Window Drawing Routines */

void drawwindow()
{
    int i, w, dw;
    XWindowAttributes attribs;

    /* if there is an active cell enter the data in it */
    closerect();

    /* now set up the window with the new dimensions */
    XGetWindowAttributes(iodisplay, iowindow, &attribs);
    bwidth = attribs.border_width;
    fullwindowWidth = attribs.width;
    fullwindowHeight = attribs.height;
    windowWidth = w = 2 * bwidth + boxw[0] + BOXW(colmin);
    nwide = 2;
    for (i = 2; i < 100; i++)
    { /* 100 on-screen columns cannot occur */
        dw = BOXW(i + colmin - 1);
        if ((w += dw) > fullwindowWidth)
        {
            nwide = i;
            windowWidth = w - dw;
            break;
        }
    }
    nhigh = (fullwindowHeight - 2 * bwidth - hwidth) / box_h;
    windowHeight = nhigh * box_h + 2 * bwidth;

    clearwindow();

    for (i = 1; i < nhigh; i++)
        drawrectangle(0, hwidth + i * box_h, boxw[0], box_h, 1, 1);
    /* so row 0 and col 0 are reserved for labels */
    colmax = colmin + (nwide - 2);
    rowmax = rowmin + (nhigh - 2);
    printlabs();
    if (inputlist != R_NilValue)
        for (i = colmin; i <= colmax; i++)
            drawcol(i);

    /* draw the quit box */

    i = textwidth("Quit", 4);
    drawrectangle(fullwindowWidth - 6 - bwidth - i, 3, i + 4, hwidth - 6, 1, 1);
    drawtext(fullwindowWidth - 4 - bwidth - i, hwidth - 7, "Quit", 4);

    highlightrect();

    Rsync();
}

void doHscroll(int oldcol)
{
    int i, w, dw;
    int oldnwide = nwide, oldwindowWidth = windowWidth;

    /* horizontal re-position */
    windowWidth = w = 2 * bwidth + boxw[0] + BOXW(colmin);
    nwide = 2;
    for (i = 2; i < 100; i++)
    {
        dw = BOXW(i + colmin - 1);
        if ((w += dw) > fullwindowWidth)
        {
            nwide = i;
            windowWidth = w - dw;
            break;
        }
    }
    colmax = colmin + (nwide - 2);
    if (oldcol < colmin)
    { /* drop oldcol...colmin-1 */
        dw = boxw[0];
        for (i = oldcol; i < colmin; i++)
            dw += BOXW(i);
        copyH(dw, boxw[0], oldwindowWidth - dw + 1);
        dw = oldwindowWidth - BOXW(oldcol) + 1;
        cleararea(dw, hwidth, fullwindowWidth - dw, fullwindowHeight);
        /* oldnwide includes the row labels */
        for (i = oldcol + oldnwide - 1; i <= colmax; i++)
            drawcol(i);
    }
    else
    {
        /* move one or more cols left */
        dw = BOXW(colmin);
        copyH(boxw[0], boxw[0] + dw, windowWidth - dw + 1);
        dw = windowWidth + 1;
        cleararea(dw, hwidth, fullwindowWidth - dw, fullwindowHeight);
        drawcol(colmin);
    }

    highlightrect();

    Rsync();
}

/* find_coords finds the coordinates of the upper left corner of the
   given cell on the screen: row and col are on-screen coords */

void find_coords(int row, int col, int *xcoord, int *ycoord)
{
    int i, w;
    w = bwidth;
    if (col > 0)
        w += boxw[0];
    for (i = 1; i < col; i++)
        w += BOXW(i + colmin - 1);
    *xcoord = w;
    *ycoord = bwidth + hwidth + box_h * row;
}

/* draw the window with the top left box at column wcol and row wrow */

void jumpwin(int wcol, int wrow)
{
    if (wcol < 0 || wrow < 0)
    {
        bell();
        return;
    }
    closerect();
    if (colmin != wcol || rowmin != wrow)
    {
        colmin = wcol;
        rowmin = wrow;
        drawwindow();
    }
    else
        highlightrect();
}

void advancerect(int which)
{

    /* if we are in the header, changing a name then only down is
       allowed */
    if (crow < 1 && which != DOWN)
    {
        bell();
        return;
    }

    closerect();

    switch (which)
    {
    case UP:
        if (crow == 1)
        {
            if (rowmin == 1)
                bell();
            else
                jumppage(UP);
        }
        else
            crow--;
        break;
    case DOWN:
        if (crow == (nhigh - 1))
            jumppage(DOWN);
        else
            crow++;
        break;
    case RIGHT:
        if (ccol == (nwide - 1))
            jumppage(RIGHT);
        else
            ccol++;
        break;
    case LEFT:
        if (ccol == 1)
        {
            if (colmin == 1)
                bell();
            else
                jumppage(LEFT);
        }
        else
            ccol--;
        break;
    default:
        UNIMPLEMENTED("advancerect");
    }

    highlightrect();
}

static char *get_col_name(int col)
{
    SEXP tmp;
    static char clab[15];
    if (col <= length(inputlist))
    {
        tmp = nthcdr(inputlist, col - 1);
        if (TAG(tmp) != R_NilValue)
            return CHAR(PRINTNAME(TAG(tmp)));
    }
    sprintf(clab, "var%d", col);
    return clab;
}

static int get_col_width(int col)
{
    int i, w = 0, w1;
    char *strp;
    SEXP tmp;
    if (col <= length(inputlist))
    {
        tmp = nthcdr(inputlist, col - 1);
        if (tmp == R_NilValue)
            return box_w;
        PrintDefaults(R_NilValue);
        if (TAG(tmp) != R_NilValue)
        {
            strp = CHAR(PRINTNAME(TAG(tmp)));
        }
        else
            strp = "var12";
        w = textwidth(strp, strlen(strp));
        tmp = CAR(tmp);
        PrintDefaults(R_NilValue);
        for (i = 0; i < (int)LEVELS(tmp); i++)
        {
            strp = EncodeElement(tmp, i, 0);
            w1 = textwidth(strp, strlen(strp));
            if (w1 > w)
                w = w1;
        }
        if (w < 0.5 * box_w)
            w = 0.5 * box_w;
        if (w < 0.8 * box_w)
            w += 0.1 * box_w;
        return w + 8;
    }
    return box_w;
}

typedef enum
{
    UNKNOWNN,
    NUMERIC,
    CHARACTER
} CellType;

static CellType get_col_type(int col)
{
    SEXP tmp;
    CellType res = UNKNOWNN;

    if (col <= length(inputlist))
    {
        tmp = CAR(nthcdr(inputlist, col - 1));
        if (TYPEOF(tmp) == REALSXP)
            res = NUMERIC;
        if (TYPEOF(tmp) == STRSXP)
            res = CHARACTER;
    }
    return res;
}

/* whichcol is absolute col no, col is position on screen */
void drawcol(int whichcol)
{
    int i, src_x, src_y, len, col = whichcol - colmin + 1, bw = BOXW(whichcol);
    char *clab;
    SEXP tmp;

    find_coords(0, col, &src_x, &src_y);
    cleararea(src_x, src_y, bw, windowHeight);
    for (i = 0; i < nhigh; i++)
        drawrectangle(src_x, hwidth + i * box_h, bw, box_h, 1, 1);

    /* now fill it in if it is active */
    clab = get_col_name(whichcol);
    printstring(clab, strlen(clab), 0, col, 0);

    if (length(inputlist) >= whichcol)
    {
        tmp = nthcdr(inputlist, whichcol - 1);
        if (CAR(tmp) != R_NilValue)
        {
            len = min(rowmax, LEVELS(CAR(tmp)));
            for (i = (rowmin - 1); i < len; i++)
                printelt(CAR(tmp), i, i - rowmin + 2, col);
        }
    }
    Rsync();
}

/* whichrow is absolute row no */
void drawrow(int whichrow)
{
    int i, src_x, src_y, lenip, row = whichrow - rowmin + 1, w;
    char rlab[15];
    SEXP tvec;

    find_coords(row, 0, &src_x, &src_y);
    cleararea(src_x, src_y, windowWidth, box_h);
    drawrectangle(src_x, src_y, boxw[0], box_h, 1, 1);

    sprintf(rlab, "%4d", whichrow);
    printstring(rlab, strlen(rlab), row, 0, 0);

    w = bwidth + boxw[0];
    for (i = colmin; i <= colmax; i++)
    {
        drawrectangle(w, src_y, BOXW(i), box_h, 1, 1);
        w += BOXW(i);
    }

    lenip = length(inputlist);
    for (i = colmin; i <= colmax; i++)
    {
        if (i > lenip)
            break;
        tvec = CAR(nthcdr(inputlist, i - 1));
        if (tvec != R_NilValue)
            if (whichrow <= (int)LEVELS(tvec))
                printelt(tvec, whichrow - 1, row, i - colmin + 1);
    }

    Rsync();
}

/* printelt: print the correct value from vector[vrow] into the
   spreadsheet in row ssrow and col sscol */

/* WARNING: This has no check that you're not beyond the end of the
   vector. Caller must check. */

void printelt(SEXP invec, int vrow, int ssrow, int sscol)
{
    char *strp;
    PrintDefaults(R_NilValue);
    if (TYPEOF(invec) == REALSXP)
    {
        if (REAL(invec)[vrow] != ssNA_REAL)
        {
            strp = EncodeElement(invec, vrow, 0);
            printstring(strp, strlen(strp), ssrow, sscol, 0);
        }
    }
    else if (TYPEOF(invec) == STRSXP)
    {
        if (!streql(CHAR(STRING(invec)[vrow]), CHAR(STRING(ssNA_STRING)[0])))
        {
            strp = EncodeElement(invec, vrow, 0);
            printstring(strp, strlen(strp), ssrow, sscol, 0);
        }
    }
    else
        error("dataentry: internal memory error");
}

static void drawelt(int whichrow, int whichcol)
{
    int i;
    char *clab;
    SEXP tmp;

    if (whichrow == 0)
    {
        clab = get_col_name(whichcol + colmin - 1);
        printstring(clab, strlen(clab), 0, whichcol, 0);
    }
    else
    {
        if (length(inputlist) >= whichcol + colmin - 1)
        {
            tmp = nthcdr(inputlist, whichcol + colmin - 2);
            if (CAR(tmp) != R_NilValue && (i = rowmin + whichrow - 2) < (int)LEVELS(CAR(tmp)))
                printelt(CAR(tmp), i, whichrow, whichcol);
        }
        else
            printstring("", 0, whichrow, whichcol, 0);
    }

    Rsync();
}

void jumppage(int dir)
{
    int i, w, oldcol, wcol;

    switch (dir)
    {
    case UP:
        rowmin--;
        rowmax--;
        copyarea(0, hwidth + box_h, 0, hwidth + 2 * box_h);
        drawrow(rowmin);
        break;
    case DOWN:
        rowmin++;
        rowmax++;
        copyarea(0, hwidth + 2 * box_h, 0, hwidth + box_h);
        drawrow(rowmax);
        break;
    case LEFT:
        colmin--;
        doHscroll(colmin + 1);
        break;
    case RIGHT:
        oldcol = colmin;
        wcol = colmin + ccol + 1; /* column to be selected */
                                  /* There may not be room to fit the next column in */
        w = fullwindowWidth - boxw[0] - BOXW(colmax + 1);
        for (i = colmax; i >= oldcol; i--)
        {
            w -= BOXW(i);
            if (w < 0)
            {
                colmin = i + 1;
                break;
            }
        }
        ccol = wcol - colmin;
        doHscroll(oldcol);
        break;
    }
}
/* draw a rectangle, used to highlight/downlight the current box */

void printrect(int lwd, int fore)
{
    int x, y;
    find_coords(crow, ccol, &x, &y);
    drawrectangle(x + lwd - 1, y + lwd - 1, BOXW(ccol + colmin - 1) - lwd + 1, box_h - lwd + 1, lwd, fore);
    Rsync();
}

void downlightrect()
{
    printrect(2, 0);
    printrect(1, 1);
}

void highlightrect()
{
    printrect(2, 1);
}

static SEXP getccol()
{
    SEXP tmp, tmp2;
    int i, len, newlen, wcol, wrow;
    SEXPTYPE type;
    char cname[10];

    wcol = ccol + colmin - 1;
    wrow = crow + rowmin - 1;
    if (length(inputlist) < wcol)
        inputlist = listAppend(inputlist, allocList(wcol - length(inputlist)));
    tmp = nthcdr(inputlist, wcol - 1);
    newcol = 0;
    if (CAR(tmp) == R_NilValue)
    {
        newcol = 1;
        xmaxused = wcol;
        len = max(100, wrow);
        CAR(tmp) = ssNewVector(REALSXP, len);
        if (TAG(tmp) == R_NilValue)
        {
            sprintf(cname, "var%d", wcol);
            TAG(tmp) = install(cname);
        }
    }
    if (!isVector(CAR(tmp)))
        error("internal type error in dataentry");
    len = LENGTH(CAR(tmp));
    type = TYPEOF(CAR(tmp));
    if (len < wrow)
    {
        for (newlen = len * 2; newlen < wrow; newlen *= 2)
            ;
        tmp2 = ssNewVector(type, newlen);
        for (i = 0; i < len; i++)
            if (type == REALSXP)
                REAL(tmp2)[i] = REAL(CAR(tmp))[i];
            else if (type == STRSXP)
                STRING(tmp2)[i] = STRING(CAR(tmp))[i];
            else
                error("internal type error in dataentry");
        LEVELS(tmp2) = LEVELS(CAR(tmp));
        CAR(tmp) = tmp2;
    }
    return (tmp);
}

/* close up the entry to a cell, put the value that has been entered
   into the correct place and as the correct type */

extern double R_strtod(char *c, char **end); /* in coerce.c */

void closerect()
{
    SEXP cvec, c0vec, tvec;
    int wcol = ccol + colmin - 1, wrow = rowmin + crow - 1, wrow0;

    *bufp = '\0';

    /* first check to see if anything has been entered */
    if (CellModified)
    {
        if (crow == 0)
        {
            if (clength != 0)
            {
                /* then we are entering a new column name */
                if (length(inputlist) < wcol)
                    inputlist = listAppend(inputlist, allocList((wcol - length(inputlist))));
                tvec = nthcdr(inputlist, wcol - 1);
                TAG(tvec) = install(buf);
                printstring(buf, strlen(buf), 0, wcol, 0);
            }
            else
            {
                sprintf(buf, "var%d", ccol);
                printstring(buf, strlen(buf), 0, wcol, 0);
            }
        }
        else
        {
            c0vec = getccol();
            cvec = CAR(c0vec);
            wrow0 = (int)LEVELS(cvec);
            if (wrow > wrow0)
                LEVELS(cvec) = wrow;
            ymaxused = max(ymaxused, wrow);
            if (clength != 0)
            {
                /* do it this way to ensure NA, Inf, ...  can get set */
                char *endp;
                double new = R_strtod(buf, &endp);
                int warn = !isBlankString(endp);
                if (TYPEOF(cvec) == STRSXP)
                {
                    tvec = allocString(strlen(buf));
                    strcpy(CHAR(tvec), buf);
                    STRING(cvec)[wrow - 1] = tvec;
                }
                else
                    REAL(cvec)[wrow - 1] = new;
                if (newcol & warn)
                {
                    /* change mode to character */
                    int levs = LEVELS(cvec);
                    cvec = CAR(c0vec) = coerceVector(cvec, STRSXP);
                    LEVELS(cvec) = levs;
                    tvec = allocString(strlen(buf));
                    strcpy(CHAR(tvec), buf);
                    STRING(cvec)[wrow - 1] = tvec;
                }
            }
            else
            {
                if (TYPEOF(cvec) == STRSXP)
                    STRING(cvec)[wrow - 1] = NA_STRING;
                else
                    REAL(cvec)[wrow - 1] = NA_REAL;
            }
            drawelt(crow, ccol); /* to get the cell scrolling right */
            if (wrow > wrow0)
                drawcol(wcol); /* to fill in NAs */
        }
    }
    CellModified = 0;

    downlightrect();

    ndecimal = 0;
    nneg = 0;
    ne = 0;
    currentexp = 0;
    clength = 0;
    bufp = buf;
}

/* print a null terminated string, check to see if it is longer than
   the print area and print it, left adjusted if necessary; clear the
   area of previous text; */

void printstring(char *ibuf, int buflen, int row, int col, int left)
{
    int i, x_pos, y_pos, bw;
    char buf[45], *pc = buf;

    find_coords(row, col, &x_pos, &y_pos);
    if (col == 0)
        bw = boxw[0];
    else
        bw = BOXW(col + colmin - 1);
    cleararea(x_pos + 2, y_pos + 2, bw - 3, box_h - 3);
    strncpy(buf, ibuf, buflen);
    if (left)
    {
        for (i = buflen; i > 1; i--)
        {
            if (textwidth(pc, i) < (bw - text_offset))
                break;
            *(++pc) = '<';
        }
    }
    else
    {
        for (i = buflen; i > 1; i--)
        {
            if (textwidth(buf, i) < (bw - text_offset))
                break;
            *(pc + i - 2) = '>';
        }
    }
    drawtext(x_pos + text_offset, y_pos + box_h - text_offset, pc, i);
    Rsync();
}

void clearrect()
{
    int x_pos, y_pos;

    find_coords(crow, ccol, &x_pos, &y_pos);
    cleararea(x_pos, y_pos, BOXW(ccol + colmin - 1), box_h);
    Rsync();
}

/* handlechar has to be able to parse decimal numbers and strings,
   depending on the current column type, only printing characters
   should get this far */

/* --- Not true! E.g. ESC ends up in here... */

void handlechar(char *text)
{
    int c = text[0];

    if (c == '\033')
    {
        CellModified = 0;
        clength = 0;
        drawelt(crow, ccol);
        return;
    }
    else
        CellModified = 1;

    if (clength == 0)
    {

        if (crow == 0) /* variable name */
            currentexp = 3;
        else
            switch (get_col_type(ccol + colmin - 1))
            {
            case NUMERIC:
                currentexp = 1;
                break;
            default:
                currentexp = 2;
            }
        clearrect();
        highlightrect();
    }

    if (currentexp == 1) /* we are parsing a number */
        switch (c)
        {
        case '-':
            if (nneg == 0)
                nneg++;
            else
                goto donehc;
            break;
        case '.':
            if (ndecimal == 0)
                ndecimal++;
            else
                goto donehc;
            break;
        case 'e':
        case 'E':
            if (ne == 0)
            {
                nneg = ndecimal = 0; /* might have decimal in exponent */
                ne++;
            }
            else
                goto donehc;
            break;
        default:
            if (!isdigit((int)text[0]))
                goto donehc;
            break;
        }
    if (currentexp == 3)
    {
        if (isspace(c))
            goto donehc;
        if (clength == 0)
        {
            if (c != '.' && !isalpha(c))
                goto donehc;
            else if (c != '.' && !isalnum(c))
                goto donehc;
        }
    }

    if (clength++ > 29)
    {
        warning("dataentry: expression too long");
        clength--;
        goto donehc;
    }

    *bufp++ = text[0];
    printstring(buf, clength, crow, ccol, 1);
    return;

donehc:
    bell();
}

void printlabs()
{
    char clab[10], *p;
    int i;

    for (i = colmin; i <= colmax; i++)
    {
        p = get_col_name(i);
        printstring(p, strlen(p), 0, i - colmin + 1, 0);
    }
    for (i = rowmin; i <= rowmax; i++)
    {
        sprintf(clab, "%4d", i);
        printstring(clab, strlen(clab), i - rowmin + 1, 0, 0);
    }
}

/* ================ X11-specific ================ */

/* find out whether the button click was in the quit box */
static int checkquit(int xw)
{
    int wi;

    wi = textwidth("Quit", 4);
    if ((xw < fullwindowWidth - bwidth - 2) && (xw > fullwindowWidth - bwidth - wi - 6))
        return 1;
    else
        return 0;
}

/* when a buttonpress event happens find the square that is being
   pointed to if the pointer is in the header we need to see if the
   quit button was pressed and if so quit. This is done by having
   findcell return an int which is zero if we should quit and one
   otherwise */

static int findcell()
{

    int xw, yw, xr, yr, wcol = 0, wrow, i, w;
    unsigned int keys;
    Window root, child;

    closerect();
    XQueryPointer(iodisplay, iowindow, &root, &child, &xr, &yr, &xw, &yw, &keys);

    if (keys & Button1Mask)
    { /* left click

/* check to see if the click was in the header */

        if (yw < hwidth + bwidth)
        {
            if (checkquit(xw))
                return 1;
            else
                return 0;
        }

        /* see if it is in the row labels */
        if (xw < bwidth + boxw[0])
        {
            bell();
            highlightrect();
            return 0;
        }
        /* translate to box coordinates */
        wrow = (yw - bwidth - hwidth) / box_h;
        w = bwidth + boxw[0];
        for (i = 1; i <= nwide; i++)
            if ((w += BOXW(i + colmin - 1)) > xw)
            {
                wcol = i;
                break;
            }

        /* next check to see if it is in the column labels */

        if (yw < hwidth + bwidth + box_h)
        {
            if (xw > bwidth + boxw[0])
                popupmenu(xr, yr, wcol, wrow);
            else
            {
                highlightrect();
                bell();
            }
        }
        else if (wcol != ccol || wrow != crow)
        {
            ccol = wcol;
            crow = wrow;
        }
    }
    if (keys & Button2Mask)
    { /* Paste: eventually */
    }
    highlightrect();
    return 0;
}

/* Event Loop Functions */

static void eventloop()
{
    int done;
    DEEvent ioevent;

    done = 0;
    while (done == 0)
    {
        if (NextEvent(&ioevent))
        {
            switch (WhichEvent(ioevent))
            {
            case activateEvt:
                drawwindow();
                break;
            case mouseDown:
                done = doMouseDown(&ioevent);
                break;
            case keyDown:
                doSpreadKey(0, &ioevent);
                break;
            case MappingNotify:
                RefreshKeyboardMapping(&ioevent);
                break;
            case ConfigureNotify:
                doConfigure(&ioevent);
                break;
            case ClientMessage:
                if (ioevent.xclient.message_type == _XA_WM_PROTOCOLS && ioevent.xclient.data.l[0] == protocol)
                {
                    /* user clicked on `close' aka `destroy' */
                    closewin();
                    done = 1;
                }
                break;
            }
        }
    }
}

int doMouseDown(DEEvent *event)
{
    return findcell();
}

static void doSpreadKey(int key, DEEvent *event)
{
    KeySym iokey;
    char text[1];

    iokey = GetKey(event);
    text[0] = GetCharP(event);

    if (CheckControl(event))
        doControl(event);
    else if ((iokey == XK_Return) || (iokey == XK_KP_Enter) || (iokey == XK_Linefeed) || (iokey == XK_Down))
        advancerect(DOWN);
    else if (iokey == XK_Left)
        advancerect(LEFT);
    else if (iokey == XK_Right)
        advancerect(RIGHT);
    else if (iokey == XK_Up)
        advancerect(UP);
    else if ((iokey == XK_BackSpace) || (iokey == XK_Delete))
    {
        if (clength > 0)
        {
            clength--;
            bufp--;
            printstring(buf, clength, crow, ccol, 1);
        }
        else
            bell();
    }
    else if (iokey == XK_Tab)
    {
        if (CheckShift(event))
            advancerect(LEFT);
        else
            advancerect(RIGHT);
    }
    else if (iokey == XK_Home)
        jumpwin(1, 1);
    else if (IsModifierKey(iokey))
    {
    }
    else
        handlechar(text);
}

static int NextEvent(DEEvent *ioevent)
{
    XNextEvent(iodisplay, ioevent);
    return 1;
}

static int WhichEvent(DEEvent ioevent)
{
    return ioevent.type;
}

static KeySym GetKey(DEEvent *event)
{
    int i;
    char text[1];
    KeySym iokey;

    i = XLookupString(event, text, 10, &iokey, 0);
    return iokey;
}

static char GetCharP(DEEvent *event)
{
    int i;
    char text[1];
    KeySym iokey;

    i = XLookupString(event, text, 1, &iokey, 0);
    return text[0];
}

static int CheckControl(DEEvent *event)
{
    return (*event).xkey.state & ControlMask;
}

static int CheckShift(DEEvent *event)
{
    return (*event).xkey.state & ShiftMask;
}

static void doControl(DEEvent *event)
{
    int i;
    char text[1];
    KeySym iokey;

    (*event).xkey.state = 0;
    i = XLookupString(event, text, 1, &iokey, 0);
    /* one row overlap when scrolling: top line <--> bottom line */
    switch (text[0])
    {
    case 'b':
        i = rowmin - nhigh + 2;
        jumpwin(colmin, max(1, i));
        break;
    case 'f':
        jumpwin(colmin, rowmax);
        break;
    case 'l':
        for (i = 1; i <= min(100, xmaxused); i++)
            boxw[i] = get_col_width(i);
        drawwindow();
        break;
    }
}

static void doConfigure(DEEvent *event)
{
    if ((fullwindowWidth != (*event).xconfigure.width) || (fullwindowHeight != (*event).xconfigure.height))
        drawwindow();
}

static void RefreshKeyboardMapping(DEEvent *event)
{
    XRefreshKeyboardMapping(event);
}

/* Initialize/Close Windows */

void closewin()
{
    XFreeGC(iodisplay, iogc);
    XDestroyWindow(iodisplay, iowindow);
    XCloseDisplay(iodisplay);
}

static int R_X11Err(Display *dsp, XErrorEvent *event)
{
    char buff[1000];
    XGetErrorText(dsp, event->error_code, buff, 1000);
    warning("X11 protocol error: %s", buff);
    return 0;
}

static int R_X11IOErr(Display *dsp)
{
    error("X11 fatal IO error: please save work and shut down R");
    return 0; /* but should never get here */
}

/* set up the window, print the grid and column/row labels */

int initwin()
{
    int i, twidth, w;
    int ioscreen;
    unsigned long iowhite, ioblack;
    char ioname[] = "R DataEntryWindow";
    char digits[] = "123456789.0";
    Window root;
    XEvent ioevent;
    XSetWindowAttributes winattr;
    XWindowAttributes attribs;

    if ((iodisplay = XOpenDisplay(NULL)) == NULL)
        return (1);
    XSetErrorHandler(R_X11Err);
    XSetIOErrorHandler(R_X11IOErr);

    /* Get Font Loaded if we can */

    font_info = XLoadQueryFont(iodisplay, font_name);
    if (font_info == NULL)
        return 1; /* ERROR */

    /* find out how wide the input boxes should be and set up the
       window size defaults */

    twidth = textwidth(digits, strlen(digits));
    box_w = twidth + 4;
    box_h = font_info->max_bounds.ascent + font_info->max_bounds.descent + 4;
    text_offset = 2 + font_info->max_bounds.descent;
    windowWidth = 0;
    windowHeight = 26 * box_h + hwidth + 2;
    boxw[0] = textwidth("1234 ", 5) + 8;
    for (i = 1; i < 100; i++)
        boxw[i] = get_col_width(i);
    w = 0;
    for (i = 0; i <= xmaxused; i++)
    {
        w += boxw[i];
        if (w > 800)
        {
            windowWidth = w - boxw[i];
            break;
        }
    }
    if (windowWidth == 0)
        windowWidth = w;
    windowWidth += 2;

    ioscreen = DefaultScreen(iodisplay);
    iowhite = WhitePixel(iodisplay, ioscreen);
    ioblack = BlackPixel(iodisplay, ioscreen);
    hand_cursor = XCreateFontCursor(iodisplay, XC_hand2);

    iohint.x = 0;
    iohint.y = 0;
    iohint.width = windowWidth;
    iohint.height = windowHeight;
    iohint.flags = PPosition | PSize;
    root = DefaultRootWindow(iodisplay);

    if ((iowindow = XCreateSimpleWindow(iodisplay, root, iohint.x, iohint.y, iohint.width, iohint.height, bwidth,
                                        ioblack, iowhite)) == 0)
        return 1;

    XSetStandardProperties(iodisplay, iowindow, ioname, ioname, None, ioname, 0, &iohint);

    winattr.backing_store = Always;
    XChangeWindowAttributes(iodisplay, iowindow, CWBackingStore, &winattr);

    /* set up protocols so that window manager sends */
    /* me an event when user "destroys" window */
    _XA_WM_PROTOCOLS = XInternAtom(iodisplay, "WM_PROTOCOLS", 0);
    protocol = XInternAtom(iodisplay, "WM_DELETE_WINDOW", 0);
    XSetWMProtocols(iodisplay, iowindow, &protocol, 1);

    iogc = XCreateGC(iodisplay, iowindow, 0, 0);
    XSetFont(iodisplay, iogc, font_info->fid);
    XSetBackground(iodisplay, iogc, iowhite);
    XSetForeground(iodisplay, iogc, BlackPixel(iodisplay, DefaultScreen(iodisplay)));
    XSetLineAttributes(iodisplay, iogc, 1, LineSolid, CapRound, JoinRound);

    XSelectInput(iodisplay, iowindow, ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask);
    XMapRaised(iodisplay, iowindow);

    /* now set up the menu-window, for now use the same text
       dimensions as above */

    menuwindow = XCreateSimpleWindow(iodisplay, root, 0, 0, twidth, 4 * box_h, 2, ioblack, iowhite);
    for (i = 0; i < 4; i++)
    {
        menupanes[i] = XCreateSimpleWindow(iodisplay, menuwindow, 0, box_h * i, twidth, box_h, 1, ioblack, iowhite);
        XSelectInput(iodisplay, menupanes[i], ButtonPressMask | ButtonReleaseMask | ExposureMask);
    }

    /* XMapSubwindows(iodisplay, menuwindow); */

    winattr.override_redirect = True;
    XChangeWindowAttributes(iodisplay, menuwindow, CWBackingStore | CWOverrideRedirect, &winattr);
    Rsync();

    /* this next sequence makes sure the window is up and ready before
       you start drawing in it */

    XNextEvent(iodisplay, &ioevent);
    if (ioevent.xany.type == Expose)
    {
        while (ioevent.xexpose.count)
            XNextEvent(iodisplay, &ioevent);
    }
    XGetWindowAttributes(iodisplay, iowindow, &attribs);
    bwidth = attribs.border_width;
    fullwindowWidth = attribs.width;
    fullwindowHeight = attribs.height;

    /* set the active rectangle to be the upper left one */
    crow = 1;
    ccol = 1;
    CellModified = 0;
    return 0;
}

/* MAC/X11 BASICS */

static void bell()
{
    XBell(iodisplay, 20);
}

static void cleararea(int xpos, int ypos, int width, int height)
{
    XClearArea(iodisplay, iowindow, xpos, ypos, width, height, 0);
}

static void clearwindow()
{
    XClearWindow(iodisplay, iowindow);
}

static void copyarea(int src_x, int src_y, int dest_x, int dest_y)
{
    int mx = max(src_x, dest_x), my = max(src_y, dest_y);
    XCopyArea(iodisplay, iowindow, iowindow, iogc, src_x, src_y, fullwindowWidth - mx, fullwindowHeight - my, dest_x,
              dest_y);
    Rsync();
}

static void copyH(int src_x, int dest_x, int width)
{
    XCopyArea(iodisplay, iowindow, iowindow, iogc, src_x + bwidth, hwidth, width, windowHeight + 1, dest_x + bwidth,
              hwidth);
}

#if 0
static void drawline(int fromx, int fromy, int tox, int toy)
{
    XDrawLine(iodisplay, iowindow, iogc, fromx, fromy, tox, toy);
}
#endif

static void drawrectangle(int xpos, int ypos, int width, int height, int lwd, int fore)
{
    if (fore == 0)
        XSetForeground(iodisplay, iogc, WhitePixel(iodisplay, DefaultScreen(iodisplay)));
    else
        XSetForeground(iodisplay, iogc, BlackPixel(iodisplay, DefaultScreen(iodisplay)));
    XSetLineAttributes(iodisplay, iogc, lwd, LineSolid, CapRound, JoinRound);
    XDrawRectangle(iodisplay, iowindow, iogc, xpos, ypos, width, height);
}

static void drawtext(int xpos, int ypos, char *text, int len)
{
    XDrawImageString(iodisplay, iowindow, iogc, xpos, ypos, text, len);
    Rsync();
}

static void Rsync()
{
    XSync(iodisplay, 0);
}

static int textwidth(char *text, int nchar)
{
    int t1;

    t1 = XTextWidth(font_info, text, nchar);
    return t1;
}

/* Menus */

void popupmenu(int x_pos, int y_pos, int col, int row)
{
    int i, button, levs;
    char name[20];
    XEvent event;
    Window selected_pane;
    SEXP tvec;

    XMapSubwindows(iodisplay, menuwindow);
    XMapRaised(iodisplay, menuwindow);
    XMoveWindow(iodisplay, menuwindow, x_pos, y_pos);

    /* now fill in the menu panes with the correct information */

    if (length(inputlist) < col + colmin - 1)
        inputlist = listAppend(inputlist, allocList(col + colmin - 1 - length(inputlist)));
    tvec = nthcdr(inputlist, col + colmin - 2);
    if (TAG(tvec) != R_NilValue)
        sprintf(name, "  %s", CHAR(PRINTNAME(TAG(tvec))));
    else
        sprintf(name, " COLUMN %d", col + colmin - 1);
    XDrawString(iodisplay, menupanes[0], iogc, 3, box_h - 3, name, strlen(name));
    for (i = 1; i < 4; i++)
        XDrawString(iodisplay, menupanes[i], iogc, 3, box_h - 3, menu_label[i - 1], strlen(menu_label[i - 1]));
    if (CAR(tvec) == R_NilValue || TYPEOF(CAR(tvec)) == REALSXP)
        XDrawString(iodisplay, menupanes[1], iogc, box_w - 20, box_h - 3, "X", 1);
    else
        XDrawString(iodisplay, menupanes[2], iogc, box_w - 20, box_h - 3, "X", 1);

    /*
      start an event loop; we're looking for a button press and a button
      release in the same window
    */

    while (1)
    {
        XNextEvent(iodisplay, &event);
        if (event.type == ButtonPress)
        {
            button = event.xbutton.button;
            selected_pane = event.xbutton.window;
            for (i = 0; selected_pane != menupanes[i]; i++)
                if (i >= 4)
                    goto done;
            while (1)
            {
                while (XCheckTypedEvent(iodisplay, ButtonPress, &event))
                    ;
                XMaskEvent(iodisplay, ButtonReleaseMask, &event);
                if (event.xbutton.button == button)
                    break;
            }
            if (selected_pane == event.xbutton.window)
            {
                for (i = 0; selected_pane != menupanes[i]; i++)
                    ;
                switch (i)
                {
                case 0:
                    bell();
                    break;
                case 1:
                    if (CAR(tvec) == R_NilValue)
                        CAR(tvec) = ssNewVector(REALSXP, 100);
                    levs = LEVELS(CAR(tvec));
                    CAR(tvec) = coerceVector(CAR(tvec), REALSXP);
                    LEVELS(CAR(tvec)) = levs;
                    goto done;
                case 2:
                    if (CAR(tvec) == R_NilValue)
                        CAR(tvec) = ssNewVector(STRSXP, 100);
                    levs = LEVELS(CAR(tvec));
                    CAR(tvec) = coerceVector(CAR(tvec), STRSXP);
                    LEVELS(CAR(tvec)) = levs;
                    goto done;
                case 3:
                    closerect();
                    ccol = col;
                    crow = 0;
                    clearrect();
                    goto done;
                }
            }
        }
        /* this doesn't work and perhaps I should move it up to the
           main control loop */
        else if (event.type == Expose)
        {
            if (event.xexpose.window == menuwindow)
            {
                XDrawString(iodisplay, menupanes[0], iogc, 3, box_h - 3, name, strlen(name));
                for (i = 1; i < 4; i++)
                    XDrawString(iodisplay, menupanes[i], iogc, 3, box_h - 3, menu_label[i - 1],
                                strlen(menu_label[i - 1]));
            }
        }
    }
done:
    popdownmenu();
    highlightrect();
}

void popdownmenu()
{
    XUnmapWindow(iodisplay, menuwindow);
    XUnmapSubwindows(iodisplay, menuwindow);
}
