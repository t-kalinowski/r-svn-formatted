/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998--2005  Robert Gentleman, Ross Ihaka and the
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

/* <UTF8>
   Used XTextWidth and XDrawText, so need to use fontsets

   Also needed input context.
*/

/* The version for R 2.1.0 is partly based on patches by
   Ei-ji Nakama <nakama@ki.rim.or.jp> for use with Japanese fonts. */

#define DPRINTS(x) printf(#x "=[%s]\n", x)
#define DPRINTX(x) printf(#x "=%x\n", x)
#define DPRINTD(x) printf(#x "=%d\n", x)

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>

/* don't use X11 function prototypes (which tend to ...): */
#define NeedFunctionPrototypes 0
/* XFree 4.4.0 forgets to define this if NeedFunctionPrototypes=0 */
#ifndef _Xconst
#define _Xconst const
#endif
#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/cursorfont.h>

#include "Defn.h"
#include "Print.h"

#ifdef SUPPORT_MBCS
/* This only uses a FontSet in a MBCS */
#define USE_FONTSET 1
/* In theory we should do this, but it works less well
# ifdef X_HAVE_UTF8_STRING
#  define HAVE_XUTF8TEXTEXTENTS 1
#  define HAVE_XUTF8DRAWSTRING 1
#  define HAVE_XUTF8DRAWIMAGESTRING 1
# endif */
#endif

#ifndef HAVE_KEYSYM
#define KeySym int
#endif

#define DEEvent XEvent

typedef enum
{
    UP,
    DOWN,
    LEFT,
    RIGHT
} DE_DIRECTION;

typedef enum
{
    UNKNOWNN,
    NUMERIC,
    CHARACTER
} CellType;

/* EXPORTS : */
SEXP RX11_dataentry(SEXP call, SEXP op, SEXP args, SEXP rho);

/* Local Function Definitions */

static void advancerect(DE_DIRECTION);
static int CheckControl(DEEvent *);
static int CheckShift(DEEvent *);
static int checkquit(int);
static void clearrect(void);
static void clearwindow(void);
static void closerect(void);
static void closewin(void);
static void copycell(void);
static void doControl(DEEvent *);
static int doMouseDown(DEEvent *);
static void eventloop(void);
static void doSpreadKey(int, DEEvent *);
static void downlightrect(void);
static void drawwindow(void);
static void drawcol(int);
static void drawrow(int);
static void find_coords(int, int, int *, int *);
static int findcell(void);
static char *GetCharP(DEEvent *);
static KeySym GetKey(DEEvent *);
static void handlechar(char *);
static void highlightrect(void);
static Rboolean initwin(void);
static void jumppage(DE_DIRECTION);
static void jumpwin(int, int);
static void pastecell(int, int);
static void popdownmenu(void);
static void popupmenu(int, int, int, int);
static void printlabs(void);
static void printrect(int, int);
static void printstring(char *, int, int, int, int);
static void printelt(SEXP, int, int, int);
static void RefreshKeyboardMapping(DEEvent *);
static void cell_cursor_init(void);

/* Functions to hide Xlib calls */
static void bell(void);
static void cleararea(int, int, int, int);
static void copyH(int, int, int);
static void copyarea(int, int, int, int);
static void doConfigure(DEEvent *ioevent);
#if 0
static void drawline(int, int, int, int);
#endif
static void drawrectangle(int, int, int, int, int, int);
static void drawtext(int, int, char *, int);
#if 0
static int  NextEvent(DEEvent *ioevent);
#endif
static void RefreshKeyboardMapping(DEEvent *ioevent);
static void Rsync(void);
static int textwidth(char *, int);
static int WhichEvent(DEEvent ioevent);

static char *get_col_name(int col);
static int get_col_width(int col);
static CellType get_col_type(int col);
#ifdef USE_FONTSET
static void calc_pre_edit_pos(void);
#endif
static int last_wchar_bytes(char *);
static SEXP work, names, lens;
static PROTECT_INDEX wpi, npi, lpi;
static SEXP ssNA_STRING;
static double ssNA_REAL;

/* Global variables needed for the graphics */

static int box_w;       /* width of a box */
static int boxw[100];   /* widths of cells */
static int box_h;       /* height of a box */
static int windowWidth; /* current width of the window */
static int fullwindowWidth;
static int windowHeight; /* current height of the window */
static int fullwindowHeight;
static int currentexp; /* boolean: whether an cell is active */
static int crow;       /* current row */
static int ccol;       /* current column */
static int nwide, nhigh;
static int colmax, colmin, rowmax, rowmin;
static int ndecimal; /* count decimal points */
static int ne;       /* count exponents */
static int nneg;     /* indicate whether its a negative */
static int clength;  /* number of characters currently entered */

#define BOOSTED_BUF_SIZE 201
static char buf[BOOSTED_BUF_SIZE]; /* boosted to allow for MBCS */
static char *bufp;
static int bwidth; /* width of the border */
static int hwidth; /* width of header  */
static int text_offset;

static SEXP ssNewVector(SEXPTYPE, int);

static Atom _XA_WM_PROTOCOLS, protocol;

static Rboolean newcol, CellModified;
static int nboxchars;
static int xmaxused, ymaxused;
static int box_coords[6];
static char copycontents[sizeof(buf) + 1];
static int labdigs = 4;
static char labform[6];

/* Xwindows Globals */

static Display *iodisplay;
static Window iowindow, menuwindow, menupanes[4];
static GC iogc;
static XFontStruct *font_info;
static char *font_name = "9x15";
#ifdef USE_FONTSET
static Status status;
static XFontSet font_set;
static XFontStruct **fs_list;
static int font_set_cnt;
static char *fontset_name = "-*-fixed-medium-r-normal--13-*-*-*-*-*-*-*";
static XIM ioim;
static XIMStyle ioim_style;
static XIMStyles *ioim_styles;

/*
 * XIM:
 * OverTheSpot XIMPreeditPosition | XIMStatusArea;
 * OffTheSpot  XIMPreeditArea     | XIMStatusArea;
 * Root        XIMPreeditNothing  | XIMStatusNothing;
 */
static XIMStyle preedit_styles[] = {
    XIMPreeditPosition, XIMPreeditArea, XIMPreeditNothing, XIMPreeditNone, (XIMStyle)NULL,
};
static XIMStyle status_styles[] = {
    XIMStatusArea,
    XIMStatusNothing,
    XIMStatusNone,
    (XIMStyle)NULL,
};
static XIC ioic;
#endif

#define mouseDown ButtonPress
#define keyDown KeyPress
#define activateEvt MapNotify
#define updateEvt Expose

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif
#define BOXW(x) (min(((x < 100 && nboxchars == 0) ? boxw[x] : box_w), fullwindowWidth - boxw[0] - 2 * bwidth - 2))

/*
  Underlying assumptions (for this version R >= 1.8.0)

  The data are stored in a list `work', with unused columns having
  NULL entries.  The names for the list are in `names', which should
  have a name for all displayable columns (up to xmaxused).
  The *used* lengths of the columns are in `lens': this needs only be
  set for non-NULL columns.

  If the list was originally length(0), that should work with
  0 pre-defined cols.  (It used to have 1 pre-defined numeric column.)

  All row and col numbers are 1-based.

  BDR May 2003
 */

/*
   The spreadsheet function returns a list of vectors. The types of
   these vectors can be specified by the user as can their names. It
   the names are specified they are set during initialization. The
   user can change these via a menu interface, they can also change
   the type.

   The vectors are created too long and if they need to be increased
   this is done by using the next higher power of 2. They start 100
   long.  Vectors are initialized to NA when they are created so that
   NA is returned for any cell that was not set by the user.  So that
   coercion back and forth maintains values of ssNA_REAL and
   ssNA_STRING I have set ssNA_STRING to be coerceVector(ssNA_REAL),
   very weird but easy.

   In Macintosh we needed to call the main event loop to get
   events. This ensures that the spreadsheet interacts well with the
   other windows. Under X windows we let the window manager handle
   those sorts of details.

 */

static char *menu_label[] = {
    " Real",
    " Character",
    "Change Name ",
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
            SET_STRING_ELT(tvec, j, STRING_ELT(ssNA_STRING, 0));
    return (tvec);
}

static void closewin_cend(void *data)
{
    closewin();
}

SEXP RX11_dataentry(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP colmodes, tnames, tvec, tvec2, work2;
    SEXPTYPE type;
    int i, j, cnt, len, nprotect;
    RCNTXT cntxt;
    char clab[25];

    nprotect = 0; /* count the PROTECT()s */
    PROTECT_WITH_INDEX(work = duplicate(CAR(args)), &wpi);
    nprotect++;
    colmodes = CADR(args);
    tnames = getAttrib(work, R_NamesSymbol);

    if (TYPEOF(work) != VECSXP || TYPEOF(colmodes) != VECSXP)
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

    /* setup work, names, lens  */
    xmaxused = length(work);
    ymaxused = 0;
    PROTECT_WITH_INDEX(lens = allocVector(INTSXP, xmaxused), &lpi);
    nprotect++;

    if (isNull(tnames))
    {
        PROTECT_WITH_INDEX(names = allocVector(STRSXP, xmaxused), &npi);
        for (i = 0; i < xmaxused; i++)
        {
            sprintf(clab, "var%d", i);
            SET_STRING_ELT(names, i, mkChar(clab));
        }
    }
    else
        PROTECT_WITH_INDEX(names = duplicate(tnames), &npi);
    nprotect++;
    for (i = 0; i < xmaxused; i++)
    {
        int len = LENGTH(VECTOR_ELT(work, i));
        INTEGER(lens)[i] = len;
        ymaxused = max(len, ymaxused);
        type = TYPEOF(VECTOR_ELT(work, i));
        if (LENGTH(colmodes) > 0 && !isNull(VECTOR_ELT(colmodes, i)))
            type = str2type(CHAR(STRING_ELT(VECTOR_ELT(colmodes, i), 0)));
        if (type != STRSXP)
            type = REALSXP;
        if (isNull(VECTOR_ELT(work, i)))
        {
            if (type == NILSXP)
                type = REALSXP;
            SET_VECTOR_ELT(work, i, ssNewVector(type, 100));
        }
        else if (!isVector(VECTOR_ELT(work, i)))
            errorcall(call, "invalid type for value");
        else
        {
            if (TYPEOF(VECTOR_ELT(work, i)) != type)
                SET_VECTOR_ELT(work, i, coerceVector(VECTOR_ELT(work, i), type));
        }
    }

    /* start up the window, more initializing in here */
    if (initwin())
        errorcall(call, "invalid device");

    /* set up a context which will close the window if there is an error */
    begincontext(&cntxt, CTXT_CCODE, R_NilValue, R_BaseEnv, R_BaseEnv, R_NilValue, R_NilValue);
    cntxt.cend = &closewin_cend;
    cntxt.cenddata = NULL;

    highlightrect();

    cell_cursor_init();

    eventloop();

    endcontext(&cntxt);
    closewin();

    /* drop out unused columns */
    for (i = 0, cnt = 0; i < xmaxused; i++)
        if (!isNull(VECTOR_ELT(work, i)))
            cnt++;
    if (cnt < xmaxused)
    {
        PROTECT(work2 = allocVector(VECSXP, cnt));
        nprotect++;
        for (i = 0, j = 0; i < xmaxused; i++)
        {
            if (!isNull(VECTOR_ELT(work, i)))
            {
                SET_VECTOR_ELT(work2, j, VECTOR_ELT(work, i));
                INTEGER(lens)[j] = INTEGER(lens)[i];
                SET_STRING_ELT(names, j, STRING_ELT(names, i));
                j++;
            }
        }
        REPROTECT(names = lengthgets(names, cnt), npi);
    }
    else
        work2 = work;

    for (i = 0; i < LENGTH(work2); i++)
    {
        len = INTEGER(lens)[i];
        tvec = VECTOR_ELT(work2, i);
        if (LENGTH(tvec) != len)
        {
            tvec2 = ssNewVector(TYPEOF(tvec), len);
            for (j = 0; j < len; j++)
            {
                if (TYPEOF(tvec) == REALSXP)
                {
                    if (REAL(tvec)[j] != ssNA_REAL)
                        REAL(tvec2)[j] = REAL(tvec)[j];
                    else
                        REAL(tvec2)[j] = NA_REAL;
                }
                else if (TYPEOF(tvec) == STRSXP)
                {
                    if (!streql(CHAR(STRING_ELT(tvec, j)), CHAR(STRING_ELT(ssNA_STRING, 0))))
                        SET_STRING_ELT(tvec2, j, STRING_ELT(tvec, j));
                    else
                        SET_STRING_ELT(tvec2, j, NA_STRING);
                }
                else
                    error("dataentry: internal memory problem");
            }
            SET_VECTOR_ELT(work2, i, tvec2);
        }
    }

    setAttrib(work2, R_NamesSymbol, names);
    UNPROTECT(nprotect);
    return work2;
}

/* Window Drawing Routines */

static void setcellwidths(void)
{
    int i, w, dw;

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
}

static void drawwindow(void)
{
    int i, st;
    XWindowAttributes attribs;

    /* if there is an active cell enter the data in it */
    /*
     * case colname input Expose not use.
     * closerect();
     */

    /* now set up the window with the new dimensions */
    XGetWindowAttributes(iodisplay, iowindow, &attribs);
    bwidth = attribs.border_width;
    fullwindowWidth = attribs.width;
    fullwindowHeight = attribs.height;
    setcellwidths();
    nhigh = (fullwindowHeight - 2 * bwidth - hwidth) / box_h;
    windowHeight = nhigh * box_h + 2 * bwidth;

    clearwindow();

    for (i = 1; i < nhigh; i++)
        drawrectangle(0, hwidth + i * box_h, boxw[0], box_h, 1, 1);
    /* so row 0 and col 0 are reserved for labels */
    colmax = colmin + (nwide - 2);
    rowmax = rowmin + (nhigh - 2);
    printlabs();
    for (i = colmin; i <= colmax; i++)
        drawcol(i);

    /* draw the quit etc boxes */

    i = textwidth("Quit", 4);
    box_coords[0] = st = fullwindowWidth - 6 - bwidth;
    box_coords[1] = st - i;
    drawrectangle(st - i, 3, i + 4, hwidth - 6, 1, 1);
    drawtext(st + 2 - i, hwidth - 7, "Quit", 4);

    box_coords[4] = st = st - 5 * i;
    i = textwidth("Paste", 5);
    box_coords[5] = st - i;
    drawrectangle(st - i, 3, i + 4, hwidth - 6, 1, 1);
    drawtext(st + 2 - i, hwidth - 7, "Paste", 5);

    box_coords[2] = st = st - 2 * i;
    i = textwidth("Copy", 4);
    box_coords[3] = st - i;
    drawrectangle(st - i, 3, i + 4, hwidth - 6, 1, 1);
    drawtext(st + 2 - i, hwidth - 7, "Copy", 4);

    highlightrect();

    Rsync();
}

static void doHscroll(int oldcol)
{
    int i, dw;
    int oldnwide = nwide, oldwindowWidth = windowWidth;

    /* horizontal re-position */
    setcellwidths();
    colmax = colmin + (nwide - 2);
    if (oldcol < colmin)
    { /* drop oldcol...colmin - 1 */
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
    cell_cursor_init();

    Rsync();
}

/* find_coords finds the coordinates of the upper left corner of the
   given cell on the screen: row and col are on-screen coords */

static void find_coords(int row, int col, int *xcoord, int *ycoord)
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

static void jumpwin(int wcol, int wrow)
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
        closerect();
        drawwindow();
    }
    else
        highlightrect();
}

static void advancerect(DE_DIRECTION which)
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

    cell_cursor_init();
}

static void cell_cursor_init(void)
{
    int i, whichrow = crow + rowmin - 1, whichcol = ccol + colmin - 1;
    SEXP tmp;

    memset(buf, 0, sizeof(buf));

    if (crow == 0)
    {
        strncpy(buf, get_col_name(whichcol), BOOSTED_BUF_SIZE - 1);
    }
    else
    {
        if (length(work) >= whichcol)
        {
            tmp = VECTOR_ELT(work, whichcol - 1);
            if (tmp != R_NilValue && (i = whichrow - 1) < LENGTH(tmp))
            {
                PrintDefaults(R_NilValue);
                if (TYPEOF(tmp) == REALSXP)
                {
                    if (REAL(tmp)[i] != ssNA_REAL)
                        strncpy(buf, EncodeElement(tmp, i, 0, '.'), BOOSTED_BUF_SIZE - 1);
                }
                else if (TYPEOF(tmp) == STRSXP)
                {
                    if (!streql(CHAR(STRING_ELT(tmp, i)), CHAR(STRING_ELT(ssNA_STRING, 0))))
                        strncpy(buf, EncodeElement(tmp, i, 0, '.'), BOOSTED_BUF_SIZE - 1);
                }
            }
        }
    }
    buf[BOOSTED_BUF_SIZE - 1] = '\0';
    clength = strlen(buf);
    bufp = buf + clength;
}

static char *get_col_name(int col)
{
    static char clab[25];
    if (col <= xmaxused)
    {
        /* don't use NA labels */
        SEXP tmp = STRING_ELT(names, col - 1);
        if (tmp != NA_STRING)
            return (CHAR(tmp));
    }
    sprintf(clab, "var%d", col);
    return clab;
}

static int get_col_width(int col)
{
    int i, w = 0, w1;
    char *strp;
    SEXP tmp, lab;

    if (nboxchars > 0)
        return box_w;
    if (col <= xmaxused)
    {
        tmp = VECTOR_ELT(work, col - 1);
        if (isNull(tmp))
            return box_w;
        /* don't use NA labels */
        lab = STRING_ELT(names, col - 1);
        if (lab != NA_STRING)
            strp = CHAR(lab);
        else
            strp = "var12";
        PrintDefaults(R_NilValue);

        w = textwidth(strp, strlen(strp));
        for (i = 0; i < INTEGER(lens)[col - 1]; i++)
        {
            strp = EncodeElement(tmp, i, 0, '.');
            w1 = textwidth(strp, strlen(strp));
            if (w1 > w)
                w = w1;
        }
        if (w < 0.5 * box_w)
            w = 0.5 * box_w;
        if (w < 0.8 * box_w)
            w += 0.1 * box_w;
        if (w > 600)
            w = 600;
        return w + 8;
    }
    return box_w;
}

static CellType get_col_type(int col)
{
    SEXP tmp;
    CellType res = UNKNOWNN;

    if (col <= xmaxused)
    {
        tmp = VECTOR_ELT(work, col - 1);
        if (TYPEOF(tmp) == REALSXP)
            res = NUMERIC;
        if (TYPEOF(tmp) == STRSXP)
            res = CHARACTER;
    }
    return res;
}

/* whichcol is absolute col no, col is position on screen */
static void drawcol(int whichcol)
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

    if (xmaxused >= whichcol)
    {
        tmp = VECTOR_ELT(work, whichcol - 1);
        if (!isNull(tmp))
        {
            len = min(rowmax, INTEGER(lens)[whichcol - 1]);
            for (i = (rowmin - 1); i < len; i++)
                printelt(tmp, i, i - rowmin + 2, col);
        }
    }
    Rsync();
}

/* whichrow is absolute row no */
static void drawrow(int whichrow)
{
    int i, src_x, src_y, row = whichrow - rowmin + 1, w;
    char rlab[15];
    SEXP tvec;

    find_coords(row, 0, &src_x, &src_y);
    cleararea(src_x, src_y, windowWidth, box_h);
    drawrectangle(src_x, src_y, boxw[0], box_h, 1, 1);

    sprintf(rlab, labform, whichrow);
    printstring(rlab, strlen(rlab), row, 0, 0);

    w = bwidth + boxw[0];
    for (i = colmin; i <= colmax; i++)
    {
        drawrectangle(w, src_y, BOXW(i), box_h, 1, 1);
        w += BOXW(i);
    }

    for (i = colmin; i <= colmax; i++)
    {
        if (i > xmaxused)
            break;
        if (!isNull(tvec = VECTOR_ELT(work, i - 1)))
            if (whichrow <= INTEGER(lens)[i - 1])
                printelt(tvec, whichrow - 1, row, i - colmin + 1);
    }

    Rsync();
}

/* printelt: print the correct value from vector[vrow] into the
   spreadsheet in row ssrow and col sscol */

/* WARNING: This has no check that you're not beyond the end of the
   vector. Caller must check. */

static void printelt(SEXP invec, int vrow, int ssrow, int sscol)
{
    char *strp;
    PrintDefaults(R_NilValue);
    if (TYPEOF(invec) == REALSXP)
    {
        if (REAL(invec)[vrow] != ssNA_REAL)
        {
            strp = EncodeElement(invec, vrow, 0, '.');
            printstring(strp, strlen(strp), ssrow, sscol, 0);
        }
    }
    else if (TYPEOF(invec) == STRSXP)
    {
        if (!streql(CHAR(STRING_ELT(invec, vrow)), CHAR(STRING_ELT(ssNA_STRING, 0))))
        {
            strp = EncodeElement(invec, vrow, 0, '.');
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
        if (xmaxused >= whichcol + colmin - 1)
        {
            tmp = VECTOR_ELT(work, whichcol + colmin - 2);
            if (!isNull(tmp) && (i = rowmin + whichrow - 2) < INTEGER(lens)[whichcol + colmin - 2])
                printelt(tmp, i, whichrow, whichcol);
        }
        else
            printstring("", 0, whichrow, whichcol, 0);
    }

    Rsync();
}

static void jumppage(DE_DIRECTION dir)
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
        if (rowmax >= 65535)
            return;
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

static void printrect(int lwd, int fore)
{
    int x, y;
    find_coords(crow, ccol, &x, &y);
    drawrectangle(x + lwd - 1, y + lwd - 1, BOXW(ccol + colmin - 1) - lwd + 1, box_h - lwd + 1, lwd, fore);
    Rsync();
}

static void downlightrect(void)
{
    printrect(2, 0);
    printrect(1, 1);
}

static void highlightrect(void)
{
    printrect(2, 1);
}

static void getccol()
{
    SEXP tmp, tmp2;
    int i, len, newlen, wcol, wrow;
    SEXPTYPE type;
    char clab[25];

    wcol = ccol + colmin - 1;
    wrow = crow + rowmin - 1;
    if (wcol > xmaxused)
    {
        /* extend work, names and lens */
        REPROTECT(work = lengthgets(work, wcol), wpi);
        REPROTECT(names = lengthgets(names, wcol), npi);
        for (i = xmaxused; i < wcol; i++)
        {
            sprintf(clab, "var%d", i + 1);
            SET_STRING_ELT(names, i, mkChar(clab));
        }
        REPROTECT(lens = lengthgets(lens, wcol), lpi);
        xmaxused = wcol;
    }
    newcol = FALSE;
    if (isNull(VECTOR_ELT(work, wcol - 1)))
    {
        newcol = TRUE;
        SET_VECTOR_ELT(work, wcol - 1, ssNewVector(REALSXP, max(100, wrow)));
        INTEGER(lens)[wcol - 1] = 0;
    }
    if (!isVector(tmp = VECTOR_ELT(work, wcol - 1)))
        error("internal type error in dataentry");
    len = INTEGER(lens)[wcol - 1];
    type = TYPEOF(tmp);
    if (len < wrow)
    {
        for (newlen = max(len * 2, 10); newlen < wrow; newlen *= 2)
            ;
        tmp2 = ssNewVector(type, newlen);
        for (i = 0; i < len; i++)
            if (type == REALSXP)
                REAL(tmp2)[i] = REAL(tmp)[i];
            else if (type == STRSXP)
                SET_STRING_ELT(tmp2, i, STRING_ELT(tmp, i));
            else
                error("internal type error in dataentry");
        SET_VECTOR_ELT(work, wcol - 1, tmp2);
    }
}

/* close up the entry to a cell, put the value that has been entered
   into the correct place and as the correct type */

static void closerect(void)
{
    SEXP cvec;
    int i, wcol = ccol + colmin - 1, wrow = rowmin + crow - 1, wrow0;
    char clab[25];

    *bufp = '\0';

    /* first check to see if anything has been entered */
    if (CellModified)
    {
        if (crow == 0)
        {
            if (clength != 0)
            {
                /* then we are entering a new column name */
                if (xmaxused < wcol)
                {
                    /* extend work, names and lens */
                    REPROTECT(work = lengthgets(work, wcol), wpi);
                    REPROTECT(names = lengthgets(names, wcol), npi);
                    for (i = xmaxused; i < wcol - 1; i++)
                    {
                        sprintf(clab, "var%d", i + 1);
                        SET_STRING_ELT(names, i, mkChar(clab));
                    }
                    REPROTECT(lens = lengthgets(lens, wcol), lpi);
                    xmaxused = wcol;
                }
                SET_STRING_ELT(names, wcol - 1, mkChar(buf));
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
            getccol();
            cvec = VECTOR_ELT(work, wcol - 1);
            wrow0 = INTEGER(lens)[wcol - 1];
            if (wrow > wrow0)
                INTEGER(lens)[wcol - 1] = wrow;
            ymaxused = max(ymaxused, wrow);
            if (clength != 0)
            {
                /* do it this way to ensure NA, Inf, ...  can get set */
                char *endp;
                double new = R_strtod(buf, &endp);
                Rboolean warn = !isBlankString(endp);
                if (TYPEOF(cvec) == STRSXP)
                    SET_STRING_ELT(cvec, wrow - 1, mkChar(buf));
                else
                    REAL(cvec)[wrow - 1] = new;
                if (newcol & warn)
                {
                    /* change mode to character */
                    SET_VECTOR_ELT(work, wcol - 1, coerceVector(cvec, STRSXP));
                    SET_STRING_ELT(VECTOR_ELT(work, wcol - 1), wrow - 1, mkChar(buf));
                }
            }
            else
            {
                if (TYPEOF(cvec) == STRSXP)
                    SET_STRING_ELT(cvec, wrow - 1, NA_STRING);
                else
                    REAL(cvec)[wrow - 1] = NA_REAL;
            }
            drawelt(crow, ccol); /* to get the cell scrolling right */
            if (wrow > wrow0)
                drawcol(wcol); /* to fill in NAs */
        }
    }
    CellModified = FALSE;

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

/* This version will only display 200 chars, but the maximum col width
   will not allow that many */
static void printstring(char *ibuf, int buflen, int row, int col, int left)
{
    int i, x_pos, y_pos, bw, bufw;
    char pbuf[BOOSTED_BUF_SIZE];
#ifdef USE_FONTSET
    int wcsbufw, j;
    wchar_t wcspbuf[BOOSTED_BUF_SIZE], *wcspc = wcspbuf;
    wchar_t wcs[BOOSTED_BUF_SIZE];
    char s[BOOSTED_BUF_SIZE];
    wchar_t *w_p;
    char *p;
    int cnt;

    find_coords(row, col, &x_pos, &y_pos);
    if (col == 0)
        bw = boxw[0];
    else
        bw = BOXW(col + colmin - 1);
    cleararea(x_pos + 2, y_pos + 2, bw - 3, box_h - 3);
    bufw = (buflen > BOOSTED_BUF_SIZE - 1) ? BOOSTED_BUF_SIZE - 1 : buflen;
    strncpy(pbuf, ibuf, bufw);
    pbuf[bufw] = '\0';

    p = pbuf;
    wcsbufw = mbsrtowcs(wcspbuf, (const char **)&p, bufw, NULL);
    wcspbuf[wcsbufw] = L'\0';
    if (left)
    {
        for (i = wcsbufw; i > 1; i--)
        {
            for (j = 0; *(wcspc + j) != L'\0'; j++)
                wcs[j] = *(wcspc + j);
            wcs[j] = L'\0';
            w_p = wcs;
            cnt = wcsrtombs(s, (const wchar_t **)&w_p, sizeof(wcs), NULL);
            s[cnt] = '\0';
            if (textwidth(s, strlen(s)) < (bw - text_offset))
                break;
            *(++wcspc) = L'<';
        }
    }
    else
    {
        for (i = wcsbufw; i > 1; i--)
        {
            for (j = 0; *(wcspc + j) != L'\0'; j++)
                wcs[j] = *(wcspc + j);
            wcs[j] = L'\0';
            w_p = wcs;
            cnt = wcsrtombs(s, (const wchar_t **)&w_p, sizeof(wcs), NULL);
            s[cnt] = '\0';
            if (textwidth(s, strlen(s)) < (bw - text_offset))
                break;
            *(wcspbuf + i - 2) = L'>';
            *(wcspbuf + i - 1) = L'\0';
        }
    }
    for (j = 0; *(wcspc + j) != L'\0'; j++)
        wcs[j] = *(wcspc + j);
    wcs[j] = L'\0';
    w_p = wcs;
    cnt = wcsrtombs(s, (const wchar_t **)&w_p, sizeof(wcs), NULL);

    drawtext(x_pos + text_offset, y_pos + box_h - text_offset, s, cnt);

#else  /* USE_FONTSET */
    char *pc = pbuf;
    find_coords(row, col, &x_pos, &y_pos);
    if (col == 0)
        bw = boxw[0];
    else
        bw = BOXW(col + colmin - 1);
    cleararea(x_pos + 2, y_pos + 2, bw - 3, box_h - 3);
    bufw = (buflen > 200) ? 200 : buflen;
    strncpy(pbuf, ibuf, bufw);
    if (left)
    {
        for (i = bufw; i > 1; i--)
        {
            if (textwidth(pc, i) < (bw - text_offset))
                break;
            *(++pc) = '<';
        }
    }
    else
    {
        for (i = bufw; i > 1; i--)
        {
            if (textwidth(pbuf, i) < (bw - text_offset))
                break;
            *(pbuf + i - 2) = '>';
        }
    }
    drawtext(x_pos + text_offset, y_pos + box_h - text_offset, pc, i);
#endif /* USE_FONTSET */

    Rsync();
}

static void clearrect(void)
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
#ifdef USE_FONTSET
#include <R_ext/rlocale.h>
#include <wchar.h>
#include <wctype.h>
#endif

/* <FIXME> This is not correct for stateful MBCSs, but that's hard to
   do as we get a char at a time */
static void handlechar(char *text)
{
    int i, c = text[0];
#ifdef USE_FONTSET
    wchar_t wcs[BOOSTED_BUF_SIZE];

    memset(wcs, 0, sizeof(wcs));
#endif

    if (c == '\033')
    { /* ESC */
        CellModified = FALSE;
        clength = 0;
        bufp = buf;
        drawelt(crow, ccol);
        cell_cursor_init();
        return;
    }
    else
        CellModified = TRUE;

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

    /* NA number? */
    if (1 == get_col_type(ccol + colmin - 1))
    {
        /* input numeric for NA of buffer , suppress NA.*/
        if (strcmp(buf, CHAR(STRING_ELT(ssNA_STRING, 0))) == 0)
        {
            buf[0] = '\0';
            clength = 0;
            bufp = buf;
        }
    }

    if (currentexp == 1)
    { /* we are parsing a number */

#ifdef USE_FONTSET
        char *mbs = text;
        int cnt = mbsrtowcs(wcs, (const char **)&mbs, strlen(text) + 1, NULL);

        for (i = 0; i < cnt; i++)
        {
            switch (wcs[i])
            {
            case L'-':
                if (nneg == 0)
                    nneg++;
                else
                    goto donehc;
                break;
            case L'.':
                if (ndecimal == 0)
                    ndecimal++;
                else
                    goto donehc;
                break;
            case L'e':
            case L'E':
                if (ne == 0)
                {
                    nneg = ndecimal = 0; /* might have decimal in exponent */
                    ne++;
                }
                else
                    goto donehc;
                break;
            default:
                if (!iswdigit(wcs[i]))
                    goto donehc;
                break;
            }
        }
#else
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
            if (!isdigit(c))
                goto donehc;
            break;
        }
#endif
    }
    if (currentexp == 3)
    {
#ifdef USE_FONTSET
        char *mbs = text;
        int cnt = mbsrtowcs(wcs, (const char **)&mbs, strlen(text) + 1, NULL);
        for (i = 0; i < cnt; i++)
        {
            if (iswspace(wcs[i]))
                goto donehc;
            if (clength == 0 && wcs[i] != L'.' && !iswalpha(wcs[i]))
                goto donehc;
            else if (wcs[i] != L'.' && !iswalnum(wcs[i]))
                goto donehc;
        }
#else
        if (isspace(c))
            goto donehc;
        if (clength == 0 && c != '.' && !isalpha(c))
            goto donehc;
        else if (c != '.' && !isalnum(c))
            goto donehc;
#endif
    }

    if (clength + strlen(text) > BOOSTED_BUF_SIZE - MB_CUR_MAX - 1)
    {
        warning("dataentry: expression too long");
        goto donehc;
    }

    for (i = 0; i < strlen(text); i++)
        *bufp++ = text[i];
    *(bufp + 1) = '\0';
    clength += strlen(text);
    printstring(buf, clength, crow, ccol, 1);
    return;

donehc:
    bell();
}

static void printlabs(void)
{
    char clab[15], *p;
    int i;

    for (i = colmin; i <= colmax; i++)
    {
        p = get_col_name(i);
        printstring(p, strlen(p), 0, i - colmin + 1, 0);
    }
    for (i = rowmin; i <= rowmax; i++)
    {
        sprintf(clab, labform, i);
        printstring(clab, strlen(clab), i - rowmin + 1, 0, 0);
    }
}

/* ================ X11-specific ================ */

/* find out whether the button click was in the quit box */
static int checkquit(int xw)
{
    if (xw > box_coords[1] && xw < box_coords[0])
        return 1;
    if (xw > box_coords[3] && xw < box_coords[2])
        return 2;
    if (xw > box_coords[5] && xw < box_coords[4])
        return 3;
    return 0;
}

/* when a buttonpress event happens find the square that is being
   pointed to if the pointer is in the header we need to see if the
   quit button was pressed and if so quit. This is done by having
   findcell return an int which is one if we should quit and zero
   otherwise */

static int findcell(void)
{

    int xw, yw, xr, yr, wcol = 0, wrow, i, w;
    unsigned int keys;
    Window root, child;

    closerect();
    XQueryPointer(iodisplay, iowindow, &root, &child, &xr, &yr, &xw, &yw, &keys);

    if (keys & Button1Mask)
    { /* left click */

        /* check to see if the click was in the header */

        if (yw < hwidth + bwidth)
        {
            i = checkquit(xw);
            if (i == 1)
                return 1;
            else if (i == 2)
                copycell();
            else if (i == 3)
                pastecell(crow, ccol);
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
        else if (wrow > nhigh - 1 || wcol > nwide - 1)
        {
            /* off the grid */
            highlightrect();
            bell();
        }
        else if (wcol != ccol || wrow != crow)
        {
            ccol = wcol;
            crow = wrow;
        }
    }
    if (keys & Button2Mask)
    { /* Paste */
        int row, col = 0;

        if (yw < hwidth + bwidth || xw < bwidth + boxw[0])
            return 0;

        /* translate to box coordinates */
        row = (yw - bwidth - hwidth) / box_h;
        w = bwidth + boxw[0];
        for (i = 1; i <= nwide; i++)
            if ((w += BOXW(i + colmin - 1)) > xw)
            {
                col = i;
                break;
            }
        pastecell(row, col);
    }
    highlightrect();
    return 0;
}

/* Event Loop Functions */

static void eventloop(void)
{
    int done;
    DEEvent ioevent;

    done = 0;
    while (done == 0)
    {
        /*
    if (NextEvent(&ioevent)) {
        */
        XNextEvent(iodisplay, &ioevent);
        {
#ifdef USE_FONTSET
            if (XFilterEvent(&ioevent, None))
            {
                if (ioic)
                {
                    XSetICFocus(ioic);
                    if (ioim_style & XIMPreeditPosition)
                        calc_pre_edit_pos();
                }
                continue;
            }
#endif

            switch (WhichEvent(ioevent))
            {
            case keyDown: /* KeyPress */
                doSpreadKey(0, &ioevent);
                break;
            case Expose:
                /*
                 * XIM on  - KeyPress - Expose
                 * XIM off - KeyPress - KeyRelease
                 * colname change XIM on mode. type Backspace.
                 */
                if (crow == 0)
                {
                    drawwindow();
                    printstring(buf, clength, crow, ccol, 1);
                }
                else
                {
                    closerect();
                    drawwindow();
                    cell_cursor_init();
                }
                break;
            case activateEvt: /* MapNotify */
                closerect();
                drawwindow();
                cell_cursor_init();
                break;
            case mouseDown: /* ButtonPress */
                done = doMouseDown(&ioevent);
                cell_cursor_init();
                break;
            case MappingNotify:
                RefreshKeyboardMapping(&ioevent);
                break;
            case ConfigureNotify:
                doConfigure(&ioevent);
                cell_cursor_init();
                break;
            case ClientMessage:
                if (ioevent.xclient.message_type == _XA_WM_PROTOCOLS && ioevent.xclient.data.l[0] == protocol)
                {
                    /* user clicked on `close' aka `destroy' */
                    done = 1;
                }
                break;
            }
        }
    }
}

static int doMouseDown(DEEvent *event)
{
    return findcell();
}

static void doSpreadKey(int key, DEEvent *event)
{
    KeySym iokey;
    char *text;

    iokey = GetKey(event);
    text = GetCharP(event);

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
#ifdef XK_Page_Up
    else if (iokey == XK_Page_Up)
    {
        int i = rowmin - nhigh + 2;
        jumpwin(colmin, max(1, i));
        cell_cursor_init();
    }
#elif defined(XK_Prior)
    else if (iokey == XK_Prior)
    {
        int i = rowmin - nhigh + 2;
        jumpwin(colmin, max(1, i));
        cell_cursor_init();
    }
#endif
#ifdef XK_Page_Down
    else if (iokey == XK_Page_Down)
    {
        jumpwin(colmin, rowmax);
        cell_cursor_init();
    }
#elif defined(XK_Next)
    else if (iokey == XK_Next)
    {
        jumpwin(colmin, rowmax);
        cell_cursor_init();
    }
#endif
    else if ((iokey == XK_BackSpace) || (iokey == XK_Delete))
    {
        if (clength > 0)
        {
            int last_w;
            last_w = last_wchar_bytes(NULL);
            clength -= last_w;
            bufp -= last_w;
            *bufp = '\0';
            CellModified = TRUE;
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
    {
        jumpwin(1, 1);
        downlightrect();
        crow = ccol = 1;
        highlightrect();
        cell_cursor_init();
    }
    else if (iokey == XK_End)
    {
        int i = ymaxused - nhigh + 2;
        jumpwin(xmaxused, max(i, 1));
        downlightrect();
        crow = ymaxused - rowmin + 1;
        ccol = 1;
        highlightrect();
        cell_cursor_init();
    }
    else if (IsModifierKey(iokey))
    {
    }
    else
    {
        handlechar(text);
    }
}

#if 0
static int NextEvent(DEEvent * ioevent)
{
    XNextEvent(iodisplay, ioevent);
    return 1;
}
#endif

static int WhichEvent(DEEvent ioevent)
{
    return ioevent.type;
}

static KeySym GetKey(DEEvent *event)
{
    char text[1];
    KeySym iokey;

    XLookupString((XKeyEvent *)event, text, 1, &iokey, NULL);
    return iokey;
}

static char *GetCharP(DEEvent *event)
{
    static char text[BOOSTED_BUF_SIZE];
    KeySym iokey;

    memset(text, 0, sizeof(text));

#ifdef USE_FONTSET
    if (mbcslocale)
    {
#ifdef HAVE_XUTF8LOOKUPSTRING
        if (utf8locale)
            Xutf8LookupString(ioic, (XKeyEvent *)event, text, sizeof(text) - clength, &iokey, &status);
        else
#endif
            XmbLookupString(ioic, (XKeyEvent *)event, text, sizeof(text) - clength, &iokey, &status);
        /* FIXME check the return code */
        if (status == XBufferOverflow)
            warning("dataentry: expression too long");
    }
    else
#endif
        XLookupString((XKeyEvent *)event, text, sizeof(text) - clength, &iokey, NULL);
    return text;
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
    XLookupString((XKeyEvent *)event, text, 1, &iokey, NULL);
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
        closerect();
        for (i = 1; i <= min(100, xmaxused); i++)
            boxw[i] = get_col_width(i);
        closerect();
        drawwindow();
        break;
    }
    cell_cursor_init();
}

static void doConfigure(DEEvent *event)
{
    if ((fullwindowWidth != (*event).xconfigure.width) || (fullwindowHeight != (*event).xconfigure.height))
    {
        closerect();
        drawwindow();
    }
}

static void RefreshKeyboardMapping(DEEvent *event)
{
    XRefreshKeyboardMapping((XMappingEvent *)event);
}

/* Initialize/Close Windows */

void closewin(void)
{
    XFreeGC(iodisplay, iogc);
#ifdef USE_FONTSET
    if (mbcslocale)
    {
        XDestroyIC(ioic);
        XCloseIM(ioim);
    }
#endif
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

static Rboolean initwin(void) /* TRUE = Error */
{
    int i, twidth, w, minwidth;
    int ioscreen;
    unsigned long iowhite, ioblack;
    /* char ioname[] = "R DataEntryWindow"; */
    char digits[] = "123456789.0";
    Window root;
    XEvent ioevent;
    XSetWindowAttributes winattr;
    XWindowAttributes attribs;
    XSizeHints *hint;
    unsigned long fevent = 0UL;
#ifdef USE_FONTSET
    int j, k;
    XVaNestedList xva_nlist;
    XPoint xpoint;
#endif

    strcpy(copycontents, "");

#ifdef USE_FONTSET
    if (!XSupportsLocale())
        warning("locale not supported by Xlib: some X ops will operate in C locale");
    if (!XSetLocaleModifiers(""))
        warning("X cannot set locale modifiers");
#endif

    if ((iodisplay = XOpenDisplay(NULL)) == NULL)
    {
        warning("unable to open display");
        return TRUE;
    }
    XSetErrorHandler(R_X11Err);
    XSetIOErrorHandler(R_X11IOErr);

    /* Get Font Loaded if we can */

#ifdef USE_FONTSET
    if (mbcslocale)
    {
        int missing_charset_count;
        char **missing_charset_list;
        char *def_string;

        char opt_fontset_name[512];

        /*
          options("X11fonts")[1] read font name
        */
        char *s = CHAR(STRING_ELT(GetOption(install("X11fonts"), R_NilValue), 0));

        if (s == NULL)
        {
            strcpy(opt_fontset_name, fontset_name);
        }
        else
        {
            sprintf(opt_fontset_name, s, "medium", "r", 12);
        }

        font_set =
            XCreateFontSet(iodisplay, opt_fontset_name, &missing_charset_list, &missing_charset_count, &def_string);
        if (missing_charset_count)
            XFreeStringList(missing_charset_list);
        if (font_set == NULL)
        {
            warning("unable to create fontset %s", fontset_name);
            return TRUE; /* ERROR */
        }
    }
    else
#endif
    {
        font_info = XLoadQueryFont(iodisplay, font_name);
        if (font_info == NULL)
        {
            warning("unable to losd font %s", font_name);
            return TRUE; /* ERROR */
        }
    }

    /* find out how wide the input boxes should be and set up the
       window size defaults */

    nboxchars = asInteger(GetOption(install("de.cellwidth"), R_GlobalEnv));
    if (nboxchars == NA_INTEGER || nboxchars < 0)
        nboxchars = 0;

    twidth = textwidth(digits, strlen(digits));

    if (nboxchars > 0)
        twidth = (twidth * nboxchars) / 10;
    box_w = twidth + 4;
#ifdef USE_FONTSET
    if (mbcslocale)
    {
        XFontSetExtents *extent = XExtentsOfFontSet(font_set);
        char **ml;
        box_h = (extent->max_logical_extent.height) + (extent->max_logical_extent.height / 5) + 4;
        font_set_cnt = XFontsOfFontSet(font_set, &fs_list, &ml);
        text_offset = 2 + fs_list[0]->max_bounds.descent;
    }
    else
#endif
    {
        box_h = font_info->max_bounds.ascent + font_info->max_bounds.descent + 4;
        text_offset = 2 + font_info->max_bounds.descent;
    }
    windowHeight = 26 * box_h + hwidth + 2;
    /* this used to presume 4 chars sufficed for row numbering */
    labdigs = max(3, 1 + floor(log10((double)ymaxused)));
    sprintf(labform, "%%%dd", labdigs);
    boxw[0] = 0.1 * labdigs * textwidth("0123456789", 10) + textwidth(" ", 1) + 8;
    for (i = 1; i < 100; i++)
        boxw[i] = get_col_width(i);

    /* try for a window width that covers all the columns, or is around
       800 pixels */
    w = windowWidth = 0;
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
    /* allow enough width for buttons */
    minwidth = 7.5 * textwidth("Paste", 5);
    if (windowWidth < minwidth)
        windowWidth = minwidth;

    ioscreen = DefaultScreen(iodisplay);
    iowhite = WhitePixel(iodisplay, ioscreen);
    ioblack = BlackPixel(iodisplay, ioscreen);

    hint = XAllocSizeHints();

    hint->x = 0;
    hint->y = 0;
    hint->width = windowWidth;
    hint->height = windowHeight;
    hint->flags = PPosition | PSize;
    /*
     * not necessary?
    hints.flags = InputHint;
    hints.input = True;
    */
    root = DefaultRootWindow(iodisplay);

    if ((iowindow = XCreateSimpleWindow(iodisplay, root, hint->x, hint->y, hint->width, hint->height, bwidth, ioblack,
                                        iowhite)) == 0)
    {
        warning("unable to open window for data editor");
        return TRUE;
    }

    /*
    XSetStandardProperties(iodisplay, iowindow, ioname, ioname, None,
               (char **)NULL, 0, iohint);
    */
    XSetWMNormalHints(iodisplay, iowindow, hint);
    XFree(hint);

    winattr.backing_store = WhenMapped;
    XChangeWindowAttributes(iodisplay, iowindow, CWBackingStore, &winattr);

    /* set up protocols so that window manager sends */
    /* me an event when user "destroys" window */
    _XA_WM_PROTOCOLS = XInternAtom(iodisplay, "WM_PROTOCOLS", 0);
    protocol = XInternAtom(iodisplay, "WM_DELETE_WINDOW", 0);
    XSetWMProtocols(iodisplay, iowindow, &protocol, 1);
    /*
     * not necessary
    XSetWMHints(iodisplay, iowindow, &hints);
     */

    iogc = XCreateGC(iodisplay, iowindow, 0, 0);

#ifdef USE_FONTSET
    if (mbcslocale)
    {
        ioim = XOpenIM(iodisplay, NULL, NULL, NULL);
        if (!ioim)
        {
            XDestroyWindow(iodisplay, iowindow);
            XCloseDisplay(iodisplay);
            warning("unable to open X Input Method");
            return TRUE;
        }

        /* search supported input style */
        XGetIMValues(ioim, XNQueryInputStyle, &ioim_styles, NULL);
        for (i = 0; i < ioim_styles->count_styles; i++)
        {
            for (j = 0; preedit_styles[j]; j++)
            {
                for (k = 0; status_styles[k]; k++)
                {
                    ioim_style = (preedit_styles[j] | status_styles[k]);
                    if (ioim_styles->supported_styles[i] == ioim_style)
                    {
                        goto loop_out;
                    }
                }
            }
        }
    loop_out:

        /* create input context */
        xpoint.x = 0;
        xpoint.y = 0;
        xva_nlist = XVaCreateNestedList(0, XNFontSet, font_set, XNSpotLocation, &xpoint, NULL);

        ioic = XCreateIC(ioim, XNInputStyle, ioim_style, XNClientWindow, iowindow, XNFocusWindow, iowindow,
                         XNPreeditAttributes, xva_nlist, XNStatusAttributes, xva_nlist, NULL);
        XFree(xva_nlist);
        if (!ioic)
        {
            XCloseIM(ioim);
            XDestroyWindow(iodisplay, iowindow);
            XCloseDisplay(iodisplay);
            warning("unable to open X Input Context");
            return TRUE;
        }

        /* get XIM processes event. */
        XGetICValues(ioic, XNFilterEvents, &fevent, NULL);
    }
#endif

#ifdef USE_FONTSET
    if (!mbcslocale)
#endif
        XSetFont(iodisplay, iogc, font_info->fid);
    XSetBackground(iodisplay, iogc, iowhite);
    XSetForeground(iodisplay, iogc, BlackPixel(iodisplay, DefaultScreen(iodisplay)));
    XSetLineAttributes(iodisplay, iogc, 1, LineSolid, CapRound, JoinRound);

    /*
    XSelectInput(iodisplay, iowindow,
         ButtonPressMask | KeyPressMask
         | ExposureMask | StructureNotifyMask | fevent);
    */

    XSelectInput(iodisplay, iowindow,
                 ButtonPressMask | KeyPressMask | StructureNotifyMask | ExposureMask | EnterWindowMask |
                     LeaveWindowMask | fevent);
    XMapRaised(iodisplay, iowindow);

    /* now set up the menu-window, for now use the same text
       dimensions as above */

    /* font size consideration */
    for (i = 0; i < (sizeof(menu_label) / sizeof(char *)); i++)
        twidth = (twidth < textwidth(menu_label[i], strlen(menu_label[i])))
                     ? textwidth(menu_label[i], strlen(menu_label[i]))
                     : twidth;

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
    CellModified = FALSE;
    return FALSE; /* success */
}

/* MAC/X11 BASICS */

static void bell(void)
{
    XBell(iodisplay, 20);
}

static void cleararea(int xpos, int ypos, int width, int height)
{
    XClearArea(iodisplay, iowindow, xpos, ypos, width, height, 0);
}

static void clearwindow(void)
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
#ifdef USE_FONTSET
    if (mbcslocale)
#ifdef HAVE_XUTF8DRAWIMAGESTRING
        if (utf8locale)
            Xutf8DrawImageString(iodisplay, iowindow, font_set, iogc, xpos, ypos, text, len);
        else
#endif
            XmbDrawImageString(iodisplay, iowindow, font_set, iogc, xpos, ypos, text, len);
    else
#endif
        XDrawImageString(iodisplay, iowindow, iogc, xpos, ypos, text, len);
    Rsync();
}

static void Rsync()
{
    XSync(iodisplay, 0);
}

static int textwidth(char *text, int nchar)
{

#ifdef USE_FONTSET
    if (mbcslocale)
    {
#ifdef HAVE_XUTF8TEXTESCAPEMENT
        if (utf8locale)
            return Xutf8TextEscapement(font_set, text, nchar);
        else
#endif
            return XmbTextEscapement(font_set, text, nchar);
    }
#endif
    return XTextWidth(font_info, text, nchar);
}

/* Menus */

void popupmenu(int x_pos, int y_pos, int col, int row)
{
    int i, button, popupcol = col + colmin - 1;
    char *name, clab[20];
    XEvent event;
    Window selected_pane;
    SEXP tvec;

    XMoveWindow(iodisplay, menuwindow, x_pos, y_pos);
    XMapSubwindows(iodisplay, menuwindow);
    XMapRaised(iodisplay, menuwindow);

    /* now fill in the menu panes with the correct information */

    if (popupcol > xmaxused)
    {
        /* extend work, names and lens */
        REPROTECT(work = lengthgets(work, popupcol), wpi);
        REPROTECT(names = lengthgets(names, popupcol), npi);
        for (i = xmaxused + 1; i < popupcol; i++)
        {
            sprintf(clab, "var%d", i + 1);
            SET_STRING_ELT(names, i, mkChar(clab));
        }
        REPROTECT(lens = lengthgets(lens, popupcol), lpi);
        xmaxused = popupcol;
    }
    tvec = VECTOR_ELT(work, popupcol - 1);
    name = CHAR(STRING_ELT(names, popupcol - 1));
#ifdef USE_FONTSET
    if (mbcslocale)
#ifdef HAVE_XUTF8DRAWSTRING
        if (utf8locale)
            Xutf8DrawString(iodisplay, menupanes[0], font_set, iogc, 3, box_h - 3, name, strlen(name));
        else
#endif
            XmbDrawString(iodisplay, menupanes[0], font_set, iogc, 3, box_h - 3, name, strlen(name));
    else
#endif
        XDrawString(iodisplay, menupanes[0], iogc, 3, box_h - 3, name, strlen(name));
    for (i = 1; i < 4; i++)
#ifdef USE_FONTSET
        if (mbcslocale)
#ifdef HAVE_XUTF8DRAWSTRING
            if (utf8locale)
                Xutf8DrawString(iodisplay, menupanes[i], font_set, iogc, 3, box_h - 3, menu_label[i - 1],
                                strlen(menu_label[i - 1]));
            else
#endif
                XmbDrawString(iodisplay, menupanes[i], font_set, iogc, 3, box_h - 3, menu_label[i - 1],
                              strlen(menu_label[i - 1]));
        else
#endif
            XDrawString(iodisplay, menupanes[i], iogc, 3, box_h - 3, menu_label[i - 1], strlen(menu_label[i - 1]));

    if (isNull(tvec) || TYPEOF(tvec) == REALSXP)
#ifdef USE_FONTSET
        if (mbcslocale)
#ifdef HAVE_XUTF8DRAWSTRING
            if (utf8locale)
                Xutf8DrawString(iodisplay, menupanes[1], font_set, iogc, 0, box_h - 3, "*", 1);
            else
#endif
                XmbDrawString(iodisplay, menupanes[1], font_set, iogc, 0, box_h - 3, "*", 1);
        else
#endif
            XDrawString(iodisplay, menupanes[1], iogc, 0, box_h - 3, "*", 1);
    else
#ifdef USE_FONTSET
        if (mbcslocale)
#ifdef HAVE_XUTF8DRAWSTRING
        if (utf8locale)
            Xutf8DrawString(iodisplay, menupanes[2], font_set, iogc, 0, box_h - 3, "*", 1);
        else
#endif
            XmbDrawString(iodisplay, menupanes[2], font_set, iogc, 0, box_h - 3, "*", 1);
    else
#endif
        XDrawString(iodisplay, menupanes[2], iogc, 0, box_h - 3, "*", 1);

    /*
      start an event loop; we're looking for a button press and a button
      release in the same window
    */

    while (1)
    {
        XNextEvent(iodisplay, &event);

        /* event is processed with input method */

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
                    if (isNull(tvec))
                        SET_VECTOR_ELT(work, popupcol - 1, ssNewVector(REALSXP, 100));
                    else
                        SET_VECTOR_ELT(work, popupcol - 1, coerceVector(tvec, REALSXP));
                    goto done;
                case 2:
                    if (isNull(tvec))
                        SET_VECTOR_ELT(work, popupcol - 1, ssNewVector(STRSXP, 100));
                    else
                        SET_VECTOR_ELT(work, popupcol - 1, coerceVector(tvec, STRSXP));
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

void popdownmenu(void)
{
    XUnmapWindow(iodisplay, menuwindow);
    XUnmapSubwindows(iodisplay, menuwindow);
}

static void copycell(void)
{
    /*
     * whichrow = crow + colmin - 1 => whichrow = crow + rowmin - 1
     *                   ^^^                             ^^^
     */
    int i, whichrow = crow + rowmin - 1, whichcol = ccol + colmin - 1;
    SEXP tmp;

    if (whichrow == 0)
    {
        /* won't have  cell here */
    }
    else
    {
        strcpy(copycontents, "");
        if (length(work) >= whichcol)
        {
            tmp = VECTOR_ELT(work, whichcol - 1);
            if (tmp != R_NilValue && (i = whichrow - 1) < LENGTH(tmp))
            {
                PrintDefaults(R_NilValue);
                if (TYPEOF(tmp) == REALSXP)
                {
                    if (REAL(tmp)[i] != ssNA_REAL)
                    {
                        strncpy(copycontents, EncodeElement(tmp, i, 0, '.'), BOOSTED_BUF_SIZE - 1);
                        copycontents[BOOSTED_BUF_SIZE - 1] = '\0';
                    }
                }
                else if (TYPEOF(tmp) == STRSXP)
                {
                    if (!streql(CHAR(STRING_ELT(tmp, i)), CHAR(STRING_ELT(ssNA_STRING, 0))))
                    {
                        strncpy(copycontents, EncodeElement(tmp, i, 0, '.'), BOOSTED_BUF_SIZE - 1);
                        copycontents[BOOSTED_BUF_SIZE - 1] = '\0';
                    }
                }
            }
        }
    }
    highlightrect();
}

static void pastecell(int row, int col)
{
    downlightrect();
    crow = row;
    ccol = col;
    if (strlen(copycontents))
    {
        strcpy(buf, copycontents);
        clength = strlen(copycontents);
        bufp = buf + clength;
        CellModified = TRUE;
    }
    closerect();
    highlightrect();
}

#ifdef USE_FONTSET
static void calc_pre_edit_pos(void)
{
    XVaNestedList xva_nlist;
    XPoint xpoint;
    int i;
    int w;

    xpoint.x = boxw[0];
    for (i = 1; i < ccol; i++)
        xpoint.x += BOXW(colmin + i - 1);
#ifdef HAVE_XUTF8TEXTESCAPEMENT
    if (utf8locale)
        w = Xutf8TextEscapement(font_set, buf, clength);
    else
#endif
        w = XmbTextEscapement(font_set, buf, clength);
    xpoint.x += (w > BOXW(colmin + ccol - 1)) ? BOXW(colmin + ccol - 1) : w;
    xpoint.x += text_offset;
    xpoint.y = hwidth + (crow + 1) * box_h - text_offset;

    /*
      <FIXME>
      I may depend on implementation of XIM, but I do not obey,
      setting value, and investigation in various implementation
      system is need.
      It is only a problem in an appearance.
    */
    xva_nlist = XVaCreateNestedList(0, XNSpotLocation, &xpoint, XNFontSet, font_set, NULL);
    XSetICValues(ioic, XNPreeditAttributes, xva_nlist, NULL);

    XFree(xva_nlist);
    return;
}
#endif

/* last character bytes */
static int last_wchar_bytes(char *str)
{
#ifdef USE_FONTSET
    wchar_t wcs[BOOSTED_BUF_SIZE];
    mbstate_t mb_st;
    int cnt;
    char last_mbs[8];
    char *mbs;

    mbs = (str == NULL) ? buf : str;

    memset(wcs, 0, sizeof(wcs));
    memset(&mb_st, 0, sizeof(mbstate_t));

    if ((size_t)-1 == (cnt = mbsrtowcs(wcs, (const char **)&mbs, strlen(mbs), &mb_st)))
    {
        return 0;
    }
    if (wcs[0] == L'\0')
        return 0;

    memset(last_mbs, 0, sizeof(last_mbs));
    wcrtomb(last_mbs, wcs[cnt - 1], &mb_st);
    return (strlen(last_mbs));
#else
    return (1);
#endif
}
