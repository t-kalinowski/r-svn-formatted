/*
 *  R : A Computer Language for Statistical Data Analysis
 *  file pager.c
 *  Copyright (C) 1998--2000  Guido Masarotto and Brian Ripley
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

#ifdef Win32
#define USE_MDI 1
#endif

#include "R_ext/Error.h" /* for warning() */
#include <windows.h>
#include "graphapp/ga.h"
#ifdef USE_MDI
#include "graphapp/stdimg.h"
#endif
#include "console.h"
#include "consolestructs.h"
#include "rui.h"
#include "Startup.h" /* for UImode */

extern UImode CharacterMode;

#define PAGERMAXKEPT 12
#define PAGERMAXTITLE 128

static int pagerActualKept = 0, pagerActualShown;
static pager pagerInstance = NULL;
static menubar pagerBar = NULL;
static xbuf pagerXbuf[PAGERMAXKEPT];
static char pagerTitles[PAGERMAXKEPT][PAGERMAXTITLE + 8];
static menuitem pagerMenus[PAGERMAXKEPT];
static int pagerRow[PAGERMAXKEPT];
static void pagerupdateview();

int pagerMultiple, haveusedapager;

/*
   To be fixed: during creation, memory is allocated two times
   (faster for small files but a big waste otherwise)
*/
static xbuf file2xbuf(char *name, int del)
{
    HANDLE f;
    DWORD rr, vv;
    char *q, *p;
    xlong dim;
    xint ms;
    xbuf xb;

    f = CreateFile(name, GENERIC_READ, FILE_SHARE_WRITE, NULL, OPEN_EXISTING, 0, NULL);
    if (f == INVALID_HANDLE_VALUE)
    {
        warning("File %s could not be opened by internal pager", name);
        return NULL;
    }
    vv = GetFileSize(f, NULL);
    p = (char *)winmalloc((size_t)vv + 1);
    if (!p)
    {
        CloseHandle(f);
        warning("Insufficient memory to display %s in internal pager", name);
        return NULL;
    }
    ReadFile(f, p, vv, &rr, NULL);
    CloseHandle(f);
    if (del)
        DeleteFile(name);
    p[rr] = '\0';
    for (q = p, ms = 1, dim = rr; *q; q++)
    {
        if (*q == '\t')
            dim += TABSIZE;
        else if (*q == '\n')
        {
            dim++;
            ms++;
        }
    }
    if ((xb = newxbuf(dim + 1, ms, 1)))
        for (q = p, ms = 0; *q; q++)
        {
            if (*q == '\n')
            {
                ms++;
                xbufaddc(xb, *q);
                /* next line interprets underlining in help files */
                if (q[1] == '_' && q[2] == '\b')
                    xb->user[ms] = -2;
            }
            else
                xbufaddc(xb, *q);
        }
    winfree(p);
    return xb;
}

static void delpager(control m)
{
    int i;

    ConsoleData p = getdata(m);
    if (!pagerMultiple)
    {
        for (i = 0; i < pagerActualKept; i++)
        {
            xbufdel(pagerXbuf[i]);
        }
        pagerActualKept = 0;
    }
    else
        xbufdel(p->lbuf);
    freeConsoleData(getdata(m));
}

static void pagerbclose(control m)
{
    show(RConsole);
    if (!pagerMultiple)
    {
        hide(pagerInstance);
        del(pagerInstance);
        pagerInstance = pagerBar = NULL;
    }
    else
    {
        hide(m);
        del(m);
    }
}

static void pagerclose(control m)
{
    pagerbclose(getdata(m));
}

static void pagerprint(control m)
{
    consoleprint(getdata(m));
}

static void pagercopy(control m)
{
    control c = getdata(m);

    if (consolecancopy(c))
        consolecopy(c);
    else
        R_ShowMessage("No selection");
}

static void pagerpaste(control m)
{
    control c = getdata(m);

    if (CharacterMode != RGui)
    {
        R_ShowMessage("No RGui console to paste to");
        return;
    }
    if (!consolecancopy(c))
    {
        R_ShowMessage("No selection");
        return;
    }
    else
    {
        consolecopy(c);
    }
    if (consolecanpaste(RConsole))
    {
        consolepaste(RConsole);
        show(RConsole);
    }
}

static void pagerselectall(control m)
{
    control c = getdata(m);

    consoleselectall(c);
}

static void pagerconsole(control m)
{
    show(RConsole);
}

static void pagerchangeview(control m)
{
    ConsoleData p = getdata(pagerInstance);
    int i = getvalue(m);

    if (i >= pagerActualKept)
        return;
    uncheck(pagerMenus[pagerActualShown]);
    /* save position of middle line of pager display */
    pagerRow[pagerActualShown] = FV + ROWS / 2;
    pagerActualShown = i;
    check(pagerMenus[i]);
    pagerupdateview();
}

static void pagerupdateview()
{
    control c = pagerInstance;
    ConsoleData p = getdata(c);

    settext(pagerInstance, &pagerTitles[pagerActualShown][4]);
    p->lbuf = pagerXbuf[pagerActualShown];
    setfirstvisible(c, pagerRow[pagerActualShown] - ROWS / 2);
    setfirstcol(c, 0);
    show(c);
}

static int pageraddfile(char *wtitle, char *filename, int deleteonexit)
{
    ConsoleData p = getdata(pagerInstance);
    int i;
    xbuf nxbuf = file2xbuf(filename, deleteonexit);

    if (!nxbuf)
    {
        /*	R_ShowMessage("File not found or memory insufficient"); */
        return 0;
    }
    if (pagerActualKept == PAGERMAXKEPT)
    {
        pagerActualKept -= 1;
        xbufdel(pagerXbuf[pagerActualKept]);
    }
    if (pagerActualKept > 0)
        pagerRow[0] = FV;
    for (i = pagerActualKept; i > 0; i--)
    {
        pagerXbuf[i] = pagerXbuf[i - 1];
        pagerRow[i] = pagerRow[i - 1];
        strcpy(&pagerTitles[i][4], &pagerTitles[i - 1][4]);
    }
    pagerXbuf[0] = nxbuf;
    pagerRow[0] = 0;
    strcpy(&pagerTitles[0][4], wtitle);
    pagerActualKept += 1;
    for (i = 0; i < pagerActualKept; i++)
    {
        enable(pagerMenus[i]);
        settext(pagerMenus[i], pagerTitles[i]);
    }
    for (i = pagerActualKept; i < PAGERMAXKEPT; i++)
        disable(pagerMenus[i]);
    uncheck(pagerMenus[pagerActualShown]);
    pagerActualShown = 0;
    check(pagerMenus[pagerActualShown]);
    return 1;
}

static MenuItem PagerPopup[] = {{"Copy", pagercopy, 0},
                                {"Paste to console", pagerpaste, 0},
                                {"Select all", pagerselectall, 0},
                                {"-", 0, 0},
                                {"Close", pagerclose, 0},
                                LASTMENUITEM};

static void pagermenuact(control m)
{
    control c = getdata(m);
    ConsoleData p = getdata(c);
    if (consolecancopy(c))
    {
        enable(p->mcopy);
        enable(p->mpopcopy);
        if (CharacterMode == RGui)
        {
            enable(p->mpaste);
            enable(p->mpoppaste);
        }
    }
    else
    {
        disable(p->mcopy);
        disable(p->mpopcopy);
        disable(p->mpaste);
        disable(p->mpoppaste);
    }
}

RECT *RgetMDIsize(); /* in rui.c */
#define MCHECK(a)                                                                                                      \
    if (!(a))                                                                                                          \
    {                                                                                                                  \
        freeConsoleData(p);                                                                                            \
        del(c);                                                                                                        \
        return NULL;                                                                                                   \
    }
static pager pagercreate()
{
    ConsoleData p;
    int w, h, i, x, y, w0, h0;
    pager c;
    menuitem m;

    p = newconsoledata((consolefn) ? consolefn : FixedFont, pagerrow, pagercol, consolefg, consoleuser, consolebg,
                       PAGER);
    if (!p)
        return NULL;

    /*    if (ismdi()) {
        x = y = w = h = 0;
        }
        else {
        w = WIDTH ;
        h = HEIGHT;
        x = (devicewidth(NULL) - w) / 2;
        y = (deviceheight(NULL) - h) / 2 ;
        } */
    w = WIDTH;
    h = HEIGHT;
    /* centre a single pager, randomly place each of multiple pagers */
#ifdef USE_MDI
    if (ismdi())
    {
        RECT *pR = RgetMDIsize();
        w0 = pR->right;
        h0 = pR->bottom;
    }
    else
    {
#endif
        w0 = devicewidth(NULL);
        h0 = deviceheight(NULL);
#ifdef USE_MDI
    }
#endif
    x = (w0 - w) / 2;
    x = x > 20 ? x : 20;
    y = (h0 - h) / 2;
    y = y > 20 ? y : 20;
    if (pagerMultiple)
    {
#ifdef Win32
        DWORD rand = GetTickCount();
#else
        int rand = 0;
#endif
        int w0 = 0.4 * x, h0 = 0.4 * y;
        w0 = w0 > 20 ? w0 : 20;
        h0 = h0 > 20 ? h0 : 20;
        x += (rand % w0) - w0 / 2;
        y += ((rand / w0) % h0) - h0 / 2;
    }
    c = (pager)newwindow("PAGER", rect(x, y, w, h),
                         Document | StandardWindow | Menubar | VScrollbar | HScrollbar | TrackMouse);
    if (!c)
    {
        freeConsoleData(p);
        return NULL;
    }
    setdata(c, p);
    if (h == 0)
        HEIGHT = getheight(c);
    if (w == 0)
        WIDTH = getwidth(c);
    COLS = WIDTH / FW - 1;
    ROWS = HEIGHT / FH - 1;
    BORDERX = (WIDTH - COLS * FW) / 2;
    BORDERY = (HEIGHT - ROWS * FH) / 2;
    gsetcursor(c, ArrowCursor);
    gchangescrollbar(c, VWINSB, 0, 0, ROWS, 0);
    gchangescrollbar(c, HWINSB, 0, COLS - 1, COLS, 1);
    setbackground(c, consolebg);
#ifdef USE_MDI
    if (ismdi() && (RguiMDI & RW_TOOLBAR))
    {
        int btsize = 24;
        rect r = rect(2, 2, btsize, btsize);
        control tb, bt;
        addto(c);
        MCHECK(tb = newtoolbar(btsize + 4));
        gsetcursor(tb, ArrowCursor);
        addto(tb);
        MCHECK(bt = newtoolbutton(copy1_image, r, pagerpaste));
        MCHECK(addtooltip(bt, "Paste to console"));
        gsetcursor(bt, ArrowCursor);
        setdata(bt, (void *)c);
        r.x += (btsize + 6);
        MCHECK(bt = newtoolbutton(print_image, r, pagerprint));
        MCHECK(addtooltip(bt, "Print"));
        gsetcursor(bt, ArrowCursor);
        setdata(bt, (void *)c);
        r.x += (btsize + 6);
        MCHECK(bt = newtoolbutton(console_image, r, pagerconsole));
        MCHECK(addtooltip(bt, "Return focus to Console"));
        gsetcursor(bt, ArrowCursor);
    }
#endif
    addto(c);
    MCHECK(m = gpopup(pagermenuact, PagerPopup));
    setdata(m, c);
    setdata(p->mpopcopy = PagerPopup[0].m, c);
    setdata(p->mpoppaste = PagerPopup[1].m, c);
    setdata(PagerPopup[2].m, c);
    setdata(PagerPopup[4].m, c);
    MCHECK(m = newmenubar(pagermenuact));
    setdata(m, c);
    MCHECK(newmenu("File"));
    MCHECK(m = newmenuitem("Print", 0, pagerprint));
    setdata(m, c);
    MCHECK(m = newmenuitem("-", 0, NULL));
    MCHECK(m = newmenuitem("Close", 0, pagerclose));
    setdata(m, c);
    MCHECK(newmenu("Edit"));
    MCHECK(p->mcopy = newmenuitem("Copy", 'C', pagercopy));
    setdata(p->mcopy, c);
    MCHECK(p->mpaste = newmenuitem("Paste to console", 'V', pagerpaste));
    setdata(p->mpaste, c);
    MCHECK(m = newmenuitem("Select all", 0, pagerselectall));
    setdata(m, c);
    if (!pagerMultiple)
    {
        MCHECK(newmenu("View"));
        for (i = 0; i < PAGERMAXKEPT; i++)
        {
            sprintf(pagerTitles[i], "&%c.  ", 'A' + i);
            MCHECK(pagerMenus[i] = newmenuitem(&pagerTitles[i][1], 0, pagerchangeview));
            setvalue(pagerMenus[i], i);
        }
    }
#ifdef USE_MDI
    if (ismdi())
        newmdimenu();
#endif
    MCHECK(BM = newbitmap(WIDTH, HEIGHT, 2));
    setdata(c, p);
    sethit(c, console_sbf);
    setresize(c, consoleresize);
    setredraw(c, drawconsole);
    setdel(c, delpager);
    setclose(c, pagerbclose);
    setkeyaction(c, console_ctrlkeyin);
    setkeydown(c, console_normalkeyin);
    setmousedrag(c, console_mousedrag);
    setmouserepeat(c, console_mouserep);
    setmousedown(c, console_mousedown);
    return (c);
}

static pager newpager1win(char *wtitle, char *filename, int deleteonexit)
{
    if (!pagerInstance && !(pagerInstance = pagercreate()))
    {
        R_ShowMessage("Unable to create pager windows");
        return NULL;
    }
    if (!pageraddfile(wtitle, filename, deleteonexit))
        return NULL;
    pagerupdateview();
    return pagerInstance;
}

static pager newpagerNwin(char *wtitle, char *filename, int deleteonexit)
{
    pager c = pagercreate();
    ConsoleData p;

    if (!c)
        return NULL;
    settext(c, wtitle);
    p = getdata(c);
    if (!(p->lbuf = file2xbuf(filename, deleteonexit)))
    {
        del(c);
        return NULL;
    }
    if (c)
    {
        gchangescrollbar(c, VWINSB, 0, NUMLINES - 1, ROWS, 0);
        show(c);
    }
    return c;
}

pager newpager(char *title, char *filename, char *header, int deleteonexit)
{
    char wtitle[PAGERMAXTITLE + 1];
    pager c;

    /*    if (ismdi()) pagerMultiple = 1;*/
    strncpy(wtitle, title, PAGERMAXTITLE);
    wtitle[PAGERMAXTITLE] = '\0';
    if (strlen(header) && ((strlen(header) + strlen(wtitle) + 4) < PAGERMAXTITLE))
    {
        if (strlen(wtitle))
            strcat(wtitle, " - ");
        strcat(wtitle, header);
    }
    if (!pagerMultiple)
        c = newpager1win(wtitle, filename, deleteonexit);
    else
        c = newpagerNwin(wtitle, filename, deleteonexit);
    haveusedapager++;
    BringToTop(c);
    return c;
}
