/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999  Guido Masarotto
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

/* Support for printer
 *  printer newprinter()  - return a printer object - to draw to the
 *                          printer use drawto(...) and the drawXXX
 *                          functions + nextpage(printer)
 *                          printers can be deleted by 'del(printer)'
 */

#include "internal.h"

/*
 *  Internal printer deletion function.
 */
static void private_delprinter(printer obj)
{
    HDC h = (HDC)obj->handle;
    if (!obj || !h || (obj->kind != PrinterObject))
        return;
    EndPage(h);
    EndDoc(h);
    DeleteDC(h);
    return;
}

/*
 *  Create/return the base printer object.
 */
static object get_printer_base(void)
{
    static object printer_base = NULL;

    if (!printer_base)
        printer_base = new_object(BaseObject, 0, NULL);
    return printer_base;
}

static HDC chooseprinter()
{
    PRINTDLG pd;
    HDC dc;
    char cwd[MAX_PATH];

    GetCurrentDirectory(MAX_PATH, cwd);

    pd.lStructSize = sizeof(PRINTDLG);
    pd.hwndOwner = NULL;
    pd.hDevMode = (HANDLE)NULL;
    pd.hDevNames = (HANDLE)NULL;
    pd.Flags = PD_RETURNDC | PD_NOSELECTION | PD_NOPAGENUMS | PD_USEDEVMODECOPIES;
    pd.nFromPage = 0;
    pd.nToPage = 0;
    pd.nMinPage = 0;
    pd.nMaxPage = 0;
    pd.nCopies = 1;
    pd.hInstance = (HINSTANCE)NULL;
    pd.lCustData = (LPARAM)0;
    pd.lpfnPrintHook = 0;
    pd.lpfnSetupHook = 0;
    pd.lpPrintTemplateName = (LPCSTR)0;
    pd.lpSetupTemplateName = (LPCSTR)0;
    pd.hPrintTemplate = (HGLOBAL)0;
    pd.hSetupTemplate = (HGLOBAL)0;

    dc = PrintDlg(&pd) ? pd.hDC : NULL;
    SetCurrentDirectory(cwd);

    return dc;
}

printer newprinter(double width, double height)
{
    DOCINFO docinfo;
    printer obj;
    HDC hDC = chooseprinter();
    double dd, AL;
    int ww, hh, x0, y0;
    if (!hDC)
        return NULL;
    obj = new_object(PrinterObject, (HANDLE)hDC, get_printer_base());
    if (!obj)
    {
        askok("Insufficient memory for new printer");
        DeleteDC(hDC);
        return NULL;
    }
    if ((width == 0.0) && (height == 0.0))
    {
        ww = GetDeviceCaps(hDC, HORZRES);
        hh = GetDeviceCaps(hDC, VERTRES);
    }
    else
    {
        if (width < 0.1)
            width = 0.1;
        if (height < 0.1)
            height = 0.1;
        dd = GetDeviceCaps(hDC, HORZSIZE) / width;
        AL = (dd < 1.0) ? dd : 1.0;
        dd = GetDeviceCaps(hDC, VERTSIZE) / height;
        AL = (dd < AL) ? dd : AL;
        ww = (AL * width) * GetDeviceCaps(hDC, LOGPIXELSX) / 25.4;
        hh = (AL * height) * GetDeviceCaps(hDC, LOGPIXELSY) / 25.4;
    }
    x0 = (GetDeviceCaps(hDC, HORZRES) - ww) / 2;
    y0 = (GetDeviceCaps(hDC, VERTRES) - hh) / 2;
    obj->rect = rect(x0, y0, ww, hh);
    obj->depth = GetDeviceCaps(hDC, BITSPIXEL) * GetDeviceCaps(hDC, PLANES);
    obj->die = private_delprinter;
    obj->drawstate = copydrawstate();
    obj->drawstate->dest = obj;

    docinfo.cbSize = sizeof(DOCINFO); /* set this size... */
    docinfo.lpszDocName = "GraphAppPrintJob";
    docinfo.lpszOutput = 0; /* no file output... */
    docinfo.lpszDatatype = 0;
    docinfo.fwType = 0;

    if (StartDoc(hDC, &docinfo) <= 0)
    {
        askok("Impossible to start the print job");
        del(obj);
        return NULL;
    }

    StartPage(hDC);

    return obj;
}

void nextpage(printer p)
{
    if (!p || (p->kind != PrinterObject))
        return;
    EndPage((HDC)p->handle);
    StartPage((HDC)p->handle);
}
