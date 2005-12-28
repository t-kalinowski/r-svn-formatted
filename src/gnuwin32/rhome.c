/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999--2000  Guido Masarotto and Brian Ripley
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
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define WIN32_LEAN_AND_MEAN 1
#include <windows.h>
#include <string.h> /* for strrchr(...) */
#include <stdio.h>
#include <ctype.h>
#include <Rversion.h>
#include <stdlib.h> /* for exit */

static char rhomebuf[MAX_PATH];

/* <MBCS-FIXME> We can't just use Rf_strchr as this is called
   from front-ends */
#define GOBACKONESLASH                                                                                                 \
    p = strrchr(rhomebuf, '\\');                                                                                       \
    if (!p)                                                                                                            \
    {                                                                                                                  \
        MessageBox(NULL, "Installation problem", "Terminating", MB_TASKMODAL | MB_ICONSTOP | MB_OK);                   \
        exit(1);                                                                                                       \
    }                                                                                                                  \
    *p = '\0'

char *getRHOME()
{
    DWORD nc;
    char *p;
    int hasspace = 0;

    nc = GetModuleFileName(NULL, rhomebuf, MAX_PATH);
    GOBACKONESLASH;
    GOBACKONESLASH;
    /* make sure no spaces in path */
    for (p = rhomebuf; *p; p++)
        if (isspace(*p))
        {
            hasspace = 1;
            break;
        }
    if (hasspace)
        GetShortPathName(rhomebuf, rhomebuf, MAX_PATH);
    return (rhomebuf);
}

static char DLLversion[25];

char *getDLLVersion()
{
    sprintf(DLLversion, "%s.%s", R_MAJOR, R_MINOR);
    return (DLLversion);
}

char *get_R_HOME()
{
    LONG rc;
    HKEY hkey;
    DWORD keytype = REG_SZ, cbData = sizeof(rhomebuf);

    if (getenv("R_HOME"))
    {
        strncpy(rhomebuf, getenv("R_HOME"), MAX_PATH);
        return (rhomebuf);
    }
    rc = RegOpenKeyEx(HKEY_LOCAL_MACHINE, "Software\\R-core\\R", 0, KEY_READ, &hkey);
    if (rc == ERROR_SUCCESS)
    {
        rc = RegQueryValueEx(hkey, "InstallPath", 0, &keytype, (LPBYTE)rhomebuf, &cbData);
        RegCloseKey(hkey);
    }
    else
        return NULL;
    if (rc != ERROR_SUCCESS)
        return NULL;
    return rhomebuf;
}
