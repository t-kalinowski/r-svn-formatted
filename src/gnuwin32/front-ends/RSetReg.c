/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001  R Development Core Team
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
#define NONAMELESSUNION
#include <windows.h>
#include <stdio.h>
#include <Rversion.h>

#define REG_KEY_NAME "Software\\R-core\\R"

extern char *getRHOME(); /* in ../rhome.c */

int main(int argc, char **argv)
{
    int status = 0;

    char *RHome, version[40];
    LONG rc;
    HKEY hkey;

    if (argc > 1)
    { /* remove the keys */
        printf("unregistering R ... ");
        if ((rc = RegOpenKeyEx(HKEY_LOCAL_MACHINE, REG_KEY_NAME, 0, KEY_SET_VALUE, &hkey)) == ERROR_SUCCESS)
        {
            RegDeleteKey(hkey, "InstallPath");
            RegDeleteKey(hkey, "Current Version");
            RegCloseKey(hkey);
            if ((rc = RegOpenKeyEx(HKEY_LOCAL_MACHINE, "Software\\R-core", 0, KEY_SET_VALUE, &hkey)) == ERROR_SUCCESS)
            {
                RegDeleteKey(hkey, "R");
                RegCloseKey(hkey);
                if ((rc = RegOpenKeyEx(HKEY_LOCAL_MACHINE, "Software", 0, KEY_SET_VALUE, &hkey)) == ERROR_SUCCESS)
                {
                    RegDeleteKey(hkey, "R-core");
                    RegCloseKey(hkey);
                }
            }
            printf("succeeded\n");
        }
        else
        {
            printf("was not registered\n");
            status = 1;
        }
    }
    else
    {
        RHome = getRHOME();
        sprintf(version, "%s.%s %s", R_MAJOR, R_MINOR, R_STATUS);
        if ((rc = RegOpenKeyEx(HKEY_LOCAL_MACHINE, REG_KEY_NAME, 0, KEY_ALL_ACCESS, &hkey)) != ERROR_SUCCESS)
        {
            /* failed to open key, so try to create it */
            rc = RegCreateKey(HKEY_LOCAL_MACHINE, REG_KEY_NAME, &hkey);
        }
        if (rc == ERROR_SUCCESS)
        {
            rc = RegSetValueEx(hkey, "InstallPath", 0, REG_SZ, (CONST BYTE *)RHome, lstrlen(RHome) + 1);
            if (rc == ERROR_SUCCESS)
                rc = RegSetValueEx(hkey, "Current Version", 0, REG_SZ, (CONST BYTE *)version, lstrlen(version) + 1);
            RegCloseKey(hkey);
        }
        else
        {
            status = 1;
        }
    }
    exit(status);
}
