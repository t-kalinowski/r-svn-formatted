/*
 *  R : A Computer Language for Statistical Data Analysis
 *  file shext.c
 *  Copyright (C) 2001  Guido Masarotto and Brian Ripley
 *                2004  R Development Core Team
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

/* 27/03/2000 win32-api needs this for ANSI compliance */
#define NONAMELESSUNION

#include <windows.h>
#include <shlobj.h>

int CALLBACK InitBrowseCallbackProc(HWND hwnd, UINT uMsg, LPARAM lParam, LPARAM lpData)
{
    if (uMsg == BFFM_INITIALIZED)
        SendMessage(hwnd, BFFM_SETSELECTION, 1, lpData);
    return (0);
}

/* browse for a folder under the Desktop, return the path in the argument */

void selectfolder(char *folder)
{
    char buf[MAX_PATH];
    LPMALLOC g_pMalloc;
    HWND hwnd = 0;
    BROWSEINFO bi;
    LPITEMIDLIST pidlBrowse;

    /* Get the shell's allocator. */
    if (!SUCCEEDED(SHGetMalloc(&g_pMalloc)))
        return;

    bi.hwndOwner = hwnd;
    bi.pidlRoot = NULL;
    bi.pszDisplayName = buf;
    bi.lpszTitle = "Choose a directory";
    bi.ulFlags = BIF_RETURNONLYFSDIRS;
    bi.lpfn = (BFFCALLBACK)InitBrowseCallbackProc;
    bi.lParam = (int)folder;

    /* Browse for a folder and return its PIDL. */
    pidlBrowse = SHBrowseForFolder(&bi);
    if (pidlBrowse != NULL)
    {
        SHGetPathFromIDList(pidlBrowse, folder);
        g_pMalloc->lpVtbl->Free(g_pMalloc, pidlBrowse);
    }
}

int ShellGetPersonalDirectory(char *folder) /* Folder is assumed to be at least MAX_PATH long */
{
    LPMALLOC g_pMalloc;
    LPITEMIDLIST pidlUser;
    int result;

    result = 0;

    /* Get the shell's allocator. */
    if (SUCCEEDED(SHGetMalloc(&g_pMalloc)))
    {

        /* Get the PIDL of the user's Directory. */
        if (SUCCEEDED(SHGetSpecialFolderLocation(0, CSIDL_PERSONAL, &pidlUser)))
        {
            if (SUCCEEDED(SHGetPathFromIDList(pidlUser, folder)))
                result = 1;
            g_pMalloc->lpVtbl->Free(g_pMalloc, pidlUser);
        }
    }
    return (result);
}
