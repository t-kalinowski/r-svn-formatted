/*
 *  R : A Computer Language for Statistical Data Analysis
 *  file run.c: a simple 'reading' pipe (and a command executor)
 *  Copyright (C) 1999  Guido Masarotto  and Brian Ripley
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

#include <windows.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "run.h"

static char RunError[256] = "";

static char *expandcmd(char *cmd)
{
    char c;
    char *s, *p, *q, *f, *dest, *src;
    char fl[MAX_PATH], fn[MAX_PATH];
    int d, ext;

    if (!(s = (char *)malloc(MAX_PATH + strlen(cmd))))
    {
        strcpy(RunError, "Insufficient memory (expandcmd)");
        return NULL;
    }
    for (p = cmd; *p && isspace(*p); p++)
        ;
    for (q = p, d = 0; *q && (d || !isspace(*q)); q++)
        if (*q == '\"')
            d = d ? 0 : 1;
    if (d)
    {
        strcpy(RunError, "A \" is missing(expandcmd)");
        return NULL;
    }
    c = *q;
    *q = '\0';

    /*
     * I resort to this since SearchPath returned FOUND also
     * for file name without extension -> explicitly set
     *  extension
     */
    for (f = p, ext = 0; *f; f++)
    {
        if ((*f == '\\') || (*f == '/'))
            ext = 0;
        else if (*f == '.')
            ext = 1;
    }
    /* SearchFile doesn't like \" */
    for (dest = fl, src = p; *src; src++)
        if (*src != '\"')
            *dest++ = *src;
    *dest = '\0';
    if (ext)
    {
        /*
         * user set extension; we don't check that it is executable;
         * it might get an error after; but maybe sometimes
         * in the future every extension will be executable
         */
        d = SearchPath(NULL, fl, NULL, MAX_PATH, fn, &f);
    }
    else
    {
        int iexts = 0;
        char *exts[] = {".exe", ".com", ".cmd", ".bat", NULL};
        while (exts[iexts])
        {
            strcpy(dest, exts[iexts]);
            if ((d = SearchPath(NULL, fl, NULL, MAX_PATH, fn, &f)))
                break;
            iexts++;
        }
    }
    if (!d)
    {
        free(s);
        strncpy(RunError, p, 200);
        strcat(RunError, " not found");
        *q = c;
        return NULL;
    }
    /*
       Paranoia : on my system switching to short names is not needed
       since SearchPath already returns 'short names'. However,
       this is not documented so I prefer to be explicit.
       Problem is that we have removed \" from the executable since
       SearchPath seems dislikes them
    */
    GetShortPathName(fn, s, MAX_PATH);
    *q = c;
    strcat(s, q);
    return s;
}

static HANDLE pcreate(char *cmd, char *finput, int newconsole, int visible, int inpipe)
{
    DWORD ret;
    SECURITY_ATTRIBUTES sa;
    PROCESS_INFORMATION pi;
    STARTUPINFO si;
    HANDLE hIN = INVALID_HANDLE_VALUE, hSAVED, hTHIS;
    char *ecmd;

    sa.nLength = sizeof(sa);
    sa.lpSecurityDescriptor = NULL;
    sa.bInheritHandle = TRUE;

    if (!(ecmd = expandcmd(cmd)))
        return NULL;
    hTHIS = GetCurrentProcess();
    if (finput && finput[0])
    {
        DuplicateHandle(hTHIS, GetStdHandle(STD_INPUT_HANDLE), hTHIS, &hSAVED, 0, FALSE, DUPLICATE_SAME_ACCESS);
        CloseHandle(hTHIS);
        hIN = CreateFile(finput, GENERIC_READ, 0, &sa, OPEN_EXISTING, 0, NULL);
        if (hIN == INVALID_HANDLE_VALUE)
        {
            free(ecmd);
            strcpy(RunError, "Impossible to redirect input");
            return NULL;
        }
        SetStdHandle(STD_INPUT_HANDLE, hIN);
    }
    si.cb = sizeof(si);
    si.lpReserved = NULL;
    si.lpReserved2 = NULL;
    si.cbReserved2 = 0;
    si.lpDesktop = NULL;
    si.lpTitle = NULL;
    if ((finput && finput[0]) || inpipe)
    {
        si.dwFlags = STARTF_USESHOWWINDOW | STARTF_USESTDHANDLES;
        DuplicateHandle(hTHIS, GetStdHandle(STD_INPUT_HANDLE), hTHIS, &si.hStdInput, 0, TRUE, DUPLICATE_SAME_ACCESS);
        DuplicateHandle(hTHIS, GetStdHandle(STD_OUTPUT_HANDLE), hTHIS, &si.hStdOutput, 0, TRUE, DUPLICATE_SAME_ACCESS);
        DuplicateHandle(hTHIS, GetStdHandle(STD_ERROR_HANDLE), hTHIS, &si.hStdError, 0, TRUE, DUPLICATE_SAME_ACCESS);
    }
    else
        si.dwFlags = STARTF_USESHOWWINDOW;
    switch (visible)
    {
    case -1:
        si.wShowWindow = SW_HIDE;
        break;
    case 0:
        si.wShowWindow = SW_SHOWMINIMIZED;
        break;
    case 1:
        si.wShowWindow = SW_SHOWDEFAULT;
        break;
    }
    ret = CreateProcess(0, (char *)ecmd, &sa, &sa, TRUE, (newconsole && (visible == 1)) ? CREATE_NEW_CONSOLE : 0, NULL,
                        NULL, &si, &pi);
    if (finput && finput[0])
    {
        SetStdHandle(STD_INPUT_HANDLE, hSAVED);
        CloseHandle(hIN);
    }
    if (si.dwFlags & STARTF_USESTDHANDLES)
    {
        CloseHandle(si.hStdInput);
        CloseHandle(si.hStdOutput);
        CloseHandle(si.hStdError);
    }
    if (!ret)
    {
        strcpy(RunError, "Impossible to run '");
        strncat(RunError, ecmd, 200);
        free(ecmd);
        return NULL;
    }
    free(ecmd);
    CloseHandle(pi.hThread);
    return pi.hProcess;
}

static int pwait(HANDLE p)
{
    DWORD ret;

    WaitForSingleObject(p, INFINITE);
    GetExitCodeProcess(p, &ret);
    return ret;
}

static DWORD CALLBACK threadedwait(LPVOID param)
{
    rpipe *p = (rpipe *)param;

    p->exitcode = pwait(p->process);
    FlushFileBuffers(p->write);
    FlushFileBuffers(p->read);
    p->active = 0;
    return 0;
}

char *runerror()
{
    return RunError;
}

int runcmd(char *cmd, int wait, int visible, char *finput)
{
    HANDLE p;
    int ret;
    /* I hope no program will use this as an error code */
    if (!(p = pcreate(cmd, finput, !wait, visible, 0)))
        return NOLAUNCH;
    if (wait)
        ret = pwait(p);
    else
        ret = 0;
    CloseHandle(p);
    return ret;
}

rpipe *rpipeOpen(char *cmd, int visible, char *finput)
{
    rpipe *r;
    HANDLE hOUT, hERR, hThread, hTHIS, hTemp;
    DWORD id;

    if (!(r = (rpipe *)malloc(sizeof(struct structRPIPE))))
    {
        strcpy(RunError, "Insufficient memory (rpipeOpen)");
        return NULL;
    }
    r->process = NULL;
    if (CreatePipe(&hTemp, &(r->write), NULL, 0) == FALSE)
    {
        rpipeClose(r);
        strcpy(RunError, "Inpossible to create pipe");
        return NULL;
    }
    hTHIS = GetCurrentProcess();
    DuplicateHandle(hTHIS, GetStdHandle(STD_OUTPUT_HANDLE), hTHIS, &hOUT, 0, FALSE, DUPLICATE_SAME_ACCESS);
    DuplicateHandle(hTHIS, GetStdHandle(STD_ERROR_HANDLE), hTHIS, &hERR, 0, FALSE, DUPLICATE_SAME_ACCESS);
    DuplicateHandle(hTHIS, hTemp, hTHIS, &r->read, 0, FALSE, DUPLICATE_SAME_ACCESS);
    CloseHandle(hTemp);
    CloseHandle(hTHIS);
    SetStdHandle(STD_OUTPUT_HANDLE, r->write);
    SetStdHandle(STD_ERROR_HANDLE, r->write);
    r->process = pcreate(cmd, finput, 0, visible, 1);
    r->active = 1;
    SetStdHandle(STD_OUTPUT_HANDLE, hOUT);
    SetStdHandle(STD_ERROR_HANDLE, hERR);
    if (!r->process)
        return NULL;
    if (!(hThread = CreateThread(NULL, 0, threadedwait, r, 0, &id)))
    {
        rpipeClose(r);
        strcpy(RunError, "Inpossible to create thread/pipe");
        return NULL;
    }
    CloseHandle(hThread);
    return r;
}

int rpipeGetc(rpipe *r)
{
    DWORD a, b;
    char c;

    if (!r)
        return NOLAUNCH;
    while (PeekNamedPipe(r->read, NULL, 0, NULL, &a, NULL))
    {
        if (!a && !r->active)
        {
            /* I got a case in which process terminated after Peek.. */
            PeekNamedPipe(r->read, NULL, 0, NULL, &a, NULL);
            if (!a)
                return NOLAUNCH; /* end of pipe */
        }
        if (a)
        {
            if (ReadFile(r->read, &c, 1, &b, NULL) == TRUE)
                return c;
            else
                return NOLAUNCH; /* error but...treated as eof */
        }
    }
    return NOLAUNCH; /* again.. */
}

char *rpipeGets(rpipe *r, char *buf, int len)
{
    int i, c;

    if ((len < 2) || !r)
        return NULL;
    for (i = 0; i < (len - 1); i++)
    {
        if ((c = rpipeGetc(r)) == NOLAUNCH)
        {
            if (i == 0)
                return NULL;
            else
            {
                buf[i] = '\0';
                return buf;
            }
        }
        buf[i] = c;
        if (c == '\n')
        {
            if ((i > 0) && (buf[i - 1] == '\r'))
            {
                buf[i - 1] = '\n';
                buf[i] = '\0';
            }
            else
                buf[i + 1] = '\0';
            return buf;
        }
    }
    buf[len - 1] = '\0';
    return buf;
}

int rpipeClose(rpipe *r)
{
    int i;

    if (!r)
        return NOLAUNCH;
    if (r->active)
        TerminateProcess(r->process, 99);
    CloseHandle(r->read);
    CloseHandle(r->write);
    CloseHandle(r->process);
    i = r->exitcode;
    free(r);
    return i;
}
