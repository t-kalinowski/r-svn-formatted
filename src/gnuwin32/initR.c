/*  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998--2006 R Development Core Team
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

#define Win32
#define WIN32_LEAN_AND_MEAN 1
#include <windows.h>
#include <stdio.h>
#include <Rversion.h>
#include <R_ext/RStartup.h>
/* for askok and askyesnocancel */
#include "graphapp/graphapp.h"

/* for signal-handling code */
#include <psignal.h>

/* one way to allow user interrupts: called in ProcessEvents */
extern int UserBreak;

/* calls into the R DLL */
extern char *getDLLVersion(), *getRUser(), *get_R_HOME();
extern void R_DefParams(Rstart), R_SetParams(Rstart), R_setStartTime();
extern void setup_term_ui(void), ProcessEvents(void);
extern int R_ReplDLLdo1();

/* simple input, simple output */

/* This version blocks all events: a real one needs to call ProcessEvents
   frequently. See rterm.c and ../system.c for one approach using
   a separate thread for input.
*/
int myReadConsole(char *prompt, char *buf, int len, int addtohistory)
{
    fputs(prompt, stdout);
    fflush(stdout);
    if (fgets(buf, len, stdin))
        return 1;
    else
        return 0;
}

void myWriteConsole(char *buf, int len)
{
    printf("%s", buf);
}

void myCallBack()
{
    /* called during i/o, eval, graphics in ProcessEvents */
}

void myBusy(int which)
{
    /* set a busy cursor ... if which = 1, unset if which = 0 */
}

static void my_onintr(int sig)
{
    UserBreak = 1;
}

int Rf_initEmbeddedR(int argc, char **argv)
{
    structRstart rp;
    Rstart Rp = &rp;
    char Rversion[25], *RHome;

    sprintf(Rversion, "%s.%s", R_MAJOR, R_MINOR);
    if (strcmp(getDLLVersion(), Rversion) != 0)
    {
        fprintf(stderr, "Error: R.DLL version does not match\n");
        exit(1);
    }

    R_setStartTime();
    R_DefParams(Rp);
    if ((RHome = get_R_HOME()) == NULL)
    {
        fprintf(stderr, "R_HOME must be set in the environment or Registry\n");
        exit(2);
    }
    Rp->rhome = RHome;
    Rp->home = getRUser();
    Rp->CharacterMode = LinkDLL;
    Rp->ReadConsole = myReadConsole;
    Rp->WriteConsole = myWriteConsole;
    Rp->CallBack = myCallBack;
    Rp->ShowMessage = askok;
    Rp->YesNoCancel = askyesnocancel;
    Rp->Busy = myBusy;

    Rp->R_Quiet = TRUE;
    Rp->R_Interactive = TRUE;
    Rp->RestoreAction = SA_RESTORE;
    Rp->SaveAction = SA_NOSAVE;
    R_SetParams(Rp); /* so R_ShowMessage is set */
    R_SizeFromEnv(Rp);
    R_SetParams(Rp);
    R_set_command_line_arguments(argc, argv);

    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));

    signal(SIGBREAK, my_onintr);
    setup_term_ui(); /* initialize graphapp, eventloop, read Rconsole */
    setup_Rmainloop();

    return 0;
}
