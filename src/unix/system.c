/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997-2000   Robert Gentleman, Ross Ihaka
 *                            and the R Development Core Team
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

/* See system.txt for a description of functions */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Defn.h"
#include "Fileio.h"
#include "Graphics.h" /* KillAllDevices() [nothing else?] */

#include "devX11.h" /* X11ConnectionNumber, etc */

#include "Startup.h"
#include "Runix.h"

/* necessary for some (older, i.e., ~ <= 1997) Linuxen, and apparently
   also some AIX systems.
   */
#ifndef FD_SET
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#endif

#include <unistd.h> /* isatty() */

void fpu_setup(int); /* in sys-unix.c */

int UsingReadline = 1;
int SaveAction = SA_SAVEASK;
int RestoreAction = SA_RESTORE;
int LoadSiteFile = True;
int LoadInitFile = True;
int DebugInitFile = False;

/* call pointers to allow interface switching */

void R_Suicide(char *s)
{
    ptr_R_Suicide(s);
}
void R_ShowMessage(char *s)
{
    ptr_R_ShowMessage(s);
}
int R_ReadConsole(char *prompt, unsigned char *buf, int len, int addtohistory)
{
    return ptr_R_ReadConsole(prompt, buf, len, addtohistory);
}
void R_WriteConsole(char *buf, int len)
{
    ptr_R_WriteConsole(buf, len);
}
void R_ResetConsole()
{
    ptr_R_ResetConsole();
}
void R_FlushConsole()
{
    ptr_R_FlushConsole();
}
void R_ClearerrConsole()
{
    ptr_R_ClearerrConsole();
}
void R_Busy(int which)
{
    ptr_R_Busy(which);
}
void R_CleanUp(int saveact, int status, int runLast)
{
    ptr_R_CleanUp(saveact, status, runLast);
}
int R_ShowFiles(int nfile, char **file, char **headers, char *wtitle, int del, char *pager)
{
    return ptr_R_ShowFiles(nfile, file, headers, wtitle, del, pager);
}
int R_ChooseFile(int new, char *buf, int len)
{
    return ptr_R_ChooseFile(new, buf, len);
}
void (*gnome_start)(int ac, char **av, Rstart Rp);

void R_setStartTime();     /* in sys-unix.c */
void R_load_X11_shlib();   /* in dynload.c */
void R_load_gnome_shlib(); /* in dynload.c */

int main(int ac, char **av)
{
    int i, ioff = 1, j, value, ierr, useX11 = 1, usegnome = 0;
    char *p, msg[1024], **avv;
    structRstart rstart;
    Rstart Rp = &rstart;

    ptr_R_Suicide = Rstd_Suicide;
    ptr_R_ShowMessage = Rstd_ShowMessage;
    ptr_R_ReadConsole = Rstd_ReadConsole;
    ptr_R_WriteConsole = Rstd_WriteConsole;
    ptr_R_ResetConsole = Rstd_ResetConsole;
    ptr_R_FlushConsole = Rstd_FlushConsole;
    ptr_R_ClearerrConsole = Rstd_ClearerrConsole;
    ptr_R_Busy = Rstd_Busy;
    ptr_R_CleanUp = Rstd_CleanUp;
    ptr_R_ShowFiles = Rstd_ShowFiles;
    ptr_R_ChooseFile = Rstd_ChooseFile;

#ifdef HAVE_TIMES
    R_setStartTime();
#endif
    R_DefParams(Rp);
    R_SizeFromEnv(Rp);
    /* Store the command line arguments before they are processed
       by the R option handler. These are stored in Rp and then moved
       to the global variable CommandLineArgs in R_SetParams.
     */
    R_set_command_line_arguments(ac, av, Rp);

    /* first task is to select the GUI */
    for (i = 0, avv = av; i < ac; i++, avv++)
    {
        if (!strncmp(*avv, "--gui", 5) || !strncmp(*avv, "-g", 2))
        {
            if (!strncmp(*avv, "--gui", 5) && strlen(*avv) >= 7)
                p = &(*avv)[6];
            else
            {
                if (i + 1 < ac)
                {
                    avv++;
                    p = *avv;
                    ioff++;
                }
                else
                {
                    sprintf(msg, "WARNING: --gui or -g without value ignored");
                    R_ShowMessage(msg);
                    p = "X11";
                }
            }
            if (!strcmp(p, "none"))
                useX11 = 0;
            else if (!strcmp(p, "gnome") || !strcmp(p, "GNOME"))
                usegnome = 1;
            else if (!strcmp(p, "X11") || !strcmp(p, "x11"))
                useX11 = 1;
            else
            {
#ifdef HAVE_X11
                sprintf(msg, "WARNING: unknown gui `%s', using X11\n", p);
#else
                sprintf(msg, "WARNING: unknown gui `%s', using none\n", p);
#endif
                R_ShowMessage(msg);
            }
            /* now remove it/them */
            for (j = i; j < ac - ioff; j++)
            {
                av[j] = av[j + ioff];
            }
            ac -= ioff;
            break;
        }
    }

    GnomeDeviceDriver = stub_GnomeDeviceDriver;
    GTKDeviceDriver = stub_GTKDeviceDriver;
    X11DeviceDriver = stub_X11DeviceDriver;
    ptr_dataentry = stub_dataentry;
#ifdef HAVE_X11
    if (useX11)
    {
        if (!usegnome)
        {
            R_load_X11_shlib();
            R_GUIType = "X11";
        }
        else
        {
#ifndef HAVE_GNOME
            R_Suicide("GNOME GUI is not available in this version");
#endif
            R_load_X11_shlib();
            R_load_gnome_shlib();
            R_GUIType = "GNOME";
            gnome_start(ac, av, Rp);
            /* this will never return, but for safety */
            return 0;
        }
    }
#endif

    R_common_command_line(&ac, av, Rp);
    while (--ac)
    {
        if (**++av == '-')
        {
            if (!strcmp(*av, "--no-readline"))
            {
                UsingReadline = 0;
            }
            else
            {
                sprintf(msg, "WARNING: unknown option %s\n", *av);
                R_ShowMessage(msg);
            }
        }
        else
        {
            sprintf(msg, "ARGUMENT '%s' __ignored__\n", *av);
            R_ShowMessage(msg);
        }
    }
    R_SetParams(Rp);

    /* On Unix the console is a file; we just use stdio to write on it */

    R_Interactive = isatty(0);
    R_Consolefile = stdout;
    R_Outputfile = stdout;
    R_Sinkfile = NULL;
    if ((R_Home = R_HomeDir()) == NULL)
    {
        R_Suicide("R home directory is not defined");
    }
    /*
     *  Since users' expectations for save/no-save will differ, we decided
     *  that they should be forced to specify in the non-interactive case.
     */
    if (!R_Interactive && SaveAction != SA_SAVE && SaveAction != SA_NOSAVE)
        R_Suicide("you must specify `--save', `--no-save' or `--vanilla'");

    if ((R_HistoryFile = getenv("R_HISTFILE")) == NULL)
        R_HistoryFile = ".Rhistory";
    R_HistorySize = 512;
    if ((p = getenv("R_HISTSIZE")))
    {
        value = Decode2Long(p, &ierr);
        if (ierr != 0 || value < 0)
            REprintf("WARNING: invalid R_HISTSIZE ignored;");
        else
            R_HistorySize = value;
    }
    Rstd_read_history(R_HistoryFile);
    fpu_setup(1);

    mainloop();
    /*++++++  in ../main/main.c */
    return 0;
}

/* Declarations to keep f77 happy */

int MAIN_(int ac, char **av)
{
    return 0;
}
int MAIN__(int ac, char **av)
{
    return 0;
}
int __main(int ac, char **av)
{
    return 0;
}
