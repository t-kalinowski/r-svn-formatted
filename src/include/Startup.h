/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999  The R Development Core Team
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

#ifndef True
#define True 1
#define False 0
#endif

#ifdef Win32
typedef int (*blah1) (char *, char *, int, int);
typedef void (*blah2) (char *, int);
typedef void (*blah3) ();
typedef void (*blah4) (char *);
typedef int (*blah5) (char *);
typedef void (*blah6) (int);
typedef enum {RGui, RTerm, LinkDLL} UImode;
#endif

typedef struct
{
    int R_Quiet;               /* > 0 to suppress messages */
    int R_Slave;
    int R_Interactive;
    int R_Verbose;
    int RestoreAction;
    int SaveAction;
    int LoadSiteFile;
    int LoadInitFile;
    int DebugInitFile;
    int vsize;
    int nsize;
#ifdef Win32
    char *rhome;               /* R_HOME */
    char *home;                /* HOME  */
    blah1 ReadConsole;
    blah2 WriteConsole;
    blah3 CallBack;
    blah4 message;
    blah5 yesnocancel;
    blah6 busy;
    int NoRenviron;
    UImode CharacterMode;
#endif
} structRstart;

typedef structRstart *Rstart;

void R_DefParams(Rstart);
void R_SetParams(Rstart);
void R_SetWin32(Rstart);
void R_SizeFromEnv(Rstart);
void R_common_command_line(int *, char **, Rstart);

