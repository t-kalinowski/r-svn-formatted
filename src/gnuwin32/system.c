/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997--2004  Robert Gentleman, Ross Ihaka and the
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

/* See ../unix/system.txt for a description of functions */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Defn.h"
#include "Fileio.h"
#include "Rdevices.h" /* KillAllDevices() [nothing else?] */
#include "graphapp/ga.h"
#include "console.h"
#include "rui.h"
#include "editor.h"
#include "getline/getline.h"
#include <windows.h> /* for CreateEvent,.. */
#include <process.h> /* for _beginthread,... */
#include <io.h>      /* for isatty, chdir */
#include "run.h"
#include "Startup.h"
#include <stdlib.h>    /* for exit */
#include "shext.h"     /* for ShellGetPersonalDirectory */
void CleanTempDir();   /* from extra.c */
void editorcleanall(); /* from editor.c */

R_size_t R_max_memory = INT_MAX;
Rboolean UseInternet2 = FALSE;

extern SA_TYPE SaveAction;      /* from ../main/startup.c */
Rboolean DebugMenuitem = FALSE; /* exported for rui.c */

__declspec(dllexport) UImode CharacterMode;
int ConsoleAcceptCmd;
void closeAllHlpFiles();
void UnLoad_Unzip_Dll();
void UnLoad_Rbitmap_Dll();
void set_workspace_name(char *fn); /* ../unix/sys-common.c */

/* used to avoid some flashing during cleaning up */
Rboolean AllDevicesKilled = FALSE;
int setupui(void);
void delui(void);
int (*R_YesNoCancel)(char *s);

static DWORD mainThreadId;

static char oldtitle[512];

__declspec(dllexport) Rboolean UserBreak = FALSE;

/* callbacks */
static void (*R_CallBackHook)();
static void R_DoNothing()
{
}
static void (*my_R_Busy)(int);

/*
 *   Called at I/O, during eval etc to process GUI events.
 */

void (*R_tcldo)();
static void tcl_do_none()
{
}

void R_ProcessEvents(void)
{
    while (peekevent())
        doevent();
    if (UserBreak)
    {
        UserBreak = FALSE;
        onintr();
    }
    R_CallBackHook();
    if (R_tcldo)
        R_tcldo();
}

/*
 *  1) FATAL MESSAGES AT STARTUP
 */

void R_Suicide(char *s)
{
    char pp[1024];

    snprintf(pp, 1024, "Fatal error: %s\n", s);
    R_ShowMessage(pp);
    R_CleanUp(SA_SUICIDE, 2, 0);
}

/*
 *  2. CONSOLE I/O
 */

/*
 * We support 4 different type of input.
 * 1) from the gui console;
 * 2) from a character mode console (interactive);
 * 3) from a pipe under --ess, i.e, interactive.
 * 4) from a file or from a pipe (not interactive)
 *
 * Hence, it is better to have a different function for every
 * situation.
 * Same, it is true for output (but in this case, 2=3=4)
 *
 * BTW, 3 and 4 are different on input  since fgets,ReadFile...
 * "blocks" =>  (e.g.) you cannot give focus to the graphics device if
 * you are wating for input. For this reason, input is got in a different
 * thread
 *
 * All works in this way:
 * R_ReadConsole calls TrueReadConsole which points to:
 * case 1: GuiReadConsole
 * case 2 and 3: ThreadedReadConsole
 * case 4: FileReadConsole
 * ThreadedReadConsole wake up our 'reader thread' and wait until
 * a new line of input is available. The 'reader thread' uses
 * InThreadReadConsole to get it. InThreadReadConsole points to:
 * case 2: CharReadConsole
 * case 3: FileReadConsole
 */

/* Global variables */
static int (*TrueReadConsole)(char *, char *, int, int);
static int (*InThreadReadConsole)(char *, char *, int, int);
static void (*TrueWriteConsole)(char *, int);
HANDLE EhiWakeUp;
static char *tprompt, *tbuf;
static int tlen, thist, lineavailable;

/* Fill a text buffer with user typed console input. */
int R_ReadConsole(char *prompt, unsigned char *buf, int len, int addtohistory)
{
    R_ProcessEvents();
    return TrueReadConsole(prompt, (char *)buf, len, addtohistory);
}

/* Write a text buffer to the console. */
/* All system output is filtered through this routine. */

void R_WriteConsole(char *buf, int len)
{
    R_ProcessEvents();
    TrueWriteConsole(buf, len);
}

/*1: from GUI console */
int R_is_running = 0;

void Rconsolesetwidth(int cols)
{
    if (R_is_running && setWidthOnResize)
        R_SetOptionWidth(cols);
}

static int GuiReadConsole(char *prompt, char *buf, int len, int addtohistory)
{
    int res;
    char *p;
    char *NormalPrompt = (char *)CHAR(STRING_ELT(GetOption(install("prompt"), R_NilValue), 0));

    if (!R_is_running)
    {
        R_is_running = 1;
        Rconsolesetwidth(consolecols(RConsole));
    }
    ConsoleAcceptCmd = !strcmp(prompt, NormalPrompt);
    res = consolereads(RConsole, prompt, buf, len, addtohistory);
    for (p = buf; *p; p++)
        if (*p == EOF)
            *p = '\001';
    ConsoleAcceptCmd = 0;
    return !res;
}

/* 2 and 3: reading in a thread */

/* 'Reader thread' main function */
static void __cdecl ReaderThread(void *unused)
{
    while (1)
    {
        WaitForSingleObject(EhiWakeUp, INFINITE);
        tlen = InThreadReadConsole(tprompt, tbuf, tlen, thist);
        lineavailable = 1;
        PostThreadMessage(mainThreadId, 0, 0, 0);
    }
}

static int ThreadedReadConsole(char *prompt, char *buf, int len, int addtohistory)
{
    sighandler_t oldint, oldbreak;
    /*
     *   SIGINT/SIGBREAK when ESS is waiting for output are a real pain:
     *   they get processed after user hit <return>.
     *   The '^C\n' in raw Rterm is nice. But, do we really need it ?
     */
    oldint = signal(SIGINT, SIG_IGN);
    oldbreak = signal(SIGBREAK, SIG_IGN);
    mainThreadId = GetCurrentThreadId();
    lineavailable = 0;
    tprompt = prompt;
    tbuf = buf;
    tlen = len;
    thist = addtohistory;
    SetEvent(EhiWakeUp);
    while (1)
    {
        WaitMessage();
        if (lineavailable)
            break;
        doevent();
        if (R_tcldo)
            R_tcldo();
    }
    lineavailable = 0;
    /* restore handler  */
    signal(SIGINT, oldint);
    signal(SIGBREAK, oldbreak);
    return tlen;
}

/*2: from character console with getline (only used as InThreadReadConsole)*/
static int CharReadConsole(char *prompt, char *buf, int len, int addtohistory)
{
    int res = getline(prompt, buf, len);
    if (addtohistory)
        gl_histadd(buf);
    return !res;
}

/*3: (as InThreadReadConsole) and 4: non-interactive */
static int FileReadConsole(char *prompt, char *buf, int len, int addhistory)
{
    int ll;
    if (!R_Slave)
    {
        fputs(prompt, stdout);
        fflush(stdout);
    }
    if (fgets(buf, len, stdin) == NULL)
        return 0;
    /* according to system.txt, should be terminated in \n, so check this
       at eof */
    ll = strlen((char *)buf);
    if (feof(stdin) && buf[ll - 1] != '\n' && ll < len)
    {
        buf[ll++] = '\n';
        buf[ll] = '\0';
    }
    if (!R_Interactive && !R_Slave)
        fputs(buf, stdout);
    return 1;
}

/* Rgui */
static void GuiWriteConsole(char *buf, int len)
{
    char *p;

    for (p = buf; *p; p++)
        if (*p == '\001')
            *p = EOF;
    consolewrites(RConsole, buf);
}

/* Rterm write */
static void TermWriteConsole(char *buf, int len)
{
    printf("%s", buf);
}

/* Indicate that input is coming from the console */

void R_ResetConsole()
{
}

/* Stdio support to ensure the console file buffer is flushed */

void R_FlushConsole()
{
    if (CharacterMode == RTerm && R_Interactive)
        fflush(stdin);
    else if (CharacterMode == RGui)
        consoleflush(RConsole);
}

/* Reset stdin if the user types EOF on the console. */

void R_ClearerrConsole()
{
    if (CharacterMode == RTerm)
        clearerr(stdin);
}

/*
 *  3) ACTIONS DURING (LONG) COMPUTATIONS
 */

void GuiBusy(int which)
{
    if (which == 1)
        gsetcursor(RConsole, WatchCursor);
    if (which == 0)
        gsetcursor(RConsole, ArrowCursor);
}

void CharBusy(int which)
{
}

void R_Busy(int which)
{
    my_R_Busy(which);
}

/*
 *  4) INITIALIZATION AND TERMINATION ACTIONS
 */

/*
   R_CleanUp is invoked at the end of the session to give the user the
   option of saving their data.
   If ask == SA_SAVEASK the user should be asked if possible (and this
   option should not occur in non-interactive use).
   If ask = SA_SAVE or SA_NOSAVE the decision is known.
   If ask = SA_DEFAULT use the SaveAction set at startup.
   In all these cases run .Last() unless quitting is cancelled.
   If ask = SA_SUICIDE, no save, no .Last, possibly other things.
 */

void R_CleanUp(SA_TYPE saveact, int status, int runLast)
{
    if (saveact == SA_DEFAULT) /* The normal case apart from R_Suicide */
        saveact = SaveAction;

    if (saveact == SA_SAVEASK)
    {
        if (R_Interactive)
        {
            switch (R_YesNoCancel("Save workspace image?"))
            {
            case YES:
                saveact = SA_SAVE;
                break;
            case NO:
                saveact = SA_NOSAVE;
                break;
            case CANCEL:
                jump_to_toplevel();
                break;
            }
        }
        else
            saveact = SaveAction;
    }

    switch (saveact)
    {
    case SA_SAVE:
        if (runLast)
            R_dot_Last();
        if (R_DirtyImage)
            R_SaveGlobalEnv();
        if (CharacterMode == RGui || (R_Interactive && CharacterMode == RTerm))
            gl_savehistory(R_HistoryFile);
        break;
    case SA_NOSAVE:
        if (runLast)
            R_dot_Last();
        break;
    case SA_SUICIDE:
    default:
        break;
    }
    R_RunExitFinalizers();
    editorcleanall();
    CleanEd();
    CleanTempDir();
    closeAllHlpFiles();
    KillAllDevices();
    AllDevicesKilled = TRUE;
    if (R_Interactive && CharacterMode == RTerm)
        SetConsoleTitle(oldtitle);
#if 0
    UnLoad_Unzip_Dll();
#endif
    UnLoad_Rbitmap_Dll();
    if (R_CollectWarnings && saveact != SA_SUICIDE && CharacterMode == RTerm)
        PrintWarnings();
    app_cleanup();
    exit(status);
}

/*
 *  7) PLATFORM DEPENDENT FUNCTIONS
 */

/*
   This function can be used to display the named files with the
   given titles and overall title.  On GUI platforms we could
   use a read-only window to display the result.  Here we just
   make up a temporary file and invoke a pager on it.
*/

/*
 *     nfile   = number of files
 *     file    = array of filenames
 *     headers = the `headers' args of file.show. Printed before each file.
 *     wtitle  = title for window: the `title' arg of file.show
 *     del     = flag for whether files should be deleted after use
 *     pager   = pager to be used.
 */

int R_ShowFiles(int nfile, char **file, char **headers, char *wtitle, Rboolean del, char *pager)
{
    int i;
    char buf[1024];

    if (nfile > 0)
    {
        if (pager == NULL || strlen(pager) == 0)
            pager = "internal";
        for (i = 0; i < nfile; i++)
        {
            if (!access(file[i], R_OK))
            {
                if (!strcmp(pager, "internal"))
                {
                    newpager(wtitle, file[i], headers[i], del);
                }
                else if (!strcmp(pager, "console"))
                {
                    size_t len;
                    FILE *f = fopen(file[i], "rt");
                    if (f)
                    {
                        while ((len = fread(buf, 1, 1023, f)))
                        {
                            buf[len] = '\0';
                            R_WriteConsole(buf, strlen(buf));
                        }
                        fclose(f);
                        if (del)
                            DeleteFile(file[i]);
                    }
                    else
                    {
                        snprintf(buf, 1024, "Unable to open file '%s'", file[i]);
                        warning(buf);
                    }
                }
                else
                {
                    /* Quote path if necessary */
                    if (pager[0] != '"' && strchr(pager, ' '))
                        snprintf(buf, 1024, "\"%s\" \"%s\"", pager, file[i]);
                    else
                        snprintf(buf, 1024, "%s \"%s\"", pager, file[i]);
                    runcmd(buf, 0, 1, "");
                }
            }
            else
            {
                snprintf(buf, 1024, "file.show(): file %s does not exist\n", file[i]);
                warning(buf);
            }
        }
        return 0;
    }
    return 1;
}

int internal_ShowFile(char *file, char *header)
{
    SEXP pager = GetOption(install("pager"), R_NilValue);
    char *files[1], *headers[1];

    files[0] = file;
    headers[0] = header;
    return R_ShowFiles(1, files, headers, "File", 0, CHAR(STRING_ELT(pager, 0)));
}

/*
   This function can be used to open the named files in text editors, with the
   given titles and overall title.
   If the file does not exist then the editor should be opened to create a new file.
*/

/*
 *     nfile   = number of files
 *     file    = array of filenames
 *     editor  = editor to be used.
 */

int R_EditFiles(int nfile, char **file, char **title, char *editor)
{
    int i;
    char buf[1024];

    if (nfile > 0)
    {
        if (editor == NULL || strlen(editor) == 0)
            editor = "internal";
        for (i = 0; i < nfile; i++)
        {
            if (!strcmp(editor, "internal"))
            {
                Rgui_Edit(file[i], title[i], 0);
            }
            else
            {
                /* Quote path if necessary */
                if (editor[0] != '"' && strchr(editor, ' '))
                    snprintf(buf, 1024, "\"%s\" \"%s\"", editor, file[i]);
                else
                    snprintf(buf, 1024, "%s \"%s\"", editor, file[i]);
                runcmd(buf, 0, 1, "");
            }
        }
        return 0;
    }
    return 1;
}

/* Prompt the user for a file name.  Return the length of */
/* the name typed.  On Gui platforms, this should bring up */
/* a dialog box so a user can choose files that way. */

extern int DialogSelectFile(char *buf, int len); /* from rui.c */

int R_ChooseFile(int new, char *buf, int len)
{
    return (DialogSelectFile(buf, len));
}

/* code for R_ShowMessage, R_YesNoCancel */

void (*pR_ShowMessage)(char *s);
void R_ShowMessage(char *s)
{
    (*pR_ShowMessage)(s);
}

static void char_message(char *s)
{
    if (!s)
        return;
    if (R_Consolefile)
    {
        /* flush out standard output in case it uses R_Consolefile */
        if (R_Outputfile)
            fflush(R_Outputfile);
        fprintf(R_Consolefile, "%s\n", s);
        fflush(R_Consolefile);
    }
    else
        R_WriteConsole(s, strlen(s));
}

static int char_YesNoCancel(char *s)
{
    char ss[128];
    unsigned char a[3];

    sprintf(ss, "%s [y/n/c]: ", s);
    R_ReadConsole(ss, a, 3, 0);
    switch (a[0])
    {
    case 'y':
    case 'Y':
        return YES;
    case 'n':
    case 'N':
        return NO;
    default:
        return CANCEL;
    }
}

/*--- Initialization Code ---*/

static char RHome[MAX_PATH + 7];
static char UserRHome[MAX_PATH + 7];
static char RUser[MAX_PATH];
char *getRHOME(); /* in rhome.c */
void R_setStartTime();

void R_SetWin32(Rstart Rp)
{
    R_Home = Rp->rhome;
    if (strlen(R_Home) >= MAX_PATH)
        R_Suicide("Invalid R_HOME");
    sprintf(RHome, "R_HOME=%s", R_Home);
    putenv(RHome);
    strcpy(UserRHome, "R_USER=");
    strcat(UserRHome, Rp->home);
    putenv(UserRHome);

    CharacterMode = Rp->CharacterMode;
    switch (CharacterMode)
    {
    case RGui:
        R_GUIType = "Rgui";
        break;
    case RTerm:
        R_GUIType = "RTerm";
        break;
    default:
        R_GUIType = "unknown";
    }
    TrueReadConsole = Rp->ReadConsole;
    TrueWriteConsole = Rp->WriteConsole;
    R_CallBackHook = Rp->CallBack;
    pR_ShowMessage = Rp->ShowMessage;
    R_YesNoCancel = Rp->YesNoCancel;
    my_R_Busy = Rp->Busy;
    /* Process R_HOME/etc/Renviron.site, then
       .Renviron or ~/.Renviron, if it exists.
       Only used here in embedded versions */
    if (!Rp->NoRenviron)
    {
        process_site_Renviron();
        process_user_Renviron();
    }
    _controlfp(_MCW_EM, _MCW_EM);
    _controlfp(_PC_64, _MCW_PC);
}

/* Remove and process NAME=VALUE command line arguments */

static void Putenv(char *str)
{
    char *buf;
    buf = (char *)malloc((strlen(str) + 1) * sizeof(char));
    if (!buf)
        R_ShowMessage("allocation failure in reading Renviron");
    strcpy(buf, str);
    putenv(buf);
    /* no free here: storage remains in use */
}

static void env_command_line(int *pac, char **argv)
{
    int ac = *pac, newac = 1; /* Remember argv[0] is process name */
    char **av = argv;

    while (--ac)
    {
        ++av;
        if (**av != '-' && strchr(*av, '='))
            Putenv(*av);
        else
            argv[newac++] = *av;
    }
    *pac = newac;
}

char *PrintUsage(void)
{
    static char msg[5000];
    char msg0[] = "Start R, a system for statistical computation and graphics, with the\nspecified options\n\nEnvVars: "
                  "Environmental variables can be set by NAME=value strings\n\nOptions:\n  -h, --help            Print "
                  "usage message and exit\n  --version             Print version info and exit\n  --save               "
                  " Do save workspace at the end of the session\n  --no-save             Don't save it\n",
         msg1[] = "  --no-environ          Don't read the site and user environment files\n  --no-site-file        "
                  "Don't read the site-wide Rprofile\n  --no-init-file        Don't read the .Rprofile or ~/.Rprofile "
                  "files\n  --restore             Do restore previously saved objects at startup\n  --no-restore-data  "
                  "   Don't restore previously saved objects\n  --no-restore-history  Don't restore the R history "
                  "file\n  --no-restore          Don't restore anything\n",
         msg2[] = "  --vanilla             Combine --no-save, --no-restore, --no-site-file,\n                          "
                  "--no-init-file and --no-environ\n  --min-vsize=N         Set vector heap min to N bytes; '4M' = 4 "
                  "MegaB\n  --max-vsize=N         Set vector heap max to N bytes;\n  --min-nsize=N         Set min "
                  "number of cons cells to N\n  --max-nsize=N         Set max number of cons cells to N\n",
         msg2b[] = "  --max-mem-size=N      Set limit for memory to be used by R\n  --max-ppsize=N        Set max size "
                   "of protect stack to N\n",
         msg3[] = "  -q, --quiet           Don't print startup message\n  --silent              Same as --quiet\n  "
                  "--slave               Make R run as quietly as possible\n  --verbose             Print more "
                  "information about progress\n  --args                Skip the rest of the command line\n",
         msg4[] = "  --ess                 Don't use getline for command-line editing\n                          and "
                  "assert interactive use";
    if (CharacterMode == RTerm)
        strcpy(msg, "Usage: Rterm [options] [< infile] [> outfile] [EnvVars]\n\n");
    else
        strcpy(msg, "Usage: Rgui [options] [EnvVars]\n\n");
    strcat(msg, msg0);
    strcat(msg, msg1);
    strcat(msg, msg2);
    strcat(msg, msg2b);
    strcat(msg, msg3);
    if (CharacterMode == RTerm)
        strcat(msg, msg4);
    return msg;
}

#include <winbase.h>

char *getRUser()
{
    /*
     * try R_USER then HOME then Windows homes then working directory
     */
    char *p, *q;

    if ((p = getenv("R_USER")))
    {
        if (strlen(p) >= MAX_PATH)
            R_Suicide("Invalid R_USER");
        strcpy(RUser, p);
    }
    else if ((p = getenv("HOME")))
    {
        if (strlen(p) >= MAX_PATH)
            R_Suicide("Invalid HOME");
        strcpy(RUser, p);
    }
    else if (ShellGetPersonalDirectory(RUser))
    {
        /* nothing to do */;
    }
    else if ((p = getenv("HOMEDRIVE")) && (q = getenv("HOMEPATH")))
    {
        if (strlen(p) >= MAX_PATH)
            R_Suicide("Invalid HOMEDRIVE");
        strcpy(RUser, p);
        if (strlen(RUser) + strlen(q) >= MAX_PATH)
            R_Suicide("Invalid HOMEDRIVE+HOMEPATH");
        strcat(RUser, q);
    }
    else
    {
        GetCurrentDirectory(MAX_PATH, RUser);
    }
    p = RUser + (strlen(RUser) - 1);
    if (*p == '/' || *p == '\\')
        *p = '\0';
    return RUser;
}

int cmdlineoptions(int ac, char **av)
{
    int i, ierr;
    R_size_t value;
    char *p;
    char s[1024];
    structRstart rstart;
    Rstart Rp = &rstart;
    MEMORYSTATUS ms;
    Rboolean usedRdata = FALSE, processing = TRUE;

    /* ensure R_Home gets set early: we are in rgui or rterm here */
    R_Home = getRHOME();

#ifdef _R_HAVE_TIMING_
    R_setStartTime();
#endif

    /* Store the command line arguments before they are processed
       by the different option handlers. We do this here so that
       we get all the name=value pairs. Otherwise these will
       have been removed by the time we get to call
       R_common_command_line().
     */
    R_set_command_line_arguments(ac, av);

    /* set defaults for R_max_memory. This is set here so that
       embedded applications get no limit */
    GlobalMemoryStatus(&ms);
    R_max_memory = min(1024 * Mega, ms.dwTotalPhys);
    /* need enough to start R: fails on a 8Mb system */
    R_max_memory = max(16 * Mega, R_max_memory);

    R_DefParams(Rp);
    Rp->CharacterMode = CharacterMode;
    for (i = 1; i < ac; i++)
        if (!strcmp(av[i], "--no-environ") || !strcmp(av[i], "--vanilla"))
            Rp->NoRenviron = TRUE;

    Rp->CallBack = R_DoNothing;
    /* Here so that --ess and similar can change */
    InThreadReadConsole = NULL;
    if (CharacterMode == RTerm)
    {
        if (isatty(0) && isatty(1))
        {
            Rp->R_Interactive = TRUE;
            Rp->ReadConsole = ThreadedReadConsole;
            InThreadReadConsole = CharReadConsole;
        }
        else
        {
            Rp->R_Interactive = FALSE;
            Rp->ReadConsole = FileReadConsole;
        }
        /* Windows 95/98/ME have a shell that cannot redirect stderr,
           so don't use that on those OSes */
        {
            OSVERSIONINFO verinfo;
            GetVersionEx(&verinfo);
            switch (verinfo.dwPlatformId)
            {
            case VER_PLATFORM_WIN32_WINDOWS:
                R_Consolefile = stdout; /* used for errors */
                break;
            default:
                R_Consolefile = stderr; /* used for errors */
            }
        }
        R_Consolefile = stderr; /* used for errors */
        R_Outputfile = stdout;  /* used for sink-able output */
        Rp->WriteConsole = TermWriteConsole;
        Rp->ShowMessage = char_message;
        Rp->YesNoCancel = char_YesNoCancel;
        Rp->Busy = CharBusy;
    }
    else
    {
        Rp->R_Interactive = TRUE;
        Rp->ReadConsole = GuiReadConsole;
        Rp->WriteConsole = GuiWriteConsole;
        Rp->ShowMessage = askok;
        Rp->YesNoCancel = askyesnocancel;
        Rp->Busy = GuiBusy;
    }

    pR_ShowMessage = Rp->ShowMessage; /* used here */
    TrueWriteConsole = Rp->WriteConsole;
    R_CallBackHook = Rp->CallBack;

    /* process environment variables
     * precedence:  command-line, .Renviron, inherited
     */
    if (!Rp->NoRenviron)
    {
        process_site_Renviron();
        process_user_Renviron();
        Rp->NoRenviron = TRUE;
    }
    env_command_line(&ac, av);
    /*    R_SizeFromEnv(Rp); */

    R_common_command_line(&ac, av, Rp);

    while (--ac)
    {
        if (processing && **++av == '-')
        {
            if (!strcmp(*av, "--help") || !strcmp(*av, "-h"))
            {
                R_ShowMessage(PrintUsage());
                exit(0);
            }
            else if (!strcmp(*av, "--no-environ"))
            {
                Rp->NoRenviron = TRUE;
            }
            else if (!strcmp(*av, "--ess"))
            {
                /* Assert that we are interactive even if input is from a file */
                Rp->R_Interactive = TRUE;
                Rp->ReadConsole = ThreadedReadConsole;
                InThreadReadConsole = FileReadConsole;
            }
            else if (!strcmp(*av, "--internet2"))
            {
                UseInternet2 = TRUE;
            }
            else if (!strcmp(*av, "--mdi"))
            {
                MDIset = 1;
            }
            else if (!strcmp(*av, "--sdi") || !strcmp(*av, "--no-mdi"))
            {
                MDIset = -1;
            }
            else if (!strncmp(*av, "--max-mem-size", 14))
            {
                if (strlen(*av) < 16)
                {
                    ac--;
                    av++;
                    p = *av;
                }
                else
                    p = &(*av)[15];
                if (p == NULL)
                {
                    R_ShowMessage("WARNING: no max-mem-size given\n");
                    break;
                }
                value = R_Decode2Long(p, &ierr);
                if (ierr)
                {
                    if (ierr < 0)
                        sprintf(s, "WARNING: --max-mem-size value is invalid: ignored\n");
                    else
                        sprintf(s, "WARNING: --max-mem-size=%lu`%c': too large and ignored\n", (unsigned long)value,
                                (ierr == 1) ? 'M' : ((ierr == 2) ? 'K' : 'k'));
                    R_ShowMessage(s);
                }
                else if (value < 10 * Mega)
                {
                    sprintf(s, "WARNING: max-mem-size =%4.1fM too small and ignored\n", value / (1024.0 * 1024.0));
                    R_ShowMessage(s);
                }
                else
                    R_max_memory = value;
            }
            else if (!strcmp(*av, "--debug"))
            {
                DebugMenuitem = TRUE;
                breaktodebugger();
            }
            else if (!strcmp(*av, "--args"))
            {
                break;
            }
            else
            {
                snprintf(s, 1024, "WARNING: unknown option %s\n", *av);
                R_ShowMessage(s);
            }
        }
        else
        {
            /* Look for *.RData, as given by drag-and-drop */
            WIN32_FIND_DATA pp;
            HANDLE res;
            char *p, *fn = pp.cFileName, path[MAX_PATH];

            res = FindFirstFile(*av, &pp); /* convert to long file name */
            if (!usedRdata && res != INVALID_HANDLE_VALUE && strlen(fn) >= 6 &&
                stricmp(fn + strlen(fn) - 6, ".RData") == 0)
            {
                set_workspace_name(fn);
                strcpy(path, *av); /* this was generated by Windows so must fit */
                for (p = path; *p; p++)
                    if (*p == '\\')
                        *p = '/';
                p = strrchr(path, '/');
                if (p)
                {
                    *p = '\0';
                    chdir(path);
                }
                usedRdata = TRUE;
                Rp->RestoreAction = SA_RESTORE;
            }
            else
            {
                snprintf(s, 1024, "ARGUMENT '%s' __ignored__\n", *av);
                R_ShowMessage(s);
            }
            if (res != INVALID_HANDLE_VALUE)
                FindClose(res);
        }
    }
    Rp->rhome = R_Home;

    R_tcldo = tcl_do_none;
    Rp->home = getRUser();
    R_SetParams(Rp);

    /*
     *  Since users' expectations for save/no-save will differ, we decided
     *  that they should be forced to specify in the non-interactive case.
     */
    if (!R_Interactive && Rp->SaveAction != SA_SAVE && Rp->SaveAction != SA_NOSAVE)
        R_Suicide("you must specify `--save', `--no-save' or `--vanilla'");

    if (InThreadReadConsole &&
        (!(EhiWakeUp = CreateEvent(NULL, FALSE, FALSE, NULL)) || (_beginthread(ReaderThread, 0, NULL) == -1)))
        R_Suicide("impossible to create 'reader thread'; you must free some system resources");

    if ((R_HistoryFile = getenv("R_HISTFILE")) == NULL)
        R_HistoryFile = ".Rhistory";
    R_HistorySize = 512;
    if ((p = getenv("R_HISTSIZE")))
    {
        int value, ierr;
        value = R_Decode2Long(p, &ierr);
        if (ierr != 0 || value < 0)
            REprintf("WARNING: invalid R_HISTSIZE ignored;");
        else
            R_HistorySize = value;
    }
    return 0;
}

void setup_term_ui()
{
    initapp(0, 0);
    R_tcldo = tcl_do_none;
    readconsolecfg();
}

void saveConsoleTitle()
{
    GetConsoleTitle(oldtitle, 512);
}
