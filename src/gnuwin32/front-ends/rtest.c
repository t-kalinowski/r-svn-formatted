#include <windows.h>
#include <stdio.h>
#include "Rversion.h"
#include "Startup.h"
#include "graphapp/graphapp.h"

/* for signal-handling code */
#define PSIGNAL
#include "psignal.h"

#define SA_NORESTORE 0
#define SA_RESTORE 1

#define SA_DEFAULT 1
#define SA_NOSAVE 2
#define SA_SAVE 3
#define SA_SAVEASK 4
#define SA_SUICIDE 5

/* one way to allow user interrupts: called in ProcessEvents */
#ifdef VisualC
__declspec(dllimport) int UserBreak;
#else
#define UserBreak (*__imp_UserBreak)
extern int UserBreak;
#endif

/* calls into the R DLL */
extern char *getDLLVersion();
extern void R_DefParams(Rstart);
extern void R_SetParams(Rstart);
extern void setup_term_ui(void);
extern void ProcessEvents(void);
extern void setup_Rmainloop(), end_Rmainloop(), R_ReplDLLinit(), R_CleanUp();
extern int R_ReplDLLdo1();
extern void run_Rmainloop(void);

/* getline-based input, simple output */
#include <string.h>
#include "../getline/getline.h"
static char LastLine[512];

int myReadConsole(char *prompt, char *buf, int len, int addtohistory)
{
    static char *gl = NULL;
    int i;

    if (!gl)
    {
        strcat(LastLine, prompt);
        gl = getline(LastLine);
        LastLine[0] = '\0';
        if (addtohistory)
            gl_histadd(gl);
    }
    for (i = 0; *gl && (*gl != '\n') && (i < len - 2); gl++, i++)
        buf[i] = *gl;
    buf[i] = '\n';
    buf[i + 1] = '\0';
    if (!*gl || (*gl == '\n'))
        gl = NULL;
    return 1;
}

void myWriteConsole(char *buf, int len)
{
    char *p = strrchr(buf, '\n');

    if (p)
        strcpy(LastLine, p + 1);
    else
        strcat(LastLine, buf);
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

int main(int argc, char **argv)
{
    structRstart rp;
    Rstart Rp = &rp;
    char Rversion[25];

    sprintf(Rversion, "%s.%s", R_MAJOR, R_MINOR);
    if (strcmp(getDLLVersion(), Rversion) != 0)
    {
        fprintf(stderr, "Error: R.DLL version does not match\n");
        exit(1);
    }

    R_DefParams(Rp);
    Rp->rhome = "c:/R/rw0650";
    Rp->home = "c:/bdr";
    Rp->CharacterMode = LinkDLL;
    Rp->ReadConsole = myReadConsole;
    Rp->WriteConsole = myWriteConsole;
    Rp->CallBack = myCallBack;
    Rp->message = askok;
    Rp->yesnocancel = askyesnocancel;
    Rp->busy = myBusy;

    Rp->R_Quiet = 1;
    Rp->R_Interactive = 0;
    Rp->SaveAction = SA_NOSAVE;
    Rp->nsize = 300000;
    Rp->vsize = 6e6;
    R_SetParams(Rp);

    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    LastLine[0] = 0;

    signal(SIGBREAK, my_onintr);
    setup_term_ui();
    setup_Rmainloop();
#ifdef SIMPLE_CASE
    run_Rmainloop();
    end_Rmainloop();
#else
    R_ReplDLLinit();
    while (R_ReplDLLdo1() > 0)
    {
        /* add user actions here if desired */
    }
    /* only get here on EOF (not q()) */
    R_CleanUp(SA_DEFAULT);
#endif
    end_Rmainloop();
    return 0;
}
