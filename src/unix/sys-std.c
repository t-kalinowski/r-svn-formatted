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
#include "Devices.h" /* for KillAllDevices */
#include "Runix.h"
#include "Startup.h"

#include "Runix.h"

#ifdef HAVE_LIBREADLINE
#ifdef HAVE_READLINE_READLINE_H
#include <readline/readline.h>
#endif
#ifdef HAVE_READLINE_HISTORY_H
#include <readline/history.h>
#endif
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h> /* for unlink */
#endif

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h> /* for struct timeval */
#endif

extern SA_TYPE SaveAction;
extern Rboolean UsingReadline;

/*
 *  1) FATAL MESSAGES AT STARTUP
 */

void Rstd_Suicide(char *s)
{
    REprintf("Fatal error: %s\n", s);
    R_CleanUp(SA_SUICIDE, 2, 0);
}

/*
 *  2. CONSOLE I/O
 */

/*--- I/O Support Code ---*/

/* These routines provide hooks for supporting console I/O.
 * Under raw Unix these routines simply provide a
 * connection to the stdio library.
 * Under a Motif interface the routines would be
 * considerably more complex.
 */

#define __SYSTEM__
#include "R_ext/eventloop.h"
#undef __SYSTEM__

/*
   This object is used for the standard input and its file descriptor
   value is reset by setSelectMask() each time to ensure that it points
   to the correct value of stdin.
 */
static InputHandler BasicInputHandler = {StdinActivity, -1, NULL};

/*
   This can be reset by the initialization routines which
   can ignore stdin, etc..
*/
InputHandler *R_InputHandlers = &BasicInputHandler;

/*
  Initialize the input source handlers used to check for input on the
  different file descriptors.
 */
InputHandler *initStdinHandler(void)
{
    InputHandler *inputs;

    inputs = addInputHandler(R_InputHandlers, fileno(stdin), NULL, StdinActivity);
    /* Defer the X11 registration until it is loaded and actually used. */

    return (inputs);
}

/*
  Creates and registers a new InputHandler with the linked list `handlers'.
  This sets the global variable InputHandlers if it is not already set.
  In the standard interactive case, this will have been set to be the
  BasicInputHandler object.
 */
InputHandler *addInputHandler(InputHandler *handlers, int fd, InputHandlerProc handler, int activity)
{
    InputHandler *input, *tmp;
    input = (InputHandler *)calloc(1, sizeof(InputHandler));

    input->activity = activity;
    input->fileDescriptor = fd;
    input->handler = handler;

    tmp = handlers;

    if (handlers == NULL)
    {
        R_InputHandlers = input;
        return (input);
    }

    /* Go to the end of the list to append the new one.  */
    while (tmp->next != NULL)
    {
        tmp = tmp->next;
    }
    tmp->next = input;

    return (handlers);
}

/*
  Removes the specified handler from the linked list.

  See getInputHandler() for first locating the target handler instance.
 */
int removeInputHandler(InputHandler **handlers, InputHandler *it)
{
    InputHandler *tmp;

    /* If the handler is the first one in the list, move the list to point
       to the second element. That's why we use the address of the first
       element as the first argument.
    */
    if (*handlers == it)
    {
        *handlers = (*handlers)->next;
        return (1);
    }

    tmp = *handlers;

    while (tmp)
    {
        if (tmp->next == it)
        {
            tmp->next = it->next;
            return (1);
        }
    }

    return (0);
}

InputHandler *getInputHandler(InputHandler *handlers, int fd)
{
    InputHandler *tmp;
    tmp = handlers;

    while (tmp != NULL)
    {
        if (tmp->fileDescriptor == fd)
            return (tmp);
        tmp = tmp->next;
    }

    return (tmp);
}

/*
 Arrange to wait until there is some activity or input pending
 on one of the file descriptors to which we are listening.

 We could make the file descriptor mask persistent across
 calls and change it only when a listener is added or deleted.
 Later.


 This replaces the previous version which looked only on stdin and the X11
 device connection.  This allows more than one X11 device to be open on a different
 connection. Also, it allows connections a la S4 to be developed on top of this
 mechanism. The return type of this routine has changed.
*/

/* A package can enable polled event handling by making R_PolledEvents
   point to a non-dummy routine and setting R_wait_usec to a suitable
   timeout value (e.g. 100000) */

static void nop(void)
{
}

void (*R_PolledEvents)(void) = nop;

int R_wait_usec = 0; /* 0 means no timeout */

static int setSelectMask(InputHandler *, fd_set *);

static InputHandler *waitForActivity()
{
    int maxfd;
    fd_set readMask;
    struct timeval tv;

    do
    {
        R_PolledEvents();
        tv.tv_sec = 0;
        tv.tv_usec = R_wait_usec;
        maxfd = setSelectMask(R_InputHandlers, &readMask);
    } while (!select(maxfd + 1, &readMask, NULL, NULL, (R_wait_usec) ? &tv : NULL));

    return (getSelectedHandler(R_InputHandlers, &readMask));
}

/*
  Create the mask representing the file descriptors select() should
  monitor and return the maximum of these file descriptors so that
  it can be passed directly to select().

  If the first element of the handlers is the standard input handler
  then we set its file descriptor to the current value of stdin - its
  file descriptor.
 */

static int setSelectMask(InputHandler *handlers, fd_set *readMask)
{
    int maxfd = -1;
    InputHandler *tmp = handlers;
    FD_ZERO(readMask);

    /* If we are dealing with BasicInputHandler always put stdin */
    if (handlers == &BasicInputHandler)
        handlers->fileDescriptor = fileno(stdin);

    while (tmp)
    {
        FD_SET(tmp->fileDescriptor, readMask);
        maxfd = maxfd < tmp->fileDescriptor ? tmp->fileDescriptor : maxfd;
        tmp = tmp->next;
    }

    return (maxfd);
}

/*
  Determine which handler was identified as having input pending
  by the call to select().
  We have a very simple version of scheduling. We skip the first one
  if it is the standard input one. This allows the others to not be `starved'.
  Change this as one wants by not skipping the first one.
 */
InputHandler *getSelectedHandler(InputHandler *handlers, fd_set *readMask)
{
    InputHandler *tmp = handlers;

    /*
      Temporarily skip the first one if a) there is another one, and
      b) this is the BasicInputHandler.
    */
    if (handlers == &BasicInputHandler && handlers->next)
        tmp = handlers->next;

    while (tmp)
    {
        if (FD_ISSET(tmp->fileDescriptor, readMask))
            return (tmp);
        tmp = tmp->next;
    }
    /* Now deal with the first one. */
    if (FD_ISSET(handlers->fileDescriptor, readMask))
        return (handlers);

    return ((InputHandler *)NULL);
}

#ifdef HAVE_LIBREADLINE
/* callback for rl_callback_read_char */

static int readline_gotaline;
static int readline_addtohistory;
static int readline_len;
static int readline_eof;
static unsigned char *readline_buf;

static void readline_handler(unsigned char *line)
{
    int l;
    rl_callback_handler_remove();
    if ((readline_eof = !line)) /* Yes, I don't mean ==...*/
        return;
    if (line[0])
    {
#ifdef HAVE_READLINE_HISTORY_H
        if (strlen((char *)line) && readline_addtohistory)
            add_history((char *)line);
#endif
        l = (((readline_len - 2) > strlen((char *)line)) ? strlen((char *)line) : (readline_len - 2));
        strncpy((char *)readline_buf, (char *)line, l);
        readline_buf[l] = '\n';
        readline_buf[l + 1] = '\0';
    }
    else
    {
        readline_buf[0] = '\n';
        readline_buf[1] = '\0';
    }
    readline_gotaline = 1;
}
#endif

/* Fill a text buffer from stdin or with user typed console input. */

int Rstd_ReadConsole(char *prompt, unsigned char *buf, int len, int addtohistory)
{
    if (!R_Interactive)
    {
        int ll;
        if (!R_Slave)
            fputs(prompt, stdout);
        if (fgets((char *)buf, len, stdin) == NULL)
            return 0;
        ll = strlen((char *)buf);
        /* remove CR in CRLF ending */
        if (buf[ll - 1] == '\n' && buf[ll - 2] == '\r')
        {
            buf[ll - 2] = '\n';
            buf[--ll] = '\0';
        }
        /* according to system.txt, should be terminated in \n, so check this
           at eof */
        if (feof(stdin) && buf[ll - 1] != '\n' && ll < len)
        {
            buf[ll++] = '\n';
            buf[ll] = '\0';
        }
        if (!R_Slave)
            fputs((char *)buf, stdout);
        return 1;
    }
    else
    {
#ifdef HAVE_LIBREADLINE
        if (UsingReadline)
        {
            readline_gotaline = 0;
            readline_buf = buf;
            readline_addtohistory = addtohistory;
            readline_len = len;
            readline_eof = 0;
            rl_callback_handler_install(prompt, readline_handler);
        }
        else
#endif
        {
            fputs(prompt, stdout);
            fflush(stdout);
        }

        if (R_InputHandlers == NULL)
            initStdinHandler();

        for (;;)
        {
            InputHandler *what = waitForActivity();
            if (what != NULL)
            {
                if (what->fileDescriptor == fileno(stdin))
                {
                    /* We could make this a regular handler, but we need to pass additional arguments. */
#ifdef HAVE_LIBREADLINE
                    if (UsingReadline)
                    {
                        rl_callback_read_char();
                        if (readline_eof)
                            return 0;
                        if (readline_gotaline)
                            return 1;
                    }
                    else
#endif
                    {
                        if (fgets((char *)buf, len, stdin) == NULL)
                            return 0;
                        else
                            return 1;
                    }
                }
                else
                    what->handler((void *)NULL);
            }
        }
    }
}

/* Write a text buffer to the console. */
/* All system output is filtered through this routine. */

void Rstd_WriteConsole(char *buf, int len)
{
    printf("%s", buf);
}

/* Indicate that input is coming from the console */

void Rstd_ResetConsole()
{
}

/* Stdio support to ensure the console file buffer is flushed */

void Rstd_FlushConsole()
{
    fflush(stdin);
}

/* Reset stdin if the user types EOF on the console. */

void Rstd_ClearerrConsole()
{
    clearerr(stdin);
}

/*
 *  3) ACTIONS DURING (LONG) COMPUTATIONS
 */

void Rstd_Busy(int which)
{
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

void R_dot_Last(void); /* in main.c */

void Rstd_CleanUp(SA_TYPE saveact, int status, int runLast)
{
    unsigned char buf[128];

    if (saveact == SA_DEFAULT) /* The normal case apart from R_Suicide */
        saveact = SaveAction;

    if (saveact == SA_SAVEASK)
    {
        if (R_Interactive)
        {
        qask:
            R_ClearerrConsole();
            R_FlushConsole();
            R_ReadConsole("Save workspace image? [y/n/c]: ", buf, 128, 0);
            switch (buf[0])
            {
            case 'y':
            case 'Y':
                saveact = SA_SAVE;
                break;
            case 'n':
            case 'N':
                saveact = SA_NOSAVE;
                break;
            case 'c':
            case 'C':
                jump_to_toplevel();
                break;
            default:
                goto qask;
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
#ifdef HAVE_LIBREADLINE
#ifdef HAVE_READLINE_HISTORY_H
        if (R_Interactive && UsingReadline)
        {
            stifle_history(R_HistorySize);
            write_history(R_HistoryFile);
        }
#endif
#endif
        break;
    case SA_NOSAVE:
        if (runLast)
            R_dot_Last();
        break;
    case SA_SUICIDE:
    default:
        break;
    }
    CleanEd();
    KillAllDevices();
    if (saveact != SA_SUICIDE && R_CollectWarnings)
        PrintWarnings(); /* from device close and .Last */
    fpu_setup(FALSE);

    exit(status);
}

/*
 *  7) PLATFORM DEPENDENT FUNCTIONS
 */

int Rstd_ShowFiles(int nfile,      /* number of files */
                   char **file,    /* array of filenames */
                   char **headers, /* the `headers' args of file.show.
                              Printed before each file. */
                   char *wtitle,   /* title for window
                              = `title' arg of file.show */
                   Rboolean del,   /* should files be deleted after use? */
                   char *pager)    /* pager to be used */

{
    /*
        This function can be used to display the named files with the
        given titles and overall title.	 On GUI platforms we could
        use a read-only window to display the result.  Here we just
        make up a temporary file and invoke a pager on it.
    */

    int c, i, res;
    char *filename;
    FILE *fp, *tfp;
    char buf[1024];

    if (nfile > 0)
    {
        if (pager == NULL || strlen(pager) == 0)
            pager = "more";
        filename = Runix_tmpnam(NULL);
        if ((tfp = fopen(filename, "w")) != NULL)
        {
            for (i = 0; i < nfile; i++)
            {
                if (headers[i] && *headers[i])
                    fprintf(tfp, "%s\n\n", headers[i]);
                if ((fp = R_fopen(R_ExpandFileName(file[i]), "r")) != NULL)
                {
                    while ((c = fgetc(fp)) != EOF)
                        fputc(c, tfp);
                    fprintf(tfp, "\n");
                    fclose(fp);
                    if (del)
                        unlink(R_ExpandFileName(file[i]));
                }
                else
                    fprintf(tfp, "NO FILE %s\n\n", file[i]);
            }
            fclose(tfp);
        }
        sprintf(buf, "%s < %s", pager, filename);
        res = system(buf);
        unlink(filename);
        return (res != 0);
    }
    return 1;
}

/*
   Prompt the user for a file name.  Return the length of
   the name typed.  On Gui platforms, this should bring up
   a dialog box so a user can choose files that way.
*/

int Rstd_ChooseFile(int new, char *buf, int len)
{
    int namelen;
    char *bufp;
    R_ReadConsole("Enter file name: ", (unsigned char *)buf, len, 0);
    namelen = strlen(buf);
    bufp = &buf[namelen - 1];
    while (bufp >= buf && isspace((int)*bufp))
        *bufp-- = '\0';
    return strlen(buf);
}

void Rstd_ShowMessage(char *s)
{
    REprintf(s);
}

void Rstd_read_history(char *s)
{
#ifdef HAVE_LIBREADLINE
#ifdef HAVE_READLINE_HISTORY_H
    if (R_Interactive && UsingReadline)
    {
        read_history(s);
    }
#endif
#endif
}

void Rstd_loadhistory(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP sfile;
    char file[PATH_MAX];

    checkArity(op, args);
    sfile = CAR(args);
    if (!isString(sfile) || LENGTH(sfile) < 1)
        errorcall(call, "invalid file argument");
    strcpy(file, R_ExpandFileName(CHAR(STRING_ELT(sfile, 0))));
#if defined(HAVE_LIBREADLINE) && defined(HAVE_READLINE_HISTORY_H)
    if (R_Interactive && UsingReadline)
    {
        clear_history();
        read_history(file);
    }
    else
        errorcall(call, "no history mechanism available");
#else
    errorcall(call, "no history mechanism available");
#endif
}

void Rstd_savehistory(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP sfile;
    char file[PATH_MAX];

    checkArity(op, args);
    sfile = CAR(args);
    if (!isString(sfile) || LENGTH(sfile) < 1)
        errorcall(call, "invalid file argument");
    strcpy(file, R_ExpandFileName(CHAR(STRING_ELT(sfile, 0))));
#if defined(HAVE_LIBREADLINE) && defined(HAVE_READLINE_HISTORY_H)
    if (R_Interactive && UsingReadline)
    {
        write_history(file);
        history_truncate_file(file, R_HistorySize);
    }
    else
        errorcall(call, "no history available to save");
#else
    errorcall(call, "no history available to save");
#endif
}

#include <setjmp.h>
static jmp_buf sleep_return;

static int OldTimeout;
static void (*OldHandler)(void);

#ifdef HAVE_TIMES
#include <time.h>
#include <sys/times.h>
#ifndef CLK_TCK
/* this is in ticks/second, generally 60 on BSD style Unix, 100? on SysV */
#ifdef HZ
#define CLK_TCK HZ
#else
#define CLK_TCK 60
#endif
#endif /* CLK_TCK */

static struct tms timeinfo;
static double timeint, start, elapsed;

static void SleepHandler(void)
{
    elapsed = (times(&timeinfo) - start) / (double)CLK_TCK;
    /*    Rprintf("elapsed %f,  R_wait_usec %d\n", elapsed, R_wait_usec); */
    if (elapsed >= timeint)
        longjmp(sleep_return, 100);
    if (timeint - elapsed < 0.5)
        R_wait_usec = 1e6 * (timeint - elapsed) + 10000;
    OldHandler();
}

SEXP do_syssleep(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    checkArity(op, args);
    timeint = asReal(CAR(args));
    if (ISNAN(timeint) || timeint < 0)
        errorcall(call, "invalid time value");
    OldHandler = R_PolledEvents;
    R_PolledEvents = SleepHandler;
    OldTimeout = R_wait_usec;
    if (OldTimeout == 0 || OldTimeout > 500000)
        R_wait_usec = 500000;

    start = times(&timeinfo);
    if (setjmp(sleep_return) != 100)
        for (;;)
        {
            InputHandler *what = waitForActivity();
            if (what != NULL)
            {
                if (what->fileDescriptor != fileno(stdin))
                    what->handler((void *)NULL);
                else
                    usleep(R_wait_usec / 2);
                /* we can't handle console read events here,
                   so just sleep for a while */
            }
        }

    R_PolledEvents = OldHandler;
    R_wait_usec = OldTimeout;
    return R_NilValue;
}

#else
SEXP do_syssleep(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    error("Sys.sleep is not implemented on this system")
}
#endif
