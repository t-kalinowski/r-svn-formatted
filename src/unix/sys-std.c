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
#include "Graphics.h" /* for devX11.h */
#include "devX11.h"
#include "Runix.h"

int SaveAction;

void fpu_setup(int); /* in sys-unix.c */

#ifdef HAVE_LIBREADLINE
#include <readline/readline.h>
#ifdef HAVE_READLINE_HISTORY_H
#include <readline/history.h>
#endif
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h> /* for unlink */
#endif

extern int UsingReadline;

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

#include "R_ext/eventloop.h"

/*
   This oject is used for the standard input and its file descriptor
   value is reset by setSelectMask() each time to ensure that it points
   to the correct value of stdin.
 */
static InputHandler BasicInputHandler = {StdinActivity, -1, NULL};

/*
   This can be reset by the initialization routines which
   can ignore stdin, etc..
*/
InputHandler *InputHandlers = &BasicInputHandler;

/*
  Initialize the input source handlers used to check for input on the
  different file descriptors.
 */
InputHandler *initStdinHandler(void)
{
    InputHandler *inputs;
    extern void R_processEvents(void);

    inputs = addInputHandler(InputHandlers, fileno(stdin), NULL, StdinActivity);
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
        InputHandlers = input;
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

static int setSelectMask(InputHandler *, fd_set *);

static InputHandler *waitForActivity()
{
    int maxfd;
    fd_set readMask;

    maxfd = setSelectMask(InputHandlers, &readMask);

    select(maxfd + 1, &readMask, NULL, NULL, NULL);

    return (getSelectedHandler(InputHandlers, &readMask));
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

    /* If we are dealing with BasicInputHandler Always put stdin */
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
      b) thi is the BasicInputHandler.
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
        if (strlen(line) && readline_addtohistory)
            add_history(line);
#endif
        l = (((readline_len - 2) > strlen(line)) ? strlen(line) : (readline_len - 2));
        strncpy(readline_buf, line, l);
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

/* Fill a text buffer with user typed console input. */

int Rstd_ReadConsole(char *prompt, unsigned char *buf, int len, int addtohistory)
{
    if (!R_Interactive)
    {
        if (!R_Slave)
            fputs(prompt, stdout);
        if (fgets(buf, len, stdin) == NULL)
            return 0;
        if (!R_Slave)
            fputs(buf, stdout);
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

        if (InputHandlers == NULL)
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
                        if (fgets(buf, len, stdin) == NULL)
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

void Rstd_CleanUp(int saveact, int status, int runLast)
{
    char buf[128];

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
    fpu_setup(0);

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

int Rstd_ShowFiles(int nfile, char **file, char **headers, char *wtitle, int del, char *pager)
{
    int c, i, res;
    char *filename;
    FILE *fp, *tfp;
    char buf[1024];

    if (nfile > 0)
    {
        if (pager == NULL || strlen(pager) == 0)
            pager = "more";
        filename = tmpnam(NULL);
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
    R_ReadConsole("Enter file name: ", buf, len, 0);
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
