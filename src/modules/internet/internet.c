/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2000, 2001   The R Development Core Team.
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

#include <Defn.h>
#include <Fileio.h>
#include <Rconnections.h>
#include <R_ext/R-ftp-http.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif

/* ------------------- internet access functions  --------------------- */

#if defined(USE_WININET_ASYNC) && !defined(USE_WININET)
#define USE_WININET 2
#endif

static Rboolean IDquiet = TRUE;

static void url_open(Rconnection con)
{
    void *ctxt;
    char *url = con->description;
    UrlScheme type = ((Rurlconn)(con->private))->type;

    if (con->mode[0] != 'r')
        error("can only open URLs for reading");

    switch (type)
    {
    case HTTPsh:
        ctxt = R_HTTPOpen(url);
        if (ctxt == NULL)
            error("cannot open URL `%s'", url);
        ((Rurlconn)(con->private))->ctxt = ctxt;
        break;
    case FTPsh:
        ctxt = R_FTPOpen(url);
        if (ctxt == NULL)
            error("cannot open URL `%s'", url);
        ((Rurlconn)(con->private))->ctxt = ctxt;
        break;
    default:
        error("unknown URL scheme");
    }

    con->isopen = TRUE;
    con->canwrite = (con->mode[0] == 'w' || con->mode[0] == 'a');
    con->canread = !con->canwrite;
    if (strlen(con->mode) >= 2 && con->mode[1] == 'b')
        con->text = FALSE;
    else
        con->text = TRUE;
    con->save = -1000;
}

static void url_close(Rconnection con)
{
    UrlScheme type = ((Rurlconn)(con->private))->type;
    switch (type)
    {
    case HTTPsh:
        R_HTTPClose(((Rurlconn)(con->private))->ctxt);
        break;
    case FTPsh:
        R_FTPClose(((Rurlconn)(con->private))->ctxt);
        break;
    }
    con->isopen = FALSE;
}

static int url_fgetc(Rconnection con)
{
    UrlScheme type = ((Rurlconn)(con->private))->type;
    void *ctxt = ((Rurlconn)(con->private))->ctxt;
    unsigned char c;
    size_t n = 0; /* -Wall */

    switch (type)
    {
    case HTTPsh:
        n = R_HTTPRead(ctxt, (char *)&c, 1);
        break;
    case FTPsh:
        n = R_FTPRead(ctxt, (char *)&c, 1);
        break;
    }
    return (n == 1) ? c : R_EOF;
}

static size_t url_read(void *ptr, size_t size, size_t nitems, Rconnection con)
{
    UrlScheme type = ((Rurlconn)(con->private))->type;
    void *ctxt = ((Rurlconn)(con->private))->ctxt;
    size_t n = 0; /* -Wall */

    switch (type)
    {
    case HTTPsh:
        n = R_HTTPRead(ctxt, ptr, size * nitems);
        break;
    case FTPsh:
        n = R_FTPRead(ctxt, ptr, size * nitems);
        break;
    }
    return n / size;
}

Rconnection R_newurl(char *description, char *mode)
{
    Rconnection new;

    new = (Rconnection)malloc(sizeof(struct Rconn));
    if (!new)
        error("allocation of url connection failed");
    new->class = (char *)malloc(strlen("file") + 1);
    if (!new->class)
    {
        free(new);
        error("allocation of url connection failed");
    }
    strcpy(new->class, "url");
    new->description = (char *)malloc(strlen(description) + 1);
    if (!new->description)
    {
        free(new->class);
        free(new);
        error("allocation of url connection failed");
    }
    init_con(new, description, mode);
    new->canwrite = FALSE;
    new->open = &url_open;
    new->close = &url_close;
    new->fgetc = &url_fgetc;
    new->read = &url_read;
    new->private = (void *)malloc(sizeof(struct urlconn));
    if (!new->private)
    {
        free(new->description);
        free(new->class);
        free(new);
        error("allocation of url connection failed");
    }

    IDquiet = TRUE;
    return new;
}

#ifndef Win32
static void putdots(int *pold, int new)
{
    int i, old = *pold;
    *pold = new;
    for (i = old; i < new; i++)
    {
        REprintf(".");
        if ((i + 1) % 50 == 0)
            REprintf("\n");
        else if ((i + 1) % 10 == 0)
            REprintf(" ");
    }
    if (R_Consolefile)
        fflush(R_Consolefile);
}
#endif

/* note, ALL the possible structures have the first two elements */
typedef struct
{
    int length;
    char *type;
    void *ctxt;
} inetconn;

#ifdef Win32
#include <graphapp/ga.h>
#endif

/* download(url, destfile, quiet) */

#define CPBUFSIZE 65536
#define IBUFSIZE 4096
SEXP do_download(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans, scmd, sfile, smode;
    char *url, *file, *mode;
    int quiet, status = 0;
#ifdef Win32
    window wprog;
    progressbar pb;
    label l_url;
#endif

    checkArity(op, args);
    scmd = CAR(args);
    if (!isString(scmd) || length(scmd) < 1)
        error("invalid `url' argument");
    if (length(scmd) > 1)
        warning("only first element of `url' argument used");
    url = CHAR(STRING_ELT(scmd, 0));
    sfile = CADR(args);
    if (!isString(sfile) || length(sfile) < 1)
        error("invalid `destfile' argument");
    if (length(sfile) > 1)
        warning("only first element of `destfile' argument used");
    file = CHAR(STRING_ELT(sfile, 0));
    IDquiet = quiet = asLogical(CADDR(args));
    if (quiet == NA_LOGICAL)
        error("invalid `quiet' argument");
    smode = CADDDR(args);
    if (!isString(smode) || length(smode) != 1)
        error("invalid `mode' argument");
    mode = CHAR(STRING_ELT(smode, 0));

    if (strncmp(url, "file://", 7) == 0)
    {
        FILE *in, *out;
        static char buf[CPBUFSIZE];
        size_t n;

        /* Use binary transfers */
        in = R_fopen(R_ExpandFileName(url + 7), (mode[2] == 'b') ? "rb" : "r");
        if (!in)
            error("cannot open URL `%s'", url);
        out = R_fopen(R_ExpandFileName(file), mode);
        if (!out)
            error("cannot open destfile `%s'", file);
        while ((n = fread(buf, 1, CPBUFSIZE, in)) > 0)
            fwrite(buf, 1, n, out);
        fclose(out);
        fclose(in);

#ifdef HAVE_INTERNET
    }
    else if (strncmp(url, "http://", 7) == 0)
    {

        FILE *out;
        void *ctxt;
        int len, total, guess, nnew, nbytes = 0;
        char buf[IBUFSIZE];
#ifndef Win32
        int ndots = 0;
#endif

        out = R_fopen(R_ExpandFileName(file), mode);
        if (!out)
            error("cannot open destfile `%s'", file);

        R_Busy(1);
        if (!quiet)
            REprintf("trying URL `%s'\n", url);
        ctxt = R_HTTPOpen(url);
        if (ctxt == NULL)
            status = 1;
        else
        {
            if (!quiet)
                REprintf("opened URL\n", url);
            guess = total = ((inetconn *)ctxt)->length;
            if (guess <= 0)
                guess = 100 * 1024;
#ifdef Win32
            wprog = newwindow("Download progress", rect(0, 0, 540, 100), Titlebar | Centered);
            setbackground(wprog, LightGrey);
            strcpy(buf, "URL: ");
            if (strlen(url) > 60)
            {
                strcat(buf, "... ");
                strcat(buf, url + (strlen(url) - 60));
            }
            else
                strcat(buf, url);
            l_url = newlabel(buf, rect(10, 15, 520, 25), AlignCenter);
            pb = newprogressbar(rect(20, 50, 500, 20), 0, guess, 1024, 1);
            show(wprog);
#endif
            while ((len = R_HTTPRead(ctxt, buf, sizeof(buf))) > 0)
            {
                fwrite(buf, 1, len, out);
                nbytes += len;
                nnew = nbytes / 1024;
#ifdef Win32
                if (nbytes > guess)
                {
                    guess *= 2;
                    setprogressbarrange(pb, 0, guess);
                }
                setprogressbar(pb, nbytes);
#else
                if (!quiet)
                    putdots(&ndots, nnew);
#endif
            }
            R_HTTPClose(ctxt);
            fclose(out);
            if (!quiet)
            {
                if (nbytes > 10240)
                    REprintf("\ndownloaded %dKb\n\n", nbytes / 1024, url);
                else
                    REprintf("\ndownloaded %d bytes\n\n", nbytes, url);
            }
#ifdef Win32
            hide(wprog);
            del(l_url);
            del(pb);
            del(wprog);
#endif
            if (total > 0 && total != nbytes)
                warning("downloaded length %d != reported length %d", nbytes, total);
        }
        R_Busy(0);
        if (status == 1)
            error("cannot open URL `%s'", url);
    }
    else if (strncmp(url, "ftp://", 6) == 0)
    {

        FILE *out;
        void *ctxt;
        int len, total, guess, nnew, nbytes = 0;
        char buf[IBUFSIZE];
#ifndef Win32
        int ndots = 0;
#endif

        out = R_fopen(R_ExpandFileName(file), mode);
        if (!out)
            error("cannot open destfile `%s'", file);

        R_Busy(1);
        if (!quiet)
            REprintf("trying URL `%s'\n", url);
        ctxt = R_FTPOpen(url);
        if (ctxt == NULL)
            status = 1;
        else
        {
            if (!quiet)
                REprintf("opened URL\n", url);
            guess = total = ((inetconn *)ctxt)->length;
            if (guess <= 0)
                guess = 100 * 1024;
#ifdef Win32
            wprog = newwindow("Download progress", rect(0, 0, 540, 100), Titlebar | Centered);
            setbackground(wprog, LightGrey);
            strcpy(buf, "URL: ");
            if (strlen(url) > 60)
            {
                strcat(buf, "... ");
                strcat(buf, url + (strlen(url) - 60));
            }
            else
                strcat(buf, url);
            l_url = newlabel(buf, rect(10, 15, 520, 25), AlignCenter);
            pb = newprogressbar(rect(20, 50, 500, 20), 0, guess, 1024, 1);
            show(wprog);
#endif
            while ((len = R_FTPRead(ctxt, buf, sizeof(buf))) > 0)
            {
                fwrite(buf, 1, len, out);
                nbytes += len;
                nnew = nbytes / 1024;
#ifdef Win32
                if (nbytes > guess)
                {
                    guess *= 2;
                    setprogressbarrange(pb, 0, guess);
                }
                setprogressbar(pb, nbytes);
#else
                if (!quiet)
                    putdots(&ndots, nnew);
#endif
            }
            R_FTPClose(ctxt);
            fclose(out);
            if (!quiet)
            {
                if (nbytes > 10240)
                    REprintf("\ndownloaded %dKb\n\n", nbytes / 1024, url);
                else
                    REprintf("\ndownloaded %d bytes\n\n", nbytes, url);
            }
#ifdef Win32
            hide(wprog);
            del(l_url);
            del(pb);
            del(wprog);
#endif
            if (total > 0 && total != nbytes)
                warning("downloaded length %d != reported length %d", nbytes, total);
        }
        R_Busy(0);
        if (status == 1)
            error("cannot open URL `%s'", url);
#endif
    }
    else
        error("unsupported URL scheme");

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = status;
    UNPROTECT(1);
    return ans;
}

#if defined(SUPPORT_LIBXML) && !defined(USE_WININET)

void *R_HTTPOpen(const char *url)
{
    inetconn *con;
    void *ctxt;
    int timeout = asInteger(GetOption(install("timeout"), R_NilValue));
    int len = -1;
    char *type = NULL;

    if (timeout == NA_INTEGER || timeout <= 0)
        timeout = 60;

    RxmlNanoHTTPTimeout(timeout);
    ctxt = RxmlNanoHTTPOpen(url, NULL);
    if (ctxt != NULL)
    {
        int rc = RxmlNanoHTTPReturnCode(ctxt);
        if (rc != 200)
        {
            RxmlNanoHTTPClose(ctxt);
            error("cannot open: HTTP status was `%d'", rc);
            return NULL;
        }
        else
        {
            type = RxmlNanoHTTPContentType(ctxt);
            len = RxmlNanoHTTPContentLength(ctxt);
            if (!IDquiet)
            {
                Rprintf("Content type `%s'", type ? type : "unknown");
                if (len >= 0)
                    Rprintf(" length %d bytes\n", len);
                else
                    Rprintf(" length unknown\n", len);
#ifdef Win32
                R_FlushConsole();
#endif
            }
        }
    }
    con = (inetconn *)malloc(sizeof(inetconn));
    if (con)
    {
        con->length = len;
        con->type = type;
        con->ctxt = ctxt;
    }
    return con;
}

int R_HTTPRead(void *ctx, char *dest, int len)
{
    return RxmlNanoHTTPRead(((inetconn *)ctx)->ctxt, dest, len);
}

void R_HTTPClose(void *ctx)
{
    if (ctx)
    {
        RxmlNanoHTTPClose(((inetconn *)ctx)->ctxt);
        free(ctx);
    }
}

void *R_FTPOpen(const char *url)
{
    inetconn *con;
    void *ctxt;
    int timeout = asInteger(GetOption(install("timeout"), R_NilValue));
    int len = 0;

    if (timeout == NA_INTEGER || timeout <= 0)
        timeout = 60;
    RxmlNanoFTPTimeout(timeout);
    ctxt = RxmlNanoFTPOpen(url);
    if (!ctxt)
        return NULL;
    if (!IDquiet)
    {
        len = RxmlNanoFTPContentLength(ctxt);
        if (len >= 0)
            Rprintf("ftp data connection made, file length %d bytes\n", len);
        else
            Rprintf("ftp data connection made, file length unknown\n");
#ifdef Win32
        R_FlushConsole();
#endif
    }
    con = (inetconn *)malloc(sizeof(inetconn));
    if (con)
    {
        con->length = len;
        con->type = NULL;
        con->ctxt = ctxt;
    }
    return con;
}

int R_FTPRead(void *ctx, char *dest, int len)
{
    return RxmlNanoFTPRead(((inetconn *)ctx)->ctxt, dest, len);
}

void R_FTPClose(void *ctx)
{
    if (ctx)
    {
        RxmlNanoFTPClose(((inetconn *)ctx)->ctxt);
        free(ctx);
    }
}
#endif /* SUPPORT_LIBXML */

#ifdef USE_WININET

#include <windows.h>
#include <wininet.h>
typedef struct wictxt
{
    int length;
    char *type;
    HINTERNET hand;
    HINTERNET session;
} wIctxt, *WIctxt;

#ifdef USE_WININET_ASYNC
static int timeout;

static int callback_status;
static LPINTERNET_ASYNC_RESULT callback_res;

static void CALLBACK InternetCallback(HINTERNET hInternet, DWORD context, DWORD Status, LPVOID lpvStatusInformation,
                                      DWORD dwStatusInformationLength)
{
    callback_status = Status;
    /* printf("callback with context %ld, code %ld\n", context, Status); */
    if (Status == INTERNET_STATUS_REQUEST_COMPLETE)
    {
        callback_res = (LPINTERNET_ASYNC_RESULT)lpvStatusInformation;
    }
}
#endif /* USE_WININET_ASYNC */

void *R_HTTPOpen(const char *url)
{
    WIctxt wictxt;
    DWORD status, d1 = 4, d2 = 0, d3 = 100;
    char buf[101];

    /*	BOOL res = InternetAttemptConnect(0);

        if (res != ERROR_SUCCESS) {
        warning("no Internet connection available");
        return NULL;
        }*/

    wictxt = (WIctxt)malloc(sizeof(wIctxt));
    wictxt->length = -1;
    wictxt->type = NULL;
    wictxt->hand = InternetOpen("R", INTERNET_OPEN_TYPE_PRECONFIG, NULL, NULL,
#ifdef USE_WININET_ASYNC
                                INTERNET_FLAG_ASYNC
#else
                                0
#endif
    );
    if (!wictxt->hand)
    {
        free(wictxt);
        error("cannot open Internet connection");
    }

#ifdef USE_WININET_ASYNC
    timeout = asInteger(GetOption(install("timeout"), R_NilValue));
    if (timeout == NA_INTEGER || timeout <= 0)
        timeout = 60;
    InternetSetStatusCallback(wictxt->hand, (INTERNET_STATUS_CALLBACK)InternetCallback);
    if (!IDquiet)
    {
        Rprintf("using Asynchronous WinInet calls, timeout %d secs\n", timeout);
        R_FlushConsole();
    }

    callback_status = 0;
    InternetOpenUrl(wictxt->hand, url, NULL, 0, INTERNET_FLAG_KEEP_CONNECTION | INTERNET_FLAG_NO_CACHE_WRITE, 17);

    {
        DWORD t1 = GetTickCount();
        while (callback_status != INTERNET_STATUS_REQUEST_COMPLETE && GetTickCount() < t1 + 1000 * timeout)
        {
            R_ProcessEvents();
            Sleep(100);
        }
        if (callback_status != INTERNET_STATUS_REQUEST_COMPLETE)
        {
            InternetCloseHandle(wictxt->hand);
            free(wictxt);
            error("InternetOpenUrl timed out");
        }
    }

    wictxt->session = (HINTERNET)callback_res->dwResult;
#else
    if (!IDquiet)
    {
        Rprintf("using Synchronous WinInet calls\n");
        R_FlushConsole();
    }
    wictxt->session =
        InternetOpenUrl(wictxt->hand, url, NULL, 0, INTERNET_FLAG_KEEP_CONNECTION | INTERNET_FLAG_NO_CACHE_WRITE, 0);
#endif /* USE_WININET_ASYNC */
    if (!wictxt->session)
    {
        DWORD err1 = GetLastError(), err2, blen = 101;
        InternetCloseHandle(wictxt->hand);
        free(wictxt);
        if (err1 == ERROR_INTERNET_EXTENDED_ERROR)
        {
            InternetGetLastResponseInfo(&err2, buf, &blen);
            error("InternetOpenUrl failed: `%s'", buf);
        }
        else
        {
            FormatMessage(FORMAT_MESSAGE_FROM_HMODULE, GetModuleHandle("wininet.dll"), err1,
                          MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), buf, 101, NULL);
            error("InternetOpenUrl failed: `%s'", buf);
        }
    }

    HttpQueryInfo(wictxt->session, HTTP_QUERY_STATUS_CODE | HTTP_QUERY_FLAG_NUMBER, &status, &d1, &d2);
    if (status != 200)
    {
        d2 = 0;
        HttpQueryInfo(wictxt->session, HTTP_QUERY_STATUS_TEXT, &buf, &d3, &d2);
        InternetCloseHandle(wictxt->session);
        InternetCloseHandle(wictxt->hand);
        free(wictxt);
        error("cannot open: HTTP status was `%d %s'", status, buf);
    }

    HttpQueryInfo(wictxt->session, HTTP_QUERY_CONTENT_TYPE, &buf, &d3, &d2);
    d2 = 0;
    HttpQueryInfo(wictxt->session, HTTP_QUERY_CONTENT_LENGTH | HTTP_QUERY_FLAG_NUMBER, &status, &d1, &d2);
    wictxt->length = status;
    wictxt->type = strdup(buf);
    if (!IDquiet)
    {
        Rprintf("Content type `%s' length %d bytes\n", buf, status);
        R_FlushConsole();
    }

    R_ProcessEvents();
    return (void *)wictxt;
}

int R_HTTPRead(void *ctx, char *dest, int len)
{
    DWORD nread;

    InternetReadFile(((WIctxt)ctx)->session, dest, len, &nread);
#ifdef USE_WININET_ASYNC
    {
        DWORD t1 = GetTickCount();
        while (callback_status != INTERNET_STATUS_REQUEST_COMPLETE && GetTickCount() < t1 + 1000 * timeout)
        {
            R_ProcessEvents();
            Sleep(100);
        }
        if (callback_status != INTERNET_STATUS_REQUEST_COMPLETE)
        {
            warning("Internet read timed out");
            nread = 0;
        }
    }
#endif
    R_ProcessEvents();
    return (int)nread;
}

void R_HTTPClose(void *ctx)
{
    InternetCloseHandle(((WIctxt)ctx)->session);
    InternetCloseHandle(((WIctxt)ctx)->hand);
    if (((WIctxt)ctx)->type)
        free(((WIctxt)ctx)->type);
    free(ctx);
}

void *R_FTPOpen(const char *url)
{
    WIctxt wictxt;

    wictxt = (WIctxt)malloc(sizeof(wIctxt));
    wictxt->length = -1;
    wictxt->type = NULL;

    wictxt->hand = InternetOpen("R", INTERNET_OPEN_TYPE_PRECONFIG, NULL, NULL,
#ifdef USE_WININET_ASYNC
                                INTERNET_FLAG_ASYNC
#else
                                0
#endif
    );
    if (!wictxt->hand)
    {
        free(wictxt);
        error("cannot open Internet connection");
    }

#ifdef USE_WININET_ASYNC
    timeout = asInteger(GetOption(install("timeout"), R_NilValue));
    if (timeout == NA_INTEGER || timeout <= 0)
        timeout = 60;
    InternetSetStatusCallback(wictxt->hand, (INTERNET_STATUS_CALLBACK)InternetCallback);
    if (!IDquiet)
    {
        Rprintf("using Asynchronous WinInet calls, timeout %d secs\n", timeout);
        R_FlushConsole();
    }

    callback_status = 0;
    InternetOpenUrl(wictxt->hand, url, NULL, 0, INTERNET_FLAG_KEEP_CONNECTION | INTERNET_FLAG_NO_CACHE_WRITE, 17);
    {
        DWORD t1 = GetTickCount();
        while (callback_status != INTERNET_STATUS_REQUEST_COMPLETE && GetTickCount() < t1 + 1000 * timeout)
        {
            R_ProcessEvents();
            Sleep(100);
        }
        if (callback_status != INTERNET_STATUS_REQUEST_COMPLETE)
        {
            InternetCloseHandle(wictxt->hand);
            free(wictxt);
            error("InternetOpenUrl timed out");
        }
    }

    wictxt->session = (HINTERNET)callback_res->dwResult;
#else
    if (!IDquiet)
    {
        Rprintf("using Synchronous WinInet calls\n");
        R_FlushConsole();
    }
    wictxt->session =
        InternetOpenUrl(wictxt->hand, url, NULL, 0, INTERNET_FLAG_KEEP_CONNECTION | INTERNET_FLAG_NO_CACHE_WRITE, 0);
    if (!wictxt->session)
    {
        char buf[256];
        DWORD err1 = GetLastError(), err2, blen = 256;
        InternetCloseHandle(wictxt->hand);
        free(wictxt);
        if (err1 == ERROR_INTERNET_EXTENDED_ERROR)
        {
            InternetGetLastResponseInfo(&err2, buf, &blen);
            error("InternetOpenUrl failed: `%s'", buf);
        }
        else
        {
            FormatMessage(FORMAT_MESSAGE_FROM_HMODULE, GetModuleHandle("wininet.dll"), err1,
                          MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), buf, 101, NULL);
            error("InternetOpenUrl failed: `%s'", buf);
        }
    }
#endif /* USE_WININET_ASYNC */
    R_ProcessEvents();
    return (void *)wictxt;
}

int R_FTPRead(void *ctx, char *dest, int len)
{
    return R_HTTPRead(ctx, dest, len);
}

void R_FTPClose(void *ctx)
{
    R_HTTPClose(ctx);
}
#endif

#ifndef HAVE_INTERNET
void *R_HTTPOpen(const char *url)
{
    return NULL;
}

int R_HTTPRead(void *ctx, char *dest, int len)
{
    return -1;
}

void R_HTTPClose(void *ctx)
{
}

void *R_FTPOpen(const char *url)
{
    return NULL;
}

int R_FTPRead(void *ctx, char *dest, int len)
{
    return -1;
}

void R_FTPClose(void *ctx)
{
}
#endif

#define MBUFSIZE 8192
void RxmlMessage(int level, const char *format, ...)
{
    int clevel;
    char buf[MBUFSIZE], *p;
    va_list(ap);

    clevel = asInteger(GetOption(install("internet.info"), R_NilValue));
    if (clevel == NA_INTEGER)
        clevel = 2;

    if (level < clevel)
        return;

    va_start(ap, format);
#ifdef HAVE_VSNPRINTF
    vsnprintf(buf, MBUFSIZE, format, ap);
    buf[MBUFSIZE - 1] = '\0';
#else
    vsprintf(buf, format, ap);
#endif
    va_end(ap);
    p = buf + strlen(buf) - 1;
    if (strlen(buf) > 0 && *p == '\n')
        *p = '\0';
    Rprintf(buf);
    Rprintf("\n");
}
