/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2000,2001   The R Development Core Team.
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
#include <R_ext/Complex.h>
/* #include <fcntl.h> not yet */

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* Win32 does have popen, but it does not work in GUI applications,
   so test that later */
#ifdef Win32
#define HAVE_POPEN
#include <Startup.h>
extern UImode CharacterMode;
#endif

#define NCONNECTIONS 50
#define NSINKS 21

static Rconnection Connections[NCONNECTIONS];

static int R_SinkNumber, R_SinkSaved;
static int SinkCons[NSINKS], SinkConsClose[NSINKS];

/* ------------- admin functions (see also at end) ----------------- */

int NextConnection()
{
    int i;
    for (i = 3; i < NCONNECTIONS; i++)
        if (!Connections[i])
            break;
    if (i > NCONNECTIONS)
        error("All connections are in use");
    return i;
}

/* internal, not the same as R function getConnection */
Rconnection getConnection(int n)
{
    Rconnection con = NULL;

    if (n < 0 || n == NA_INTEGER || !(con = Connections[n]))
        error("invalid connection");
    return con;
}

void Rconn_setEncoding(Rconnection con, SEXP enc)
{
    int i;

    if (!isInteger(enc) || length(enc) != 256)
        error("invalid `enc' argument");
    for (i = 0; i < 256; i++)
        con->encoding[i] = (unsigned char)INTEGER(enc)[i];
}

/* ------------------- null connection functions --------------------- */

static void null_open(Rconnection con)
{
    error("open/close not enabled for this connection");
}

static int null_vfprintf(Rconnection con, const char *format, va_list ap)
{
    error("printing not enabled for this connection");
    return 0; /* -Wall */
}

static int null_fgetc(Rconnection con)
{
    error("getc not enabled for this connection");
    return 0; /* -Wall */
}

static long null_seek(Rconnection con, int where, int origin)
{
    error("seek not enabled for this connection");
    return 0; /* -Wall */
}

static void null_truncate(Rconnection con)
{
    error("truncate not enabled for this connection");
}

static int null_fflush(Rconnection con)
{
    return 0;
}

static size_t null_read(void *ptr, size_t size, size_t nitems, Rconnection con)
{
    error("read not enabled for this connection");
    return 0; /* -Wall */
}

static size_t null_write(const void *ptr, size_t size, size_t nitems, Rconnection con)
{
    error("write not enabled for this connection");
    return 0; /* -Wall */
}

/* ------------------- file connections --------------------- */

static void file_open(Rconnection con)
{
    FILE *fp;
    /*    int fd, flags; */ /* fcntl does not exist on Windows */

    fp = R_fopen(R_ExpandFileName(con->description), con->mode);
    if (!fp)
        error("cannot open file `%s'", R_ExpandFileName(con->description));
    ((Rfileconn)(con->private))->fp = fp;
    con->isopen = TRUE;
    con->canwrite = (con->mode[0] == 'w' || con->mode[0] == 'a');
    con->canread = !con->canwrite;
    if (strlen(con->mode) >= 2 && con->mode[1] == 'b')
        con->text = FALSE;
    else
        con->text = TRUE;
    con->save = -1000;

    /*    if(!con->blocking) {
        fd = fileno(fp);
        flags = fcntl(fd, F_GETFL);
        flags |= O_NONBLOCK;
        fcntl(fd, F_SETFL, flags);
        }*/
}

static void file_close(Rconnection con)
{
    fclose(((Rfileconn)(con->private))->fp);
    con->isopen = FALSE;
}

static void file_destroy(Rconnection con)
{
    free(con->private);
}

static int file_vfprintf(Rconnection con, const char *format, va_list ap)
{
    FILE *fp = ((Rfileconn)(con->private))->fp;
    return vfprintf(fp, format, ap);
}

static int file_fgetc(Rconnection con)
{
    FILE *fp = ((Rfileconn)(con->private))->fp;
    int c = con->encoding[fgetc(fp)];
    return feof(fp) ? R_EOF : c;
}

static long file_seek(Rconnection con, int where, int origin)
{
    FILE *fp = ((Rfileconn)(con->private))->fp;
    long pos = ftell(fp);
    int whence = SEEK_SET;

    switch (origin)
    {
    case 2:
        whence = SEEK_CUR;
    case 3:
        whence = SEEK_END;
    default:
        whence = SEEK_SET;
    }
    if (where >= 0)
        fseek(fp, where, whence);
    return pos;
}

static void file_truncate(Rconnection con)
{
    FILE *fp = ((Rfileconn)(con->private))->fp;
    int fd = fileno(fp);
    int size = lseek(fd, 0, SEEK_CUR);

    if (!con->isopen || !con->canwrite)
        error("can only truncate connections open for writing");

#ifdef HAVE_FTRUNCATE
    if (ftruncate(fd, size))
        error("file truncation failed");
#else
#ifdef Win32
    if (chsize(fd, size))
        error("file truncation failed");
#else
    error("Unavailable on this platform");
#endif
#endif
}

static int file_fflush(Rconnection con)
{
    FILE *fp = ((Rfileconn)(con->private))->fp;

    return fflush(fp);
}

static size_t file_read(void *ptr, size_t size, size_t nitems, Rconnection con)
{
    FILE *fp = ((Rfileconn)(con->private))->fp;
    return fread(ptr, size, nitems, fp);
}

static size_t file_write(const void *ptr, size_t size, size_t nitems, Rconnection con)
{
    FILE *fp = ((Rfileconn)(con->private))->fp;
    return fwrite(ptr, size, nitems, fp);
}

static Rconnection newfile(char *description, char *mode)
{
    Rconnection new;
    new = (Rconnection)malloc(sizeof(struct Rconn));
    if (!new)
        error("allocation of file connection failed");
    new->class = (char *)malloc(strlen("file") + 1);
    if (!new->class)
    {
        free(new);
        error("allocation of file connection failed");
    }
    strcpy(new->class, "file");
    new->description = (char *)malloc(strlen(description) + 1);
    if (!new->description)
    {
        free(new->class);
        free(new);
        error("allocation of file connection failed");
    }
    strcpy(new->description, description);
    strncpy(new->mode, mode, 4);
    new->mode[4] = '\0';
    new->isopen = new->incomplete = FALSE;
    new->canread = new->canwrite = TRUE; /* in principle */
    new->canseek = TRUE;
    new->text = TRUE;
    new->open = &file_open;
    new->close = &file_close;
    new->destroy = &file_destroy;
    new->vfprintf = &file_vfprintf;
    new->fgetc = &file_fgetc;
    new->seek = &file_seek;
    new->truncate = &file_truncate;
    new->fflush = &file_fflush;
    new->read = &file_read;
    new->write = &file_write;
    new->nPushBack = 0;
    new->private = (void *)malloc(sizeof(struct fileconn));
    if (!new->private)
    {
        free(new->description);
        free(new->class);
        free(new);
        error("allocation of file connection failed");
    }
    return new;
}

SEXP do_file(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP sfile, sopen, ans, class, enc;
    char *file, *open;
    int i, ncon, block;
    Rconnection con = NULL;

    checkArity(op, args);
    sfile = CAR(args);
    if (!isString(sfile) || length(sfile) < 1)
        errorcall(call, "invalid `description' argument");
    if (length(sfile) > 1)
        warning("only first element of `description' argument used");
    file = CHAR(STRING_ELT(sfile, 0));
    sopen = CADR(args);
    if (!isString(sopen) || length(sopen) != 1)
        error("invalid `open' argument");
    block = asLogical(CADDR(args));
    if (block == NA_LOGICAL)
        error("invalid `block' argument");
    enc = CADDDR(args);
    if (!isInteger(enc) || length(enc) != 256)
        error("invalid `enc' argument");
    open = CHAR(STRING_ELT(sopen, 0));
    ncon = NextConnection();
    con = Connections[ncon] = newfile(file, strlen(open) ? open : "r");
    for (i = 0; i < 256; i++)
        con->encoding[i] = (unsigned char)INTEGER(enc)[i];

    /* open it if desired */
    if (strlen(open))
        con->open(con);

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = ncon;
    PROTECT(class = allocVector(STRSXP, 2));
    SET_STRING_ELT(class, 0, mkChar("file"));
    SET_STRING_ELT(class, 1, mkChar("connection"));
    classgets(ans, class);
    UNPROTECT(2);

    return ans;
}

/* ------------------- pipe connections --------------------- */

#ifdef HAVE_POPEN
static void pipe_open(Rconnection con)
{
    FILE *fp;

    fp = popen(con->description, con->mode);
    if (!fp)
        error("cannot open cmd `%s'", con->description);
    ((Rfileconn)(con->private))->fp = fp;
    con->isopen = TRUE;
    con->canwrite = (con->mode[0] == 'w');
    con->canread = !con->canwrite;
    if (strlen(con->mode) >= 2 && con->mode[1] == 'b')
        con->text = FALSE;
    else
        con->text = TRUE;
    con->save = -1000;
}

static void pipe_close(Rconnection con)
{
    pclose(((Rfileconn)(con->private))->fp);
    con->isopen = FALSE;
}

static void pipe_destroy(Rconnection con)
{
    free(con->private);
}

static Rconnection newpipe(char *description, char *mode)
{
    Rconnection new;
    new = (Rconnection)malloc(sizeof(struct Rconn));
    if (!new)
        error("allocation of pipe connection failed");
    new->class = (char *)malloc(strlen("pipe") + 1);
    if (!new->class)
    {
        free(new);
        error("allocation of pipe connection failed");
    }
    strcpy(new->class, "pipe");
    new->description = (char *)malloc(strlen(description) + 1);
    if (!new->description)
    {
        free(new->class);
        free(new);
        error("allocation of pipe connection failed");
    }
    strcpy(new->description, description);
    strncpy(new->mode, mode, 4);
    new->mode[4] = '\0';
    new->isopen = new->incomplete = FALSE;
    new->canread = new->canwrite = TRUE; /* in principle */
    new->canseek = FALSE;
    new->text = TRUE;
    new->open = &pipe_open;
    new->close = &pipe_close;
    new->destroy = &pipe_destroy;
    new->vfprintf = &file_vfprintf;
    new->fgetc = &file_fgetc;
    new->seek = &null_seek;
    new->truncate = &null_truncate;
    new->fflush = &file_fflush;
    new->read = &file_read;
    new->write = &file_write;
    new->nPushBack = 0;
    new->private = (void *)malloc(sizeof(struct fileconn));
    if (!new->private)
    {
        free(new->description);
        free(new->class);
        free(new);
        error("allocation of pipe connection failed");
    }
    return new;
}
#endif

#ifdef Win32
extern Rconnection newWpipe(char *description, char *mode);
#endif

SEXP do_pipe(SEXP call, SEXP op, SEXP args, SEXP env)
{
#ifdef HAVE_POPEN
    SEXP scmd, sopen, ans, class, enc;
    char *file, *open;
    int i, ncon;
    Rconnection con = NULL;

    checkArity(op, args);
    scmd = CAR(args);
    if (!isString(scmd) || length(scmd) < 1)
        error("invalid `description' argument");
    if (length(scmd) > 1)
        warning("only first element of `description' argument used");
    file = CHAR(STRING_ELT(scmd, 0));
    sopen = CADR(args);
    if (!isString(sopen) || length(sopen) != 1)
        error("invalid `open' argument");
    open = CHAR(STRING_ELT(sopen, 0));
    enc = CADDR(args);
    if (!isInteger(enc) || length(enc) != 256)
        error("invalid `enc' argument");

    ncon = NextConnection();
#ifdef Win32
    if (CharacterMode != RTerm)
        con = newWpipe(file, strlen(open) ? open : "r");
    else
        con = newpipe(file, strlen(open) ? open : "r");
#else
    con = newpipe(file, strlen(open) ? open : "r");
#endif
    Connections[ncon] = con;
    for (i = 0; i < 256; i++)
        con->encoding[i] = (unsigned char)INTEGER(enc)[i];

    /* open it if desired */
    if (strlen(open))
        con->open(con);

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = ncon;
    PROTECT(class = allocVector(STRSXP, 2));
    SET_STRING_ELT(class, 0, mkChar("pipe"));
#ifdef Win32
    if (CharacterMode != RTerm)
        SET_STRING_ELT(class, 0, mkChar("pipeWin32"));
#endif
    SET_STRING_ELT(class, 1, mkChar("connection"));
    classgets(ans, class);
    UNPROTECT(2);

    return ans;
#else
    error("pipe connections are not available on this system");
    return R_NilValue; /* -Wall */
#endif
}

/* ------------------- gzipped file connections --------------------- */

#ifdef HAVE_LIBZ
#include <zlib.h>

static void gzfile_open(Rconnection con)
{
    gzFile fp;

    fp = gzopen(R_ExpandFileName(con->description), con->mode);
    if (!fp)
        error("cannot open compressed file `%s'", R_ExpandFileName(con->description));
    ((Rgzfileconn)(con->private))->fp = fp;
    con->isopen = TRUE;
    con->canwrite = (con->mode[0] == 'w' || con->mode[0] == 'a');
    con->canread = !con->canwrite;
    con->save = -1000;
}

static void gzfile_close(Rconnection con)
{
    gzclose(((Rgzfileconn)(con->private))->fp);
    con->isopen = FALSE;
}

static void gzfile_destroy(Rconnection con)
{
    free(con->private);
}

#define BUFSIZE 1000
static int gzfile_vfprintf(Rconnection con, const char *format, va_list ap)
{
    gzFile fp = ((Rgzfileconn)(con->private))->fp;
    char buf[BUFSIZE], *b = buf, *vmax = vmaxget();
    int res, usedRalloc = FALSE;

#ifdef HAVE_VSNPRINTF
    res = vsnprintf(buf, BUFSIZE, format, ap);
    if (res >= BUFSIZE)
    { /* res is the desired output length */
        usedRalloc = TRUE;
        b = R_alloc(res + 1, sizeof(char));
        vsprintf(buf, format, ap);
    }
    else if (res < 0)
    { /* just a failure indication */
        usedRalloc = TRUE;
        b = R_alloc(10 * BUFSIZE, sizeof(char));
        res = vsnprintf(buf, 10 * BUFSIZE, format, ap);
        if (res < 0)
        {
            *(b + 10 * BUFSIZE) = '\0';
            warning("printing of extremely long output is truncated");
            res = 10 * BUFSIZE;
        }
    }
    gzwrite(fp, b, res);
#else
    /* allocate a large buffer and hope */
    b = R_alloc(10 * BUFSIZE, sizeof(char));
    res = vsprintf(b, format, ap);
    gzwrite(fp, b, res);
#endif
    if (usedRalloc)
        vmaxset(vmax);
    return res;
}

static int gzfile_fgetc(Rconnection con)
{
    gzFile fp = ((Rgzfileconn)(con->private))->fp;

    /* Looks like eof is signalled one char early */
    if (gzeof(fp))
        return R_EOF;
    return con->encoding[gzgetc(fp)];
}

static long gzfile_seek(Rconnection con, int where, int origin)
{
    gzFile fp = ((Rgzfileconn)(con->private))->fp;
    long pos = gztell(fp);
    int whence = SEEK_SET;

    switch (origin)
    {
    case 2:
        whence = SEEK_CUR;
    case 3:
        whence = SEEK_END;
    default:
        whence = SEEK_SET;
    }
    if (where >= 0)
        gzseek(fp, where, whence);
    return pos;
}

static int gzfile_fflush(Rconnection con)
{
    gzFile fp = ((Rgzfileconn)(con->private))->fp;

    return gzflush(fp, Z_SYNC_FLUSH);
}

static size_t gzfile_read(void *ptr, size_t size, size_t nitems, Rconnection con)
{
    gzFile fp = ((Rgzfileconn)(con->private))->fp;
    return gzread(fp, ptr, size * nitems) / size;
}

static size_t gzfile_write(const void *ptr, size_t size, size_t nitems, Rconnection con)
{
    gzFile fp = ((Rgzfileconn)(con->private))->fp;
    return gzwrite(fp, (const voidp)ptr, size * nitems) / size;
}

static Rconnection newgzfile(char *description, char *mode, int compress)
{
    Rconnection new;
    new = (Rconnection)malloc(sizeof(struct Rconn));
    if (!new)
        error("allocation of file connection failed");
    new->class = (char *)malloc(strlen("gzfile") + 1);
    if (!new->class)
    {
        free(new);
        error("allocation of gzfile connection failed");
    }
    strcpy(new->class, "gzfile");
    new->description = (char *)malloc(strlen(description) + 1);
    if (!new->description)
    {
        free(new->class);
        free(new);
        error("allocation of gzfile connection failed");
    }
    strcpy(new->description, description);
    strncpy(new->mode, mode, 1);
    sprintf(new->mode + 1, "b%1d", compress);

    new->isopen = new->incomplete = FALSE;
    new->canread = new->canwrite = TRUE; /* in principle */
    new->canseek = TRUE;
    new->text = FALSE;
    new->open = &gzfile_open;
    new->close = &gzfile_close;
    new->destroy = &gzfile_destroy;
    new->vfprintf = &gzfile_vfprintf;
    new->fgetc = &gzfile_fgetc;
    new->seek = &gzfile_seek;
    new->truncate = &null_truncate;
    new->fflush = &gzfile_fflush;
    new->read = &gzfile_read;
    new->write = &gzfile_write;
    new->nPushBack = 0;
    new->private = (void *)malloc(sizeof(struct fileconn));
    if (!new->private)
    {
        free(new->description);
        free(new->class);
        free(new);
        error("allocation of gzfile connection failed");
    }
    return new;
}

SEXP do_gzfile(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP sfile, sopen, ans, class, enc;
    char *file, *open;
    int i, ncon, compress;
    Rconnection con = NULL;

    checkArity(op, args);
    sfile = CAR(args);
    if (!isString(sfile) || length(sfile) < 1)
        errorcall(call, "invalid `description' argument");
    if (length(sfile) > 1)
        warning("only first element of `description' argument used");
    file = CHAR(STRING_ELT(sfile, 0));
    sopen = CADR(args);
    if (!isString(sopen) || length(sopen) != 1)
        error("invalid `open' argument");
    enc = CADDR(args);
    if (!isInteger(enc) || length(enc) != 256)
        error("invalid `enc' argument");
    compress = asInteger(CADDDR(args));
    if (compress == NA_LOGICAL || compress < 0 || compress > 9)
        error("invalid `compress' argument");
    open = CHAR(STRING_ELT(sopen, 0));
    ncon = NextConnection();
    con = Connections[ncon] = newgzfile(file, strlen(open) ? open : "r", compress);

    for (i = 0; i < 256; i++)
        con->encoding[i] = (unsigned char)INTEGER(enc)[i];

    /* open it if desired */
    if (strlen(open))
        con->open(con);

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = ncon;
    PROTECT(class = allocVector(STRSXP, 2));
    SET_STRING_ELT(class, 0, mkChar("file"));
    SET_STRING_ELT(class, 1, mkChar("connection"));
    classgets(ans, class);
    UNPROTECT(2);

    return ans;
}
#else
SEXP do_gzfile(SEXP call, SEXP op, SEXP args, SEXP env)
{
    error("zlib is not available on this system");
    return R_NilValue; /* -Wall */
}
#endif

/* ------------------- terminal connections --------------------- */

/* The size of the console buffer */
#define CONSOLE_BUFFER_SIZE 1024

static unsigned char ConsoleBuf[CONSOLE_BUFFER_SIZE];
static unsigned char *ConsoleBufp;
static int ConsoleBufCnt;

static int ConsoleGetchar()
{
    if (--ConsoleBufCnt < 0)
    {
        if (R_ReadConsole("", ConsoleBuf, CONSOLE_BUFFER_SIZE, 0) == 0)
        {
            R_ClearerrConsole();
            return R_EOF;
        }
        R_ParseCnt++;
        ConsoleBufp = ConsoleBuf;
        ConsoleBufCnt = strlen((char *)ConsoleBuf);
        ConsoleBufCnt--;
    }
    return *ConsoleBufp++;
}

static int stdin_fgetc(Rconnection con)
{
    return ConsoleGetchar();
}

static int stdout_vfprintf(Rconnection con, const char *format, va_list ap)
{
    if (R_Outputfile)
        vfprintf(R_Outputfile, format, ap);
    else
        Rcons_vprintf(format, ap);
    return 0;
}

static int stdout_fflush(Rconnection con)
{
    if (R_Outputfile)
        return fflush(R_Outputfile);
    return 0;
}

static int stderr_vfprintf(Rconnection con, const char *format, va_list ap)
{
    REvprintf(format, ap);
    return 0;
}

static int stderr_fflush(Rconnection con)
{
    if (R_Consolefile)
        return fflush(R_Consolefile);
    return 0;
}

static Rconnection newterminal(char *description, char *mode)
{
    Rconnection new;
    new = (Rconnection)malloc(sizeof(struct Rconn));
    if (!new)
        error("allocation of terminal connection failed");
    new->class = (char *)malloc(strlen("terminal") + 1);
    if (!new->class)
    {
        free(new);
        error("allocation of terminal connection failed");
    }
    strcpy(new->class, "terminal");
    new->description = (char *)malloc(strlen(description) + 1);
    if (!new->description)
    {
        free(new->class);
        free(new);
        error("allocation of terminal connection failed");
    }
    strcpy(new->description, description);
    strncpy(new->mode, mode, 4);
    new->mode[4] = '\0';
    new->isopen = TRUE;
    new->incomplete = FALSE;
    new->text = TRUE;
    new->canread = (strcmp(mode, "r") == 0);
    new->canwrite = (strcmp(mode, "w") == 0);
    new->canseek = FALSE;
    new->open = &null_open;
    new->close = &null_open;
    new->destroy = &null_open;
    new->vfprintf = &null_vfprintf;
    new->fgetc = &null_fgetc;
    new->seek = &null_seek;
    new->truncate = &null_truncate;
    new->fflush = &null_fflush;
    new->read = &null_read;
    new->write = &null_write;
    new->nPushBack = 0;
    new->save = -1000;
    new->private = NULL;
    return new;
}

SEXP do_stdin(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans, class;
    Rconnection con = getConnection(0);

    checkArity(op, args);
    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = 0;
    PROTECT(class = allocVector(STRSXP, 2));
    SET_STRING_ELT(class, 0, mkChar(con->class));
    SET_STRING_ELT(class, 1, mkChar("connection"));
    classgets(ans, class);
    UNPROTECT(2);
    return ans;
}

SEXP do_stdout(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans, class;
    Rconnection con = getConnection(R_OutputCon);

    checkArity(op, args);
    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = R_OutputCon;
    PROTECT(class = allocVector(STRSXP, 2));
    SET_STRING_ELT(class, 0, mkChar(con->class));
    SET_STRING_ELT(class, 1, mkChar("connection"));
    classgets(ans, class);
    UNPROTECT(2);
    return ans;
}

SEXP do_stderr(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans, class;
    Rconnection con = getConnection(2);

    checkArity(op, args);
    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = 2;
    PROTECT(class = allocVector(STRSXP, 2));
    SET_STRING_ELT(class, 0, mkChar(con->class));
    SET_STRING_ELT(class, 1, mkChar("connection"));
    classgets(ans, class);
    UNPROTECT(2);
    return ans;
}

/* ------------------- text connections --------------------- */

/* read a R character vector into a buffer */
static void text_init(Rconnection con, SEXP text)
{
    int i, nlines = length(text), nchars = 0;
    Rtextconn this = (Rtextconn)con->private;

    for (i = 0; i < nlines; i++)
        nchars += strlen(CHAR(STRING_ELT(text, i))) + 1;
    this->data = (char *)malloc(nchars + 1);
    if (!this->data)
    {
        free(this);
        free(con->description);
        free(con->class);
        free(con);
        error("cannot allocate memory for text connection");
    }
    *(this->data) = '\0';
    for (i = 0; i < nlines; i++)
    {
        strcat(this->data, CHAR(STRING_ELT(text, i)));
        strcat(this->data, "\n");
    }
    this->nchars = nchars;
    this->cur = this->save = 0;
}

static void text_open(Rconnection con)
{
    con->save = -1000;
}

static void text_close(Rconnection con)
{
}

static void text_destroy(Rconnection con)
{
    Rtextconn this = (Rtextconn)con->private;

    free(this->data);
    this->cur = this->nchars = 0;
}

static int text_fgetc(Rconnection con)
{
    Rtextconn this = (Rtextconn)con->private;
    if (this->save)
    {
        int c;
        c = this->save;
        this->save = 0;
        return c;
    }
    if (this->cur >= this->nchars)
        return R_EOF;
    else
        return (int)(this->data[this->cur++]);
}

static long text_seek(Rconnection con, int where, int origin)
{
    if (where >= 0)
        error("seek is not relevant for text connection");
    return 0; /* if just asking, always at the beginning */
}

static Rconnection newtext(char *description, SEXP text)
{
    Rconnection new;
    new = (Rconnection)malloc(sizeof(struct Rconn));
    if (!new)
        error("allocation of text connection failed");
    new->class = (char *)malloc(strlen("textConnection") + 1);
    if (!new->class)
    {
        free(new);
        error("allocation of text connection failed");
    }
    strcpy(new->class, "textConnection");
    new->description = (char *)malloc(strlen(description) + 1);
    if (!new->description)
    {
        free(new->class);
        free(new);
        error("allocation of text connection failed");
    }
    strcpy(new->description, description);
    strcpy(new->mode, "r");
    new->isopen = new->text = TRUE;
    new->incomplete = FALSE;
    new->canread = TRUE;
    new->canwrite = FALSE;
    new->canseek = FALSE;
    new->open = &text_open;
    new->close = &text_close;
    new->destroy = &text_destroy;
    new->vfprintf = &null_vfprintf;
    new->fgetc = &text_fgetc;
    new->seek = &text_seek;
    new->truncate = &null_truncate;
    new->fflush = &null_fflush;
    new->read = &null_read;
    new->write = &null_write;
    new->nPushBack = 0;
    new->private = (void *)malloc(sizeof(struct textconn));
    new->save = -1000;
    if (!new->private)
    {
        free(new->description);
        free(new->class);
        free(new);
        error("allocation of text connection failed");
    }
    text_init(new, text);
    return new;
}

static void outtext_close(Rconnection con)
{
    Routtextconn this = (Routtextconn)con->private;
    SEXP tmp;

    if (strlen(this->lastline) > 0)
    {
        PROTECT(tmp = lengthgets(this->data, ++this->len));
        SET_STRING_ELT(tmp, this->len - 1, mkChar(this->lastline));
        defineVar(this->namesymbol, tmp, R_GlobalEnv);
        this->data = tmp;
        UNPROTECT(1);
    }
}

static void outtext_destroy(Rconnection con)
{
}

#define BUFSIZE 1000
static int text_vfprintf(Rconnection con, const char *format, va_list ap)
{
    Routtextconn this = (Routtextconn)con->private;
    char buf[BUFSIZE], *b = buf, *p, *q, *vmax = vmaxget();
    int res = 0, usedRalloc = FALSE, buffree, already = strlen(this->lastline);
    SEXP tmp;

    strcpy(b, this->lastline);
    p = b + already;
    buffree = BUFSIZE - already;

#ifdef HAVE_VSNPRINTF
    res = vsnprintf(p, buffree, format, ap);
    if (res >= buffree)
    { /* res is the desired output length */
        usedRalloc = TRUE;
        b = R_alloc(res + already + 1, sizeof(char));
        strcpy(b, this->lastline);
        p = b + already;
        vsprintf(p, format, ap);
    }
    else if (res < 0)
    { /* just a failure indication */
        usedRalloc = TRUE;
        b = R_alloc(10 * BUFSIZE, sizeof(char));
        strcpy(b, this->lastline);
        p = b + already;
        res = vsnprintf(p, 10 * BUFSIZE - already, format, ap);
        if (res < 0)
        {
            *(b + 10 * BUFSIZE) = '\0';
            warning("printing of extremely long output is truncated");
        }
    }
#else
    /* allocate a large buffer and hope */
    b = R_alloc(10 * BUFSIZE, sizeof(char));
    strcpy(b, this->lastline);
    p = b + already;
    res = vsprintf(p, format, ap);
#endif

    /* copy buf line-by-line to object */
    for (p = buf;; p = q + 1)
    {
        q = strchr(p, '\n');
        if (q)
        {
            *q = '\0';
            PROTECT(tmp = lengthgets(this->data, ++this->len));
            SET_STRING_ELT(tmp, this->len - 1, mkChar(p));
            defineVar(this->namesymbol, tmp, R_GlobalEnv);
            this->data = tmp;
            UNPROTECT(1);
        }
        else
        {
            /* retain the last line */
            if (strlen(this->lastline) < LAST_LINE_LEN)
            {
                strcpy(this->lastline, p);
            }
            else
            {
                strncpy(this->lastline, p, LAST_LINE_LEN - 1);
                this->lastline[LAST_LINE_LEN - 1] = '\0';
                warning("line truncated in output text connection");
            }
            con->incomplete = strlen(this->lastline) > 0;
            break;
        }
    }
    if (usedRalloc)
        vmaxset(vmax);
    return res;
}

static void outtext_init(Rconnection con, char *mode)
{
    Routtextconn this = (Routtextconn)con->private;
    SEXP val;

    this->namesymbol = install(con->description);
    if (strcmp(mode, "w") == 0)
    {
        /* create variable pointed to by con->description */
        PROTECT(val = allocVector(STRSXP, 0));
        defineVar(this->namesymbol, val, R_GlobalEnv);
        UNPROTECT(1);
    }
    else
    {
        /* take over existing variable */
        val = findVar1(this->namesymbol, R_GlobalEnv, STRSXP, FALSE);
        if (val == R_UnboundValue)
        {
            warning("text connection: appending to a non-existent char vector");
            PROTECT(val = allocVector(STRSXP, 0));
            defineVar(this->namesymbol, val, R_GlobalEnv);
            UNPROTECT(1);
        }
    }
    this->len = LENGTH(val);
    this->data = val;
    this->lastline[0] = '\0';
}

static Rconnection newouttext(char *description, SEXP sfile, char *mode)
{
    Rconnection new;
    new = (Rconnection)malloc(sizeof(struct Rconn));
    if (!new)
        error("allocation of text connection failed");
    new->class = (char *)malloc(strlen("textConnection") + 1);
    if (!new->class)
    {
        free(new);
        error("allocation of text connection failed");
    }
    strcpy(new->class, "textConnection");
    new->description = (char *)malloc(strlen(description) + 1);
    if (!new->description)
    {
        free(new->class);
        free(new);
        error("allocation of text connection failed");
    }
    strcpy(new->description, description);
    strcpy(new->mode, mode); /* must be "w" or "a" at this point */
    new->isopen = new->text = TRUE;
    new->incomplete = FALSE;
    new->canread = FALSE;
    new->canwrite = TRUE;
    new->canseek = FALSE;
    new->open = &text_open;
    new->close = &outtext_close;
    new->destroy = &outtext_destroy;
    new->vfprintf = &text_vfprintf;
    new->fgetc = &null_fgetc;
    new->seek = &text_seek;
    new->truncate = &null_truncate;
    new->fflush = &null_fflush;
    new->read = &null_read;
    new->write = &null_write;
    new->nPushBack = 0;
    new->private = (void *)malloc(sizeof(struct outtextconn));
    new->save = -1000;
    if (!new->private)
    {
        free(new->description);
        free(new->class);
        free(new);
        error("allocation of text connection failed");
    }
    outtext_init(new, mode);
    return new;
}

SEXP do_textconnection(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP sfile, stext, sopen, ans, class;
    char *desc, *open;
    int ncon;
    Rconnection con = NULL;

    checkArity(op, args);
    sfile = CAR(args);
    if (!isString(sfile) || length(sfile) != 1)
        error("invalid `description' argument");
    desc = CHAR(STRING_ELT(sfile, 0));
    stext = CADR(args);
    if (!isString(stext))
        error("invalid `text' argument");
    sopen = CADDR(args);
    if (!isString(sopen) || length(sopen) != 1)
        error("invalid `open' argument");
    open = CHAR(STRING_ELT(sopen, 0));
    ncon = NextConnection();
    if (!strlen(open) || strncmp(open, "r", 1) == 0)
        con = Connections[ncon] = newtext(desc, stext);
    else if (strncmp(open, "w", 1) == 0 || strncmp(open, "a", 1) == 0)
        con = Connections[ncon] = newouttext(CHAR(STRING_ELT(stext, 0)), sfile, open);
    else
        errorcall(call, "unsupported mode");
    /* already opened */

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = ncon;
    PROTECT(class = allocVector(STRSXP, 2));
    SET_STRING_ELT(class, 0, mkChar("textConnection"));
    SET_STRING_ELT(class, 1, mkChar("connection"));
    classgets(ans, class);
    UNPROTECT(2);
    return ans;
}

/* ------------------- open, close, seek --------------------- */

SEXP do_open(SEXP call, SEXP op, SEXP args, SEXP env)
{
    int i, block;
    Rconnection con = NULL;
    SEXP sopen;
    char *open;

    checkArity(op, args);
    i = asInteger(CAR(args));
    con = getConnection(i);
    if (i < 3)
        error("cannot open standard connections");
    sopen = CADR(args);
    if (!isString(sopen) || length(sopen) != 1)
        error("invalid `open' argument");
    block = asLogical(CADDR(args));
    if (block == NA_LOGICAL)
        error("invalid `block' argument");
    open = CHAR(STRING_ELT(sopen, 0));
    if (strlen(open) > 0)
        strcpy(con->mode, open);
    con->blocking = block;
    con->open(con);
    return R_NilValue;
}

SEXP do_isopen(SEXP call, SEXP op, SEXP args, SEXP env)
{
    Rconnection con;
    SEXP ans;

    checkArity(op, args);
    con = getConnection(asInteger(CAR(args)));
    PROTECT(ans = allocVector(LGLSXP, 1));
    LOGICAL(ans)[0] = con->isopen != FALSE;
    UNPROTECT(1);
    return ans;
}

SEXP do_isincomplete(SEXP call, SEXP op, SEXP args, SEXP env)
{
    Rconnection con;
    SEXP ans;

    checkArity(op, args);
    con = getConnection(asInteger(CAR(args)));
    PROTECT(ans = allocVector(LGLSXP, 1));
    LOGICAL(ans)[0] = con->incomplete != FALSE;
    UNPROTECT(1);
    return ans;
}

SEXP do_isseekable(SEXP call, SEXP op, SEXP args, SEXP env)
{
    Rconnection con;
    SEXP ans;

    checkArity(op, args);
    con = getConnection(asInteger(CAR(args)));
    PROTECT(ans = allocVector(LGLSXP, 1));
    LOGICAL(ans)[0] = con->canseek != FALSE;
    UNPROTECT(1);
    return ans;
}

void con_close(int i)
{
    int j;
    Rconnection con = NULL;

    con = getConnection(i);
    if (con->isopen)
        con->close(con);
    con->destroy(con);
    free(con->class);
    free(con->description);
    /* clear the pushBack */
    if (con->nPushBack > 0)
    {
        for (j = 0; j < con->nPushBack; j++)
            free(con->PushBack[j]);
        free(con->PushBack);
    }
    free(Connections[i]);
    Connections[i] = NULL;
}

SEXP do_close(SEXP call, SEXP op, SEXP args, SEXP env)
{
    int i, j;

    checkArity(op, args);
    i = asInteger(CAR(args));
    if (i < 3)
        error("cannot close standard connections");
    for (j = 0; j < R_SinkNumber; j++)
        if (i == SinkCons[j])
            error("cannot close output sink connection");
    if (i == R_ErrorCon)
        error("cannot close messages sink connection");
    con_close(i);
    return R_NilValue;
}

/* seek(con, where = numeric(), origin = "start") */
SEXP do_seek(SEXP call, SEXP op, SEXP args, SEXP env)
{
    int where, origin;
    SEXP ans;
    Rconnection con = NULL;

    checkArity(op, args);
    con = getConnection(asInteger(CAR(args)));
    where = asInteger(CADR(args));
    origin = asInteger(CADDR(args));
    if (where == NA_INTEGER || where < 0)
        where = -1;
    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = con->seek(con, where, origin);
    UNPROTECT(1);
    return ans;
}

/* truncate(con) */
SEXP do_truncate(SEXP call, SEXP op, SEXP args, SEXP env)
{
    Rconnection con = NULL;

    checkArity(op, args);
    con = getConnection(asInteger(CAR(args)));
    con->truncate(con);
    return R_NilValue;
}

/* ------------------- read, write  text --------------------- */

int Rconn_fgetc(Rconnection con)
{
    char *curLine;
    int c;

    if (con->nPushBack <= 0)
    {
        /* map CR or CRLF to LF */
        if (con->save != -1000)
        {
            c = con->save;
            con->save = con->save2;
            con->save2 = -1000;
            return c;
        }
        c = con->fgetc(con);
        if (c == '\r')
        {
            c = con->fgetc(con);
            if (c != '\n')
            {
                con->save2 = con->save;
                con->save = (c != '\r') ? c : '\n';
                return ('\n');
            }
        }
        return c;
    }
    curLine = con->PushBack[con->nPushBack - 1];
    c = curLine[con->posPushBack++];
    if (con->posPushBack >= strlen(curLine))
    {
        /* last character on a line, so pop the line */
        free(curLine);
        con->nPushBack--;
        con->posPushBack = 0;
        if (con->nPushBack == 0)
            free(con->PushBack);
    }
    return c;
}

int Rconn_ungetc(int c, Rconnection con)
{
    con->save2 = con->save;
    con->save = c;
    return c;
}

int Rconn_printf(Rconnection con, const char *format, ...)
{
    int res;
    va_list(ap);

    va_start(ap, format);
    res = con->vfprintf(con, format, ap);
    va_end(ap);
    return res;
}

/* readLines(con = stdin(), n = 1, ok = TRUE) */
#define BUF_SIZE 1000
SEXP do_readLines(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans = R_NilValue, ans2;
    int i, n, nn, nnn, ok, nread, c, nbuf, buf_size = BUF_SIZE;
    Rconnection con = NULL;
    Rboolean wasopen;
    char *buf;

    checkArity(op, args);
    con = getConnection(asInteger(CAR(args)));
    n = asInteger(CADR(args));
    if (n == NA_INTEGER)
        errorcall(call, "invalid value for `n'");
    ok = asLogical(CADDR(args));
    if (ok == NA_LOGICAL)
        errorcall(call, "invalid value for `ok'");
    if (!con->canread)
        errorcall(call, "cannot read from this connection");
    wasopen = con->isopen;
    if (!wasopen)
        con->open(con);
    buf = (char *)malloc(buf_size);
    if (!buf)
        error("cannot allocate buffer in readLines");
    nn = (n < 0) ? 1000 : n; /* initially allocate space for 1000 lines */
    nnn = (n < 0) ? INT_MAX : n;
    PROTECT(ans = allocVector(STRSXP, nn));
    for (nread = 0; nread < nnn; nread++)
    {
        if (nread >= nn)
        {
            ans2 = allocVector(STRSXP, 2 * nn);
            for (i = 0; i < nn; i++)
                SET_STRING_ELT(ans2, i, STRING_ELT(ans, i));
            nn *= 2;
            UNPROTECT(1); /* old ans */
            PROTECT(ans = ans2);
        }
        nbuf = 0;
        while ((c = Rconn_fgetc(con)) != R_EOF)
        {
            if (nbuf == buf_size)
            {
                buf_size *= 2;
                buf = (char *)realloc(buf, buf_size);
                if (!buf)
                    error("cannot allocate buffer in readLines");
            }
            if (c != '\n')
                buf[nbuf++] = c;
            else
                break;
        }
        buf[nbuf] = '\0';
        SET_STRING_ELT(ans, nread, mkChar(buf));
        if (c == R_EOF)
            goto no_more_lines;
    }
    UNPROTECT(1);
    free(buf);
    if (!wasopen)
        con->close(con);
    return ans;
no_more_lines:
    free(buf);
    if (!wasopen)
        con->close(con);
    if (nbuf > 0)
    { /* incomplete last line */
        nread++;
        warningcall(call, "incomplete final line");
    }
    if (n < nnn && !ok)
        errorcall(call, "too few lines read");
    PROTECT(ans2 = allocVector(STRSXP, nread));
    for (i = 0; i < nread; i++)
        SET_STRING_ELT(ans2, i, STRING_ELT(ans, i));
    UNPROTECT(2);
    return ans2;
}

static void writecon(Rconnection con, char *format, ...)
{
    va_list(ap);
    va_start(ap, format);
    con->vfprintf(con, format, ap);
    va_end(ap);
}

/* writelines(text, con = stdout(), sep = "\n") */
SEXP do_writelines(SEXP call, SEXP op, SEXP args, SEXP env)
{
    int i;
    Rboolean wasopen;
    Rconnection con = NULL;
    SEXP text, sep;

    checkArity(op, args);
    text = CAR(args);
    if (!isString(text))
        error("invalid `text' argument");
    con = getConnection(asInteger(CADR(args)));
    sep = CADDR(args);
    if (!isString(sep))
        error("invalid `sep' argument");
    if (!con->canwrite)
        error("cannot write to this connection");
    wasopen = con->isopen;
    if (!wasopen)
        con->open(con);
    for (i = 0; i < length(text); i++)
        writecon(con, "%s%s", CHAR(STRING_ELT(text, i)), CHAR(STRING_ELT(sep, 0)));
    if (!wasopen)
        con->close(con);
    return R_NilValue;
}

/* ------------------- read, write  binary --------------------- */

static void swapb(void *result, int size)
{
    int i;
    char *p = result, tmp;

    if (size == 1)
        return;
    for (i = 0; i < size / 2; i++)
    {
        tmp = p[i];
        p[i] = p[size - i - 1];
        p[size - i - 1] = tmp;
    }
}

static SEXP readOneString(Rconnection con)
{
    char *buf, *p, *new;
    int ibfs = 500, pos, m;

    buf = (char *)R_alloc(ibfs, sizeof(char));
    for (pos = 0; pos < 10000; pos++)
    {
        p = buf + pos;
        m = con->read(p, sizeof(char), 1, con);
        if (!m)
            return R_NilValue;
        if (*p == '\0')
            break;
        if (pos >= ibfs - 1)
        {
            new = (char *)R_alloc(2 * ibfs, sizeof(char));
            memcpy(new, buf, pos + 1);
            buf = new;
            ibfs *= 2;
        }
    }
    return mkChar(buf);
}

/* readBin(con, what, n, swap) */
SEXP do_readbin(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans = R_NilValue, swhat;
    int i, size, swap, n, m = 0, sizedef = 4, mode = 1;
    char *what;
    void *p = NULL;
    Rboolean wasopen;
    Rconnection con = NULL;

    checkArity(op, args);
    i = asInteger(CAR(args));
    args = CDR(args);
    if (i == NA_INTEGER || !(con = Connections[i]))
        error("invalid connection");
    if (con->text)
        error("can only read from a binary connection");
    swhat = CAR(args);
    args = CDR(args);
    if (!isString(swhat) || length(swhat) != 1)
        error("invalid value of `what'");
    what = CHAR(STRING_ELT(swhat, 0));
    n = asInteger(CAR(args));
    args = CDR(args);
    if (n == NA_INTEGER || n < 0)
        error("invalid value of `n'");
    size = asInteger(CAR(args));
    args = CDR(args);
    swap = asLogical(CAR(args));
    if (swap == NA_LOGICAL)
        error("invalid value of `swap'");
    if (!con->canread)
        error("cannot read from this connection");

    wasopen = con->isopen;
    if (!wasopen)
        con->open(con);

    if (!strcmp(what, "character"))
    {
        SEXP onechar;
        PROTECT(ans = allocVector(STRSXP, n));
        for (i = 0, m = i + 1; i < n; i++)
        {
            onechar = readOneString(con);
            if (onechar != R_NilValue)
            {
                SET_STRING_ELT(ans, i, onechar);
                m++;
            }
            else
                break;
        }
    }
    else if (!strcmp(what, "complex"))
    {
        if (size == NA_INTEGER)
            size = sizeof(Rcomplex);
        if (size != sizeof(Rcomplex))
            error("size changing is not supported for complex vectors");
        PROTECT(ans = allocVector(CPLXSXP, n));
        p = (void *)COMPLEX(ans);
        m = con->read(p, size, n, con);
        if (swap)
            for (i = 0; i < m; i++)
            {
                swapb(&(COMPLEX(ans)[i].r), sizeof(double));
                swapb(&(COMPLEX(ans)[i].i), sizeof(double));
            }
    }
    else
    {
        if (!strcmp(what, "integer") || !strcmp(what, "int"))
        {
            sizedef = sizeof(int);
            mode = 1;
            if (size == NA_INTEGER)
                size = sizedef;
            switch (size)
            {
            case sizeof(signed char):
            case sizeof(short):
            case sizeof(int):
#if SIZEOF_LONG == 8
            case sizeof(long):
#elif SIZEOF_LONG_LONG == 8
            case sizeof(long long):
#endif
                break;
            default:
                error("That size is unknown on this machine");
            }
            PROTECT(ans = allocVector(INTSXP, n));
            p = (void *)INTEGER(ans);
        }
        else if (!strcmp(what, "logical"))
        {
            sizedef = sizeof(int);
            mode = 1;
            if (size == NA_INTEGER)
                size = sizedef;
            switch (size)
            {
            case sizeof(signed char):
            case sizeof(short):
            case sizeof(int):
#if SIZEOF_LONG == 8
            case sizeof(long):
#elif SIZEOF_LONG_LONG == 8
            case sizeof(long long):
#endif
                break;
            default:
                error("That size is unknown on this machine");
            }
            PROTECT(ans = allocVector(LGLSXP, n));
            p = (void *)LOGICAL(ans);
        }
        else if (!strcmp(what, "numeric") || !strcmp(what, "double"))
        {
            sizedef = sizeof(double);
            mode = 2;
            if (size == NA_INTEGER)
                size = sizedef;
            switch (size)
            {
            case sizeof(double):
            case sizeof(float):
#if SIZEOF_LONG_DOUBLE > 8
            case sizeof(long double):
#endif
                break;
            default:
                error("That size is unknown on this machine");
            }
            PROTECT(ans = allocVector(REALSXP, n));
            p = (void *)REAL(ans);
        }
        if (size == sizedef)
        {
            m = con->read(p, size, n, con);
            if (swap && size > 1)
                for (i = 0; i < m; i++)
                    swapb((char *)p + i * size, size);
        }
        else
        {
            char *buf = R_alloc(1, size);
            int s;
            if (mode == 1)
            {
                for (i = 0, m = 0; i < n; i++)
                {
                    s = con->read(buf, size, 1, con);
                    if (s)
                        m++;
                    else
                        break;
                    if (swap && size > 1)
                        swapb(buf, size);
                    switch (size)
                    {
                    case sizeof(signed char):
                        INTEGER(ans)[i] = (int)*((signed char *)buf);
                        break;
                    case sizeof(short):
                        INTEGER(ans)[i] = (int)*((short *)buf);
                        break;
#if SIZEOF_LONG == 8
                    case sizeof(long):
                        INTEGER(ans)[i] = (int)*((long *)buf);
                        break;
#elif SIZEOF_LONG_LONG == 8
                    case sizeof(long long):
                        INTEGER(ans)[i] = (int)*((long long *)buf);
                        break;
#endif
                    }
                }
            }
            else if (mode == 2)
            {
                for (i = 0, m = 0; i < n; i++)
                {
                    s = con->read(buf, size, 1, con);
                    if (s)
                        m++;
                    else
                        break;
                    if (swap && size > 1)
                        swapb(buf, size);
                    switch (size)
                    {
                    case sizeof(float):
                        REAL(ans)[i] = (double)*((float *)buf);
                        break;
#if SIZEOF_LONG_DOUBLE > 8
                    case sizeof(long double):
                        REAL(ans)[i] = (double)*((long double *)buf);
                        break;
#endif
                    }
                }
            }
        }
    }
    if (!wasopen)
        con->close(con);
    if (m < n)
    {
        PROTECT(ans = lengthgets(ans, m));
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return ans;
}

/* writeBin(object, con, swap) */
SEXP do_writebin(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP object;
    int i, j, size, swap, len, n = 0;
    char *s, *buf;
    Rboolean wasopen;
    Rconnection con = NULL;

    checkArity(op, args);
    object = CAR(args);
    if (!isVectorAtomic(object))
        error("`x' is not an atomic vector type");
    i = asInteger(CADR(args));
    if (i == NA_INTEGER || !(con = Connections[i]))
        error("invalid connection");
    if (con->text)
        error("can only write to a binary connection");
    size = asInteger(CADDR(args));
    swap = asLogical(CADDDR(args));
    if (swap == NA_LOGICAL)
        error("invalid value of `swap'");
    if (!con->canwrite)
        error("cannot write to this connection");
    len = LENGTH(object);
    if (len == 0)
        return R_NilValue;

    wasopen = con->isopen;
    if (!wasopen)
        con->open(con);

    if (TYPEOF(object) == STRSXP)
    {
        for (i = 0; i < len; i++)
        {
            s = CHAR(STRING_ELT(object, i));
            n = con->write(s, sizeof(char), strlen(s) + 1, con);
            if (!n)
            {
                warning("problem writing to connection");
                break;
            }
        }
    }
    else
    {
        switch (TYPEOF(object))
        {
        case LGLSXP:
        case INTSXP:
            if (size == NA_INTEGER)
                size = sizeof(int);
            switch (size)
            {
            case sizeof(signed char):
            case sizeof(short):
            case sizeof(int):
#if SIZEOF_LONG == 8
            case sizeof(long):
#elif SIZEOF_LONG_LONG == 8
            case sizeof(long long):
#endif
                break;
            default:
                error("That size is unknown on this machine");
            }
            break;
        case REALSXP:
            if (size == NA_INTEGER)
                size = sizeof(double);
            switch (size)
            {
            case sizeof(double):
            case sizeof(float):
#if SIZEOF_LONG_DOUBLE > 8
            case sizeof(long double):
#endif
                break;
            default:
                error("That size is unknown on this machine");
            }
            break;
        case CPLXSXP:
            if (size == NA_INTEGER)
                size = sizeof(Rcomplex);
            if (size != sizeof(Rcomplex))
                error("size changing is not supported for complex vectors");
            break;
        default:
            error("That type is unimplemented");
        }
        buf = R_alloc(len, size);
        switch (TYPEOF(object))
        {
        case LGLSXP:
        case INTSXP:
            switch (size)
            {
            case sizeof(int):
                memcpy(buf, INTEGER(object), size * len);
                break;
#if SIZEOF_LONG == 8
            case sizeof(long): {
                long l1;
                for (i = 0, j = 0; i < len; i++, j += size)
                {
                    l1 = (long)INTEGER(object)[i];
                    memcpy(buf + j, &l1, size);
                }
                break;
            }
#elif SIZEOF_LONG_LONG == 8
            case sizeof(long long): {
                long long ll1;
                for (i = 0, j = 0; i < len; i++, j += size)
                {
                    ll1 = (long long)INTEGER(object)[i];
                    memcpy(buf + j, &ll1, size);
                }
                break;
            }
#endif
            case 2: {
                short s1;
                for (i = 0, j = 0; i < len; i++, j += size)
                {
                    s1 = (short)INTEGER(object)[i];
                    memcpy(buf + j, &s1, size);
                }
                break;
            }
            case 1:
                for (i = 0; i < len; i++)
                    buf[i] = (signed char)INTEGER(object)[i];
                break;
            }
            break;
        case REALSXP:
            switch (size)
            {
            case sizeof(double):
                memcpy(buf, REAL(object), size * len);
                break;
            case sizeof(float): {
                float f1;
                for (i = 0, j = 0; i < len; i++, j += size)
                {
                    f1 = (float)REAL(object)[i];
                    memcpy(buf + j, &f1, size);
                }
                break;
            }
#if SIZEOF_LONG_DOUBLE > 8
            case sizeof(long double): {
                long double ld1;
                for (i = 0, j = 0; i < len; i++, j += size)
                {
                    ld1 = (long double)REAL(object)[i];
                    memcpy(buf + j, &ld1, size);
                }
                break;
            }
#endif
            }
            break;
        case CPLXSXP:
            memcpy(buf, COMPLEX(object), size * len);
            break;
        }

        if (swap && size > 1)
            for (i = 0; i < len; i++)
                swapb(buf + size * i, size);

        /* write it now */
        n = con->write(buf, size, len, con);
        if (n < len)
            warning("problem writing to connection");
    }

    if (!wasopen)
        con->close(con);
    return R_NilValue;
}

static SEXP readFixedString(Rconnection con, int len)
{
    char *buf, *p;
    int pos, m;

    buf = (char *)R_alloc(len + 1, sizeof(char));
    buf[len] = '\0';
    for (pos = 0; pos < len; pos++)
    {
        p = buf + pos;
        m = con->read(p, sizeof(char), 1, con);
        if (!m)
            return R_NilValue;
    }
    return mkChar(buf);
}

/* readChar(con, nchars) */
SEXP do_readchar(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans = R_NilValue, onechar, nchars;
    int i, len, n, m = 0;
    Rboolean wasopen;
    Rconnection con = NULL;

    checkArity(op, args);
    i = asInteger(CAR(args));
    if (i == NA_INTEGER || !(con = Connections[i]))
        error("invalid connection");
    if (!con->canread)
        error("cannot read from this connection");
    nchars = CADR(args);
    n = LENGTH(nchars);
    if (n == 0)
        return allocVector(STRSXP, 0);

    wasopen = con->isopen;
    if (!wasopen)
        con->open(con);

    PROTECT(ans = allocVector(STRSXP, n));
    for (i = 0, m = i + 1; i < n; i++)
    {
        len = INTEGER(nchars)[i];
        if (len == NA_INTEGER || len < 0)
            error("supplied length is invalid");
        onechar = readFixedString(con, len);
        if (onechar != R_NilValue)
        {
            SET_STRING_ELT(ans, i, onechar);
            m++;
        }
        else
            break;
    }
    if (!wasopen)
        con->close(con);
    if (m < n)
    {
        PROTECT(ans = lengthgets(ans, m));
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return ans;
}

/* writeChar(object, con, nchars, sep) */
SEXP do_writechar(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP object, nchars, sep;
    int i, len, n, nwrite = 0, slen, tlen;
    char *s, *buf, *ssep = "";
    Rboolean wasopen, usesep;
    Rconnection con = NULL;

    checkArity(op, args);
    object = CAR(args);
    i = asInteger(CADR(args));
    if (i == NA_INTEGER || !(con = Connections[i]))
        error("invalid connection");
    if (!con->canwrite)
        error("cannot write to this connection");

    nchars = CADDR(args);
    sep = CADDDR(args);
    if (isNull(sep))
    {
        usesep = FALSE;
        slen = 0;
    }
    else
    {
        usesep = TRUE;
        if (!isString(sep) || length(sep) != 1)
            error("invalid value of `sep'");
        ssep = CHAR(STRING_ELT(sep, 0));
        slen = strlen(ssep) + 1;
    }
    n = LENGTH(nchars);
    if (n == 0)
        return R_NilValue;

    len = 0;
    for (i = 0; i < n; i++)
    {
        tlen = strlen(CHAR(STRING_ELT(object, i)));
        if (tlen > len)
            len = tlen;
    }
    buf = (char *)R_alloc(len + slen, sizeof(char));

    wasopen = con->isopen;
    if (!wasopen)
        con->open(con);

    if (TYPEOF(object) == STRSXP)
    {
        for (i = 0; i < n; i++)
        {
            len = INTEGER(nchars)[i];
            s = CHAR(STRING_ELT(object, i));
            memset(buf, '\0', len + slen);
            strncpy(buf, s, len);
            if (usesep)
            {
                strcat(buf, ssep);
                len += slen;
            }
            nwrite = con->write(buf, sizeof(char), len, con);
            if (!nwrite)
            {
                warning("problem writing to connection");
                break;
            }
        }
    }
    if (!wasopen)
        con->close(con);
    return R_NilValue;
}

/* ------------------- push back text  --------------------- */

SEXP do_pushback(SEXP call, SEXP op, SEXP args, SEXP env)
{
    int i, n, nexists, newLine;
    Rconnection con = NULL;
    SEXP stext;
    char *p, **q;

    checkArity(op, args);

    stext = CAR(args);
    if (!isString(stext))
        error("invalid `data' argument");
    i = asInteger(CADR(args));
    if (i == NA_INTEGER || !(con = Connections[i]))
        error("invalid connection");
    newLine = asLogical(CADDR(args));
    if (newLine == NA_LOGICAL)
        error("invalid `newLine' argument");
    if (!con->canread && !con->isopen)
        error("can only push back on open readable connections");
    if (!con->text)
        error("can only push back on text-mode connections");
    nexists = con->nPushBack;
    if ((n = length(stext)) > 0)
    {
        if (nexists > 0)
        {
            q = con->PushBack = (char **)realloc(con->PushBack, (n + nexists) * sizeof(char *));
        }
        else
        {
            q = con->PushBack = (char **)malloc(n * sizeof(char *));
        }
        if (!q)
            error("could not allocate space for pushBack");
        for (i = 0; i < n; i++)
        {
            p = CHAR(STRING_ELT(stext, n - i - 1));
            q += nexists + i;
            *q = (char *)malloc(strlen(p) + 1 + newLine);
            if (!(*q))
                error("could not allocate space for pushBack");
            strcpy(*q, p);
            if (newLine)
                strcat(*q, "\n");
        }
        con->posPushBack = 0;
        con->nPushBack += n;
    }
    return R_NilValue;
}

SEXP do_pushbacklength(SEXP call, SEXP op, SEXP args, SEXP env)
{
    int i;
    Rconnection con = NULL;
    SEXP ans;

    i = asInteger(CAR(args));
    if (i == NA_INTEGER || !(con = Connections[i]))
        error("invalid connection");
    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = con->nPushBack;
    UNPROTECT(1);
    return ans;
}

/* ------------------- sink functions  --------------------- */

/* Switch output to connection number icon, or popd stack if icon < 0
 */
Rboolean switch_stdout(int icon, int closeOnExit)
{
    int toclose;

    if (icon == R_OutputCon)
        return FALSE;

    if (icon >= 0 && R_SinkNumber >= NSINKS - 1)
        error("sink stack is full");

    if (icon == 0)
        error("cannot switch output to stdin");
    else if (icon == 1)
    {
        R_OutputCon = SinkCons[++R_SinkNumber] = 1;
        SinkConsClose[R_SinkNumber] = 0;
    }
    else if (icon == 2)
        error("cannot switch output to stderr");
    else if (icon >= 3)
    {
        Rconnection con = getConnection(icon); /* checks validity */
        toclose = 2 * closeOnExit;
        if (!con->isopen)
        {
            con->open(con);
            toclose = 1;
        }
        R_OutputCon = SinkCons[++R_SinkNumber] = icon;
        SinkConsClose[R_SinkNumber] = toclose;
    }
    else
    { /* removing a sink */
        if (R_SinkNumber <= 0)
        {
            warning("no sink to remove");
            return FALSE;
        }
        else
        {
            R_SinkNumber--;
            if ((icon = SinkCons[R_SinkNumber + 1]) >= 3)
            {
                Rconnection con = getConnection(icon);
                if (SinkConsClose[R_SinkNumber + 1] == 1) /* close it */
                    con->close(con);
                else if (SinkConsClose[R_SinkNumber + 1] == 2) /* destroy it */
                    con_close(icon);
            }
            R_OutputCon = SinkCons[R_SinkNumber];
        }
    }
    return TRUE;
}

SEXP do_sink(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    int icon, closeOnExit, errcon;

    checkArity(op, args);
    icon = asInteger(CAR(args));
    closeOnExit = asLogical(CADR(args));
    if (closeOnExit == NA_LOGICAL)
        error("invalid value for closeOnExit");
    errcon = asLogical(CADDR(args));
    if (errcon == NA_LOGICAL)
        error("invalid value for type");

    if (!errcon)
    {
        /* allow space for cat() to use sink() */
        if (icon >= 0 && R_SinkNumber >= NSINKS - 2)
            error("sink stack is full");
        switch_stdout(icon, closeOnExit);
        R_SinkSaved = R_SinkNumber;
    }
    else
    {
        if (icon < 0)
            R_ErrorCon = 2;
        else
        {
            getConnection(icon); /* check validity */
            R_ErrorCon = icon;
        }
    }

    return R_NilValue;
}

SEXP do_sinknumber(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP ans;
    int errcon;
    checkArity(op, args);

    errcon = asLogical(CAR(args));
    if (errcon == NA_LOGICAL)
        error("invalid value for type");
    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = errcon ? R_SinkNumber : R_ErrorCon;
    UNPROTECT(1);
    return ans;
}

void R_SinkReset()
{
    R_SinkNumber = R_SinkSaved;
}

/* ------------------- admin functions  --------------------- */

void InitConnections()
{
    int i;
    Connections[0] = newterminal("stdin", "r");
    Connections[0]->fgetc = stdin_fgetc;
    Connections[1] = newterminal("stdout", "w");
    Connections[1]->vfprintf = stdout_vfprintf;
    Connections[1]->fflush = stdout_fflush;
    Connections[2] = newterminal("stderr", "w");
    Connections[2]->vfprintf = stderr_vfprintf;
    Connections[2]->fflush = stderr_fflush;
    for (i = 3; i < NCONNECTIONS; i++)
        Connections[i] = NULL;
    R_OutputCon = 1;
    R_SinkSaved = R_SinkNumber = 0;
    SinkCons[0] = 1;
    R_ErrorCon = 2;
}

SEXP do_getallconnections(SEXP call, SEXP op, SEXP args, SEXP env)
{
    int i, j = 0, n = 0;
    SEXP ans;
    checkArity(op, args);
    for (i = 0; i < NCONNECTIONS; i++)
        if (Connections[i])
            n++;
    PROTECT(ans = allocVector(INTSXP, n));
    for (i = 0; i < NCONNECTIONS; i++)
        if (Connections[i])
            INTEGER(ans)[j++] = i;
    UNPROTECT(1);
    return ans;
}

SEXP do_sumconnection(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans, names;
    Rconnection Rcon;

    checkArity(op, args);
    Rcon = getConnection(asInteger(CAR(args)));
    PROTECT(ans = allocVector(VECSXP, 7));
    PROTECT(names = allocVector(STRSXP, 7));
    SET_STRING_ELT(names, 0, mkChar("description"));
    SET_VECTOR_ELT(ans, 0, mkString(Rcon->description));
    SET_STRING_ELT(names, 1, mkChar("class"));
    SET_VECTOR_ELT(ans, 1, mkString(Rcon->class));
    SET_STRING_ELT(names, 2, mkChar("mode"));
    SET_VECTOR_ELT(ans, 2, mkString(Rcon->mode));
    SET_STRING_ELT(names, 3, mkChar("text"));
    SET_VECTOR_ELT(ans, 3, mkString(Rcon->text ? "text" : "binary"));
    SET_STRING_ELT(names, 4, mkChar("opened"));
    SET_VECTOR_ELT(ans, 4, mkString(Rcon->isopen ? "opened" : "closed"));
    SET_STRING_ELT(names, 5, mkChar("can read"));
    SET_VECTOR_ELT(ans, 5, mkString(Rcon->canread ? "yes" : "no"));
    SET_STRING_ELT(names, 6, mkChar("can write"));
    SET_VECTOR_ELT(ans, 6, mkString(Rcon->canwrite ? "yes" : "no"));
    setAttrib(ans, R_NamesSymbol, names);
    UNPROTECT(2);
    return ans;
}

/* ------------------- internet access functions  --------------------- */

/* TODO  set timeout */
#ifdef HAVE_LIBXML
void *xmlNanoHTTPOpen(const char *URL, char **contentType);
int xmlNanoHTTPRead(void *ctx, void *dest, int len);
void xmlNanoHTTPClose(void *ctx);
int xmlNanoHTTPReturnCode(void *ctx);

void *xmlNanoFTPOpen(const char *URL);
int xmlNanoFTPRead(void *ctx, void *dest, int len);
void xmlNanoFTPClose(void *ctx);
#endif

static void url_open(Rconnection con)
{
    void *ctxt;
    char *url = con->description;
    UrlScheme type = ((Rurlconn)(con->private))->type;
    int rc;

    if (con->mode[0] != 'r')
        error("can only open URLs for reading");

    switch (type)
    {
#ifdef HAVE_LIBXML
    case HTTPsh:
        /* xmlNanoHTTPInit(); */
        ctxt = xmlNanoHTTPOpen(url, NULL);
        if (ctxt == NULL)
            error("cannot open URL `%s'", url);
        rc = xmlNanoHTTPReturnCode(ctxt);
        if (rc != 200)
        {
            xmlNanoHTTPClose(ctxt);
            error("cannot open URL `%s'", url);
        }
        ((Rurlconn)(con->private))->ctxt = ctxt;
        break;
    case FTPsh:
        /* xmlNanoFTPInit(); */
        ctxt = xmlNanoFTPOpen(url);
        if (ctxt == NULL)
            error("cannot open URL `%s'", url);
        ((Rurlconn)(con->private))->ctxt = ctxt;
        break;
#endif
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
        xmlNanoHTTPClose(((Rurlconn)(con->private))->ctxt);
        break;
    case FTPsh:
        xmlNanoFTPClose(((Rurlconn)(con->private))->ctxt);
        break;
    }
    con->isopen = FALSE;
}

static void url_destroy(Rconnection con)
{
    free(con->private);
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
        n = xmlNanoHTTPRead(ctxt, &c, 1);
        break;
    case FTPsh:
        n = xmlNanoFTPRead(ctxt, &c, 1);
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
        n = xmlNanoHTTPRead(ctxt, ptr, size * nitems);
        break;
    case FTPsh:
        n = xmlNanoFTPRead(ctxt, ptr, size * nitems);
        break;
    }
    return n / size;
}

static Rconnection newurl(char *description, char *mode)
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
    strcpy(new->description, description);
    strncpy(new->mode, mode, 4);
    new->mode[4] = '\0';
    new->isopen = new->incomplete = FALSE;
    new->canread = new->canwrite = TRUE; /* in principle */
    new->canseek = FALSE;
    new->text = TRUE;
    new->open = &url_open;
    new->close = &url_close;
    new->destroy = &url_destroy;
    new->vfprintf = &null_vfprintf;
    new->fgetc = &url_fgetc;
    new->seek = &null_seek;
    new->truncate = &null_truncate;
    new->fflush = &null_fflush;
    new->read = &url_read;
    new->write = &null_write;
    new->nPushBack = 0;
    new->private = (void *)malloc(sizeof(struct urlconn));
    if (!new->private)
    {
        free(new->description);
        free(new->class);
        free(new);
        error("allocation of url connection failed");
    }
    return new;
}

/* url(description, open, encoding) */
SEXP do_url(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP scmd, sopen, ans, class, enc;
    char *url, *open, *class2 = "url";
    int i, ncon;
    Rconnection con = NULL;
    UrlScheme type = HTTPsh; /* -Wall */

    checkArity(op, args);
    scmd = CAR(args);
    if (!isString(scmd) || length(scmd) < 1)
        error("invalid `description' argument");
    if (length(scmd) > 1)
        warning("only first element of `description' argument used");
    url = CHAR(STRING_ELT(scmd, 0));
#ifdef HAVE_LIBXML
    if (strncmp(url, "http://", 7) == 0)
    {
        type = HTTPsh;
    }
    else if (strncmp(url, "ftp://", 6) == 0)
    {
        type = FTPsh;
    }
    else
#endif
        error("unsupported URL scheme");

    sopen = CADR(args);
    if (!isString(sopen) || length(sopen) != 1)
        error("invalid `open' argument");
    open = CHAR(STRING_ELT(sopen, 0));
    enc = CADDR(args);
    if (!isInteger(enc) || length(enc) != 256)
        error("invalid `enc' argument");

    ncon = NextConnection();
    if (strncmp(url, "file://", 7) == 0)
    {
        con = newfile(url + 7, strlen(open) ? open : "r");
        class2 = "file";
    }
    else if (strncmp(url, "http://", 7) == 0 || strncmp(url, "ftp://", 6) == 0)
    {
        con = newurl(url, strlen(open) ? open : "r");
        ((Rurlconn)con->private)->type = type;
    }
    else
        error("unsupported URL schema");

    Connections[ncon] = con;
    for (i = 0; i < 256; i++)
        con->encoding[i] = (unsigned char)INTEGER(enc)[i];

    /* open it if desired */
    if (strlen(open))
        con->open(con);

    PROTECT(ans = allocVector(INTSXP, 1));
    INTEGER(ans)[0] = ncon;
    PROTECT(class = allocVector(STRSXP, 2));
    SET_STRING_ELT(class, 0, mkChar(class2));
    SET_STRING_ELT(class, 1, mkChar("connection"));
    classgets(ans, class);
    UNPROTECT(2);

    return ans;
}

void putdots(int *pold, int new)
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
    fflush(stderr);
}

/* TODO select file mode based on ContentType ? */

/* download(url, destfile, quiet) */

#define CPBUFSIZE 65536
#define IBUFSIZE 4096
SEXP do_download(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP ans, scmd, sfile, smode;
    char *url, *file, *mode;
    int quiet, status = 0;

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
    quiet = asLogical(CADDR(args));
    if (quiet == NA_LOGICAL)
        error("invalid `quiet' argument");
    smode = CADDDR(args);
    if (!isString(smode) || length(smode) != 1)
        error("invalid `mode' argument");
    mode = CHAR(STRING_ELT(smode, 0));

    if (strncmp(url, "file://", 7) == 0)
    {
        FILE *in, *out;
        char *buf[CPBUFSIZE];
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

#ifdef HAVE_LIBXML
    }
    else if (strncmp(url, "http://", 7) == 0)
    {

        FILE *out;
        void *ctxt;
        int len, ndots = 0, nnew, nbytes = 0;
        char buf[IBUFSIZE];

        out = R_fopen(R_ExpandFileName(file), mode);
        if (!out)
            error("cannot open destfile `%s'", file);

        /* xmlNanoHTTPInit(); */
        R_Busy(1);
        if (!quiet)
            REprintf("trying URL `%s'\n", url);
        ctxt = xmlNanoHTTPOpen(url, NULL);
        if (ctxt == NULL)
            status = 1;
        else
        {
            int rc = xmlNanoHTTPReturnCode(ctxt);
            if (rc != 200)
            {
                xmlNanoHTTPClose(ctxt);
                status = 1;
                if (!quiet)
                    REprintf("return code %d\n", rc);
            }
            else
            {
                if (!quiet)
                    REprintf("opened URL\n", url);
                while ((len = xmlNanoHTTPRead(ctxt, buf, sizeof(buf))) > 0)
                {
                    fwrite(buf, 1, len, out);
                    nbytes += len;
                    nnew = nbytes / 1024;
                    if (!quiet)
                        putdots(&ndots, nnew);
#ifdef Win32
                    R_ProcessEvents();
#endif
                }
                xmlNanoHTTPClose(ctxt);
                fclose(out);
                R_Busy(0);
                if (!quiet)
                {
                    if (nbytes > 10240)
                        REprintf("\ndownloaded %dKb\n\n", nbytes / 1024, url);
                    else
                        REprintf("\ndownloaded %d bytes\n\n", nbytes, url);
                }
            }
        }
        /* xmlNanoHTTPCleanup();  might have a url/http connection open */
        if (status == 1)
            error("cannot open URL `%s'", url);
    }
    else if (strncmp(url, "ftp://", 6) == 0)
    {

        FILE *out;
        void *ctxt;
        int len, ndots = 0, nnew, nbytes = 0;
        char buf[IBUFSIZE];

        out = R_fopen(R_ExpandFileName(file), mode);
        if (!out)
            error("cannot open destfile `%s'", file);

        /* xmlNanoFTPInit(); */
        R_Busy(1);
        if (!quiet)
            REprintf("trying URL `%s'\n", url);
        ctxt = xmlNanoFTPOpen(url);
        if (ctxt == NULL)
            status = 1;
        else
        {
            if (!quiet)
                REprintf("opened URL\n", url);
            while ((len = xmlNanoFTPRead(ctxt, buf, sizeof(buf))) > 0)
            {
                fwrite(buf, 1, len, out);
                nbytes += len;
                nnew = nbytes / 1024;
                if (!quiet)
                    putdots(&ndots, nnew);
#ifdef Win32
                R_ProcessEvents();
#endif
            }
            xmlNanoFTPClose(ctxt);
            fclose(out);
            R_Busy(0);
            if (!quiet)
            {
                if (nbytes > 10240)
                    REprintf("\ndownloaded %dKb\n\n", nbytes / 1024, url);
                else
                    REprintf("\ndownloaded %d bytes\n\n", nbytes, url);
            }
        }
        /* xmlNanoFTPCleanup();  might have a url/hhtp connection open */
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
