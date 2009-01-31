/*
 *  R : A Computer Language for Statistical Data Analysis
 *  file dounzip.c
 *  first part Copyright (C) 2002-8  the R Development Core Team
 *  second part Copyright (C) 1998 Gilles Vollant
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
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Defn.h>
#include <Fileio.h> /* for R_fopen */
#include "unzip.h"
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif
#include <errno.h>

#ifdef Win32
#include <io.h> /* for mkdir */
#endif

/* cf do_dircreate in platform.c */
static int R_mkdir(char *path)
{
#ifdef Win32
    char local[PATH_MAX];
    strcpy(local, path);
    /* need DOS paths on Win 9x */
    R_fixbackslash(local);
    return mkdir(local);
#endif
#ifdef Unix
    return mkdir(path, 0777);
#endif
}

#define BUF_SIZE 4096
static int extract_one(unzFile uf, const char *const dest, const char *const filename, SEXP names, int *nnames)
{
    int err = UNZ_OK;
    FILE *fout;
    char outname[PATH_MAX], dirs[PATH_MAX], buf[BUF_SIZE], *p, *pp;

    err = unzOpenCurrentFile(uf);
    if (err != UNZ_OK)
        return err;
    if (strlen(dest) > PATH_MAX - 1)
        return 1;
    strcpy(outname, dest);
    strcat(outname, FILESEP);
    if (filename)
    {
        if (strlen(dest) + strlen(filename) > PATH_MAX - 2)
            return 1;
        strcat(outname, filename);
    }
    else
    {
        unz_file_info file_info;
        char filename_inzip[PATH_MAX];
        err = unzGetCurrentFileInfo(uf, &file_info, filename_inzip, sizeof(filename_inzip), NULL, 0, NULL, 0);
        strcat(outname, filename_inzip);
    }
#ifdef Win32
    R_fixslash(outname);
#endif
    p = outname + strlen(outname) - 1;
    if (*p == '/')
    { /* Don't know how these are stored in Mac zip files */
        *p = '\0';
        if (!R_FileExists(outname))
            err = R_mkdir(outname);
    }
    else
    {
        /* make parents as required: have already checked dest exists */
        pp = outname + strlen(dest) + 1;
        while ((p = Rf_strrchr(pp, '/')))
        {
            strcpy(dirs, outname);
            dirs[p - outname] = '\0';
            /* Rprintf("dirs is %s\n", dirs); */
            if (!R_FileExists(dirs))
                R_mkdir(dirs);
            pp = p + 1;
        }
        /* Rprintf("extracting %s\n", outname); */
        fout = R_fopen(outname, "wb");
        if (!fout)
        {
            unzCloseCurrentFile(uf);
            error(_("cannot open file '%s': %s"), outname, strerror(errno));
            return 3; /* not reached */
        }
        while (1)
        {
            err = unzReadCurrentFile(uf, buf, BUF_SIZE);
            /* Rprintf("read %d bytes\n", err); */
            if (err <= 0)
                break;
            if (fwrite(buf, err, 1, fout) != 1)
            {
                err = -200;
                break;
            }
            if (err < BUF_SIZE)
            {
                err = 0;
                break;
            }
        }
        fclose(fout);
        SET_STRING_ELT(names, (*nnames)++, mkChar(outname));
    }
    unzCloseCurrentFile(uf);
    return err;
}

static int do_unzip(const char *zipname, const char *dest, int nfiles, const char **files, SEXP *pnames, int *nnames)
{
    int i, err = UNZ_OK;
    unzFile uf;
    SEXP names = *pnames;

    uf = unzOpen(zipname);
    if (!uf)
        return 1;
    if (nfiles == 0)
    { /* all files */
        unz_global_info gi;
        unzGetGlobalInfo(uf, &gi);
        for (i = 0; i < gi.number_entry; i++)
        {
            if (i > 0)
                if ((err = unzGoToNextFile(uf)) != UNZ_OK)
                    break;
            if (*nnames + 1 >= LENGTH(names))
            {
                SEXP onames = names;
                names = allocVector(STRSXP, 2 * LENGTH(names));
                UNPROTECT(1);
                PROTECT(names);
                copyVector(names, onames);
            }
            if ((err = extract_one(uf, dest, NULL, names, nnames)) != UNZ_OK)
                break;
#ifdef Win32
            R_ProcessEvents();
#else
            R_CheckUserInterrupt();
#endif
        }
    }
    else
    {
        for (i = 0; i < nfiles; i++)
        {
            if ((err = unzLocateFile(uf, files[i], 1)) != UNZ_OK)
                break;
            if ((err = extract_one(uf, dest, files[i], names, nnames)) != UNZ_OK)
                break;
#ifdef Win32
            R_ProcessEvents();
#else
            R_CheckUserInterrupt();
#endif
        }
    }
    *pnames = names;
    unzClose(uf);
    return err;
}

SEXP attribute_hidden do_int_unzip(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP fn, ans, names = R_NilValue;
    char zipname[PATH_MAX], dest[PATH_MAX];
    const char *p, *topics[500];
    int i, ntopics, rc, nnames = 0;

    if (!isString(CAR(args)) || LENGTH(CAR(args)) != 1)
        error(_("invalid zip name argument"));
    p = translateChar(STRING_ELT(CAR(args), 0));
    if (strlen(p) > PATH_MAX - 1)
        error(_("zip path is too long"));
    strcpy(zipname, p);
    args = CDR(args);
    fn = CAR(args);
    ntopics = length(fn);
    if (ntopics > 0)
    {
        if (!isString(fn) || ntopics > 500)
            error(_("invalid '%s' argument"), "topics");
        for (i = 0; i < ntopics; i++)
            topics[i] = translateChar(STRING_ELT(fn, i));
    }
    args = CDR(args);
    if (!isString(CAR(args)) || LENGTH(CAR(args)) != 1)
        error(_("invalid '%s' argument"), "destination");
    p = R_ExpandFileName(translateChar(STRING_ELT(CAR(args), 0)));
    if (strlen(p) > PATH_MAX - 1)
        error(_("'destination' is too long"));
    strcpy(dest, p);
    if (!R_FileExists(dest))
        error(_("'destination' does not exist"));

    if (ntopics > 0)
        PROTECT(names = allocVector(STRSXP, ntopics));
    else
        PROTECT(names = allocVector(STRSXP, 5000));
    rc = do_unzip(zipname, dest, ntopics, topics, &names, &nnames);
    if (rc != UNZ_OK)
        switch (rc)
        {
        case UNZ_END_OF_LIST_OF_FILE:
            warning(_("requested file not found in the zip file"));
            break;
        case UNZ_BADZIPFILE:
            warning(_("zip file is corrupt"));
            break;
        case UNZ_CRCERROR:
            warning(_("CRC error in zip file"));
            break;
        case UNZ_PARAMERROR:
        case UNZ_INTERNALERROR:
            warning(_("internal error in unz code"));
            break;
        case -200:
            warning(_("write error in extracting from zip file"));
            break;
        default:
            warning(_("error %d in extracting from zip file"), rc);
        }
    PROTECT(ans = ScalarLogical(rc));
    PROTECT(names = lengthgets(names, nnames));
    setAttrib(ans, install("extracted"), names);
    UNPROTECT(3);
    return ans;
}

/* ------------------- unz connections --------------------- */

#include <Rconnections.h>

static Rboolean unz_open(Rconnection con)
{
    unzFile uf;
    char path[2 * PATH_MAX], *p;
    const char *tmp;

    if (con->mode[0] != 'r')
    {
        warning(_("unz connections can only be opened for reading"));
        return FALSE;
    }
    tmp = R_ExpandFileName(con->description);
    if (strlen(tmp) > PATH_MAX - 1)
    {
        warning(_("zip path is too long"));
        return FALSE;
    }
    strcpy(path, tmp);
    p = Rf_strrchr(path, ':');
    if (!p)
    {
        warning(_("invalid description of unz connection"));
        return FALSE;
    }
    *p = '\0';
    uf = unzOpen(path);
    if (!uf)
    {
        warning(_("cannot open zip file '%s'"), path);
        return FALSE;
    }
    if (unzLocateFile(uf, p + 1, 1) != UNZ_OK)
    {
        warning(_("cannot locate file '%s' in zip file '%s'"), p + 1, path);
        unzClose(uf);
        return FALSE;
    }
    unzOpenCurrentFile(uf);
    ((Runzconn)(con->private))->uf = uf;
    con->isopen = TRUE;
    con->canwrite = FALSE;
    con->canread = TRUE;
    if (strlen(con->mode) >= 2 && con->mode[1] == 'b')
        con->text = FALSE;
    else
        con->text = TRUE;
    /* set_iconv(); not yet */
    con->save = -1000;
    return TRUE;
}

static void unz_close(Rconnection con)
{
    unzFile uf = ((Runzconn)(con->private))->uf;
    unzCloseCurrentFile(uf);
    unzClose(uf);
    con->isopen = FALSE;
}

static int unz_fgetc_internal(Rconnection con)
{
    unzFile uf = ((Runzconn)(con->private))->uf;
    char buf[1];
    int err, p;

    err = unzReadCurrentFile(uf, buf, 1);
    p = buf[0] % 256;
    return (err < 1) ? R_EOF : p;
}

static size_t unz_read(void *ptr, size_t size, size_t nitems, Rconnection con)
{
    unzFile uf = ((Runzconn)(con->private))->uf;
    return unzReadCurrentFile(uf, ptr, size * nitems) / size;
}

static int null_vfprintf(Rconnection con, const char *format, va_list ap)
{
    error(_("printing not enabled for this connection"));
    return 0; /* -Wall */
}

static size_t null_write(const void *ptr, size_t size, size_t nitems, Rconnection con)
{
    error(_("write not enabled for this connection"));
    return 0; /* -Wall */
}

static double null_seek(Rconnection con, double where, int origin, int rw)
{
    error(_("seek not enabled for this connection"));
    return 0; /* -Wall */
}

static int null_fflush(Rconnection con)
{
    return 0;
}

Rconnection attribute_hidden R_newunz(const char *description, const char *const mode)
{
    Rconnection new;
    new = (Rconnection)malloc(sizeof(struct Rconn));
    if (!new)
        error(_("allocation of unz connection failed"));
    new->class = (char *)malloc(strlen("unz") + 1);
    if (!new->class)
    {
        free(new);
        error(_("allocation of unz connection failed"));
    }
    strcpy(new->class, "unz");
    new->description = (char *)malloc(strlen(description) + 1);
    if (!new->description)
    {
        free(new->class);
        free(new);
        error(_("allocation of unz connection failed"));
    }
    init_con(new, description, CE_NATIVE, mode);

    new->canseek = TRUE;
    new->open = &unz_open;
    new->close = &unz_close;
    new->vfprintf = &null_vfprintf;
    new->fgetc_internal = &unz_fgetc_internal;
    new->fgetc = &dummy_fgetc;
    new->seek = &null_seek;
    new->fflush = &null_fflush;
    new->read = &unz_read;
    new->write = &null_write;
    new->private = (void *)malloc(sizeof(struct fileconn));
    if (!new->private)
    {
        free(new->description);
        free(new->class);
        free(new);
        error(_("allocation of unz connection failed"));
    }
    return new;
}

/* =================== second part ====================== */

/* From minizip contribution to zlib 1.2.3, */

/* unzip.c -- IO for uncompress .zip files using zlib
   Version 1.01e, February 12th, 2005

   Copyright (C) 1998-2005 Gilles Vollant

   Read unzip.h for more info
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zlib.h"
#include "unzip.h"
#ifdef HAVE_ERRNO_H
#include <errno.h>
#else
extern int errno;
#endif

#define local static
#define NOUNCRYPT

#ifndef CASESENSITIVITYDEFAULT_NO
#if !defined(unix) && !defined(CASESENSITIVITYDEFAULT_YES)
#define CASESENSITIVITYDEFAULT_NO
#endif
#endif

#ifndef UNZ_BUFSIZE
#define UNZ_BUFSIZE (16384)
#endif

#ifndef UNZ_MAXFILENAMEINZIP
#define UNZ_MAXFILENAMEINZIP (256)
#endif

#ifndef ALLOC
#define ALLOC(size) (malloc(size))
#endif
#ifndef TRYFREE
#define TRYFREE(p)                                                                                                     \
    {                                                                                                                  \
        if (p)                                                                                                         \
            free(p);                                                                                                   \
    }
#endif

#define SIZECENTRALDIRITEM (0x2e)
#define SIZEZIPLOCALHEADER (0x1e)

static const char unz_copyright[] = " unzip 1.01 Copyright 1998-2004 Gilles Vollant - http://www.winimage.com/zLibDll";

/* unz_file_info_interntal contain internal info about a file in zipfile*/
typedef struct unz_file_info_internal_s
{
    uLong offset_curfile; /* relative offset of local header 4 bytes */
} unz_file_info_internal;

/* file_in_zip_read_info_s contain internal information about a file in zipfile,
    when reading and decompress it */
typedef struct
{
    char *read_buffer; /* internal buffer for compressed data */
    z_stream stream;   /* zLib stream structure for inflate */

    uLong pos_in_zipfile;     /* position in byte on the zipfile, for fseek*/
    uLong stream_initialised; /* flag set if stream structure is initialised*/

    uLong offset_local_extrafield; /* offset of the local extra field */
    uInt size_local_extrafield;    /* size of the local extra field */
    uLong pos_local_extrafield;    /* position in the local extra field in read*/

    uLong crc32;                  /* crc32 of all data uncompressed */
    uLong crc32_wait;             /* crc32 we must obtain after decompress all */
    uLong rest_read_compressed;   /* number of byte to be decompressed */
    uLong rest_read_uncompressed; /*number of byte to be obtained after decomp*/
    zlib_filefunc_def z_filefunc;
    voidpf filestream;             /* io structore of the zipfile */
    uLong compression_method;      /* compression method (0==store) */
    uLong byte_before_the_zipfile; /* byte before the zipfile, (>0 for sfx)*/
    int raw;
} file_in_zip_read_info_s;

/* unz_s contain internal information about the zipfile
 */
typedef struct
{
    zlib_filefunc_def z_filefunc;
    voidpf filestream;             /* io structore of the zipfile */
    unz_global_info gi;            /* public global information */
    uLong byte_before_the_zipfile; /* byte before the zipfile, (>0 for sfx)*/
    uLong num_file;                /* number of the current file in the zipfile*/
    uLong pos_in_central_dir;      /* pos of the current file in the central dir*/
    uLong current_file_ok;         /* flag about the usability of the current file*/
    uLong central_pos;             /* position of the beginning of the central dir*/

    uLong size_central_dir;   /* size of the central directory  */
    uLong offset_central_dir; /* offset of start of central directory with
                                 respect to the starting disk number */

    unz_file_info cur_file_info;                   /* public info about the current file in zip*/
    unz_file_info_internal cur_file_info_internal; /* private info about it*/
    file_in_zip_read_info_s *pfile_in_zip_read;    /* structure about the current
                                           file if we are decompressing it */
    int encrypted;
#ifndef NOUNCRYPT
    unsigned long keys[3]; /* keys defining the pseudo-random sequence */
    const unsigned long *pcrc_32_tab;
#endif
} unz_s;

/* ===========================================================================
     Read a byte from a gz_stream; update next_in and avail_in. Return EOF
   for end of file.
   IN assertion: the stream s has been sucessfully opened for reading.
*/

local int unzlocal_getByte OF((const zlib_filefunc_def *pzlib_filefunc_def, voidpf filestream, int *pi));

local int unzlocal_getByte(pzlib_filefunc_def, filestream, pi) const zlib_filefunc_def *pzlib_filefunc_def;
voidpf filestream;
int *pi;
{
    unsigned char c;
    int err = (int)ZREAD(*pzlib_filefunc_def, filestream, &c, 1);
    if (err == 1)
    {
        *pi = (int)c;
        return UNZ_OK;
    }
    else
    {
        if (ZERROR(*pzlib_filefunc_def, filestream))
            return UNZ_ERRNO;
        else
            return UNZ_EOF;
    }
}

/* ===========================================================================
   Reads a long in LSB order from the given gz_stream. Sets
*/
local int unzlocal_getShort OF((const zlib_filefunc_def *pzlib_filefunc_def, voidpf filestream, uLong *pX));

local int unzlocal_getShort(pzlib_filefunc_def, filestream, pX) const zlib_filefunc_def *pzlib_filefunc_def;
voidpf filestream;
uLong *pX;
{
    uLong x;
    int i;
    int err;

    err = unzlocal_getByte(pzlib_filefunc_def, filestream, &i);
    x = (uLong)i;

    if (err == UNZ_OK)
        err = unzlocal_getByte(pzlib_filefunc_def, filestream, &i);
    x += ((uLong)i) << 8;

    if (err == UNZ_OK)
        *pX = x;
    else
        *pX = 0;
    return err;
}

local int unzlocal_getLong OF((const zlib_filefunc_def *pzlib_filefunc_def, voidpf filestream, uLong *pX));

local int unzlocal_getLong(pzlib_filefunc_def, filestream, pX) const zlib_filefunc_def *pzlib_filefunc_def;
voidpf filestream;
uLong *pX;
{
    uLong x;
    int i;
    int err;

    err = unzlocal_getByte(pzlib_filefunc_def, filestream, &i);
    x = (uLong)i;

    if (err == UNZ_OK)
        err = unzlocal_getByte(pzlib_filefunc_def, filestream, &i);
    x += ((uLong)i) << 8;

    if (err == UNZ_OK)
        err = unzlocal_getByte(pzlib_filefunc_def, filestream, &i);
    x += ((uLong)i) << 16;

    if (err == UNZ_OK)
        err = unzlocal_getByte(pzlib_filefunc_def, filestream, &i);
    x += ((uLong)i) << 24;

    if (err == UNZ_OK)
        *pX = x;
    else
        *pX = 0;
    return err;
}

/* My own strcmpi / strcasecmp */
local int strcmpcasenosensitive_internal(fileName1, fileName2) const char *fileName1;
const char *fileName2;
{
    for (;;)
    {
        char c1 = *(fileName1++);
        char c2 = *(fileName2++);
        if ((c1 >= 'a') && (c1 <= 'z'))
            c1 -= 0x20;
        if ((c2 >= 'a') && (c2 <= 'z'))
            c2 -= 0x20;
        if (c1 == '\0')
            return ((c2 == '\0') ? 0 : -1);
        if (c2 == '\0')
            return 1;
        if (c1 < c2)
            return -1;
        if (c1 > c2)
            return 1;
    }
}

#ifdef CASESENSITIVITYDEFAULT_NO
#define CASESENSITIVITYDEFAULTVALUE 2
#else
#define CASESENSITIVITYDEFAULTVALUE 1
#endif

#ifndef STRCMPCASENOSENTIVEFUNCTION
#define STRCMPCASENOSENTIVEFUNCTION strcmpcasenosensitive_internal
#endif

/*
   Compare two filename (fileName1,fileName2).
   If iCaseSenisivity = 1, comparision is case sensitivity (like strcmp)
   If iCaseSenisivity = 2, comparision is not case sensitivity (like strcmpi
                                                                or strcasecmp)
   If iCaseSenisivity = 0, case sensitivity is defaut of your operating system
        (like 1 on Unix, 2 on Windows)

*/
extern int ZEXPORT unzStringFileNameCompare(fileName1, fileName2, iCaseSensitivity) const char *fileName1;
const char *fileName2;
int iCaseSensitivity;
{
    if (iCaseSensitivity == 0)
        iCaseSensitivity = CASESENSITIVITYDEFAULTVALUE;

    if (iCaseSensitivity == 1)
        return strcmp(fileName1, fileName2);

    return STRCMPCASENOSENTIVEFUNCTION(fileName1, fileName2);
}

#ifndef BUFREADCOMMENT
#define BUFREADCOMMENT (0x400)
#endif

/*
  Locate the Central directory of a zipfile (at the end, just before
    the global comment)
*/
local uLong unzlocal_SearchCentralDir OF((const zlib_filefunc_def *pzlib_filefunc_def, voidpf filestream));

local uLong unzlocal_SearchCentralDir(pzlib_filefunc_def, filestream) const zlib_filefunc_def *pzlib_filefunc_def;
voidpf filestream;
{
    unsigned char *buf;
    uLong uSizeFile;
    uLong uBackRead;
    uLong uMaxBack = 0xffff; /* maximum size of global comment */
    uLong uPosFound = 0;

    if (ZSEEK(*pzlib_filefunc_def, filestream, 0, ZLIB_FILEFUNC_SEEK_END) != 0)
        return 0;

    uSizeFile = ZTELL(*pzlib_filefunc_def, filestream);

    if (uMaxBack > uSizeFile)
        uMaxBack = uSizeFile;

    buf = (unsigned char *)ALLOC(BUFREADCOMMENT + 4);
    if (buf == NULL)
        return 0;

    uBackRead = 4;
    while (uBackRead < uMaxBack)
    {
        uLong uReadSize, uReadPos;
        int i;
        if (uBackRead + BUFREADCOMMENT > uMaxBack)
            uBackRead = uMaxBack;
        else
            uBackRead += BUFREADCOMMENT;
        uReadPos = uSizeFile - uBackRead;

        uReadSize = ((BUFREADCOMMENT + 4) < (uSizeFile - uReadPos)) ? (BUFREADCOMMENT + 4) : (uSizeFile - uReadPos);
        if (ZSEEK(*pzlib_filefunc_def, filestream, uReadPos, ZLIB_FILEFUNC_SEEK_SET) != 0)
            break;

        if (ZREAD(*pzlib_filefunc_def, filestream, buf, uReadSize) != uReadSize)
            break;

        for (i = (int)uReadSize - 3; (i--) > 0;)
            if (((*(buf + i)) == 0x50) && ((*(buf + i + 1)) == 0x4b) && ((*(buf + i + 2)) == 0x05) &&
                ((*(buf + i + 3)) == 0x06))
            {
                uPosFound = uReadPos + i;
                break;
            }

        if (uPosFound != 0)
            break;
    }
    TRYFREE(buf);
    return uPosFound;
}

/*
  Open a Zip file. path contain the full pathname (by example,
     on a Windows NT computer "c:\\test\\zlib114.zip" or on an Unix computer
     "zlib/zlib114.zip".
     If the zipfile cannot be opened (file doesn't exist or in not valid), the
       return value is NULL.
     Else, the return value is a unzFile Handle, usable with other function
       of this unzip package.
*/
extern unzFile ZEXPORT unzOpen2(path, pzlib_filefunc_def) const char *path;
zlib_filefunc_def *pzlib_filefunc_def;
{
    unz_s us;
    unz_s *s;
    uLong central_pos, uL;

    uLong number_disk;         /* number of the current dist, used for
                                  spaning ZIP, unsupported, always 0*/
    uLong number_disk_with_CD; /* number the the disk with central dir, used
                                  for spaning ZIP, unsupported, always 0*/
    uLong number_entry_CD;     /* total number of entries in
                                  the central dir
                                  (same than number_entry on nospan) */

    int err = UNZ_OK;

    if (unz_copyright[0] != ' ')
        return NULL;

    if (pzlib_filefunc_def == NULL)
        fill_fopen_filefunc(&us.z_filefunc);
    else
        us.z_filefunc = *pzlib_filefunc_def;

    us.filestream = (*(us.z_filefunc.zopen_file))(us.z_filefunc.opaque, path,
                                                  ZLIB_FILEFUNC_MODE_READ | ZLIB_FILEFUNC_MODE_EXISTING);
    if (us.filestream == NULL)
        return NULL;

    central_pos = unzlocal_SearchCentralDir(&us.z_filefunc, us.filestream);
    if (central_pos == 0)
        err = UNZ_ERRNO;

    if (ZSEEK(us.z_filefunc, us.filestream, central_pos, ZLIB_FILEFUNC_SEEK_SET) != 0)
        err = UNZ_ERRNO;

    /* the signature, already checked */
    if (unzlocal_getLong(&us.z_filefunc, us.filestream, &uL) != UNZ_OK)
        err = UNZ_ERRNO;

    /* number of this disk */
    if (unzlocal_getShort(&us.z_filefunc, us.filestream, &number_disk) != UNZ_OK)
        err = UNZ_ERRNO;

    /* number of the disk with the start of the central directory */
    if (unzlocal_getShort(&us.z_filefunc, us.filestream, &number_disk_with_CD) != UNZ_OK)
        err = UNZ_ERRNO;

    /* total number of entries in the central dir on this disk */
    if (unzlocal_getShort(&us.z_filefunc, us.filestream, &us.gi.number_entry) != UNZ_OK)
        err = UNZ_ERRNO;

    /* total number of entries in the central dir */
    if (unzlocal_getShort(&us.z_filefunc, us.filestream, &number_entry_CD) != UNZ_OK)
        err = UNZ_ERRNO;

    if ((number_entry_CD != us.gi.number_entry) || (number_disk_with_CD != 0) || (number_disk != 0))
        err = UNZ_BADZIPFILE;

    /* size of the central directory */
    if (unzlocal_getLong(&us.z_filefunc, us.filestream, &us.size_central_dir) != UNZ_OK)
        err = UNZ_ERRNO;

    /* offset of start of central directory with respect to the
          starting disk number */
    if (unzlocal_getLong(&us.z_filefunc, us.filestream, &us.offset_central_dir) != UNZ_OK)
        err = UNZ_ERRNO;

    /* zipfile comment length */
    if (unzlocal_getShort(&us.z_filefunc, us.filestream, &us.gi.size_comment) != UNZ_OK)
        err = UNZ_ERRNO;

    if ((central_pos < us.offset_central_dir + us.size_central_dir) && (err == UNZ_OK))
        err = UNZ_BADZIPFILE;

    if (err != UNZ_OK)
    {
        ZCLOSE(us.z_filefunc, us.filestream);
        return NULL;
    }

    us.byte_before_the_zipfile = central_pos - (us.offset_central_dir + us.size_central_dir);
    us.central_pos = central_pos;
    us.pfile_in_zip_read = NULL;
    us.encrypted = 0;

    s = (unz_s *)ALLOC(sizeof(unz_s));
    *s = us;
    unzGoToFirstFile((unzFile)s);
    return (unzFile)s;
}

extern unzFile ZEXPORT unzOpen(path) const char *path;
{
    return unzOpen2(path, NULL);
}

/*
  Close a ZipFile opened with unzipOpen.
  If there is files inside the .Zip opened with unzipOpenCurrentFile (see later),
    these files MUST be closed with unzipCloseCurrentFile before call unzipClose.
  return UNZ_OK if there is no problem. */
extern int ZEXPORT unzClose(file) unzFile file;
{
    unz_s *s;
    if (file == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;

    if (s->pfile_in_zip_read != NULL)
        unzCloseCurrentFile(file);

    ZCLOSE(s->z_filefunc, s->filestream);
    TRYFREE(s);
    return UNZ_OK;
}

/*
  Write info about the ZipFile in the *pglobal_info structure.
  No preparation of the structure is needed
  return UNZ_OK if there is no problem. */
extern int ZEXPORT unzGetGlobalInfo(file, pglobal_info) unzFile file;
unz_global_info *pglobal_info;
{
    unz_s *s;
    if (file == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;
    *pglobal_info = s->gi;
    return UNZ_OK;
}

/*
   Translate date/time from Dos format to tm_unz (readable more easilty)
*/
local void unzlocal_DosDateToTmuDate(ulDosDate, ptm) uLong ulDosDate;
tm_unz *ptm;
{
    uLong uDate;
    uDate = (uLong)(ulDosDate >> 16);
    ptm->tm_mday = (uInt)(uDate & 0x1f);
    ptm->tm_mon = (uInt)((((uDate)&0x1E0) / 0x20) - 1);
    ptm->tm_year = (uInt)(((uDate & 0x0FE00) / 0x0200) + 1980);

    ptm->tm_hour = (uInt)((ulDosDate & 0xF800) / 0x800);
    ptm->tm_min = (uInt)((ulDosDate & 0x7E0) / 0x20);
    ptm->tm_sec = (uInt)(2 * (ulDosDate & 0x1f));
}

/*
  Get Info about the current file in the zipfile, with internal only info
*/
local int unzlocal_GetCurrentFileInfoInternal OF(
    (unzFile file, unz_file_info *pfile_info, unz_file_info_internal *pfile_info_internal, char *szFileName,
     uLong fileNameBufferSize, void *extraField, uLong extraFieldBufferSize, char *szComment, uLong commentBufferSize));

local int unzlocal_GetCurrentFileInfoInternal(file, pfile_info, pfile_info_internal, szFileName, fileNameBufferSize,
                                              extraField, extraFieldBufferSize, szComment,
                                              commentBufferSize) unzFile file;
unz_file_info *pfile_info;
unz_file_info_internal *pfile_info_internal;
char *szFileName;
uLong fileNameBufferSize;
void *extraField;
uLong extraFieldBufferSize;
char *szComment;
uLong commentBufferSize;
{
    unz_s *s;
    unz_file_info file_info;
    unz_file_info_internal file_info_internal;
    int err = UNZ_OK;
    uLong uMagic;
    long lSeek = 0;

    if (file == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;
    if (ZSEEK(s->z_filefunc, s->filestream, s->pos_in_central_dir + s->byte_before_the_zipfile,
              ZLIB_FILEFUNC_SEEK_SET) != 0)
        err = UNZ_ERRNO;

    /* we check the magic */
    if (err == UNZ_OK)
        if (unzlocal_getLong(&s->z_filefunc, s->filestream, &uMagic) != UNZ_OK)
            err = UNZ_ERRNO;
        else if (uMagic != 0x02014b50)
            err = UNZ_BADZIPFILE;

    if (unzlocal_getShort(&s->z_filefunc, s->filestream, &file_info.version) != UNZ_OK)
        err = UNZ_ERRNO;

    if (unzlocal_getShort(&s->z_filefunc, s->filestream, &file_info.version_needed) != UNZ_OK)
        err = UNZ_ERRNO;

    if (unzlocal_getShort(&s->z_filefunc, s->filestream, &file_info.flag) != UNZ_OK)
        err = UNZ_ERRNO;

    if (unzlocal_getShort(&s->z_filefunc, s->filestream, &file_info.compression_method) != UNZ_OK)
        err = UNZ_ERRNO;

    if (unzlocal_getLong(&s->z_filefunc, s->filestream, &file_info.dosDate) != UNZ_OK)
        err = UNZ_ERRNO;

    unzlocal_DosDateToTmuDate(file_info.dosDate, &file_info.tmu_date);

    if (unzlocal_getLong(&s->z_filefunc, s->filestream, &file_info.crc) != UNZ_OK)
        err = UNZ_ERRNO;

    if (unzlocal_getLong(&s->z_filefunc, s->filestream, &file_info.compressed_size) != UNZ_OK)
        err = UNZ_ERRNO;

    if (unzlocal_getLong(&s->z_filefunc, s->filestream, &file_info.uncompressed_size) != UNZ_OK)
        err = UNZ_ERRNO;

    if (unzlocal_getShort(&s->z_filefunc, s->filestream, &file_info.size_filename) != UNZ_OK)
        err = UNZ_ERRNO;

    if (unzlocal_getShort(&s->z_filefunc, s->filestream, &file_info.size_file_extra) != UNZ_OK)
        err = UNZ_ERRNO;

    if (unzlocal_getShort(&s->z_filefunc, s->filestream, &file_info.size_file_comment) != UNZ_OK)
        err = UNZ_ERRNO;

    if (unzlocal_getShort(&s->z_filefunc, s->filestream, &file_info.disk_num_start) != UNZ_OK)
        err = UNZ_ERRNO;

    if (unzlocal_getShort(&s->z_filefunc, s->filestream, &file_info.internal_fa) != UNZ_OK)
        err = UNZ_ERRNO;

    if (unzlocal_getLong(&s->z_filefunc, s->filestream, &file_info.external_fa) != UNZ_OK)
        err = UNZ_ERRNO;

    if (unzlocal_getLong(&s->z_filefunc, s->filestream, &file_info_internal.offset_curfile) != UNZ_OK)
        err = UNZ_ERRNO;

    lSeek += file_info.size_filename;
    if ((err == UNZ_OK) && (szFileName != NULL))
    {
        uLong uSizeRead;
        if (file_info.size_filename < fileNameBufferSize)
        {
            *(szFileName + file_info.size_filename) = '\0';
            uSizeRead = file_info.size_filename;
        }
        else
            uSizeRead = fileNameBufferSize;

        if ((file_info.size_filename > 0) && (fileNameBufferSize > 0))
            if (ZREAD(s->z_filefunc, s->filestream, szFileName, uSizeRead) != uSizeRead)
                err = UNZ_ERRNO;
        lSeek -= uSizeRead;
    }

    if ((err == UNZ_OK) && (extraField != NULL))
    {
        uLong uSizeRead;
        if (file_info.size_file_extra < extraFieldBufferSize)
            uSizeRead = file_info.size_file_extra;
        else
            uSizeRead = extraFieldBufferSize;

        if (lSeek != 0)
            if (ZSEEK(s->z_filefunc, s->filestream, lSeek, ZLIB_FILEFUNC_SEEK_CUR) == 0)
                lSeek = 0;
            else
                err = UNZ_ERRNO;
        if ((file_info.size_file_extra > 0) && (extraFieldBufferSize > 0))
            if (ZREAD(s->z_filefunc, s->filestream, extraField, uSizeRead) != uSizeRead)
                err = UNZ_ERRNO;
        lSeek += file_info.size_file_extra - uSizeRead;
    }
    else
        lSeek += file_info.size_file_extra;

    if ((err == UNZ_OK) && (szComment != NULL))
    {
        uLong uSizeRead;
        if (file_info.size_file_comment < commentBufferSize)
        {
            *(szComment + file_info.size_file_comment) = '\0';
            uSizeRead = file_info.size_file_comment;
        }
        else
            uSizeRead = commentBufferSize;

        if (lSeek != 0)
            if (ZSEEK(s->z_filefunc, s->filestream, lSeek, ZLIB_FILEFUNC_SEEK_CUR) == 0)
                lSeek = 0;
            else
                err = UNZ_ERRNO;
        if ((file_info.size_file_comment > 0) && (commentBufferSize > 0))
            if (ZREAD(s->z_filefunc, s->filestream, szComment, uSizeRead) != uSizeRead)
                err = UNZ_ERRNO;
        lSeek += file_info.size_file_comment - uSizeRead;
    }
    else
        lSeek += file_info.size_file_comment;

    if ((err == UNZ_OK) && (pfile_info != NULL))
        *pfile_info = file_info;

    if ((err == UNZ_OK) && (pfile_info_internal != NULL))
        *pfile_info_internal = file_info_internal;

    return err;
}

/*
  Write info about the ZipFile in the *pglobal_info structure.
  No preparation of the structure is needed
  return UNZ_OK if there is no problem.
*/
extern int ZEXPORT unzGetCurrentFileInfo(file, pfile_info, szFileName, fileNameBufferSize, extraField,
                                         extraFieldBufferSize, szComment, commentBufferSize) unzFile file;
unz_file_info *pfile_info;
char *szFileName;
uLong fileNameBufferSize;
void *extraField;
uLong extraFieldBufferSize;
char *szComment;
uLong commentBufferSize;
{
    return unzlocal_GetCurrentFileInfoInternal(file, pfile_info, NULL, szFileName, fileNameBufferSize, extraField,
                                               extraFieldBufferSize, szComment, commentBufferSize);
}

/*
  Set the current file of the zipfile to the first file.
  return UNZ_OK if there is no problem
*/
extern int ZEXPORT unzGoToFirstFile(file) unzFile file;
{
    int err = UNZ_OK;
    unz_s *s;
    if (file == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;
    s->pos_in_central_dir = s->offset_central_dir;
    s->num_file = 0;
    err = unzlocal_GetCurrentFileInfoInternal(file, &s->cur_file_info, &s->cur_file_info_internal, NULL, 0, NULL, 0,
                                              NULL, 0);
    s->current_file_ok = (err == UNZ_OK);
    return err;
}

/*
  Set the current file of the zipfile to the next file.
  return UNZ_OK if there is no problem
  return UNZ_END_OF_LIST_OF_FILE if the actual file was the latest.
*/
extern int ZEXPORT unzGoToNextFile(file) unzFile file;
{
    unz_s *s;
    int err;

    if (file == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;
    if (!s->current_file_ok)
        return UNZ_END_OF_LIST_OF_FILE;
    if (s->gi.number_entry != 0xffff) /* 2^16 files overflow hack */
        if (s->num_file + 1 == s->gi.number_entry)
            return UNZ_END_OF_LIST_OF_FILE;

    s->pos_in_central_dir += SIZECENTRALDIRITEM + s->cur_file_info.size_filename + s->cur_file_info.size_file_extra +
                             s->cur_file_info.size_file_comment;
    s->num_file++;
    err = unzlocal_GetCurrentFileInfoInternal(file, &s->cur_file_info, &s->cur_file_info_internal, NULL, 0, NULL, 0,
                                              NULL, 0);
    s->current_file_ok = (err == UNZ_OK);
    return err;
}

/*
  Try locate the file szFileName in the zipfile.
  For the iCaseSensitivity signification, see unzipStringFileNameCompare

  return value :
  UNZ_OK if the file is found. It becomes the current file.
  UNZ_END_OF_LIST_OF_FILE if the file is not found
*/
extern int ZEXPORT unzLocateFile(file, szFileName, iCaseSensitivity) unzFile file;
const char *szFileName;
int iCaseSensitivity;
{
    unz_s *s;
    int err;

    /* We remember the 'current' position in the file so that we can jump
     * back there if we fail.
     */
    unz_file_info cur_file_infoSaved;
    unz_file_info_internal cur_file_info_internalSaved;
    uLong num_fileSaved;
    uLong pos_in_central_dirSaved;

    if (file == NULL)
        return UNZ_PARAMERROR;

    if (strlen(szFileName) >= UNZ_MAXFILENAMEINZIP)
        return UNZ_PARAMERROR;

    s = (unz_s *)file;
    if (!s->current_file_ok)
        return UNZ_END_OF_LIST_OF_FILE;

    /* Save the current state */
    num_fileSaved = s->num_file;
    pos_in_central_dirSaved = s->pos_in_central_dir;
    cur_file_infoSaved = s->cur_file_info;
    cur_file_info_internalSaved = s->cur_file_info_internal;

    err = unzGoToFirstFile(file);

    while (err == UNZ_OK)
    {
        char szCurrentFileName[UNZ_MAXFILENAMEINZIP + 1];
        err = unzGetCurrentFileInfo(file, NULL, szCurrentFileName, sizeof(szCurrentFileName) - 1, NULL, 0, NULL, 0);
        if (err == UNZ_OK)
        {
            if (unzStringFileNameCompare(szCurrentFileName, szFileName, iCaseSensitivity) == 0)
                return UNZ_OK;
            err = unzGoToNextFile(file);
        }
    }

    /* We failed, so restore the state of the 'current file' to where we
     * were.
     */
    s->num_file = num_fileSaved;
    s->pos_in_central_dir = pos_in_central_dirSaved;
    s->cur_file_info = cur_file_infoSaved;
    s->cur_file_info_internal = cur_file_info_internalSaved;
    return err;
}

/*
///////////////////////////////////////////
// Contributed by Ryan Haksi (mailto://cryogen@infoserve.net)
// I need random access
//
// Further optimization could be realized by adding an ability
// to cache the directory in memory. The goal being a single
// comprehensive file read to put the file I need in a memory.
*/

/*
typedef struct unz_file_pos_s
{
    uLong pos_in_zip_directory;   // offset in file
    uLong num_of_file;            // # of file
} unz_file_pos;
*/

extern int ZEXPORT unzGetFilePos(file, file_pos) unzFile file;
unz_file_pos *file_pos;
{
    unz_s *s;

    if (file == NULL || file_pos == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;
    if (!s->current_file_ok)
        return UNZ_END_OF_LIST_OF_FILE;

    file_pos->pos_in_zip_directory = s->pos_in_central_dir;
    file_pos->num_of_file = s->num_file;

    return UNZ_OK;
}

extern int ZEXPORT unzGoToFilePos(file, file_pos) unzFile file;
unz_file_pos *file_pos;
{
    unz_s *s;
    int err;

    if (file == NULL || file_pos == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;

    /* jump to the right spot */
    s->pos_in_central_dir = file_pos->pos_in_zip_directory;
    s->num_file = file_pos->num_of_file;

    /* set the current file */
    err = unzlocal_GetCurrentFileInfoInternal(file, &s->cur_file_info, &s->cur_file_info_internal, NULL, 0, NULL, 0,
                                              NULL, 0);
    /* return results */
    s->current_file_ok = (err == UNZ_OK);
    return err;
}

/*
// Unzip Helper Functions - should be here?
///////////////////////////////////////////
*/

/*
  Read the local header of the current zipfile
  Check the coherency of the local header and info in the end of central
        directory about this file
  store in *piSizeVar the size of extra info in local header
        (filename and size of extra field data)
*/
local int unzlocal_CheckCurrentFileCoherencyHeader(s, piSizeVar, poffset_local_extrafield,
                                                   psize_local_extrafield) unz_s *s;
uInt *piSizeVar;
uLong *poffset_local_extrafield;
uInt *psize_local_extrafield;
{
    uLong uMagic, uData, uFlags;
    uLong size_filename;
    uLong size_extra_field;
    int err = UNZ_OK;

    *piSizeVar = 0;
    *poffset_local_extrafield = 0;
    *psize_local_extrafield = 0;

    if (ZSEEK(s->z_filefunc, s->filestream, s->cur_file_info_internal.offset_curfile + s->byte_before_the_zipfile,
              ZLIB_FILEFUNC_SEEK_SET) != 0)
        return UNZ_ERRNO;

    if (err == UNZ_OK)
        if (unzlocal_getLong(&s->z_filefunc, s->filestream, &uMagic) != UNZ_OK)
            err = UNZ_ERRNO;
        else if (uMagic != 0x04034b50)
            err = UNZ_BADZIPFILE;

    if (unzlocal_getShort(&s->z_filefunc, s->filestream, &uData) != UNZ_OK)
        err = UNZ_ERRNO;
    /*
        else if ((err==UNZ_OK) && (uData!=s->cur_file_info.wVersion))
            err=UNZ_BADZIPFILE;
    */
    if (unzlocal_getShort(&s->z_filefunc, s->filestream, &uFlags) != UNZ_OK)
        err = UNZ_ERRNO;

    if (unzlocal_getShort(&s->z_filefunc, s->filestream, &uData) != UNZ_OK)
        err = UNZ_ERRNO;
    else if ((err == UNZ_OK) && (uData != s->cur_file_info.compression_method))
        err = UNZ_BADZIPFILE;

    if ((err == UNZ_OK) && (s->cur_file_info.compression_method != 0) &&
        (s->cur_file_info.compression_method != Z_DEFLATED))
        err = UNZ_BADZIPFILE;

    if (unzlocal_getLong(&s->z_filefunc, s->filestream, &uData) != UNZ_OK) /* date/time */
        err = UNZ_ERRNO;

    if (unzlocal_getLong(&s->z_filefunc, s->filestream, &uData) != UNZ_OK) /* crc */
        err = UNZ_ERRNO;
    else if ((err == UNZ_OK) && (uData != s->cur_file_info.crc) && ((uFlags & 8) == 0))
        err = UNZ_BADZIPFILE;

    if (unzlocal_getLong(&s->z_filefunc, s->filestream, &uData) != UNZ_OK) /* size compr */
        err = UNZ_ERRNO;
    else if ((err == UNZ_OK) && (uData != s->cur_file_info.compressed_size) && ((uFlags & 8) == 0))
        err = UNZ_BADZIPFILE;

    if (unzlocal_getLong(&s->z_filefunc, s->filestream, &uData) != UNZ_OK) /* size uncompr */
        err = UNZ_ERRNO;
    else if ((err == UNZ_OK) && (uData != s->cur_file_info.uncompressed_size) && ((uFlags & 8) == 0))
        err = UNZ_BADZIPFILE;

    if (unzlocal_getShort(&s->z_filefunc, s->filestream, &size_filename) != UNZ_OK)
        err = UNZ_ERRNO;
    else if ((err == UNZ_OK) && (size_filename != s->cur_file_info.size_filename))
        err = UNZ_BADZIPFILE;

    *piSizeVar += (uInt)size_filename;

    if (unzlocal_getShort(&s->z_filefunc, s->filestream, &size_extra_field) != UNZ_OK)
        err = UNZ_ERRNO;
    *poffset_local_extrafield = s->cur_file_info_internal.offset_curfile + SIZEZIPLOCALHEADER + size_filename;
    *psize_local_extrafield = (uInt)size_extra_field;

    *piSizeVar += (uInt)size_extra_field;

    return err;
}

/*
  Open for reading data the current file in the zipfile.
  If there is no error and the file is opened, the return value is UNZ_OK.
*/
extern int ZEXPORT unzOpenCurrentFile3(file, method, level, raw, password) unzFile file;
int *method;
int *level;
int raw;
const char *password;
{
    int err = UNZ_OK;
    uInt iSizeVar;
    unz_s *s;
    file_in_zip_read_info_s *pfile_in_zip_read_info;
    uLong offset_local_extrafield; /* offset of the local extra field */
    uInt size_local_extrafield;    /* size of the local extra field */
#ifndef NOUNCRYPT
    char source[12];
#else
    if (password != NULL)
        return UNZ_PARAMERROR;
#endif

    if (file == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;
    if (!s->current_file_ok)
        return UNZ_PARAMERROR;

    if (s->pfile_in_zip_read != NULL)
        unzCloseCurrentFile(file);

    if (unzlocal_CheckCurrentFileCoherencyHeader(s, &iSizeVar, &offset_local_extrafield, &size_local_extrafield) !=
        UNZ_OK)
        return UNZ_BADZIPFILE;

    pfile_in_zip_read_info = (file_in_zip_read_info_s *)ALLOC(sizeof(file_in_zip_read_info_s));
    if (pfile_in_zip_read_info == NULL)
        return UNZ_INTERNALERROR;

    pfile_in_zip_read_info->read_buffer = (char *)ALLOC(UNZ_BUFSIZE);
    pfile_in_zip_read_info->offset_local_extrafield = offset_local_extrafield;
    pfile_in_zip_read_info->size_local_extrafield = size_local_extrafield;
    pfile_in_zip_read_info->pos_local_extrafield = 0;
    pfile_in_zip_read_info->raw = raw;

    if (pfile_in_zip_read_info->read_buffer == NULL)
    {
        TRYFREE(pfile_in_zip_read_info);
        return UNZ_INTERNALERROR;
    }

    pfile_in_zip_read_info->stream_initialised = 0;

    if (method != NULL)
        *method = (int)s->cur_file_info.compression_method;

    if (level != NULL)
    {
        *level = 6;
        switch (s->cur_file_info.flag & 0x06)
        {
        case 6:
            *level = 1;
            break;
        case 4:
            *level = 2;
            break;
        case 2:
            *level = 9;
            break;
        }
    }

    if ((s->cur_file_info.compression_method != 0) && (s->cur_file_info.compression_method != Z_DEFLATED))
        err = UNZ_BADZIPFILE;

    pfile_in_zip_read_info->crc32_wait = s->cur_file_info.crc;
    pfile_in_zip_read_info->crc32 = 0;
    pfile_in_zip_read_info->compression_method = s->cur_file_info.compression_method;
    pfile_in_zip_read_info->filestream = s->filestream;
    pfile_in_zip_read_info->z_filefunc = s->z_filefunc;
    pfile_in_zip_read_info->byte_before_the_zipfile = s->byte_before_the_zipfile;

    pfile_in_zip_read_info->stream.total_out = 0;

    if ((s->cur_file_info.compression_method == Z_DEFLATED) && (!raw))
    {
        pfile_in_zip_read_info->stream.zalloc = (alloc_func)0;
        pfile_in_zip_read_info->stream.zfree = (free_func)0;
        pfile_in_zip_read_info->stream.opaque = (voidpf)0;
        pfile_in_zip_read_info->stream.next_in = (voidpf)0;
        pfile_in_zip_read_info->stream.avail_in = 0;

        err = inflateInit2(&pfile_in_zip_read_info->stream, -MAX_WBITS);
        if (err == Z_OK)
            pfile_in_zip_read_info->stream_initialised = 1;
        else
        {
            TRYFREE(pfile_in_zip_read_info);
            return err;
        }
        /* windowBits is passed < 0 to tell that there is no zlib header.
         * Note that in this case inflate *requires* an extra "dummy" byte
         * after the compressed stream in order to complete decompression and
         * return Z_STREAM_END.
         * In unzip, i don't wait absolutely Z_STREAM_END because I known the
         * size of both compressed and uncompressed data
         */
    }
    pfile_in_zip_read_info->rest_read_compressed = s->cur_file_info.compressed_size;
    pfile_in_zip_read_info->rest_read_uncompressed = s->cur_file_info.uncompressed_size;

    pfile_in_zip_read_info->pos_in_zipfile = s->cur_file_info_internal.offset_curfile + SIZEZIPLOCALHEADER + iSizeVar;

    pfile_in_zip_read_info->stream.avail_in = (uInt)0;

    s->pfile_in_zip_read = pfile_in_zip_read_info;

#ifndef NOUNCRYPT
    if (password != NULL)
    {
        int i;
        s->pcrc_32_tab = get_crc_table();
        init_keys(password, s->keys, s->pcrc_32_tab);
        if (ZSEEK(s->z_filefunc, s->filestream,
                  s->pfile_in_zip_read->pos_in_zipfile + s->pfile_in_zip_read->byte_before_the_zipfile, SEEK_SET) != 0)
            return UNZ_INTERNALERROR;
        if (ZREAD(s->z_filefunc, s->filestream, source, 12) < 12)
            return UNZ_INTERNALERROR;

        for (i = 0; i < 12; i++)
            zdecode(s->keys, s->pcrc_32_tab, source[i]);

        s->pfile_in_zip_read->pos_in_zipfile += 12;
        s->encrypted = 1;
    }
#endif

    return UNZ_OK;
}

extern int ZEXPORT unzOpenCurrentFile(file) unzFile file;
{
    return unzOpenCurrentFile3(file, NULL, NULL, 0, NULL);
}

extern int ZEXPORT unzOpenCurrentFilePassword(file, password) unzFile file;
const char *password;
{
    return unzOpenCurrentFile3(file, NULL, NULL, 0, password);
}

extern int ZEXPORT unzOpenCurrentFile2(file, method, level, raw) unzFile file;
int *method;
int *level;
int raw;
{
    return unzOpenCurrentFile3(file, method, level, raw, NULL);
}

/*
  Read bytes from the current file.
  buf contain buffer where data must be copied
  len the size of buf.

  return the number of byte copied if somes bytes are copied
  return 0 if the end of file was reached
  return <0 with error code if there is an error
    (UNZ_ERRNO for IO error, or zLib error for uncompress error)
*/
extern int ZEXPORT unzReadCurrentFile(file, buf, len) unzFile file;
voidp buf;
unsigned len;
{
    int err = UNZ_OK;
    uInt iRead = 0;
    unz_s *s;
    file_in_zip_read_info_s *pfile_in_zip_read_info;
    if (file == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;
    pfile_in_zip_read_info = s->pfile_in_zip_read;

    if (pfile_in_zip_read_info == NULL)
        return UNZ_PARAMERROR;

    if ((pfile_in_zip_read_info->read_buffer == NULL))
        return UNZ_END_OF_LIST_OF_FILE;
    if (len == 0)
        return 0;

    pfile_in_zip_read_info->stream.next_out = (Bytef *)buf;

    pfile_in_zip_read_info->stream.avail_out = (uInt)len;

    if ((len > pfile_in_zip_read_info->rest_read_uncompressed) && (!(pfile_in_zip_read_info->raw)))
        pfile_in_zip_read_info->stream.avail_out = (uInt)pfile_in_zip_read_info->rest_read_uncompressed;

    if ((len > pfile_in_zip_read_info->rest_read_compressed + pfile_in_zip_read_info->stream.avail_in) &&
        (pfile_in_zip_read_info->raw))
        pfile_in_zip_read_info->stream.avail_out =
            (uInt)pfile_in_zip_read_info->rest_read_compressed + pfile_in_zip_read_info->stream.avail_in;

    while (pfile_in_zip_read_info->stream.avail_out > 0)
    {
        if ((pfile_in_zip_read_info->stream.avail_in == 0) && (pfile_in_zip_read_info->rest_read_compressed > 0))
        {
            uInt uReadThis = UNZ_BUFSIZE;
            if (pfile_in_zip_read_info->rest_read_compressed < uReadThis)
                uReadThis = (uInt)pfile_in_zip_read_info->rest_read_compressed;
            if (uReadThis == 0)
                return UNZ_EOF;
            if (ZSEEK(pfile_in_zip_read_info->z_filefunc, pfile_in_zip_read_info->filestream,
                      pfile_in_zip_read_info->pos_in_zipfile + pfile_in_zip_read_info->byte_before_the_zipfile,
                      ZLIB_FILEFUNC_SEEK_SET) != 0)
                return UNZ_ERRNO;
            if (ZREAD(pfile_in_zip_read_info->z_filefunc, pfile_in_zip_read_info->filestream,
                      pfile_in_zip_read_info->read_buffer, uReadThis) != uReadThis)
                return UNZ_ERRNO;

#ifndef NOUNCRYPT
            if (s->encrypted)
            {
                uInt i;
                for (i = 0; i < uReadThis; i++)
                    pfile_in_zip_read_info->read_buffer[i] =
                        zdecode(s->keys, s->pcrc_32_tab, pfile_in_zip_read_info->read_buffer[i]);
            }
#endif

            pfile_in_zip_read_info->pos_in_zipfile += uReadThis;

            pfile_in_zip_read_info->rest_read_compressed -= uReadThis;

            pfile_in_zip_read_info->stream.next_in = (Bytef *)pfile_in_zip_read_info->read_buffer;
            pfile_in_zip_read_info->stream.avail_in = (uInt)uReadThis;
        }

        if ((pfile_in_zip_read_info->compression_method == 0) || (pfile_in_zip_read_info->raw))
        {
            uInt uDoCopy, i;

            if ((pfile_in_zip_read_info->stream.avail_in == 0) && (pfile_in_zip_read_info->rest_read_compressed == 0))
                return (iRead == 0) ? UNZ_EOF : iRead;

            if (pfile_in_zip_read_info->stream.avail_out < pfile_in_zip_read_info->stream.avail_in)
                uDoCopy = pfile_in_zip_read_info->stream.avail_out;
            else
                uDoCopy = pfile_in_zip_read_info->stream.avail_in;

            for (i = 0; i < uDoCopy; i++)
                *(pfile_in_zip_read_info->stream.next_out + i) = *(pfile_in_zip_read_info->stream.next_in + i);

            pfile_in_zip_read_info->crc32 =
                crc32(pfile_in_zip_read_info->crc32, pfile_in_zip_read_info->stream.next_out, uDoCopy);
            pfile_in_zip_read_info->rest_read_uncompressed -= uDoCopy;
            pfile_in_zip_read_info->stream.avail_in -= uDoCopy;
            pfile_in_zip_read_info->stream.avail_out -= uDoCopy;
            pfile_in_zip_read_info->stream.next_out += uDoCopy;
            pfile_in_zip_read_info->stream.next_in += uDoCopy;
            pfile_in_zip_read_info->stream.total_out += uDoCopy;
            iRead += uDoCopy;
        }
        else
        {
            uLong uTotalOutBefore, uTotalOutAfter;
            const Bytef *bufBefore;
            uLong uOutThis;
            int flush = Z_SYNC_FLUSH;

            uTotalOutBefore = pfile_in_zip_read_info->stream.total_out;
            bufBefore = pfile_in_zip_read_info->stream.next_out;

            /*
            if ((pfile_in_zip_read_info->rest_read_uncompressed ==
                     pfile_in_zip_read_info->stream.avail_out) &&
                (pfile_in_zip_read_info->rest_read_compressed == 0))
                flush = Z_FINISH;
            */
            err = inflate(&pfile_in_zip_read_info->stream, flush);

            if ((err >= 0) && (pfile_in_zip_read_info->stream.msg != NULL))
                err = Z_DATA_ERROR;

            uTotalOutAfter = pfile_in_zip_read_info->stream.total_out;
            uOutThis = uTotalOutAfter - uTotalOutBefore;

            pfile_in_zip_read_info->crc32 = crc32(pfile_in_zip_read_info->crc32, bufBefore, (uInt)(uOutThis));

            pfile_in_zip_read_info->rest_read_uncompressed -= uOutThis;

            iRead += (uInt)(uTotalOutAfter - uTotalOutBefore);

            if (err == Z_STREAM_END)
                return (iRead == 0) ? UNZ_EOF : iRead;
            if (err != Z_OK)
                break;
        }
    }

    if (err == Z_OK)
        return iRead;
    return err;
}

/*
  Give the current position in uncompressed data
*/
extern z_off_t ZEXPORT unztell(file) unzFile file;
{
    unz_s *s;
    file_in_zip_read_info_s *pfile_in_zip_read_info;
    if (file == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;
    pfile_in_zip_read_info = s->pfile_in_zip_read;

    if (pfile_in_zip_read_info == NULL)
        return UNZ_PARAMERROR;

    return (z_off_t)pfile_in_zip_read_info->stream.total_out;
}

/*
  return 1 if the end of file was reached, 0 elsewhere
*/
extern int ZEXPORT unzeof(file) unzFile file;
{
    unz_s *s;
    file_in_zip_read_info_s *pfile_in_zip_read_info;
    if (file == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;
    pfile_in_zip_read_info = s->pfile_in_zip_read;

    if (pfile_in_zip_read_info == NULL)
        return UNZ_PARAMERROR;

    if (pfile_in_zip_read_info->rest_read_uncompressed == 0)
        return 1;
    else
        return 0;
}

/*
  Read extra field from the current file (opened by unzOpenCurrentFile)
  This is the local-header version of the extra field (sometimes, there is
    more info in the local-header version than in the central-header)

  if buf==NULL, it return the size of the local extra field that can be read

  if buf!=NULL, len is the size of the buffer, the extra header is copied in
    buf.
  the return value is the number of bytes copied in buf, or (if <0)
    the error code
*/
extern int ZEXPORT unzGetLocalExtrafield(file, buf, len) unzFile file;
voidp buf;
unsigned len;
{
    unz_s *s;
    file_in_zip_read_info_s *pfile_in_zip_read_info;
    uInt read_now;
    uLong size_to_read;

    if (file == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;
    pfile_in_zip_read_info = s->pfile_in_zip_read;

    if (pfile_in_zip_read_info == NULL)
        return UNZ_PARAMERROR;

    size_to_read = (pfile_in_zip_read_info->size_local_extrafield - pfile_in_zip_read_info->pos_local_extrafield);

    if (buf == NULL)
        return (int)size_to_read;

    if (len > size_to_read)
        read_now = (uInt)size_to_read;
    else
        read_now = (uInt)len;

    if (read_now == 0)
        return 0;

    if (ZSEEK(pfile_in_zip_read_info->z_filefunc, pfile_in_zip_read_info->filestream,
              pfile_in_zip_read_info->offset_local_extrafield + pfile_in_zip_read_info->pos_local_extrafield,
              ZLIB_FILEFUNC_SEEK_SET) != 0)
        return UNZ_ERRNO;

    if (ZREAD(pfile_in_zip_read_info->z_filefunc, pfile_in_zip_read_info->filestream, buf, read_now) != read_now)
        return UNZ_ERRNO;

    return (int)read_now;
}

/*
  Close the file in zip opened with unzipOpenCurrentFile
  Return UNZ_CRCERROR if all the file was read but the CRC is not good
*/
extern int ZEXPORT unzCloseCurrentFile(file) unzFile file;
{
    int err = UNZ_OK;

    unz_s *s;
    file_in_zip_read_info_s *pfile_in_zip_read_info;
    if (file == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;
    pfile_in_zip_read_info = s->pfile_in_zip_read;

    if (pfile_in_zip_read_info == NULL)
        return UNZ_PARAMERROR;

    if ((pfile_in_zip_read_info->rest_read_uncompressed == 0) && (!pfile_in_zip_read_info->raw))
    {
        if (pfile_in_zip_read_info->crc32 != pfile_in_zip_read_info->crc32_wait)
            err = UNZ_CRCERROR;
    }

    TRYFREE(pfile_in_zip_read_info->read_buffer);
    pfile_in_zip_read_info->read_buffer = NULL;
    if (pfile_in_zip_read_info->stream_initialised)
        inflateEnd(&pfile_in_zip_read_info->stream);

    pfile_in_zip_read_info->stream_initialised = 0;
    TRYFREE(pfile_in_zip_read_info);

    s->pfile_in_zip_read = NULL;

    return err;
}

/*
  Get the global comment string of the ZipFile, in the szComment buffer.
  uSizeBuf is the size of the szComment buffer.
  return the number of byte copied or an error code <0
*/
extern int ZEXPORT unzGetGlobalComment(file, szComment, uSizeBuf) unzFile file;
char *szComment;
uLong uSizeBuf;
{
    int err = UNZ_OK;
    unz_s *s;
    uLong uReadThis;
    if (file == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;

    uReadThis = uSizeBuf;
    if (uReadThis > s->gi.size_comment)
        uReadThis = s->gi.size_comment;

    if (ZSEEK(s->z_filefunc, s->filestream, s->central_pos + 22, ZLIB_FILEFUNC_SEEK_SET) != 0)
        return UNZ_ERRNO;

    if (uReadThis > 0)
    {
        *szComment = '\0';
        if (ZREAD(s->z_filefunc, s->filestream, szComment, uReadThis) != uReadThis)
            return UNZ_ERRNO;
    }

    if ((szComment != NULL) && (uSizeBuf > s->gi.size_comment))
        *(szComment + s->gi.size_comment) = '\0';
    return (int)uReadThis;
}

/* Additions by RX '2004 */
extern uLong ZEXPORT unzGetOffset(file) unzFile file;
{
    unz_s *s;

    if (file == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;
    if (!s->current_file_ok)
        return 0;
    if (s->gi.number_entry != 0 && s->gi.number_entry != 0xffff)
        if (s->num_file == s->gi.number_entry)
            return 0;
    return s->pos_in_central_dir;
}

extern int ZEXPORT unzSetOffset(file, pos) unzFile file;
uLong pos;
{
    unz_s *s;
    int err;

    if (file == NULL)
        return UNZ_PARAMERROR;
    s = (unz_s *)file;

    s->pos_in_central_dir = pos;
    s->num_file = s->gi.number_entry; /* hack */
    err = unzlocal_GetCurrentFileInfoInternal(file, &s->cur_file_info, &s->cur_file_info_internal, NULL, 0, NULL, 0,
                                              NULL, 0);
    s->current_file_ok = (err == UNZ_OK);
    return err;
}
