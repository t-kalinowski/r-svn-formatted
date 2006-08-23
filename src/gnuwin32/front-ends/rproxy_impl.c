/*******************************************************************************
 *  RProxy: Connector implementation between application and R language
 *  Copyright (C) 1999--2005 Thomas Baier
 *  Copyright 2006 R Development Core Team
 *
 *  R_Proxy_init based on rtest.c,  Copyright (C) 1998--2000
 *                                  R Development Core Team
 *
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Library General Public
 *  License as published by the Free Software Foundation; either
 *  version 2 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Library General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public
 *  License along with this library; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA 02110-1301, USA.
 *
 ******************************************************************************/

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>

#define NEW
#include <config.h>

#ifdef NEW
#include <Rinternals.h>
#else
#include <Defn.h>
#endif
#include <Rversion.h>
#include <Rembedded.h>
#include <R_ext/RStartup.h>
#include <R_ext/GraphicsDevice.h>
#include <graphapp.h>

#include "bdx_SEXP.h"
#include "SC_proxy.h"
#include "rproxy.h"
#include "rproxy_impl.h"

#ifdef NEW
#include <R_ext/Parse.h>
#else
/* <FIXME> Thees are private header files */
#include <IOStuff.h>
#include <Parse.h>
#endif

struct _R_Proxy_init_parameters
{
    int vsize;
    int vsize_valid;
    int nsize;
    int nsize_valid;
};

/* calls into the R DLL */
extern char *getRHOME();

int R_Proxy_Graphics_Driver(NewDevDesc *pDD, char *pDisplay, double pWidth, double pHeight, double pPointSize);

extern SC_CharacterDevice *__output_device;

/* trace to DebugView */

int R_Proxy_printf(char const *pFormat, ...)
{
    static char __tracebuf[2048];

    va_list lArgs;
    va_start(lArgs, pFormat);
    vsnprintf(__tracebuf, 2048, pFormat, lArgs);
    OutputDebugString(__tracebuf);
    return 0;
}

#ifndef NEW
static int s_EvalInProgress = 0;
#endif

static void R_Proxy_askok(char *pMsg)
{
    askok(pMsg);
    return;
}

static int R_Proxy_askyesnocancel(char *pMsg)
{
    return YES;
}

static int R_Proxy_ReadConsole(char *prompt, char *buf, int len, int addtohistory)
{
    return 0;
}

static void R_Proxy_WriteConsole(char *buf, int len)
{
    if (__output_device)
        __output_device->vtbl->write_string(__output_device, buf);
}

static void R_Proxy_CallBack()
{
    /* called during i/o, eval, graphics in ProcessEvents */
}

static void R_Proxy_Busy(int which)
{
    /* set a busy cursor ... in which = 1, unset if which = 0 */
}

#ifdef UNUSED
/* <FIXME> unused here */
int R_Proxy_parse_parameters(char const *pParameterString, struct _R_Proxy_init_parameters *pParameterStruct)
{
    /*
     * parameter string is of the form name1=value1;name2=value2;...
     *
     * currently recognized parameter names (case-sensitive):
     *
     *   NSIZE ... number of cons cells, (unsigned int) parameter
     *   VSIZE ... size of vector heap, (unsigned int) parameter
     */
    return 0;
}
#endif

/* 00-02-18 | baier | R_Proxy_init() now takes parameter string, parse it */
/* 03-06-01 | baier | now we add %R_HOME%\bin to %PATH% */
int R_Proxy_init(char const *pParameterString)
{
    structRstart rp;
    Rstart Rp = &rp;
    char Rversion[25];
    static char RHome[MAX_PATH];

    snprintf(Rversion, 25, "%s.%s", R_MAJOR, R_MINOR);
    if (strncmp(getDLLVersion(), Rversion, 25) != 0)
    {
        fprintf(stderr, "Error: R.DLL version does not match\n");
        return SC_PROXY_ERR_UNKNOWN;
    }

    R_DefParams(Rp);

    /* <FIXME> the documented interface is get_R_HOME() */

    /* first, try process-local environment space (CRT) */
    if (getenv("R_HOME"))
    {
        strcpy(RHome, getenv("R_HOME"));
    }
    else
    {
        /* get variable from process-local environment space (Windows API) */
        if (GetEnvironmentVariable("R_HOME", RHome, sizeof(RHome)) == 0)
        {
            /* not found, fall back to getRHOME() */
            strcpy(RHome, getRHOME());
        }
    }

    /* now we add %R_HOME%\bin to %PATH% (for dynamically loaded modules there) */
    {
        char buf[2048];
        snprintf(buf, 2048, "PATH=%s\\bin;%s", RHome, getenv("PATH"));
        putenv(buf);
    }

    Rp->rhome = RHome;
    Rp->home = getRUser();
    Rp->CharacterMode = LinkDLL;
    Rp->ReadConsole = R_Proxy_ReadConsole;
    Rp->WriteConsole = R_Proxy_WriteConsole;
    Rp->CallBack = R_Proxy_CallBack;
    Rp->ShowMessage = R_Proxy_askok;
    Rp->YesNoCancel = R_Proxy_askyesnocancel;
    Rp->Busy = R_Proxy_Busy;
    Rp->R_Quiet = 1;
    Rp->RestoreAction = SA_NORESTORE;
    Rp->SaveAction = SA_NOSAVE; /* had 2, with comment 'no save' which is 3 */

    R_SetParams(Rp);
    R_set_command_line_arguments(0, NULL);

    GA_initapp(0, 0);
    readconsolecfg();
    setup_Rmainloop();
    R_ReplDLLinit();

    return SC_PROXY_OK;
}

#ifdef NEW
int R_Proxy_evaluate(char const *pCmd, BDX_Data **pData)
{
    SEXP lSexp;
    int lRc = SC_PROXY_OK, evalError = 0;
    ParseStatus lStatus;
    SEXP lResult;

    lSexp = R_ParseVector(mkString(pCmd), 1, &lStatus);
    /* This is an EXPRSXP: we assume just one expression */

    switch (lStatus)
    {
    case PARSE_OK:
        PROTECT(lSexp);
        lResult = R_tryEval(VECTOR_ELT(lSexp, 0), R_GlobalEnv, &evalError);
        UNPROTECT(1);
        if (evalError)
            lRc = SC_PROXY_ERR_EVALUATE_STOP;
        else
            lRc = SEXP2BDX(lResult, pData);
        break;
    case PARSE_INCOMPLETE:
        lRc = SC_PROXY_ERR_PARSE_INCOMPLETE;
        break;
    default:
        lRc = SC_PROXY_ERR_PARSE_INVALID;
        break;
    }
    return lRc;
}

#else

/* 01-06-05 | baier | SETJMP and fatal error handling around eval() */
/* 04-08-01 | baier | ref-counting in case of error */
/* 04-10-11 | baier | restore original ref-counting */
/* 05-05-15 | baier | rework SETJMP code (store/restore jmp_buf) */
int R_Proxy_evaluate(char const *pCmd, BDX_Data **pData)
{
    SEXP rho = R_GlobalEnv;
    IoBuffer lBuffer;
    SEXP lSexp;
    int lRc;
    ParseStatus lStatus;
    SEXP lResult;

    /* for SETJMP/LONGJMP */
    s_EvalInProgress = 0;

    R_IoBufferInit(&lBuffer);
    R_IoBufferPuts((char *)pCmd, &lBuffer);
    R_IoBufferPuts("\n", &lBuffer);

    R_IoBufferReadReset(&lBuffer);
    lSexp = R_Parse1Buffer(&lBuffer, 1, &lStatus);
    PrintValue(lSexp);

    switch (lStatus)
    {
    case PARSE_OK:
        R_Visible = 0; /* Not printing, so not used */
        R_EvalDepth = 0;
        PROTECT(lSexp);
        {
            JMP_BUF lJmpBuf;
            memcpy(lJmpBuf, R_Toplevel.cjmpbuf, sizeof(lJmpBuf));
            SETJMP(R_Toplevel.cjmpbuf);
            R_GlobalContext = R_ToplevelContext = &R_Toplevel;

            if (!s_EvalInProgress)
            {
                /* <FIXME> This does not set .Last.value, does not
               print result and does not print warnings */
                s_EvalInProgress = 1;
                lResult = eval(lSexp, rho);
                memcpy(R_Toplevel.cjmpbuf, lJmpBuf, sizeof(lJmpBuf));
                s_EvalInProgress = 0;
            }
            else
            {
                memcpy(R_Toplevel.cjmpbuf, lJmpBuf, sizeof(lJmpBuf));
                return SC_PROXY_ERR_EVALUATE_STOP;
            }
        }
        lRc = SEXP2BDX(lResult, pData);
        /* no last value */
        UNPROTECT(1);
        break;
    case PARSE_INCOMPLETE:
        lRc = SC_PROXY_ERR_PARSE_INCOMPLETE;
        break;
    default:
        lRc = SC_PROXY_ERR_PARSE_INVALID;
        break;
    }

    return lRc;
}
#endif

#ifdef NEW
int R_Proxy_evaluate_noreturn(char const *pCmd)
{
    SEXP lSexp;
    int lRc = SC_PROXY_OK, evalError = 0;
    ParseStatus lStatus;
    SEXP lResult;

    lSexp = R_ParseVector(mkString(pCmd), 1, &lStatus);
    /* It would make sense to allow multiple expressions here */

    switch (lStatus)
    {
    case PARSE_OK:
        PROTECT(lSexp);
        lResult = R_tryEval(VECTOR_ELT(lSexp, 0), R_GlobalEnv, &evalError);
        UNPROTECT(1);
        if (evalError)
            lRc = SC_PROXY_ERR_EVALUATE_STOP;
        else
            lRc = SC_PROXY_OK;
        break;
    case PARSE_INCOMPLETE:
        lRc = SC_PROXY_ERR_PARSE_INCOMPLETE;
        break;
    default:
        lRc = SC_PROXY_ERR_PARSE_INVALID;
        break;
    }
    return lRc;
}

#else

/* 01-06-05 | baier | SETJMP and fatal error handling around eval() */
/* 04-08-01 | baier | ref-counting in case of error */
/* 04-10-11 | baier | restore original ref-counting */
/* 05-05-15 | baier | rework SETJMP code (store/restore jmp_buf) */
int R_Proxy_evaluate_noreturn(char const *pCmd)
{
    SEXP rho = R_GlobalEnv;
    IoBuffer lBuffer;
    SEXP lSexp;
    int lRc;
    ParseStatus lStatus;

    /* for SETJMP/LONGJMP */
    s_EvalInProgress = 0;

    R_IoBufferInit(&lBuffer);
    R_IoBufferPuts((char *)pCmd, &lBuffer);
    R_IoBufferPuts("\n", &lBuffer);

    R_IoBufferReadReset(&lBuffer);
    lSexp = R_Parse1Buffer(&lBuffer, 1, &lStatus);
    PrintValue(lSexp);

    switch (lStatus)
    {
    case PARSE_OK:
        R_Visible = 0;
        R_EvalDepth = 0;
        PROTECT(lSexp);
        {
            JMP_BUF lJmpBuf;
            memcpy(lJmpBuf, R_Toplevel.cjmpbuf, sizeof(lJmpBuf));
            SETJMP(R_Toplevel.cjmpbuf);
            R_GlobalContext = R_ToplevelContext = &R_Toplevel;

            if (!s_EvalInProgress)
            {
                s_EvalInProgress = 1;
                eval(lSexp, rho);
                memcpy(R_Toplevel.cjmpbuf, lJmpBuf, sizeof(lJmpBuf));
                s_EvalInProgress = 0;
            }
            else
            {
                memcpy(R_Toplevel.cjmpbuf, lJmpBuf, sizeof(lJmpBuf));
                return SC_PROXY_ERR_EVALUATE_STOP;
            }
        }
        /* no last value */
        UNPROTECT(1);
        lRc = SC_PROXY_OK;
        break;
    case PARSE_INCOMPLETE:
        lRc = SC_PROXY_ERR_PARSE_INCOMPLETE;
        break;
    default:
        lRc = SC_PROXY_ERR_PARSE_INVALID;
        break;
    }

    return lRc;
}
#endif

#ifdef NEW
int R_Proxy_get_symbol(char const *pSymbol, BDX_Data **pData)
{
    SEXP lVar = findVar(install((char *)pSymbol), R_GlobalEnv);

    if (lVar == R_UnboundValue)
    {
        RPROXY_TRACE(printf(">> %s is an unbound value\n", pSymbol));
        return SC_PROXY_ERR_INVALIDSYMBOL;
    }
    else if (SEXP2BDX(lVar, pData) == 0)
        return SC_PROXY_OK;
    else
        return SC_PROXY_ERR_UNSUPPORTEDTYPE;
}

#else

int R_Proxy_get_symbol(char const *pSymbol, BDX_Data **pData)
{
    IoBuffer lBuffer;
    SEXP lSexp;
    SEXP lVar;
    ParseStatus lStatus;

    R_IoBufferInit(&lBuffer);
    R_IoBufferPuts((char *)pSymbol, &lBuffer);
    R_IoBufferPuts("\n", &lBuffer);

    /* don't generate code, just a try */
    R_IoBufferReadReset(&lBuffer);
    lSexp = R_Parse1Buffer(&lBuffer, 0, &lStatus);

    if (lStatus == PARSE_OK)
    {
        /* now generate code */
        R_IoBufferReadReset(&lBuffer);
        lSexp = R_Parse1Buffer(&lBuffer, 1, &lStatus);
        R_Visible = 0;
        R_EvalDepth = 0;
        PROTECT(lSexp);

        /* check for valid symbol... */
        if (TYPEOF(lSexp) != SYMSXP)
        {
            RPROXY_TRACE(printf(">> %s is not a symbol\n", pSymbol));
            UNPROTECT(1);
            return SC_PROXY_ERR_INVALIDSYMBOL;
        }

        lVar = findVar(lSexp, R_GlobalEnv);

        if (lVar == R_UnboundValue)
        {
            RPROXY_TRACE(printf(">> %s is an unbound value\n", pSymbol));
            UNPROTECT(1);
            return SC_PROXY_ERR_INVALIDSYMBOL;
        }
        {
            int lRc = SEXP2BDX(lVar, pData);
            UNPROTECT(1);

            if (lRc == 0)
            {
                return SC_PROXY_OK;
            }
            else
            {
                return SC_PROXY_ERR_UNSUPPORTEDTYPE;
            }
        }
    }
    return SC_PROXY_OK; /* Really? - gets here on invalid input! */
}
#endif

/* 04-02-19 | baier | don't PROTECT strings in a vector, new data structs */
/* 04-03-02 | baier | removed traces */
/* 04-10-15 | baier | no more BDX_VECTOR (only BDX_ARRAY) */
/* 05-05-16 | baier | use BDX2SEXP, clean-up */
int R_Proxy_set_symbol(char const *pSymbol, BDX_Data const *pData)
{
    SEXP lSymbol = 0;
    SEXP lData = 0;

    /*  RPROXY_TRACE(printf("calling BDX2SEXP\n")); */
    if (BDX2SEXP(pData, &lData) != 0)
    {
        /*    RPROXY_TRACE(printf("error BDX2SEXP\n")); */
        return SC_PROXY_ERR_UNSUPPORTEDTYPE;
    }
    /*  RPROXY_TRACE(printf("ok BDX2SEXP\n")); */

    /* install a new symbol or get the existing symbol */
    lSymbol = install((char *)pSymbol);

    /* and set the data to the symbol */
    setVar(lSymbol, lData, R_GlobalEnv);

    return SC_PROXY_OK;
}

int R_Proxy_term()
{
    /* end_Rmainloop(); note, this never returns */
    Rf_endEmbeddedR(0);

    return SC_PROXY_OK;
}
