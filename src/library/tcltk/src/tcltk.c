#include <stdio.h>
#include <string.h>
#include <tcl.h>
#include <tk.h>

/*
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
*/

#include <Rinternals.h>
#include "R_ext/PrtUtil.h"
#ifndef Win32
#include "R_ext/eventloop.h"
#endif

/* From Parse.h -- must find better solution: */
#define PARSE_NULL 0
#define PARSE_OK 1
#define PARSE_INCOMPLETE 2
#define PARSE_ERROR 3
#define PARSE_EOF 4

SEXP R_ParseVector(SEXP, int, int *);

#include "tcltk.h" /* declarations of our `public' interface */
extern int (*R_timeout_handler)();
extern long R_timeout_val;
extern int (*ptr_gnome_start)();

static Tcl_Interp *RTcl_interp; /* Interpreter for this application. */

static void RTcl_dec_refcount(SEXP R_tclobj)
{
    Tcl_DecrRefCount((Tcl_Obj *)R_ExternalPtrAddr(R_tclobj));
}

static SEXP makeRTclObject(Tcl_Obj *tclobj)
{
    SEXP obj;

    obj = R_MakeExternalPtr(tclobj, R_NilValue, R_NilValue);
    Tcl_IncrRefCount(tclobj);
    R_RegisterCFinalizer(obj, RTcl_dec_refcount);
    return obj;
}

static int R_eval(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[])
{
    int status, i;
    SEXP text, expr, ans;

    text = PROTECT(allocVector(STRSXP, argc - 1));
    for (i = 1; i < argc; i++)
        SET_STRING_ELT(text, i - 1, mkChar(argv[i]));

    expr = PROTECT(R_ParseVector(text, -1, &status));
    if (status != PARSE_OK)
    {
        UNPROTECT(2);
        Tcl_SetResult(interp, "parse error in R expression", TCL_STATIC);
        return TCL_ERROR;
    }

    /* Note that expr becomes an EXPRSXP and hence we need the loop
       below (a straight eval(expr, R_GlobalEnv) won't work) */
    {
        int n = length(expr);
        for (i = 0; i < n; i++)
            ans = eval(VECTOR_ELT(expr, i), R_GlobalEnv);
    }
    UNPROTECT(2);
    return TCL_OK;
}

/* Functions to evaluate callbacks. Notice that these have no error
   checks, since the calls are assumed to be generated internally, and
   not by users. There are two forms:

   R_call 0xnnnnnnnn arg1 arg2 arg3...
   R_call_lang 0xnnnnnnnn 0xmmmmmmmm

   In the former, the hex address is assumed to correspond to a
   function closure which gets called with the arguments given.  The
   latter assumes that the first hex address is a LANGSXP and the
   second is an ENVSXP. */

static int R_call(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[])
{
    int i;
    SEXP expr, fun, alist;

    alist = R_NilValue;
    for (i = argc - 1; i > 1; i--)
    {
        PROTECT(alist);
        alist = LCONS(mkString(argv[i]), alist);
        UNPROTECT(1);
    }

    fun = (SEXP)strtoul(argv[1], NULL, 16);

    expr = LCONS(fun, alist);
    expr = LCONS(install("try"), LCONS(expr, R_NilValue));

    eval(expr, R_GlobalEnv);

    return TCL_OK;
}

static int R_call_lang(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[])
{
    SEXP expr, env;

    expr = (SEXP)strtoul(argv[1], NULL, 16);
    env = (SEXP)strtoul(argv[2], NULL, 16);

    expr = LCONS(install("try"), LCONS(expr, R_NilValue));

    eval(expr, env);

    return TCL_OK;
}

Tcl_Obj *tk_eval(char *cmd)
{
    if (Tcl_Eval(RTcl_interp, cmd) == TCL_ERROR)
    {
        char p[512];
        if (strlen(RTcl_interp->result) > 500)
            strcpy(p, "tcl error.\n");
        else
            sprintf(p, "[tcl] %s.\n", RTcl_interp->result);
        error(p);
    }
    return Tcl_GetObjResult(RTcl_interp);
}

/* FIXME get rid of fixed size buffers in do_Tclcallback et al. (low
   priority, since there is no meaningful way to override the current
   buffer length) */

SEXP dotTcl(SEXP args)
{
    SEXP ans;
    char *cmd;
    Tcl_Obj *val;
    if (!isValidString(CADR(args)))
        error("invalid argument");
    cmd = CHAR(STRING_ELT(CADR(args), 0));
    val = tk_eval(cmd);
    ans = makeRTclObject(val);
    return ans;
}

SEXP RTcl_ObjFromVar(SEXP args)
{
    Tcl_Obj *tclobj;

    tclobj = Tcl_GetVar2Ex(RTcl_interp, CHAR(STRING_ELT(CADR(args), 0)), NULL, 0);
    return makeRTclObject(tclobj);
}

SEXP RTcl_AssignObjToVar(SEXP args)
{
    Tcl_Obj *tclobj;

    tclobj =
        Tcl_SetVar2Ex(RTcl_interp, CHAR(STRING_ELT(CADR(args), 0)), NULL, (Tcl_Obj *)R_ExternalPtrAddr(CADDR(args)), 0);
    return R_NilValue;
}

SEXP RTcl_StringFromObj(SEXP args)
{
    char *str;

    str = Tcl_GetString((Tcl_Obj *)R_ExternalPtrAddr(CADR(args)));
    return mkString(str);
}

SEXP RTcl_ObjAsCharVector(SEXP args)
{
    int count;
    Tcl_Obj **elem;
    int ret, i;
    SEXP ans;

    ret = Tcl_ListObjGetElements(RTcl_interp, (Tcl_Obj *)R_ExternalPtrAddr(CADR(args)), &count, &elem);
    if (ret != TCL_OK)
        return RTcl_StringFromObj(args);

    PROTECT(ans = allocVector(STRSXP, count));
    for (i = 0; i < count; i++)
        SET_STRING_ELT(ans, i, mkChar(Tcl_GetString(elem[i])));
    UNPROTECT(1);
    return ans;
}

SEXP RTcl_ObjFromCharVector(SEXP args)
{
    int count;
    Tcl_Obj *tclobj, *elem;
    int ret, i;
    SEXP obj, val;

    val = CADR(args);

    tclobj = Tcl_NewObj();

    count = length(val);
    for (i = 0; i < count; i++)
    {
        elem = Tcl_NewObj();
        Tcl_SetStringObj(elem, CHAR(STRING_ELT(val, i)), -1);
        Tcl_ListObjAppendElement(RTcl_interp, tclobj, elem);
    }

    return makeRTclObject(tclobj);
}

SEXP RTcl_ObjAsDoubleVector(SEXP args)
{
    int count;
    Tcl_Obj **elem, *obj;
    int ret, i;
    double x;
    SEXP ans;

    obj = (Tcl_Obj *)R_ExternalPtrAddr(CADR(args));

    /* First try for single value */
    ret = Tcl_GetDoubleFromObj(RTcl_interp, obj, &x);
    if (ret == TCL_OK)
    {
        ans = allocVector(REALSXP, 1);
        REAL(ans)[0] = x;
        return ans;
    }

    /* Then try as list */
    ret = Tcl_ListObjGetElements(RTcl_interp, obj, &count, &elem);
    if (ret != TCL_OK) /* didn't work, return NULL */
        return R_NilValue;

    ans = allocVector(REALSXP, count);
    for (i = 0; i < count; i++)
    {
        ret = Tcl_GetDoubleFromObj(RTcl_interp, elem[i], &x);
        if (ret != TCL_OK)
            x = NA_REAL;
        REAL(ans)[i] = x;
    }
    return ans;
}

SEXP RTcl_ObjFromDoubleVector(SEXP args)
{
    int count;
    Tcl_Obj *tclobj, *elem;
    int ret, i;
    SEXP obj, val;

    val = CADR(args);

    tclobj = Tcl_NewObj();

    count = length(val);
    for (i = 0; i < count; i++)
    {
        elem = Tcl_NewDoubleObj(REAL(val)[i]);
        Tcl_ListObjAppendElement(RTcl_interp, tclobj, elem);
    }

    return makeRTclObject(tclobj);
}

SEXP RTcl_ObjAsIntVector(SEXP args)
{
    int count;
    Tcl_Obj **elem, *obj;
    int ret, i;
    int x;
    SEXP ans;

    obj = (Tcl_Obj *)R_ExternalPtrAddr(CADR(args));

    /* First try for single value */
    ret = Tcl_GetIntFromObj(RTcl_interp, obj, &x);
    if (ret == TCL_OK)
    {
        ans = allocVector(INTSXP, 1);
        INTEGER(ans)[0] = x;
        return ans;
    }

    /* Then try as list */
    ret = Tcl_ListObjGetElements(RTcl_interp, obj, &count, &elem);
    if (ret != TCL_OK) /* didn't work, return NULL */
        return R_NilValue;

    ans = allocVector(INTSXP, count);
    for (i = 0; i < count; i++)
    {
        ret = Tcl_GetIntFromObj(RTcl_interp, elem[i], &x);
        if (ret != TCL_OK)
            x = NA_REAL;
        INTEGER(ans)[i] = x;
    }
    return ans;
}

SEXP RTcl_ObjFromIntVector(SEXP args)
{
    int count;
    Tcl_Obj *tclobj, *elem;
    int ret, i;
    SEXP obj, val;

    val = CADR(args);

    tclobj = Tcl_NewObj();

    count = length(val);
    for (i = 0; i < count; i++)
    {
        elem = Tcl_NewIntObj(INTEGER(val)[i]);
        Tcl_ListObjAppendElement(RTcl_interp, tclobj, elem);
    }

    return makeRTclObject(tclobj);
}

/* Warning: These two functions return a pointer to internal static
   data. Copy immediately. */

static char *callback_closure(SEXP closure)
{
    static char buf[256], tmp[20];
    SEXP formals;

    formals = FORMALS(closure);

    sprintf(buf, "R_call 0x%lx", (unsigned long)closure);

    while (formals != R_NilValue)
    {
        if (TAG(formals) == R_DotsSymbol)
            break;
        sprintf(tmp, " %%%s", CHAR(PRINTNAME(TAG(formals))));
        strcat(buf, tmp);
        formals = CDR(formals);
    }
    return buf;
}

static char *callback_lang(SEXP call, SEXP env)
{
    static char buf[256];

    sprintf(buf, "R_call_lang 0x%lx 0x%lx", (unsigned long)call, (unsigned long)env);

    return buf;
}

/* Setup to invoke callback from Tcl. Notice that something needs to
   ensure that the callback is protected from GC as long as anything
   is around to invoke it. This is handled in interpreted code by
   assigning into the environment of the window with which the
   callback is associated */

SEXP dotTclcallback(SEXP args)
{
    SEXP ans, callback = CADR(args), env;

    char *rval = NULL; /* -Wall */

    if (isFunction(callback))
        rval = callback_closure(callback);
    else if (isLanguage(callback))
    {
        env = CADDR(args);
        rval = callback_lang(callback, env);
    }
    else
        error("argument is not of correct type");

    ans = mkString(rval);
    return ans;
}

#ifndef Win32
/* Add/delete Tcl/Tk event handler */

static void (*OldHandler)(void);
static int OldTimeout;
static int Tcl_loaded = 0;

void TclHandler(void)
{
    while (Tcl_DoOneEvent(TCL_DONT_WAIT))
        ;
    OldHandler();
}

int Gtk_TclHandler(void)
{
    while (Tcl_DoOneEvent(TCL_DONT_WAIT))
        ;
    return 1;
}

void addTcl(void)
{
    if (Tcl_loaded)
        error("Tcl already loaded");
    Tcl_loaded = 1;
    if (ptr_gnome_start != NULL)
    {
        R_timeout_handler = Gtk_TclHandler;
        R_timeout_val = 500;
    }
    else
    {
        OldHandler = R_PolledEvents;
        OldTimeout = R_wait_usec;
        R_PolledEvents = TclHandler;
        if (R_wait_usec > 10000 || R_wait_usec == 0)
            R_wait_usec = 10000;
    }
}

void delTcl(void)
{
    if (!Tcl_loaded)
        error("Tcl is not loaded");
    if (ptr_gnome_start != NULL)
    {
        R_timeout_handler = NULL;
        R_timeout_val = 0;
    }
    else
    {
        if (R_PolledEvents != TclHandler)
            error("Tcl is not last loaded handler");
        R_PolledEvents = OldHandler;
        R_wait_usec = OldTimeout;
    }
    Tcl_loaded = 0;
}
#endif

void tcltk_init(void)
{
    int code;

    RTcl_interp = Tcl_CreateInterp();
    code = Tcl_Init(RTcl_interp); /* Undocumented... If omitted, you
                    get the windows but no event
                    handling. */
    if (code != TCL_OK)
        error(RTcl_interp->result);

    code = Tk_Init(RTcl_interp); /* Load Tk into interpreter */
    if (code != TCL_OK)
        error(RTcl_interp->result);

    Tcl_StaticPackage(RTcl_interp, "Tk", Tk_Init, Tk_SafeInit);

    code = Tcl_Eval(RTcl_interp, "wm withdraw ."); /* Hide window */
    if (code != TCL_OK)
        error(RTcl_interp->result);

    Tcl_CreateCommand(RTcl_interp, "R_eval", R_eval, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

    Tcl_CreateCommand(RTcl_interp, "R_call", R_call, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

    Tcl_CreateCommand(RTcl_interp, "R_call_lang", R_call_lang, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

#ifdef Win32
    Tcl_SetServiceMode(TCL_SERVICE_ALL);
#else
    addTcl();
#endif

    /*** We may want to revive this at some point ***/

#if 0
  code = Tcl_EvalFile(RTcl_interp, "init.tcl");
  if (code != TCL_OK)
    error("%s\n", RTcl_interp->result);
#endif
}
