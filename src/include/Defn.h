/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998--2000  The R Development Core Team.
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

#ifndef DEFN_H_
#define DEFN_H_

#define COUNTING

#include "config.h"
#include "Rinternals.h"		/*-> Arith.h, Complex.h, Error.h, Memory.h
				  PrtUtil.h, Utils.h */
#include "Errormsg.h"

/* PSIGNAL may be defined on Win32 in config.h */
#ifdef PSIGNAL
#include <psignal.h>
#else
#include <signal.h>
#include <setjmp.h>
#endif
#include <time.h>
#ifdef HAVE_LOCALE_H
#include <locale.h>
#endif


/*  Heap and Pointer Protection Stack Sizes.  */

#define Mega 1048576. /* 1 Mega Byte := 2^20 (= 1048576) Bytes */

/*	R_PPSSIZE  The pointer protection stack size  */
/*	R_NSIZE	   The number of cons cells	 */
/*	R_VSIZE	   The vector heap size in bytes */
/*  These values are minima and can be overriden in config.h	*/

#ifndef R_PPSSIZE
#define	R_PPSSIZE	10000L
#endif
#ifndef R_NSIZE
#define	R_NSIZE		250000L
#endif
#ifndef R_VSIZE
#define	R_VSIZE		6291456L
#endif

#include <math.h>
#ifdef Macintosh
#define PosixArith
#define QUICKDRAW_GRAPHICS
#endif

/* all these are in Rinternals.h
#include <errno.h>
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
*/

/* Formerly in Arith.h */
#ifdef IEEE_754
# define MATH_CHECK(call)	(call)
#else
# ifdef __MAIN__
    double R_tmp;
# else
    extern double R_tmp;
# endif
#  define MATH_CHECK(call)	(errno=0,R_tmp=call,(errno==0)?R_tmp:R_NaN)
#endif


/* Getting the working directory */
#if defined(HAVE_GETCWD)
#define R_GETCWD(x, y) getcwd(x, y)
#elif defined(Win32)
#define R_GETCWD(x, y) GetCurrentDirectory(y, x)
#else
#undef R_GETCWD
#endif

/* Maximal length of an entire file name */
#if !defined(PATH_MAX)
#  if defined(HAVE_SYS_PARAM_H)
#    include <sys/param.h>
#  endif
#  if defined(MAXPATHLEN) && !defined(PATH_MAX)
#    define PATH_MAX MAXPATHLEN
#  elif defined(Win32)
#    define PATH_MAX 260
#  else
#    define PATH_MAX 255
#  endif
#endif

#ifdef HAVE_POSIX_SETJMP
# define JMP_BUF sigjmp_buf
# define SETJMP(x) sigsetjmp(x,1)
# define LONGJMP(x,i) siglongjmp(x,i)
#else
# define JMP_BUF jmp_buf
# define SETJMP(x) setjmp(x)
# define LONGJMP(x,i) longjmp(x,i)
#endif

#define HSIZE	   4119	/* The size of the hash table for symbols */
#define MAXELTSIZE 8192 /* The largest string size */
#define MAXIDSIZE   256	/* Largest symbol size possible */

/* The type of the do_xxxx functions. */
/* These are the built-in R functions. */
typedef SEXP (*CCODE)();


/* The type definitions for the table of built-in functions. */
/* This table can be found in ../main/names.c */
typedef struct {
    char   *name;    /* print name */
    CCODE  cfun;     /* c-code address */
    int	   code;     /* offset within c-code */
    int	   eval;     /* evaluate args? */
    int	   arity;    /* function arity */
    int	   gram;     /* pretty-print info */
} FUNTAB;

/* General Cons Cell Attributes */
#define ATTRIB(x)	((x)->attrib)
#define OBJECT(x)	((x)->sxpinfo.obj)
#define MARK(x)		((x)->sxpinfo.mark)
#define TYPEOF(x)	((x)->sxpinfo.type)
#define NAMED(x)	((x)->sxpinfo.named)

/* Primitive Access Macros */
#define PRIMOFFSET(x)	((x)->u.primsxp.offset)
#define PRIMFUN(x)	(R_FunTab[(x)->u.primsxp.offset].cfun)
#define PRIMNAME(x)	(R_FunTab[(x)->u.primsxp.offset].name)
#define PRIMVAL(x)	(R_FunTab[(x)->u.primsxp.offset].code)
#define PRIMARITY(x)	(R_FunTab[(x)->u.primsxp.offset].arity)
#define PPINFO(x)	(R_FunTab[(x)->u.primsxp.offset].gram)
#define PRIMPRINT(x)	(((R_FunTab[(x)->u.primsxp.offset].eval)/100)%10)

/* Promise Access Macros */
#define PREXPR(x)	((x)->u.promsxp.expr)
#define PRENV(x)	((x)->u.promsxp.env)
#define PRVALUE(x)	((x)->u.promsxp.value)
#define PRSEEN(x)	((x)->sxpinfo.gp)


/* Vector Heap Structure */
typedef struct {
	union {
		SEXP		backpointer;
		double		align;
	} u;
} VECREC, *VECP;

/* Vector Heap Macros */
#define BACKPOINTER(v)	((v).u.backpointer)
#define BYTE2VEC(n)	(((n)>0)?(((n)-1)/sizeof(VECREC)+1):0)
#define INT2VEC(n)	(((n)>0)?(((n)*sizeof(int)-1)/sizeof(VECREC)+1):0)
#define FLOAT2VEC(n)	(((n)>0)?(((n)*sizeof(double)-1)/sizeof(VECREC)+1):0)
#define COMPLEX2VEC(n)	(((n)>0)?(((n)*sizeof(Rcomplex)-1)/sizeof(VECREC)+1):0)
#define PTR2VEC(n)	(((n)>0)?(((n)*sizeof(SEXP)-1)/sizeof(VECREC)+1):0)


/* Evaluation Context Structure */
typedef struct RCNTXT {
    struct RCNTXT *nextcontext;	/* The next context up the chain */
    int callflag;		/* The context "type" */
    JMP_BUF cjmpbuf;		/* C stack and register information */
    int cstacktop;		/* Top of the pointer protection stack */
    SEXP promargs;		/* Promises supplied to closure */
    SEXP sysparent;		/* environment the closure was called from */
    SEXP call;			/* The call that effected this context*/
    SEXP cloenv;		/* The environment */
    SEXP conexit;		/* Interpreted "on.exit" code */
    void (*cend)();		/* C "on.exit" thunk */
} RCNTXT, *context;

/* The Various Context Types.

 * In general the type is a bitwise OR of the values below.
 * Note that CTXT_LOOP is already the or of CTXT_NEXT and CTXT_BREAK.
 * Only functions should have the third bit turned on;
 * this allows us to move up the context stack easily
 * with either RETURN's or GENERIC's or RESTART's.
 * If you add a new context type for functions make sure
 *   CTXT_NEWTYPE & CTXT_FUNCTION > 0
 */
enum {
    CTXT_TOPLEVEL = 0,
    CTXT_NEXT	  = 1,
    CTXT_BREAK	  = 2,
    CTXT_LOOP	  = 3,	/* break OR next target */
    CTXT_FUNCTION = 4,
    CTXT_CCODE	  = 8,
    CTXT_RETURN	  = 12,
    CTXT_BROWSER  = 16,
    CTXT_GENERIC  = 20,
    CTXT_RESTART  = 28
};

/*
TOP   0 0 0 0 0 0  = 0
NEX   1 0 0 0 0 0  = 1
BRE   0 1 0 0 0 0  = 2
LOO   1 1 0 0 0 0  = 3
FUN   0 0 1 0 0 0  = 4
CCO   0 0 0 1 0 0  = 8
BRO   0 0 0 0 1 0  = 16
RET   0 0 1 1 0 0  = 12
GEN   0 0 1 0 1 0  = 20
RES   0 0 1 1 1 0  = 28
*/
/* Miscellaneous Definitions */
#define streql(s, t)	(!strcmp((s), (t)))

/* Arithmetic and Relation Operators */
#define	PLUSOP	1
#define	MINUSOP	2
#define	TIMESOP	3
#define	DIVOP	4
#define	POWOP	5
#define	MODOP	6
#define IDIVOP	7

#define	EQOP	1
#define	NEOP	2
#define	LTOP	3
#define	LEOP	4
#define	GEOP	5
#define	GTOP	6

/* File Handling */
/*
#define R_EOF	65535
*/
#define R_EOF	-1

/* MAGIC Numbers for files */
#define R_MAGIC_BINARY 1975
#define R_MAGIC_ASCII  1976
#define R_MAGIC_XDR    1977

#define R_MAGIC_BINARY_VERSION16 1971
#define R_MAGIC_ASCII_VERSION16	 1972

/* Startup Actions */

#define SA_NORESTORE 0
#define SA_RESTORE   1

#define SA_DEFAULT   1
#define SA_NOSAVE    2
#define SA_SAVE	     3
#define SA_SAVEASK   4
#define SA_SUICIDE   5


/*--- Global Variables ---------------------------------------------------- */

/* Defined and initialized in names.c (not main.c) :*/
#ifndef __R_Names__
extern
#endif
FUNTAB	R_FunTab[];	    /* Built in functions */


#ifdef __MAIN__
#define extern
#define INI_as(v) = v
#else
#define INI_as(v)
#endif

/* extern int	errno; already have errno.h ! */
extern int	gc_inhibit_torture INI_as(1);

/* R Home Directory */
extern char*	R_Home;		    /* Root of the R tree */

/* Memory Management */
extern int	R_NSize		INI_as(R_NSIZE);/* Size of cons cell heap */
extern int	R_VSize		INI_as(R_VSIZE);/* Size of the vector heap */
extern SEXP	R_NHeap;	    /* Start of the cons cell heap */
extern SEXP	R_FreeSEXP;	    /* Cons cell free list */
extern VECREC*	R_VHeap;	    /* Base of the vector heap */
extern VECREC*	R_VTop;		    /* Current top of the vector heap */
extern VECREC*	R_VMax;		    /* bottom of R_alloc'ed heap */
extern long	R_Collected;	    /* Number of free cons cells (after gc) */
extern SEXP	R_PreciousList;	    /* List of Persistent Objects */

/* The Pointer Protection Stack */
extern int	R_PPStackSize	INI_as(R_PPSSIZE); /* The stack size (elements) */
extern int	R_PPStackTop;	    /* The top of the stack */
extern SEXP*	R_PPStack;	    /* The pointer protection stack */

/* Evaluation Environment */
extern SEXP	R_Call;		    /* The current call */
extern SEXP	R_CurrentExpr;	    /* Currently evaluating expression */
extern SEXP	R_ReturnedValue;    /* Slot for return-ing values */
extern SEXP*	R_SymbolTable;	    /* The symbol table */
extern RCNTXT	R_Toplevel;	    /* Storage for the toplevel environment */
extern RCNTXT*	R_ToplevelContext;  /* The toplevel environment */
extern RCNTXT*	R_GlobalContext;    /* The global environment */
extern int	R_Visible;	    /* Value visibility flag */
extern int	R_EvalDepth	INI_as(0);	/* Evaluation recursion depth */
extern int	R_EvalCount	INI_as(0);	/* Evaluation count */
extern int	R_BrowseLevel	INI_as(0);	/* how deep the browser is */

/* File Input/Output */
extern int	R_Interactive	INI_as(1);	/* Non-zero during interactive use */
extern int	R_Quiet		INI_as(0);	/* Be as quiet as possible */
extern int	R_Slave		INI_as(0);	/* Run as a slave process */
extern int	R_Verbose	INI_as(0);	/* Be verbose */
/* extern int	R_Console; */	    /* Console active flag */
/* IoBuffer R_ConsoleIob; : --> ./IOStuff.h */
extern FILE*	R_Inputfile	INI_as(NULL);	/* Current input flag */
extern FILE*	R_Consolefile	INI_as(NULL);	/* Console output file */
extern FILE*	R_Outputfile	INI_as(NULL);	/* Output file */
extern FILE*	R_Sinkfile	INI_as(NULL);	/* Sink file */

/* Objects Used In Parsing  */
extern SEXP	R_CommentSxp;	    /* Comments accumulate here */
extern SEXP	R_ParseText;	    /* Text to be parsed */
extern int	R_ParseCnt;	    /* Count of lines of text to be parsed */
extern int	R_ParseError	INI_as(0); /* Line where parse error occured */

/* Image Dump/Restore */
extern char	R_ImageName[256];   /* Default image name */
extern int	R_Unnamed	INI_as(1);	/* Use default name? */
extern int	R_DirtyImage	INI_as(0);	/* Current image dirty */
extern int	R_Init		INI_as(0);	/* Do we have an image loaded */
/* extern FILE*	R_FileRef;	    the environment file pointer  */

/* History */
extern char*	R_HistoryFile;	/* Name of the history file */
extern int	R_HistorySize;	/* Size of the history file */

/* Warnings/Errors */
extern int	R_CollectWarnings INI_as(0);	/* the number of warnings */
extern SEXP	R_Warnings;	    /* the warnings and their calls */



extern char ** CommandLineArgs	  INI_as(NULL); /* permanent copy of the command line arguments passed to the application. */
extern int     NumCommandLineArgs INI_as(0); /* the number of command line arguments. */

#ifdef __MAIN__
#undef extern
#endif
#undef INI_as


/*--- FUNCTIONS ------------------------------------------------------ */

#ifdef R_REMAP
#define applyRelOp		Rf_applyRelOp
#define begincontext		Rf_begincontext
#define checkArity		Rf_checkArity
#define CheckFormals		Rf_CheckFormals
#define classgets		Rf_classgets
#define CleanEd			Rf_CleanEd
#define compactPhase		Rf_compactPhase
#define DataFrameClass		Rf_DataFrameClass
#define ddfindVar		Rf_ddfindVar
#define deparse1		Rf_deparse1
#define deparse1line		Rf_deparse1line
#define DispatchGroup		Rf_DispatchGroup
#define DispatchOrEval		Rf_DispatchOrEval
#define DropDims		Rf_DropDims
#define duplicated		Rf_duplicated
#define dynamicfindVar		Rf_dynamicfindVar
#define endcontext		Rf_endcontext
#define factorsConform		Rf_factorsConform
#define findcontext		Rf_findcontext
#define findVar1		Rf_findVar1
#define findVarInFrame		Rf_findVarInFrame
#define findVarLocInFrame	Rf_findVarLocInFrame
#define FrameClassFix		Rf_FrameClassFix
#define framedepth		Rf_framedepth
#define frameSubscript		Rf_frameSubscript
#define get1index		Rf_get1index
#define getVar			Rf_getVar
#define getVarInFrame		Rf_getVarInFrame
#define hashpjw			Rf_hashpjw
#define InitArithmetic		Rf_InitArithmetic
#define InitColors		Rf_InitColors
#define InitEd			Rf_InitEd
#define InitFunctionHashing	Rf_InitFunctionHashing
#define InitGlobalEnv		Rf_InitGlobalEnv
#define InitMemory		Rf_InitMemory
#define InitNames		Rf_InitNames
#define InitOptions		Rf_InitOptions
#define initStack		Rf_initStack
#define internalTypeCheck	Rf_internalTypeCheck
#define isValidName		Rf_isValidName
#define jump_to_toplevel	Rf_jump_to_toplevel
#define levelsgets		Rf_levelsgets
#define mainloop		Rf_mainloop
#define markPhase		Rf_markPhase
#define markSExp		Rf_markSExp
#define mat2indsub		Rf_mat2indsub
#define match			Rf_match
#define mkCLOSXP		Rf_mkCLOSXP
#define mkFalse			Rf_mkFalse
#define mkPROMISE		Rf_mkPROMISE
#define mkQUOTE			Rf_mkQUOTE
#define mkSYMSXP		Rf_mkSYMSXP
#define NewEnvironment		Rf_NewEnvironment
#define OneIndex		Rf_OneIndex
#define onintr			Rf_onintr
#define parse			Rf_parse
#define PrintGreeting		Rf_PrintGreeting
#define PrintVersion		Rf_PrintVersion
#define PrintWarnings		Rf_PrintWarnings
#define promiseArgs		Rf_promiseArgs
#define scanPhase		Rf_scanPhase
#define setVarInFrame		Rf_setVarInFrame
#define sortVector		Rf_sortVector
#define ssort			Rf_ssort
#define str2type		Rf_str2type
#define StrToInternal		Rf_StrToInternal
#define substituteList		Rf_substituteList
#define tsConform		Rf_tsConform
#define tspgets			Rf_tspgets
#define type2str		Rf_type2str
#define unbindVar		Rf_unbindVar
#define unmarkPhase		Rf_unmarkPhase
#define usemethod		Rf_usemethod
#define warningcall		Rf_warningcall
#define WarningMessage		Rf_WarningMessage
#define yyerror			Rf_yyerror
#define yyinit			Rf_yyinit
#define yylex			Rf_yylex
#define yyparse			Rf_yyparse
#define yyprompt		Rf_yyprompt
#define yywrap			Rf_yywrap
#endif

/* Platform Dependent Gui Hooks */

#define	R_CONSOLE	1
#define	R_FILE		2
#define R_TEXT		3

int	R_ReadConsole(char*, unsigned char*, int, int);
void	R_WriteConsole(char*, int);
void	R_ResetConsole(void);
void	R_FlushConsole(void);
void	R_ClearerrConsole(void);
void	R_Busy(int);
void	R_CleanUp(int, int, int);
void	R_StartUp(void);
int	R_ShowFile(char*, char*);
int	R_ShowFiles(int, char **, char **, char *, int, char *);
int	R_ChooseFile(int, char*, int);
char*	R_HomeDir(void);
int	R_HiddenFile(char*);
char*	R_Date(void);
int	R_FileExists(char*);

/* Other Internally Used Functions */

SEXP applyRelOp(int, int, int);
void begincontext(RCNTXT*, int, SEXP, SEXP, SEXP, SEXP);
void checkArity(SEXP, SEXP);
void CheckFormals(SEXP);
SEXP classgets(SEXP, SEXP);
void CleanEd(void);
#ifdef Macintosh
	void CleanUpMemory( void );
#endif
void compactPhase(void);
void DataFrameClass(SEXP);
SEXP ddfindVar(SEXP, SEXP);
SEXP deparse1(SEXP,int);
SEXP deparse1line(SEXP,int);
int DispatchOrEval(SEXP, SEXP, SEXP, SEXP, SEXP*, int);
int DispatchGroup(char*, SEXP,SEXP,SEXP,SEXP,SEXP*);
SEXP DropDims(SEXP);
SEXP duplicated(SEXP);
SEXP dynamicfindVar(SEXP, RCNTXT*);
void endcontext(RCNTXT*);
int factorsConform(SEXP, SEXP);
void findcontext(int, SEXP, SEXP);
SEXP findVar1(SEXP, SEXP, SEXPTYPE, int);
SEXP findVarInFrame(SEXP, SEXP);
SEXP findVarLocInFrame(SEXP, SEXP);
void FrameClassFix(SEXP);
int framedepth(RCNTXT*);
SEXP frameSubscript(int, SEXP, SEXP);
int get1index(SEXP, SEXP, int, int);
SEXP getVar(SEXP, SEXP);
SEXP getVarInFrame(SEXP, SEXP);
int hashpjw(char*);
void InitArithmetic(void);
void InitColors(void);
void InitEd(void);
void InitFunctionHashing(void);
void InitGlobalEnv(void);
void R_InitialData(void);
void InitMemory(void);
void InitNames(void);
void InitOptions(void);
void initStack(void);
void internalTypeCheck(SEXP, SEXP, SEXPTYPE);
int isValidName(char *);
void jump_to_toplevel(void);
SEXP levelsgets(SEXP, SEXP);
void mainloop(void);
void markPhase(void);
void markSExp(SEXP);
SEXP mat2indsub(SEXP, SEXP);
SEXP match(SEXP, SEXP, int);
SEXP mkCLOSXP(SEXP, SEXP, SEXP);
/* SEXP mkEnv(SEXP, SEXP, SEXP); */
SEXP mkPRIMSXP (int, int);
SEXP mkPROMISE(SEXP, SEXP);
SEXP mkQUOTE(SEXP);
SEXP mkSYMSXP(SEXP, SEXP);
SEXP mkFalse(void);
SEXP NewEnvironment(SEXP, SEXP, SEXP);
void onintr();
int OneIndex(SEXP, SEXP, int, int, SEXP*);
SEXP parse(FILE*, int);
void PrintGreeting(void);
void PrintVersion(char *);
void PrintWarnings(void);
SEXP promiseArgs(SEXP, SEXP);
SEXP R_LoadFromFile(FILE*, int);
FILE* R_OpenLibraryFile(char *);
void R_PreserveObject(SEXP);
void R_ReleaseObject(SEXP);
void R_RestoreGlobalEnv(void);
void R_SaveGlobalEnv(void);
void R_SaveToFile(SEXP, FILE*, int);
int R_SetOptionWarn(int);
int R_SetOptionWidth(int);
void R_Suicide(char*);
void scanPhase(void);
SEXP setVarInFrame(SEXP, SEXP, SEXP);
void sortVector(SEXP);
void ssort(SEXP*,int);
SEXPTYPE str2type(char*);
int StrToInternal(char*);
SEXP substituteList(SEXP, SEXP);
SEXP R_syscall(int,RCNTXT*);
int R_sysparent(int,RCNTXT*);
SEXP R_sysframe(int,RCNTXT*);
SEXP R_sysfunction(int,RCNTXT*);
int tsConform(SEXP,SEXP);
SEXP tspgets(SEXP, SEXP);
SEXP type2str(SEXPTYPE);
void unbindVar(SEXP, SEXP);
void unmarkPhase(void);
int usemethod(char*, SEXP, SEXP, SEXP, SEXP, SEXP*);
void warningcall(SEXP, char*,...);
void WarningMessage(SEXP, int, ...);

/* gram.y & gram.c : */
void yyerror(char *);
void yyinit(void);
int yylex();
int yyparse(void);
void yyprompt(char *format, ...);
int yywrap(void);

#endif
/*
 *- Local Variables:
 *- page-delimiter: "^/\\*---"
 *- End:
 */
