/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
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
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *
 *      MEMORY MANAGEMENT
 *
 *	Separate areas are maintained for fixed and variable sized
 *	objects.  The first of these is allocated as an array of
 *	SEXPRECs and the second as an array of VECRECs.  The fixed
 *	sized objects are assembled into a free list and cons cells
 *	are allocated from it.  When the list is exhausted, a
 *	mark-sweep garbarge collection takes place and the free list
 *	is rebuilt.  Variable size objects are allocated in the VECREC
 *	area.  During a garbage collection, these are compacted to the
 *	beginning of the VECREC array.
 *
 *	The top end of the VECREC array is also used by R_alloc to
 *	maintain a stack of non-relocatable memory blocks.  These are
 *	used in calls to .Fortran and .C (and for other temporary
 *	purposes).  They are freed using a stack discipline.
 *
 *	+---------------------------------------------------+
 *	| allocated vectors |	free	| R-alloc'ed blocks |
 *	+---------------------------------------------------+
 *	                    ^           ^
 *	                    |           |
 *	                    R_VTop      R_VMax
 *
 *	If a piece of code R-allocs some blocks, it is required to
 *	reset the R_VMax pointer back to its original value before it
 *	exits.  This can be done with the functions getvmax and
 *	setvmax.
 */

#include "Memory.h"
#include "Defn.h"
#include "Graphics.h"

static int gc_reporting = 0;
static int gc_count = 0;
int gc_inhibit_torture = 1; /* gets set to zero after initialisations */

/*
#define GC_TORTURE
*/

#ifdef GC_TORTURE
#define FORCE_GC !gc_inhibit_torture
#else
#define FORCE_GC 0
#endif

#define GC_PROT(X)                                                                                                     \
    {                                                                                                                  \
        int __t = gc_inhibit_torture;                                                                                  \
        gc_inhibit_torture = 1;                                                                                        \
        X;                                                                                                             \
        gc_inhibit_torture = __t;                                                                                      \
    }

void installIntVector(SEXP, int, FILE *);

SEXP do_gcinfo(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    int i;
    SEXP old = allocVector(LGLSXP, 1);

    checkArity(op, args);
    i = asLogical(CAR(args));
    LOGICAL(old)[0] = gc_reporting;
    if (i != NA_LOGICAL)
        gc_reporting = i;
    return old;
}

SEXP do_gc(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP value;
    int ogc;
    checkArity(op, args);
    ogc = gc_reporting;
    gc_reporting = asLogical(CAR(args));
    gc();
    gc_reporting = ogc;
    /*- now return the [free , total ] for cells and heap */
    PROTECT(value = allocVector(INTSXP, 4));
    INTEGER(value)[0] = R_Collected;
    INTEGER(value)[1] = (int)(R_VSize - (R_VTop - R_VHeap));
    INTEGER(value)[2] = R_NSize;
    INTEGER(value)[3] = R_VSize;
    UNPROTECT(1);
    return value;
}

void mem_err_heap(long size)
{
    error("heap memory (%ld Kb) exhausted [needed %ld Kb more]\n", (R_VSize * sizeof(VECREC)) / 1024,
          (size * sizeof(VECREC)) / 1024);
}

void mem_err_cons()
{
    error("cons memory (%ld cells) exhausted\n", R_NSize);
}

#ifdef Macintosh
Handle gStackH;
Handle gNHeapH;
Handle gVHeapH;

/* CleanUpRMemory : This routine releases the memory that R has */
/* allocated.  This is only needed for the Mac because the memory */
/* is in system memory so not naturally cleaned up at the end of */
/* the application execution. */

void CleanUpMemory(void)
{
    OSErr result;
    if (gStackH != nil)
        TempDisposeHandle(gStackH, &result);
    if (gNHeapH != nil)
        TempDisposeHandle(gNHeapH, &result);
    if (gVHeapH != nil)
        TempDisposeHandle(gVHeapH, &result);
}

#endif

/* InitMemory : Initialise the memory to be used in R. */
/* This includes: stack space, node space and vector space */
/* This is a ghastly mess and the Mac code needs to be separated. */

void InitMemory()
{
    int i;

    gc_reporting = R_Verbose;
#ifdef Macintosh
    OSErr result;

    gStackH = TempNewHandle(R_PPStackSize * sizeof(SEXP), &result);
    if ((gStackH == NULL) || (result != noErr))
        R_Suicide("couldn't allocate system memory for pointer stack");
    TempHLock(gStackH, &result);
    R_PPStack = (SEXP *)*gStackH;
#else
    if (!(R_PPStack = (SEXP *)malloc(R_PPStackSize * sizeof(SEXP))))
        R_Suicide("couldn't allocate memory for pointer stack");
#endif

    R_PPStackTop = 0;

#ifdef Macintosh
    gNHeapH = TempNewHandle(R_NSize * sizeof(SEXPREC), &result);
    if ((gNHeapH == NULL) || (result != noErr))
        R_Suicide("couldn't allocate system memory for node heap");
    TempHLock(gNHeapH, &result);
    R_NHeap = (SEXPREC *)*gNHeapH;
#else
    if (!(R_NHeap = (SEXPREC *)malloc(R_NSize * sizeof(SEXPREC))))
        R_Suicide("couldn't allocate memory for node heap");
#endif

    R_VSize = (((R_VSize + 1) / sizeof(VECREC)));

#ifdef Macintosh
    gVHeapH = TempNewHandle(R_VSize * sizeof(VECREC), &result);
    if ((gVHeapH == NULL) || (result != noErr))
        R_Suicide("couldn't allocate system memory for vector heap");
    TempHLock(gVHeapH, &result);
    R_VHeap = (VECREC *)*gVHeapH;
#else
#ifdef DEBUGGING
    printf("R_VSize = %d malloc-ed\n", R_VSize * sizeof(VECREC));
#endif
    if (!(R_VHeap = (VECREC *)malloc(R_VSize * sizeof(VECREC))))
        R_Suicide("couldn't allocate memory for vector heap");
#endif

    R_VTop = &R_VHeap[0];
    R_VMax = &R_VHeap[R_VSize - 1];

    for (i = 0; i < R_NSize - 1; i++)
        CDR(&R_NHeap[i]) = &R_NHeap[i + 1];
    CDR(&R_NHeap[R_NSize - 1]) = NULL;
    R_FreeSEXP = &R_NHeap[0];
}

char *vmaxget(void)
{
    return (char *)R_VMax;
}

void vmaxset(char *ovmax)
{
    if (ovmax)
        R_VMax = (VECREC *)ovmax;
    else
        R_VMax = &R_VHeap[R_VSize - 1];
}

char *R_alloc(long nelem, int eltsize)
{
    unsigned int size = BYTE2VEC(nelem * eltsize);
    if (size != 0)
    {
        if (FORCE_GC || R_VMax - R_VTop < size)
        {
            gc();
            if (R_VMax - R_VTop < size)
                mem_err_heap(size);
        }
        R_VMax -= size;
    }
    return (char *)R_VMax;
}

/* S COMPATIBILITY */

char *S_alloc(long nelem, int eltsize)
{
    unsigned int i, size = nelem * eltsize;
    char *p = R_alloc(nelem, eltsize);
    for (i = 0; i < size; i++)
        p[i] = 0;
    return p;
}

char *S_realloc(char *p, long new, long old, int size)
{
    int i, nold;
    char *q;
    /* shrinking is a no-op */
    if (new <= old)
        return p;
    q = R_alloc(new, size);
    nold = old * size;
    for (i = 0; i < nold; i++)
        q[i] = p[i];
    return q;
}

/* "allocSExp" allocate a SEXPREC from free list */
/* call gc if necessary */

SEXP allocSExp(SEXPTYPE t)
{
    SEXP s;
    if (FORCE_GC || R_FreeSEXP == NULL)
    {
        gc();
        if (R_FreeSEXP == NULL)
            mem_err_cons();
    }
    s = R_FreeSEXP;
    R_FreeSEXP = CDR(s);
    CAR(s) = R_NilValue;
    CDR(s) = R_NilValue;
    TAG(s) = R_NilValue;
    ATTRIB(s) = R_NilValue;
    *(int *)(&(s)->sxpinfo) = 0;
    TYPEOF(s) = t;
    return s;
}

/* "allocString" allocate a string on the (vector) heap. */
/* All vector objects  must be a multiple of sizeof(ALIGN) */
/* bytes so that alignment is preserved for all objects */

SEXP allocString(int length)
{
    SEXP s;
    long size;
    /* number of vector cells to allocate */
    size = 1 + BYTE2VEC(length + 1);
    /* we need to do the gc here so allocSExp doesn't! */
    if (FORCE_GC || R_FreeSEXP == NULL || R_VMax - R_VTop < size)
    {
        gc();
        if (R_FreeSEXP == NULL)
            mem_err_cons();
        if (R_VMax - R_VTop < size)
            mem_err_heap(size);
    }

    GC_PROT(s = allocSExp(CHARSXP));

    CHAR(s) = (char *)(R_VTop + 1);
    LENGTH(s) = length;
    BACKPOINTER(*R_VTop) = s;
    R_VTop += size;
    return s;
}

/* Allocate a vector object on the heap */

SEXP allocVector(SEXPTYPE type, int length)
{
    SEXP s;
    int i;
    long size = 0;
    if (length < 0)
        errorcall(R_GlobalContext->call, "negative length vectors are not allowed\n");
    /* number of vector cells to allocate */
    switch (type)
    {
    case NILSXP:
        return R_NilValue;
    case CHARSXP:
        size = 1 + BYTE2VEC(length + 1);
        break;
    case LGLSXP:
    case INTSXP:
        if (length <= 0)
            size = 0;
        else
            size = 1 + INT2VEC(length);
        break;
    case REALSXP:
        if (length <= 0)
            size = 0;
        else
            size = 1 + FLOAT2VEC(length);
        break;
    case CPLXSXP:
        if (length <= 0)
            size = 0;
        else
            size = 1 + COMPLEX2VEC(length);
        break;
    case STRSXP:
    case EXPRSXP:
    case VECSXP:
        if (length <= 0)
            size = 0;
        else
            size = 1 + PTR2VEC(length);
        break;
    case LANGSXP:
        if (length == 0)
            return R_NilValue;
        s = allocList(length);
        TYPEOF(s) = LANGSXP;
        return s;
    case LISTSXP:
        return allocList(length);
    default:
        error("invalid type/length (%d/%d) in vector allocation\n", type, length);
    }
    /* we need to do the gc here so allocSExp doesn't! */
    if (FORCE_GC || R_FreeSEXP == NULL || R_VMax - R_VTop < size)
    {
        gc();
        if (R_FreeSEXP == NULL)
            mem_err_cons();
        if (R_VMax - R_VTop < size)
            mem_err_heap(size);
    }
    GC_PROT(s = allocSExp(type));

    LENGTH(s) = length;
    NAMED(s) = 0;
    ATTRIB(s) = R_NilValue;
    if (size > 0)
    {
        CHAR(s) = (char *)(R_VTop + 1);
        BACKPOINTER(*R_VTop) = s;
        R_VTop += size;
    }
    else
        CHAR(s) = (char *)0;
    /* The following prevents disaster in the case */
    /* that an uninitialised string vector is marked */
    if (type == STRSXP || type == EXPRSXP || type == VECSXP)
    {
        for (i = 0; i < length; i++)
            STRING(s)[i] = R_NilValue;
    }
    return s;
}

SEXP allocList(int n)
{
    int i;
    SEXP result;
    result = R_NilValue;
    for (i = 0; i < n; i++)
    {
        result = CONS(R_NilValue, result);
    }
    return result;
}

/* "gc" a mark-sweep garbage collector */

void gc(void)
{
    sigset_t mask, omask;
    int vcells, vfrac;

    gc_count++;
    if (gc_reporting)
        REprintf("Garbage collection [nr. %d]...", gc_count);
    sigemptyset(&mask);
    sigaddset(&mask, SIGINT);
    sigprocmask(SIG_BLOCK, &mask, &omask);
    unmarkPhase();
    markPhase();
    compactPhase();
    scanPhase();
    sigprocmask(SIG_SETMASK, &omask, &mask);
    if (gc_reporting)
    {
        REprintf("\n%ld cons cells free (%ld%%)\n", R_Collected, (100 * R_Collected / R_NSize));
        vcells = R_VSize - (R_VTop - R_VHeap);
        vfrac = 100 * vcells / R_VSize;
        REprintf("%ld Kbytes of heap free (%ld%%)\n", vcells * sizeof(VECREC) / 1024, vfrac);
    }
}

/* "unmarkPhase" reset mark in ALL cons cells */

void unmarkPhase(void)
{
    int i;

    for (i = 0; i < R_NSize; i++)
        MARK(&R_NHeap[i]) = 0;
}

/* "markPhase" set mark in all accessible cons cells */

void markPhase(void)
{
    int i;
    DevDesc *dd;

    markSExp(R_NilValue); /* Builtin constants */
    markSExp(NA_STRING);
    markSExp(R_BlankString);
    markSExp(R_UnboundValue);
    markSExp(R_MissingArg);
    markSExp(R_CommentSxp);

    markSExp(R_GlobalEnv); /* Global environent */

    for (i = 0; i < HSIZE; i++) /* Symbol table */
        markSExp(R_SymbolTable[i]);

    if (R_CurrentExpr != NULL) /* Current expression */
        markSExp(R_CurrentExpr);

    for (i = 0; i < R_MaxDevices; i++)
    { /* Device display lists */
        dd = GetDevice(i);
        if (dd)
            markSExp(dd->displayList);
    }

    for (i = 0; i < R_PPStackTop; i++) /* Protected pointers */
        markSExp(R_PPStack[i]);
}

/* "markSExp" set mark in s and all cells accessible from it */

void markSExp(SEXP s)
{
    int i;

    if (s && !MARK(s))
    {
        MARK(s) = 1;
        if (ATTRIB(s) != R_NilValue)
            markSExp(ATTRIB(s));
        switch (TYPEOF(s))
        {
        case NILSXP:
        case BUILTINSXP:
        case SPECIALSXP:
        case CHARSXP:
        case LGLSXP:
        case INTSXP:
        case REALSXP:
        case CPLXSXP:
            break;
        case STRSXP:
        case EXPRSXP:
        case VECSXP:
            for (i = 0; i < LENGTH(s); i++)
                markSExp(STRING(s)[i]);
            break;
        case ENVSXP:
            markSExp(FRAME(s));
            markSExp(ENCLOS(s));
            break;
        case CLOSXP:
        case PROMSXP:
        case LISTSXP:
        case LANGSXP:
        case DOTSXP:
        case SYMSXP:
            markSExp(TAG(s));
            markSExp(CAR(s));
            markSExp(CDR(s));
            break;
        default:
            abort();
        }
    }
}

/* "compactPhase" compact the vector heap */

void compactPhase(void)
{
    VECREC *vto, *vfrom;
    SEXP s;
    int i, size;

    vto = vfrom = R_VHeap;

    while (vfrom < R_VTop)
    {
        s = BACKPOINTER(*vfrom);
        switch (TYPEOF(s))
        { /* get size in bytes */
        case CHARSXP:
            size = LENGTH(s) + 1;
            break;
        case LGLSXP:
        case INTSXP:
            size = LENGTH(s) * sizeof(int);
            break;
        case REALSXP:
            size = LENGTH(s) * sizeof(double);
            break;
        case CPLXSXP:
            size = LENGTH(s) * sizeof(complex);
            break;
        case STRSXP:
        case EXPRSXP:
        case VECSXP:
            size = LENGTH(s) * sizeof(SEXP);
            break;
        default:
            abort();
        }
        size = 1 + BYTE2VEC(size);
        if (MARK(s))
        {
            if (vfrom != vto)
            {
                for (i = 0; i < size; i++)
                    vto[i] = vfrom[i];
            }
            CHAR(BACKPOINTER(*vto)) = (char *)(vto + 1);
            vto += size;
            vfrom += size;
        }
        else
        {
            vfrom += size;
        }
    }
    R_VTop = vto;
}

/* "scanPhase" reconstruct free list from cells not marked */

void scanPhase(void)
{
    int i;

    R_FreeSEXP = NULL;
    R_Collected = 0;
    for (i = 0; i < R_NSize; i++)
    {
        if (!MARK(&R_NHeap[i]))
        {
            CDR(&R_NHeap[i]) = R_FreeSEXP;
            R_FreeSEXP = &R_NHeap[i];
            R_Collected++;
        }
    }
}

/* "protect" push a single argument onto R_PPStack */

void protect(SEXP s)
{
    if (R_PPStackTop >= R_PPStackSize)
        error("protect(): stack overflow\n");
    R_PPStack[R_PPStackTop] = s;
    R_PPStackTop++;
}

/* "unprotect" pop argument list from top of R_PPStack */

void unprotect(int l)
{
    if (R_PPStackTop > 0)
        R_PPStackTop = R_PPStackTop - l;
    else
        error("unprotect(): stack imbalance\n");
}

/* "unprotect_ptr" remove pointer from somewhere in R_PPStack */

void unprotect_ptr(SEXP s)
{
    int i = R_PPStackTop;

    /* go look for  s  in  R_PPStack */
    /* (should be among the top few items) */
    do
    {
        if (i == 0)
            error("unprotect_ptr: pointer not found\n");
    } while (R_PPStack[--i] != s);

    /* OK, got it, and  i  is indexing its location */
    /* Now drop stack above it */

    do
    {
        R_PPStack[i] = R_PPStack[i + 1];
    } while (i++ < R_PPStackTop);

    R_PPStackTop--;
}

/* "initStack" initialize environment stack */
void initStack(void)
{
    R_PPStackTop = 0;
}

/* Wrappers for malloc/alloc/free */
/* These allow automatic freeing of malloc-ed */
/* blocks during error recovery. */

#define MAXPOINTERS 100
static char *C_Pointers[MAXPOINTERS];

void Init_C_alloc()
{
    int i;
    for (i = 0; i < MAXPOINTERS; i++)
        C_Pointers[i] = NULL;
}

void Reset_C_alloc()
{
    int i;
    for (i = 0; i < MAXPOINTERS; i++)
    {
        if (C_Pointers[i] != NULL)
            free(C_Pointers[i]);
        C_Pointers[i] = NULL;
    }
}

char *C_alloc(long nelem, int eltsize)
{
    int i;
    for (i = 0; i < MAXPOINTERS; i++)
    {
        if (C_Pointers[i] == NULL)
        {
            C_Pointers[i] = malloc(nelem * eltsize);
            if (C_Pointers[i] == NULL)
                error("C_alloc(): unable to malloc memory\n");
            else
                return C_Pointers[i];
        }
    }
    error("C_alloc(): all pointers in use (sorry)\n");
    /*-Wall:*/ return C_Pointers[0];
}

void C_free(char *p)
{
    int i;
    for (i = 0; i < MAXPOINTERS; i++)
    {
        if (C_Pointers[i] == p)
        {
            free(p);
            C_Pointers[i] = NULL;
            return;
        }
    }
    error("C_free(): attempt to free pointer not allocated by C_alloc()\n");
}
