/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 1

/* Tokens.  */
#ifndef YYTOKENTYPE
#define YYTOKENTYPE
/* Put the tokens into the symbol table, so that GDB and other debuggers
   know about them.  */
enum yytokentype
{
    END_OF_INPUT = 258,
    ERROR = 259,
    SECTIONHEADER = 260,
    RSECTIONHEADER = 261,
    VSECTIONHEADER = 262,
    SECTIONHEADER2 = 263,
    RCODEMACRO = 264,
    LATEXMACRO = 265,
    VERBMACRO = 266,
    OPTMACRO = 267,
    ESCAPE = 268,
    LISTSECTION = 269,
    ITEMIZE = 270,
    DESCRIPTION = 271,
    NOITEM = 272,
    RCODEMACRO2 = 273,
    LATEXMACRO2 = 274,
    VERBMACRO2 = 275,
    IFDEF = 276,
    ENDIF = 277,
    TEXT = 278,
    RCODE = 279,
    VERB = 280,
    COMMENT = 281,
    UNKNOWN = 282
};
#endif
/* Tokens.  */
#define END_OF_INPUT 258
#define ERROR 259
#define SECTIONHEADER 260
#define RSECTIONHEADER 261
#define VSECTIONHEADER 262
#define SECTIONHEADER2 263
#define RCODEMACRO 264
#define LATEXMACRO 265
#define VERBMACRO 266
#define OPTMACRO 267
#define ESCAPE 268
#define LISTSECTION 269
#define ITEMIZE 270
#define DESCRIPTION 271
#define NOITEM 272
#define RCODEMACRO2 273
#define LATEXMACRO2 274
#define VERBMACRO2 275
#define IFDEF 276
#define ENDIF 277
#define TEXT 278
#define RCODE 279
#define VERB 280
#define COMMENT 281
#define UNKNOWN 282

/* Copy the first part of user declarations.  */

/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996, 1997  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997--2008  Robert Gentleman, Ross Ihaka and the
 *                            R Development Core Team
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
#include "Parse.h"

#define DEBUGVALS 0 /* 1 causes detailed internal state output to R console */
#define DEBUGMODE 0 /* 1 causes Bison output of parse state, to stdout or stderr */

#define YYERROR_VERBOSE 1

static void yyerror(char *);
static int yylex();
static int yyparse(void);

#define yyconst const

typedef struct yyltype
{
    int first_line;
    int first_column;

    int last_line;
    int last_column;
} yyltype;

#define YYLTYPE yyltype

/* Useful defines so editors don't get confused ... */

#define LBRACE '{'
#define RBRACE '}'

/* Functions used in the parsing process */

static SEXP GrowList(SEXP, SEXP);
static int KeywordLookup(const char *);
static SEXP NewList(void);
static SEXP makeSrcref(YYLTYPE *, SEXP);

/* Internal lexer / parser state variables */

static int xxinRString, xxQuoteLine, xxQuoteCol;
static int xxgetc();
static int xxungetc(int);
static int xxlineno, xxcolno;
static int xxlastlinelen;
static int xxmode, xxitemType, xxbraceDepth; /* context for lexer */
static int xxDebugTokens;                    /* non-zero causes debug output to R console */
static SEXP Value;

#define RLIKE 1 /* Includes R strings; xxinRString holds the opening quote char, or 0 outside a string */
#define LATEXLIKE 2
#define VERBATIM 3
#define INOPTION 4

static SEXP SrcFile = NULL;

/* Routines used to build the parse tree */

static SEXP xxpushMode(int, int);
static void xxpopMode(SEXP);
static SEXP xxnewlist(SEXP);
static SEXP xxlist(SEXP, SEXP);
static SEXP xxmarkup(SEXP, SEXP, YYLTYPE *);
static SEXP xxmarkup2(SEXP, SEXP, SEXP, YYLTYPE *);
static SEXP xxOptionmarkup(SEXP, SEXP, SEXP, YYLTYPE *);
static SEXP xxtag(SEXP, int, YYLTYPE *);
static void xxsavevalue(SEXP, YYLTYPE *);

static int mkMarkup(int);
static int mkIfdef(int);
static int mkCode(int);
static int mkText(int);
static int mkVerb(int);
static int mkComment(int);

#define YYSTYPE SEXP

/* Enabling traces.  */
#ifndef YYDEBUG
#define YYDEBUG 1
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
#undef YYERROR_VERBOSE
#define YYERROR_VERBOSE 1
#else
#define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
#define YYTOKEN_TABLE 0
#endif

#if !defined YYSTYPE && !defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
#define yystype YYSTYPE /* obsolescent; will be withdrawn */
#define YYSTYPE_IS_DECLARED 1
#define YYSTYPE_IS_TRIVIAL 1
#endif

#if !defined YYLTYPE && !defined YYLTYPE_IS_DECLARED
typedef struct YYLTYPE
{
    int first_line;
    int first_column;
    int last_line;
    int last_column;
} YYLTYPE;
#define yyltype YYLTYPE /* obsolescent; will be withdrawn */
#define YYLTYPE_IS_DECLARED 1
#define YYLTYPE_IS_TRIVIAL 1
#endif

/* Copy the second part of user declarations.  */

/* Line 216 of yacc.c.  */

#ifdef short
#undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
#ifdef __SIZE_TYPE__
#define YYSIZE_T __SIZE_TYPE__
#elif defined size_t
#define YYSIZE_T size_t
#elif !defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
#include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#define YYSIZE_T size_t
#else
#define YYSIZE_T unsigned int
#endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T)-1)

#ifndef YY_
#if defined YYENABLE_NLS && YYENABLE_NLS
#if ENABLE_NLS
#include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#define YY_(msgid) dgettext("bison-runtime", msgid)
#endif
#endif
#ifndef YY_
#define YY_(msgid) msgid
#endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if !defined lint || defined __GNUC__
#define YYUSE(e) ((void)(e))
#else
#define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
#define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
static int YYID(int i)
#else
static int YYID(i) int i;
#endif
{
    return i;
}
#endif

#if !defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

#ifdef YYSTACK_USE_ALLOCA
#if YYSTACK_USE_ALLOCA
#ifdef __GNUC__
#define YYSTACK_ALLOC __builtin_alloca
#elif defined __BUILTIN_VA_ARG_INCR
#include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#elif defined _AIX
#define YYSTACK_ALLOC __alloca
#elif defined _MSC_VER
#include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#define alloca _alloca
#else
#define YYSTACK_ALLOC alloca
#if !defined _ALLOCA_H && !defined _STDLIB_H &&                                                                        \
    (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
#include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#ifndef _STDLIB_H
#define _STDLIB_H 1
#endif
#endif
#endif
#endif
#endif

#ifdef YYSTACK_ALLOC
/* Pacify GCC's `empty if-body' warning.  */
#define YYSTACK_FREE(Ptr)                                                                                              \
    do                                                                                                                 \
    { /* empty */                                                                                                      \
        ;                                                                                                              \
    } while (YYID(0))
#ifndef YYSTACK_ALLOC_MAXIMUM
/* The OS might guarantee only one guard page at the bottom of the stack,
   and a page size can be as small as 4096 bytes.  So we cannot safely
   invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
   to allow for a few compiler-allocated temporary stack slots.  */
#define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#endif
#else
#define YYSTACK_ALLOC YYMALLOC
#define YYSTACK_FREE YYFREE
#ifndef YYSTACK_ALLOC_MAXIMUM
#define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#endif
#if (defined __cplusplus && !defined _STDLIB_H &&                                                                      \
     !((defined YYMALLOC || defined malloc) && (defined YYFREE || defined free)))
#include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#ifndef _STDLIB_H
#define _STDLIB_H 1
#endif
#endif
#ifndef YYMALLOC
#define YYMALLOC malloc
#if !defined malloc && !defined _STDLIB_H &&                                                                           \
    (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
void *malloc(YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#endif
#endif
#ifndef YYFREE
#define YYFREE free
#if !defined free && !defined _STDLIB_H &&                                                                             \
    (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
void free(void *);      /* INFRINGES ON USER NAME SPACE */
#endif
#endif
#endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */

#if (!defined yyoverflow && (!defined __cplusplus || (defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL &&              \
                                                      defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc {
    yytype_int16 yyss;
    YYSTYPE yyvs;
    YYLTYPE yyls;
};

/* The size of the maximum gap between one aligned stack and the next.  */
#define YYSTACK_GAP_MAXIMUM (sizeof(union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
#define YYSTACK_BYTES(N) ((N) * (sizeof(yytype_int16) + sizeof(YYSTYPE) + sizeof(YYLTYPE)) + 2 * YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
#ifndef YYCOPY
#if defined __GNUC__ && 1 < __GNUC__
#define YYCOPY(To, From, Count) __builtin_memcpy(To, From, (Count) * sizeof(*(From)))
#else
#define YYCOPY(To, From, Count)                                                                                        \
    do                                                                                                                 \
    {                                                                                                                  \
        YYSIZE_T yyi;                                                                                                  \
        for (yyi = 0; yyi < (Count); yyi++)                                                                            \
            (To)[yyi] = (From)[yyi];                                                                                   \
    } while (YYID(0))
#endif
#endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
#define YYSTACK_RELOCATE(Stack)                                                                                        \
    do                                                                                                                 \
    {                                                                                                                  \
        YYSIZE_T yynewbytes;                                                                                           \
        YYCOPY(&yyptr->Stack, Stack, yysize);                                                                          \
        Stack = &yyptr->Stack;                                                                                         \
        yynewbytes = yystacksize * sizeof(*Stack) + YYSTACK_GAP_MAXIMUM;                                               \
        yyptr += yynewbytes / sizeof(*yyptr);                                                                          \
    } while (YYID(0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL 24
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST 247

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS 32
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS 23
/* YYNRULES -- Number of rules.  */
#define YYNRULES 53
/* YYNRULES -- Number of states.  */
#define YYNSTATES 87

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK 2
#define YYMAXUTOK 282

#define YYTRANSLATE(YYX) ((unsigned int)(YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] = {
    0, 2, 2, 2, 2, 2, 2, 2, 2, 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 30, 2, 31, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 28, 2, 29, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2, 2,  2, 2,  2, 2,
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] = {0,   0,   3,   6,   8,   10,  13,  16,  19,  22,  25,  29,  34,  36,
                                      38,  40,  43,  45,  47,  49,  51,  53,  55,  57,  60,  64,  67,  70,
                                      74,  79,  82,  86,  89,  92,  96,  98,  103, 106, 109, 112, 115, 118,
                                      123, 127, 130, 131, 132, 133, 134, 135, 136, 137, 141, 144};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] = {
    33, 0,  -1, 34, 3,  -1, 1,  -1, 35, -1, 34, 35, -1, 7,  43, -1, 6,  42, -1, 5,  39, -1, 14, 41, -1,
    8,  39, 39, -1, 21, 45, 34, 22, -1, 26, -1, 23, -1, 37, -1, 36, 37, -1, 23, -1, 24, -1, 25, -1, 26,
    -1, 27, -1, 53, -1, 38, -1, 10, 39, -1, 19, 39, 39, -1, 15, 40, -1, 16, 41, -1, 12, 48, 39, -1, 12,
    48, 54, 39, -1, 9,  42, -1, 18, 42, 42, -1, 11, 43, -1, 20, 43, -1, 20, 43, 44, -1, 13, -1, 21, 45,
    36, 22, -1, 46, 53, -1, 51, 53, -1, 52, 53, -1, 47, 53, -1, 49, 53, -1, 28, 50, 36, 29, -1, 28, 50,
    29, -1, 46, 23, -1, -1, -1, -1, -1, -1, -1, -1, 28, 36, 29, -1, 28, 29, -1, 30, 37, 31, -1};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] = {0,   119, 119, 120, 123, 124, 126, 127, 128, 129, 130, 131, 132, 133,
                                       135, 136, 138, 139, 140, 141, 142, 143, 144, 146, 147, 148, 149, 150,
                                       151, 152, 153, 154, 155, 156, 157, 158, 160, 162, 164, 166, 168, 172,
                                       173, 175, 178, 180, 182, 184, 186, 188, 190, 192, 193, 195};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] = {"$end",
                                      "error",
                                      "$undefined",
                                      "END_OF_INPUT",
                                      "ERROR",
                                      "SECTIONHEADER",
                                      "RSECTIONHEADER",
                                      "VSECTIONHEADER",
                                      "SECTIONHEADER2",
                                      "RCODEMACRO",
                                      "LATEXMACRO",
                                      "VERBMACRO",
                                      "OPTMACRO",
                                      "ESCAPE",
                                      "LISTSECTION",
                                      "ITEMIZE",
                                      "DESCRIPTION",
                                      "NOITEM",
                                      "RCODEMACRO2",
                                      "LATEXMACRO2",
                                      "VERBMACRO2",
                                      "IFDEF",
                                      "ENDIF",
                                      "TEXT",
                                      "RCODE",
                                      "VERB",
                                      "COMMENT",
                                      "UNKNOWN",
                                      "'{'",
                                      "'}'",
                                      "'['",
                                      "']'",
                                      "$accept",
                                      "RdFile",
                                      "SectionList",
                                      "Section",
                                      "ArgItems",
                                      "Item",
                                      "Markup",
                                      "LatexArg",
                                      "Item0Arg",
                                      "Item2Arg",
                                      "RLikeArg",
                                      "VerbatimArg",
                                      "VerbatimArg2",
                                      "IfDefTarget",
                                      "goLatexLike",
                                      "goRLike",
                                      "goOption",
                                      "goVerbatim",
                                      "goVerbatim2",
                                      "goItem0",
                                      "goItem2",
                                      "Arg",
                                      "Option",
                                      0};
#endif

#ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] = {0,   256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270,
                                         271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 123, 125, 91,  93};
#endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] = {0,  32, 33, 33, 34, 34, 35, 35, 35, 35, 35, 35, 35, 35, 36, 36, 37, 37,
                                    37, 37, 37, 37, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38,
                                    39, 40, 41, 42, 43, 44, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 53, 54};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] = {0, 2, 2, 1, 1, 2, 2, 2, 2, 2, 3, 4, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 2,
                                    3, 4, 2, 3, 2, 2, 3, 1, 4, 2, 2, 2, 2, 2, 4, 3, 2, 0, 0, 0, 0, 0, 0, 0, 3, 2, 3};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] = {
    0,  3,  44, 45, 47, 44, 50, 44, 13, 12, 0,  0,  4,  8,  0,  7,  0,  6,  0,  44, 9,  0,  0,  0, 1,  2,  5,  0,  36,
    39, 40, 10, 38, 0,  43, 45, 44, 47, 46, 34, 49, 50, 45, 44, 47, 44, 16, 17, 18, 19, 20, 52, 0, 14, 22, 21, 11, 29,
    23, 31, 44, 25, 0,  26, 45, 44, 32, 0,  51, 15, 0,  27, 44, 37, 30, 24, 48, 33, 0,  0,  28, 0, 35, 53, 42, 0,  41};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] = {-1, 10, 11, 12, 52, 53, 54, 13, 61, 20, 15, 17,
                                        77, 22, 14, 16, 60, 18, 81, 62, 21, 55, 72};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -59
static const yytype_int16 yypact[] = {
    40,  -59, -59, -59, -59, -59, -59, -59, -59, -59, 11,  65,  -59, -59, -16, -59, -16, -59, -16, -59, -59, -16,
    97,  -8,  -59, -59, -59, 115, -59, -59, -59, -59, -59, 75,  -59, -59, -59, -59, -59, -59, -59, -59, -59, -59,
    -59, -59, -59, -59, -59, -59, -59, -59, 136, -59, -59, -59, -59, -59, -59, -59, -14, -59, -16, -59, -59, -59,
    -10, 219, -59, -59, 219, -59, -59, -59, -59, -59, -59, -59, 199, -12, -59, 157, -59, -59, -59, 178, -59};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] = {-59, -59, -1, -4,  -58, -50, -59, -5,  -59, -19, -25, -31,
                                      -59, -21, -3, -59, -59, -59, -59, -59, -59, -13, -59};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] = {
    19, 28, 69, 29, 23, 30, 59, 26, 32, 78, 57, 24, 27, 66, 31, 34, 70, 64, 76, 83, 79, 33, 63, 85, 67, 0,  0,  0,
    69, 26, 0,  58, 0,  0,  0,  69, 0,  0,  65, 74, 0,  1,  23, 0,  0,  2,  3,  4,  5,  73, 0,  0,  0,  0,  6,  71,
    0,  0,  0,  0,  75, 7,  0,  8,  0,  0,  9,  80, 25, 0,  2,  3,  4,  5,  0,  0,  0,  0,  0,  6,  2,  3,  4,  5,
    0,  0,  7,  0,  8,  6,  0,  9,  0,  0,  0,  0,  7,  56, 8,  0,  0,  9,  2,  3,  4,  5,  0,  0,  0,  0,  0,  6,
    0,  0,  0,  0,  0,  0,  7,  0,  8,  0,  0,  9,  35, 36, 37, 38, 39, 0,  40, 41, 0,  42, 43, 44, 45, 0,  46, 47,
    48, 49, 50, 27, 51, 35, 36, 37, 38, 39, 0,  40, 41, 0,  42, 43, 44, 45, 0,  46, 47, 48, 49, 50, 27, 68, 35, 36,
    37, 38, 39, 0,  40, 41, 0,  42, 43, 44, 45, 0,  46, 47, 48, 49, 50, 27, 84, 35, 36, 37, 38, 39, 0,  40, 41, 0,
    42, 43, 44, 45, 0,  46, 47, 48, 49, 50, 27, 86, 35, 36, 37, 38, 39, 0,  40, 41, 0,  42, 43, 44, 45, 82, 46, 47,
    48, 49, 50, 27, 35, 36, 37, 38, 39, 0,  40, 41, 0,  42, 43, 44, 45, 0,  46, 47, 48, 49, 50, 27};

static const yytype_int8 yycheck[] = {
    5,  14, 52, 16, 7,  18, 37, 11, 21, 67, 35, 0,  28, 44, 19, 23, 30, 42, 28, 31, 70, 22, 41, 81, 45, -1, -1, -1,
    78, 33, -1, 36, -1, -1, -1, 85, -1, -1, 43, 64, -1, 1,  45, -1, -1, 5,  6,  7,  8,  62, -1, -1, -1, -1, 14, 60,
    -1, -1, -1, -1, 65, 21, -1, 23, -1, -1, 26, 72, 3,  -1, 5,  6,  7,  8,  -1, -1, -1, -1, -1, 14, 5,  6,  7,  8,
    -1, -1, 21, -1, 23, 14, -1, 26, -1, -1, -1, -1, 21, 22, 23, -1, -1, 26, 5,  6,  7,  8,  -1, -1, -1, -1, -1, 14,
    -1, -1, -1, -1, -1, -1, 21, -1, 23, -1, -1, 26, 9,  10, 11, 12, 13, -1, 15, 16, -1, 18, 19, 20, 21, -1, 23, 24,
    25, 26, 27, 28, 29, 9,  10, 11, 12, 13, -1, 15, 16, -1, 18, 19, 20, 21, -1, 23, 24, 25, 26, 27, 28, 29, 9,  10,
    11, 12, 13, -1, 15, 16, -1, 18, 19, 20, 21, -1, 23, 24, 25, 26, 27, 28, 29, 9,  10, 11, 12, 13, -1, 15, 16, -1,
    18, 19, 20, 21, -1, 23, 24, 25, 26, 27, 28, 29, 9,  10, 11, 12, 13, -1, 15, 16, -1, 18, 19, 20, 21, 22, 23, 24,
    25, 26, 27, 28, 9,  10, 11, 12, 13, -1, 15, 16, -1, 18, 19, 20, 21, -1, 23, 24, 25, 26, 27, 28};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] = {
    0,  1,  5,  6,  7,  8,  14, 21, 23, 26, 33, 34, 35, 39, 46, 42, 47, 43, 49, 39, 41, 52, 45, 46, 0,  3,  35, 28, 53,
    53, 53, 39, 53, 34, 23, 9,  10, 11, 12, 13, 15, 16, 18, 19, 20, 21, 23, 24, 25, 26, 27, 29, 36, 37, 38, 53, 22, 42,
    39, 43, 48, 40, 51, 41, 42, 39, 43, 45, 29, 37, 30, 39, 54, 53, 42, 39, 28, 44, 36, 37, 39, 50, 22, 31, 29, 36, 29};

#define yyerrok (yyerrstatus = 0)
#define yyclearin (yychar = YYEMPTY)
#define YYEMPTY (-2)
#define YYEOF 0

#define YYACCEPT goto yyacceptlab
#define YYABORT goto yyabortlab
#define YYERROR goto yyerrorlab

/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL goto yyerrlab

#define YYRECOVERING() (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                                                                         \
    do                                                                                                                 \
        if (yychar == YYEMPTY && yylen == 1)                                                                           \
        {                                                                                                              \
            yychar = (Token);                                                                                          \
            yylval = (Value);                                                                                          \
            yytoken = YYTRANSLATE(yychar);                                                                             \
            YYPOPSTACK(1);                                                                                             \
            goto yybackup;                                                                                             \
        }                                                                                                              \
        else                                                                                                           \
        {                                                                                                              \
            yyerror(YY_("syntax error: cannot back up"));                                                              \
            YYERROR;                                                                                                   \
        }                                                                                                              \
    while (YYID(0))

#define YYTERROR 1
#define YYERRCODE 256

/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
#define YYLLOC_DEFAULT(Current, Rhs, N)                                                                                \
    do                                                                                                                 \
        if (YYID(N))                                                                                                   \
        {                                                                                                              \
            (Current).first_line = YYRHSLOC(Rhs, 1).first_line;                                                        \
            (Current).first_column = YYRHSLOC(Rhs, 1).first_column;                                                    \
            (Current).last_line = YYRHSLOC(Rhs, N).last_line;                                                          \
            (Current).last_column = YYRHSLOC(Rhs, N).last_column;                                                      \
        }                                                                                                              \
        else                                                                                                           \
        {                                                                                                              \
            (Current).first_line = (Current).last_line = YYRHSLOC(Rhs, 0).last_line;                                   \
            (Current).first_column = (Current).last_column = YYRHSLOC(Rhs, 0).last_column;                             \
        }                                                                                                              \
    while (YYID(0))
#endif

/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
#if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
#define YY_LOCATION_PRINT(File, Loc)                                                                                   \
    fprintf(File, "%d.%d-%d.%d", (Loc).first_line, (Loc).first_column, (Loc).last_line, (Loc).last_column)
#else
#define YY_LOCATION_PRINT(File, Loc) ((void)0)
#endif
#endif

/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
#define YYLEX yylex(YYLEX_PARAM)
#else
#define YYLEX yylex()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

#ifndef YYFPRINTF
#include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#define YYFPRINTF fprintf
#endif

#define YYDPRINTF(Args)                                                                                                \
    do                                                                                                                 \
    {                                                                                                                  \
        if (yydebug)                                                                                                   \
            YYFPRINTF Args;                                                                                            \
    } while (YYID(0))

#define YY_SYMBOL_PRINT(Title, Type, Value, Location)                                                                  \
    do                                                                                                                 \
    {                                                                                                                  \
        if (yydebug)                                                                                                   \
        {                                                                                                              \
            YYFPRINTF(stderr, "%s ", Title);                                                                           \
            yy_symbol_print(stderr, Type, Value, Location);                                                            \
            YYFPRINTF(stderr, "\n");                                                                                   \
        }                                                                                                              \
    } while (YYID(0))

/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
static void yy_symbol_value_print(FILE *yyoutput, int yytype, YYSTYPE const *const yyvaluep,
                                  YYLTYPE const *const yylocationp)
#else
static void yy_symbol_value_print(yyoutput, yytype, yyvaluep, yylocationp) FILE *yyoutput;
int yytype;
YYSTYPE const *const yyvaluep;
YYLTYPE const *const yylocationp;
#endif
{
    if (!yyvaluep)
        return;
    YYUSE(yylocationp);
#ifdef YYPRINT
    if (yytype < YYNTOKENS)
        YYPRINT(yyoutput, yytoknum[yytype], *yyvaluep);
#else
    YYUSE(yyoutput);
#endif
    switch (yytype)
    {
    default:
        break;
    }
}

/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
static void yy_symbol_print(FILE *yyoutput, int yytype, YYSTYPE const *const yyvaluep, YYLTYPE const *const yylocationp)
#else
static void yy_symbol_print(yyoutput, yytype, yyvaluep, yylocationp) FILE *yyoutput;
int yytype;
YYSTYPE const *const yyvaluep;
YYLTYPE const *const yylocationp;
#endif
{
    if (yytype < YYNTOKENS)
        YYFPRINTF(yyoutput, "token %s (", yytname[yytype]);
    else
        YYFPRINTF(yyoutput, "nterm %s (", yytname[yytype]);

    YY_LOCATION_PRINT(yyoutput, *yylocationp);
    YYFPRINTF(yyoutput, ": ");
    yy_symbol_value_print(yyoutput, yytype, yyvaluep, yylocationp);
    YYFPRINTF(yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
static void yy_stack_print(yytype_int16 *bottom, yytype_int16 *top)
#else
static void yy_stack_print(bottom, top) yytype_int16 *bottom;
yytype_int16 *top;
#endif
{
    YYFPRINTF(stderr, "Stack now");
    for (; bottom <= top; ++bottom)
        YYFPRINTF(stderr, " %d", *bottom);
    YYFPRINTF(stderr, "\n");
}

#define YY_STACK_PRINT(Bottom, Top)                                                                                    \
    do                                                                                                                 \
    {                                                                                                                  \
        if (yydebug)                                                                                                   \
            yy_stack_print((Bottom), (Top));                                                                           \
    } while (YYID(0))

/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
static void yy_reduce_print(YYSTYPE *yyvsp, YYLTYPE *yylsp, int yyrule)
#else
static void yy_reduce_print(yyvsp, yylsp, yyrule) YYSTYPE *yyvsp;
YYLTYPE *yylsp;
int yyrule;
#endif
{
    int yynrhs = yyr2[yyrule];
    int yyi;
    unsigned long int yylno = yyrline[yyrule];
    YYFPRINTF(stderr, "Reducing stack by rule %d (line %lu):\n", yyrule - 1, yylno);
    /* The symbols being reduced.  */
    for (yyi = 0; yyi < yynrhs; yyi++)
    {
        fprintf(stderr, "   $%d = ", yyi + 1);
        yy_symbol_print(stderr, yyrhs[yyprhs[yyrule] + yyi], &(yyvsp[(yyi + 1) - (yynrhs)]),
                        &(yylsp[(yyi + 1) - (yynrhs)]));
        fprintf(stderr, "\n");
    }
}

#define YY_REDUCE_PRINT(Rule)                                                                                          \
    do                                                                                                                 \
    {                                                                                                                  \
        if (yydebug)                                                                                                   \
            yy_reduce_print(yyvsp, yylsp, Rule);                                                                       \
    } while (YYID(0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
#define YYDPRINTF(Args)
#define YY_SYMBOL_PRINT(Title, Type, Value, Location)
#define YY_STACK_PRINT(Bottom, Top)
#define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */

/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
#define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
#define YYMAXDEPTH 10000
#endif

#if YYERROR_VERBOSE

#ifndef yystrlen
#if defined __GLIBC__ && defined _STRING_H
#define yystrlen strlen
#else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T yystrlen(const char *yystr)
#else
static YYSIZE_T yystrlen(yystr) const char *yystr;
#endif
{
    YYSIZE_T yylen;
    for (yylen = 0; yystr[yylen]; yylen++)
        continue;
    return yylen;
}
#endif
#endif

#ifndef yystpcpy
#if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#define yystpcpy stpcpy
#else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
static char *yystpcpy(char *yydest, const char *yysrc)
#else
static char *yystpcpy(yydest, yysrc) char *yydest;
const char *yysrc;
#endif
{
    char *yyd = yydest;
    const char *yys = yysrc;

    while ((*yyd++ = *yys++) != '\0')
        continue;

    return yyd - 1;
}
#endif
#endif

#ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T yytnamerr(char *yyres, const char *yystr)
{
    if (*yystr == '"')
    {
        YYSIZE_T yyn = 0;
        char const *yyp = yystr;

        for (;;)
            switch (*++yyp)
            {
            case '\'':
            case ',':
                goto do_not_strip_quotes;

            case '\\':
                if (*++yyp != '\\')
                    goto do_not_strip_quotes;
                /* Fall through.  */
            default:
                if (yyres)
                    yyres[yyn] = *yyp;
                yyn++;
                break;

            case '"':
                if (yyres)
                    yyres[yyn] = '\0';
                return yyn;
            }
    do_not_strip_quotes:;
    }

    if (!yyres)
        return yystrlen(yystr);

    return yystpcpy(yyres, yystr) - yyres;
}
#endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T yysyntax_error(char *yyresult, int yystate, int yychar)
{
    int yyn = yypact[yystate];

    if (!(YYPACT_NINF < yyn && yyn <= YYLAST))
        return 0;
    else
    {
        int yytype = YYTRANSLATE(yychar);
        YYSIZE_T yysize0 = yytnamerr(0, yytname[yytype]);
        YYSIZE_T yysize = yysize0;
        YYSIZE_T yysize1;
        int yysize_overflow = 0;
        enum
        {
            YYERROR_VERBOSE_ARGS_MAXIMUM = 5
        };
        char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
        int yyx;

#if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
#endif
        char *yyfmt;
        char const *yyf;
        static char const yyunexpected[] = "syntax error, unexpected %s";
        static char const yyexpecting[] = ", expecting %s";
        static char const yyor[] = " or %s";
        char yyformat[sizeof yyunexpected + sizeof yyexpecting - 1 +
                      ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2) * (sizeof yyor - 1))];
        char const *yyprefix = yyexpecting;

        /* Start YYX at -YYN if negative to avoid negative indexes in
       YYCHECK.  */
        int yyxbegin = yyn < 0 ? -yyn : 0;

        /* Stay within bounds of both yycheck and yytname.  */
        int yychecklim = YYLAST - yyn + 1;
        int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
        int yycount = 1;

        yyarg[0] = yytname[yytype];
        yyfmt = yystpcpy(yyformat, yyunexpected);

        for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
            {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                {
                    yycount = 1;
                    yysize = yysize0;
                    yyformat[sizeof yyunexpected - 1] = '\0';
                    break;
                }
                yyarg[yycount++] = yytname[yyx];
                yysize1 = yysize + yytnamerr(0, yytname[yyx]);
                yysize_overflow |= (yysize1 < yysize);
                yysize = yysize1;
                yyfmt = yystpcpy(yyfmt, yyprefix);
                yyprefix = yyor;
            }

        yyf = YY_(yyformat);
        yysize1 = yysize + yystrlen(yyf);
        yysize_overflow |= (yysize1 < yysize);
        yysize = yysize1;

        if (yysize_overflow)
            return YYSIZE_MAXIMUM;

        if (yyresult)
        {
            /* Avoid sprintf, as that infringes on the user's name space.
               Don't have undefined behavior even if the translation
               produced a string with the wrong number of "%s"s.  */
            char *yyp = yyresult;
            int yyi = 0;
            while ((*yyp = *yyf) != '\0')
            {
                if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
                {
                    yyp += yytnamerr(yyp, yyarg[yyi++]);
                    yyf += 2;
                }
                else
                {
                    yyp++;
                    yyf++;
                }
            }
        }
        return yysize;
    }
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
static void yydestruct(const char *yymsg, int yytype, YYSTYPE *yyvaluep, YYLTYPE *yylocationp)
#else
static void yydestruct(yymsg, yytype, yyvaluep, yylocationp) const char *yymsg;
int yytype;
YYSTYPE *yyvaluep;
YYLTYPE *yylocationp;
#endif
{
    YYUSE(yyvaluep);
    YYUSE(yylocationp);

    if (!yymsg)
        yymsg = "Deleting";
    YY_SYMBOL_PRINT(yymsg, yytype, yyvaluep, yylocationp);

    switch (yytype)
    {

    default:
        break;
    }
}

/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse(void *YYPARSE_PARAM);
#else
int yyparse();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse(void);
#else
int yyparse();
#endif
#endif /* ! YYPARSE_PARAM */

/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;
/* Location data for the look-ahead symbol.  */
YYLTYPE yylloc;

/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
int yyparse(void *YYPARSE_PARAM)
#else
int yyparse(YYPARSE_PARAM) void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ || defined __cplusplus || defined _MSC_VER)
int yyparse(void)
#else
int yyparse()

#endif
#endif
{

    int yystate;
    int yyn;
    int yyresult;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;
    /* Look-ahead token as an internal (translated) token number.  */
    int yytoken = 0;
#if YYERROR_VERBOSE
    /* Buffer for error messages, and its allocated size.  */
    char yymsgbuf[128];
    char *yymsg = yymsgbuf;
    YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

    /* Three stacks and their tools:
       `yyss': related to states,
       `yyvs': related to semantic values,
       `yyls': related to locations.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss = yyssa;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp;

    /* The location stack.  */
    YYLTYPE yylsa[YYINITDEPTH];
    YYLTYPE *yyls = yylsa;
    YYLTYPE *yylsp;
    /* The locations where the error started and ended.  */
    YYLTYPE yyerror_range[2];

#define YYPOPSTACK(N) (yyvsp -= (N), yyssp -= (N), yylsp -= (N))

    YYSIZE_T yystacksize = YYINITDEPTH;

    /* The variables used to return semantic value and location from the
       action routines.  */
    YYSTYPE yyval;
    YYLTYPE yyloc;

    /* The number of symbols on the RHS of the reduced rule.
       Keep to zero when no symbol should be popped.  */
    int yylen = 0;

    YYDPRINTF((stderr, "Starting parse\n"));

    yystate = 0;
    yyerrstatus = 0;
    yynerrs = 0;
    yychar = YYEMPTY; /* Cause a token to be read.  */

    /* Initialize stack pointers.
       Waste one element of value and location stack
       so that they stay on the same level as the state stack.
       The wasted elements are never initialized.  */

    yyssp = yyss;
    yyvsp = yyvs;
    yylsp = yyls;
#if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
    /* Initialize the default location before parsing starts.  */
    yylloc.first_line = yylloc.last_line = 1;
    yylloc.first_column = yylloc.last_column = 0;
#endif

    goto yysetstate;

    /*------------------------------------------------------------.
    | yynewstate -- Push a new state, which is found in yystate.  |
    `------------------------------------------------------------*/
yynewstate:
    /* In all cases, when you get here, the value and location stacks
       have just been pushed.  So pushing a state here evens the stacks.  */
    yyssp++;

yysetstate:
    *yyssp = yystate;

    if (yyss + yystacksize - 1 <= yyssp)
    {
        /* Get the current used size of the three stacks, in elements.  */
        YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
        {
            /* Give user a chance to reallocate the stack.  Use copies of
               these so that the &'s don't force the real ones into
               memory.  */
            YYSTYPE *yyvs1 = yyvs;
            yytype_int16 *yyss1 = yyss;
            YYLTYPE *yyls1 = yyls;

            /* Each stack pointer address is followed by the size of the
               data in use in that stack, in bytes.  This used to be a
               conditional around just the two extra args, but that might
               be undefined if yyoverflow is a macro.  */
            yyoverflow(YY_("memory exhausted"), &yyss1, yysize * sizeof(*yyssp), &yyvs1, yysize * sizeof(*yyvsp),
                       &yyls1, yysize * sizeof(*yylsp), &yystacksize);
            yyls = yyls1;
            yyss = yyss1;
            yyvs = yyvs1;
        }
#else /* no yyoverflow */
#ifndef YYSTACK_RELOCATE
        goto yyexhaustedlab;
#else
        /* Extend the stack our own way.  */
        if (YYMAXDEPTH <= yystacksize)
            goto yyexhaustedlab;
        yystacksize *= 2;
        if (YYMAXDEPTH < yystacksize)
            yystacksize = YYMAXDEPTH;

        {
            yytype_int16 *yyss1 = yyss;
            union yyalloc *yyptr = (union yyalloc *)YYSTACK_ALLOC(YYSTACK_BYTES(yystacksize));
            if (!yyptr)
                goto yyexhaustedlab;
            YYSTACK_RELOCATE(yyss);
            YYSTACK_RELOCATE(yyvs);
            YYSTACK_RELOCATE(yyls);
#undef YYSTACK_RELOCATE
            if (yyss1 != yyssa)
                YYSTACK_FREE(yyss1);
        }
#endif
#endif /* no yyoverflow */

        yyssp = yyss + yysize - 1;
        yyvsp = yyvs + yysize - 1;
        yylsp = yyls + yysize - 1;

        YYDPRINTF((stderr, "Stack size increased to %lu\n", (unsigned long int)yystacksize));

        if (yyss + yystacksize - 1 <= yyssp)
            YYABORT;
    }

    YYDPRINTF((stderr, "Entering state %d\n", yystate));

    goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

    /* Do appropriate processing given the current state.  Read a
       look-ahead token if we need one and don't already have one.  */

    /* First try to decide what to do without reference to look-ahead token.  */
    yyn = yypact[yystate];
    if (yyn == YYPACT_NINF)
        goto yydefault;

    /* Not known => get a look-ahead token if don't already have one.  */

    /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
    if (yychar == YYEMPTY)
    {
        YYDPRINTF((stderr, "Reading a token: "));
        yychar = YYLEX;
    }

    if (yychar <= YYEOF)
    {
        yychar = yytoken = YYEOF;
        YYDPRINTF((stderr, "Now at end of input.\n"));
    }
    else
    {
        yytoken = YYTRANSLATE(yychar);
        YY_SYMBOL_PRINT("Next token is", yytoken, &yylval, &yylloc);
    }

    /* If the proper action on seeing token YYTOKEN is to reduce or to
       detect an error, take that action.  */
    yyn += yytoken;
    if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
        goto yydefault;
    yyn = yytable[yyn];
    if (yyn <= 0)
    {
        if (yyn == 0 || yyn == YYTABLE_NINF)
            goto yyerrlab;
        yyn = -yyn;
        goto yyreduce;
    }

    if (yyn == YYFINAL)
        YYACCEPT;

    /* Count tokens shifted since error; after three, turn off error
       status.  */
    if (yyerrstatus)
        yyerrstatus--;

    /* Shift the look-ahead token.  */
    YY_SYMBOL_PRINT("Shifting", yytoken, &yylval, &yylloc);

    /* Discard the shifted token unless it is eof.  */
    if (yychar != YYEOF)
        yychar = YYEMPTY;

    yystate = yyn;
    *++yyvsp = yylval;
    *++yylsp = yylloc;
    goto yynewstate;

/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
    yyn = yydefact[yystate];
    if (yyn == 0)
        goto yyerrlab;
    goto yyreduce;

/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
    /* yyn is the number of a rule to reduce with.  */
    yylen = yyr2[yyn];

    /* If YYLEN is nonzero, implement the default value of the action:
       `$$ = $1'.

       Otherwise, the following line sets YYVAL to garbage.
       This behavior is undocumented and Bison
       users should not rely upon it.  Assigning to YYVAL
       unconditionally makes the parser a bit smaller, and it avoids a
       GCC warning that YYVAL may be used uninitialized.  */
    yyval = yyvsp[1 - yylen];

    /* Default location.  */
    YYLLOC_DEFAULT(yyloc, (yylsp - yylen), yylen);
    YY_REDUCE_PRINT(yyn);
    switch (yyn)
    {
    case 2:

    {
        xxsavevalue((yyvsp[(1) - (2)]), &(yyloc));
        return 0;
    }
    break;

    case 3:

    {
        PROTECT(Value = R_NilValue);
        YYABORT;
    }
    break;

    case 4:

    {
        (yyval) = xxnewlist((yyvsp[(1) - (1)]));
    }
    break;

    case 5:

    {
        (yyval) = xxlist((yyvsp[(1) - (2)]), (yyvsp[(2) - (2)]));
    }
    break;

    case 6:

    {
        (yyval) = xxmarkup((yyvsp[(1) - (2)]), (yyvsp[(2) - (2)]), &(yyloc));
    }
    break;

    case 7:

    {
        (yyval) = xxmarkup((yyvsp[(1) - (2)]), (yyvsp[(2) - (2)]), &(yyloc));
    }
    break;

    case 8:

    {
        (yyval) = xxmarkup((yyvsp[(1) - (2)]), (yyvsp[(2) - (2)]), &(yyloc));
    }
    break;

    case 9:

    {
        (yyval) = xxmarkup((yyvsp[(1) - (2)]), (yyvsp[(2) - (2)]), &(yyloc));
    }
    break;

    case 10:

    {
        (yyval) = xxmarkup2((yyvsp[(1) - (3)]), (yyvsp[(2) - (3)]), (yyvsp[(3) - (3)]), &(yyloc));
    }
    break;

    case 11:

    {
        (yyval) = xxmarkup2((yyvsp[(1) - (4)]), (yyvsp[(2) - (4)]), (yyvsp[(3) - (4)]), &(yyloc));
        UNPROTECT_PTR((yyvsp[(4) - (4)]));
    }
    break;

    case 12:

    {
        (yyval) = xxtag((yyvsp[(1) - (1)]), COMMENT, &(yyloc));
    }
    break;

    case 13:

    {
        (yyval) = xxtag((yyvsp[(1) - (1)]), TEXT, &(yyloc));
    }
    break;

    case 14:

    {
        (yyval) = xxnewlist((yyvsp[(1) - (1)]));
    }
    break;

    case 15:

    {
        (yyval) = xxlist((yyvsp[(1) - (2)]), (yyvsp[(2) - (2)]));
    }
    break;

    case 16:

    {
        (yyval) = xxtag((yyvsp[(1) - (1)]), TEXT, &(yyloc));
    }
    break;

    case 17:

    {
        (yyval) = xxtag((yyvsp[(1) - (1)]), RCODE, &(yyloc));
    }
    break;

    case 18:

    {
        (yyval) = xxtag((yyvsp[(1) - (1)]), VERB, &(yyloc));
    }
    break;

    case 19:

    {
        (yyval) = xxtag((yyvsp[(1) - (1)]), COMMENT, &(yyloc));
    }
    break;

    case 20:

    {
        (yyval) = xxtag((yyvsp[(1) - (1)]), UNKNOWN, &(yyloc));
    }
    break;

    case 21:

    {
        (yyval) = xxmarkup(R_NilValue, (yyvsp[(1) - (1)]), &(yyloc));
    }
    break;

    case 22:

    {
        (yyval) = (yyvsp[(1) - (1)]);
    }
    break;

    case 23:

    {
        (yyval) = xxmarkup((yyvsp[(1) - (2)]), (yyvsp[(2) - (2)]), &(yyloc));
    }
    break;

    case 24:

    {
        (yyval) = xxmarkup2((yyvsp[(1) - (3)]), (yyvsp[(2) - (3)]), (yyvsp[(3) - (3)]), &(yyloc));
    }
    break;

    case 25:

    {
        (yyval) = xxmarkup((yyvsp[(1) - (2)]), (yyvsp[(2) - (2)]), &(yyloc));
    }
    break;

    case 26:

    {
        (yyval) = xxmarkup((yyvsp[(1) - (2)]), (yyvsp[(2) - (2)]), &(yyloc));
    }
    break;

    case 27:

    {
        (yyval) = xxmarkup((yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]), &(yyloc));
        xxpopMode((yyvsp[(2) - (3)]));
    }
    break;

    case 28:

    {
        (yyval) = xxOptionmarkup((yyvsp[(1) - (4)]), (yyvsp[(3) - (4)]), (yyvsp[(4) - (4)]), &(yyloc));
        xxpopMode((yyvsp[(2) - (4)]));
    }
    break;

    case 29:

    {
        (yyval) = xxmarkup((yyvsp[(1) - (2)]), (yyvsp[(2) - (2)]), &(yyloc));
    }
    break;

    case 30:

    {
        (yyval) = xxmarkup2((yyvsp[(1) - (3)]), (yyvsp[(2) - (3)]), (yyvsp[(2) - (3)]), &(yyloc));
    }
    break;

    case 31:

    {
        (yyval) = xxmarkup((yyvsp[(1) - (2)]), (yyvsp[(2) - (2)]), &(yyloc));
    }
    break;

    case 32:

    {
        (yyval) = xxmarkup((yyvsp[(1) - (2)]), (yyvsp[(2) - (2)]), &(yyloc));
    }
    break;

    case 33:

    {
        (yyval) = xxmarkup2((yyvsp[(1) - (3)]), (yyvsp[(2) - (3)]), (yyvsp[(3) - (3)]), &(yyloc));
    }
    break;

    case 34:

    {
        (yyval) = xxmarkup((yyvsp[(1) - (1)]), R_NilValue, &(yyloc));
    }
    break;

    case 35:

    {
        (yyval) = xxmarkup2((yyvsp[(1) - (4)]), (yyvsp[(2) - (4)]), (yyvsp[(3) - (4)]), &(yyloc));
        UNPROTECT_PTR((yyvsp[(4) - (4)]));
    }
    break;

    case 36:

    {
        xxpopMode((yyvsp[(1) - (2)]));
        (yyval) = (yyvsp[(2) - (2)]);
    }
    break;

    case 37:

    {
        xxpopMode((yyvsp[(1) - (2)]));
        (yyval) = (yyvsp[(2) - (2)]);
    }
    break;

    case 38:

    {
        xxpopMode((yyvsp[(1) - (2)]));
        (yyval) = (yyvsp[(2) - (2)]);
    }
    break;

    case 39:

    {
        xxpopMode((yyvsp[(1) - (2)]));
        (yyval) = (yyvsp[(2) - (2)]);
    }
    break;

    case 40:

    {
        xxpopMode((yyvsp[(1) - (2)]));
        (yyval) = (yyvsp[(2) - (2)]);
    }
    break;

    case 41:

    {
        xxpopMode((yyvsp[(2) - (4)]));
        (yyval) = (yyvsp[(3) - (4)]);
    }
    break;

    case 42:

    {
        xxpopMode((yyvsp[(2) - (3)]));
        (yyval) = xxnewlist(NULL);
    }
    break;

    case 43:

    {
        xxpopMode((yyvsp[(1) - (2)]));
        (yyval) = xxnewlist((yyvsp[(2) - (2)]));
    }
    break;

    case 44:

    {
        (yyval) = xxpushMode(LATEXLIKE, UNKNOWN);
    }
    break;

    case 45:

    {
        (yyval) = xxpushMode(RLIKE, UNKNOWN);
    }
    break;

    case 46:

    {
        (yyval) = xxpushMode(INOPTION, UNKNOWN);
    }
    break;

    case 47:

    {
        (yyval) = xxpushMode(VERBATIM, UNKNOWN);
    }
    break;

    case 48:

    {
        xxbraceDepth--;
        (yyval) = xxpushMode(VERBATIM, UNKNOWN);
        xxbraceDepth++;
    }
    break;

    case 49:

    {
        (yyval) = xxpushMode(LATEXLIKE, ESCAPE);
    }
    break;

    case 50:

    {
        (yyval) = xxpushMode(LATEXLIKE, LATEXMACRO2);
    }
    break;

    case 51:

    {
        (yyval) = (yyvsp[(2) - (3)]);
    }
    break;

    case 52:

    {
        (yyval) = xxnewlist(NULL);
    }
    break;

    case 53:

    {
        (yyval) = (yyvsp[(2) - (3)]);
    }
    break;

        /* Line 1267 of yacc.c.  */

    default:
        break;
    }
    YY_SYMBOL_PRINT("-> $$ =", yyr1[yyn], &yyval, &yyloc);

    YYPOPSTACK(yylen);
    yylen = 0;
    YY_STACK_PRINT(yyss, yyssp);

    *++yyvsp = yyval;
    *++yylsp = yyloc;

    /* Now `shift' the result of the reduction.  Determine what state
       that goes to, based on the state we popped back to and the rule
       number reduced by.  */

    yyn = yyr1[yyn];

    yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
    if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
        yystate = yytable[yystate];
    else
        yystate = yydefgoto[yyn - YYNTOKENS];

    goto yynewstate;

/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
    /* If not already recovering from an error, report this error.  */
    if (!yyerrstatus)
    {
        ++yynerrs;
#if !YYERROR_VERBOSE
        yyerror(YY_("syntax error"));
#else
        {
            YYSIZE_T yysize = yysyntax_error(0, yystate, yychar);
            if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
            {
                YYSIZE_T yyalloc = 2 * yysize;
                if (!(yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
                    yyalloc = YYSTACK_ALLOC_MAXIMUM;
                if (yymsg != yymsgbuf)
                    YYSTACK_FREE(yymsg);
                yymsg = (char *)YYSTACK_ALLOC(yyalloc);
                if (yymsg)
                    yymsg_alloc = yyalloc;
                else
                {
                    yymsg = yymsgbuf;
                    yymsg_alloc = sizeof yymsgbuf;
                }
            }

            if (0 < yysize && yysize <= yymsg_alloc)
            {
                (void)yysyntax_error(yymsg, yystate, yychar);
                yyerror(yymsg);
            }
            else
            {
                yyerror(YY_("syntax error"));
                if (yysize != 0)
                    goto yyexhaustedlab;
            }
        }
#endif
    }

    yyerror_range[0] = yylloc;

    if (yyerrstatus == 3)
    {
        /* If just tried and failed to reuse look-ahead token after an
       error, discard it.  */

        if (yychar <= YYEOF)
        {
            /* Return failure if at end of input.  */
            if (yychar == YYEOF)
                YYABORT;
        }
        else
        {
            yydestruct("Error: discarding", yytoken, &yylval, &yylloc);
            yychar = YYEMPTY;
        }
    }

    /* Else will try to reuse look-ahead token after shifting the error
       token.  */
    goto yyerrlab1;

/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

    /* Pacify compilers like GCC when the user code never invokes
       YYERROR and the label yyerrorlab therefore never appears in user
       code.  */
    if (/*CONSTCOND*/ 0)
        goto yyerrorlab;

    yyerror_range[0] = yylsp[1 - yylen];
    /* Do not reclaim the symbols of the rule which action triggered
       this YYERROR.  */
    YYPOPSTACK(yylen);
    yylen = 0;
    YY_STACK_PRINT(yyss, yyssp);
    yystate = *yyssp;
    goto yyerrlab1;

/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
    yyerrstatus = 3; /* Each real token shifted decrements this.  */

    for (;;)
    {
        yyn = yypact[yystate];
        if (yyn != YYPACT_NINF)
        {
            yyn += YYTERROR;
            if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
                yyn = yytable[yyn];
                if (0 < yyn)
                    break;
            }
        }

        /* Pop the current state because it cannot handle the error token.  */
        if (yyssp == yyss)
            YYABORT;

        yyerror_range[0] = *yylsp;
        yydestruct("Error: popping", yystos[yystate], yyvsp, yylsp);
        YYPOPSTACK(1);
        yystate = *yyssp;
        YY_STACK_PRINT(yyss, yyssp);
    }

    if (yyn == YYFINAL)
        YYACCEPT;

    *++yyvsp = yylval;

    yyerror_range[1] = yylloc;
    /* Using YYLLOC is tempting, but would change the location of
       the look-ahead.  YYLOC is available though.  */
    YYLLOC_DEFAULT(yyloc, (yyerror_range - 1), 2);
    *++yylsp = yyloc;

    /* Shift the error token.  */
    YY_SYMBOL_PRINT("Shifting", yystos[yyn], yyvsp, yylsp);

    yystate = yyn;
    goto yynewstate;

/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
    yyresult = 0;
    goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
    yyresult = 1;
    goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
    yyerror(YY_("memory exhausted"));
    yyresult = 2;
    /* Fall through.  */
#endif

yyreturn:
    if (yychar != YYEOF && yychar != YYEMPTY)
        yydestruct("Cleanup: discarding lookahead", yytoken, &yylval, &yylloc);
    /* Do not reclaim the symbols of the rule which action triggered
       this YYABORT or YYACCEPT.  */
    YYPOPSTACK(yylen);
    YY_STACK_PRINT(yyss, yyssp);
    while (yyssp != yyss)
    {
        yydestruct("Cleanup: popping", yystos[*yyssp], yyvsp, yylsp);
        YYPOPSTACK(1);
    }
#ifndef yyoverflow
    if (yyss != yyssa)
        YYSTACK_FREE(yyss);
#endif
#if YYERROR_VERBOSE
    if (yymsg != yymsgbuf)
        YYSTACK_FREE(yymsg);
#endif
    /* Make sure YYID is used.  */
    return YYID(yyresult);
}

static SEXP xxpushMode(int newmode, int newitem)
{
    SEXP ans;
    PROTECT(ans = allocVector(INTSXP, 6));

    INTEGER(ans)[0] = xxmode;       /* Lexer mode */
    INTEGER(ans)[1] = xxitemType;   /* What is \item? */
    INTEGER(ans)[2] = xxbraceDepth; /* Brace depth used in RCODE and VERBATIM */
    INTEGER(ans)[3] = xxinRString;  /* Quote char that started a string */
    INTEGER(ans)[4] = xxQuoteLine;  /* Where the quote was */
    INTEGER(ans)[5] = xxQuoteCol;   /*           "         */

#if DEBUGMODE
    Rprintf("xxpushMode(%d, %s) pushes %d, %s, %d\n", newmode, yytname[YYTRANSLATE(newitem)], xxmode,
            yytname[YYTRANSLATE(xxitemType)], xxbraceDepth);
#endif
    xxmode = newmode;
    xxitemType = newitem;
    xxbraceDepth = 0;
    xxinRString = 0;

    return ans;
}

static void xxpopMode(SEXP oldmode)
{
#if DEBUGVALS
    Rprintf("xxpopMode(%d, %s, %d) replaces %d, %s, %d\n", INTEGER(oldmode)[0],
            yytname[YYTRANSLATE(INTEGER(oldmode)[1])], INTEGER(oldmode)[2], xxmode, yytname[YYTRANSLATE(xxitemType)],
            xxbraceDepth);
#endif
    xxmode = INTEGER(oldmode)[0];
    xxitemType = INTEGER(oldmode)[1];
    xxbraceDepth = INTEGER(oldmode)[2];
    xxinRString = INTEGER(oldmode)[3];
    xxQuoteLine = INTEGER(oldmode)[4];
    xxQuoteCol = INTEGER(oldmode)[5];

    UNPROTECT_PTR(oldmode);
}

static SEXP xxnewlist(SEXP item)
{
    SEXP ans, tmp;
#if DEBUGVALS
    Rprintf("xxnewlist(item=%p)", item);
#endif
    PROTECT(tmp = NewList());
    if (item)
    {
        PROTECT(ans = GrowList(tmp, item));
        UNPROTECT_PTR(tmp);
        UNPROTECT_PTR(item);
    }
    else
        ans = tmp;
#if DEBUGVALS
    Rprintf(" result: %p is length %d\n", ans, length(ans));
#endif
    return ans;
}

static SEXP xxlist(SEXP oldlist, SEXP item)
{
    SEXP ans;
#if DEBUGVALS
    Rprintf("xxlist(oldlist=%p, item=%p)", oldlist, item);
#endif
    PROTECT(ans = GrowList(oldlist, item));
    UNPROTECT_PTR(item);
    UNPROTECT_PTR(oldlist);
#if DEBUGVALS
    Rprintf(" result: %p is length %d\n", ans, length(ans));
#endif
    return ans;
}

static SEXP xxmarkup(SEXP header, SEXP body, YYLTYPE *lloc)
{
    SEXP ans;
#if DEBUGVALS
    Rprintf("xxmarkup(header=%p, body=%p)", header, body);
#endif
    if (isNull(body))
        PROTECT(ans = allocVector(VECSXP, 0));
    else
    {
        PROTECT(ans = PairToVectorList(CDR(body)));
        UNPROTECT_PTR(body);
    }
    if (isNull(header))
        PROTECT(header = mkString("LIST"));

    setAttrib(ans, install("Rd_tag"), header);
    if (SrcFile)
        setAttrib(ans, R_SrcrefSymbol, makeSrcref(lloc, SrcFile));
    UNPROTECT_PTR(header);
#if DEBUGVALS
    Rprintf(" result: %p\n", ans);
#endif
    return ans;
}

static SEXP xxOptionmarkup(SEXP header, SEXP option, SEXP body, YYLTYPE *lloc)
{
    SEXP ans;
#if DEBUGVALS
    Rprintf("xxOptionmarkup(header=%p, option=%p, body=%p)", header, option, body);
#endif
    PROTECT(ans = PairToVectorList(CDR(body)));
    UNPROTECT_PTR(body);
    setAttrib(ans, install("Rd_tag"), header);
    UNPROTECT_PTR(header);
    setAttrib(ans, install("Rd_option"), option);
    UNPROTECT_PTR(option);
    if (SrcFile)
        setAttrib(ans, R_SrcrefSymbol, makeSrcref(lloc, SrcFile));
#if DEBUGVALS
    Rprintf(" result: %p\n", ans);
#endif
    return ans;
}

static SEXP xxmarkup2(SEXP header, SEXP body1, SEXP body2, YYLTYPE *lloc)
{
    SEXP ans;
#if DEBUGVALS
    Rprintf("xxmarkup2(header=%p, body1=%p, body2=%p)", header, body1, body2);
#endif
    PROTECT(ans = allocVector(VECSXP, 2));
    if (!isNull(body1))
    {
        SET_VECTOR_ELT(ans, 0, PairToVectorList(CDR(body1)));
        UNPROTECT_PTR(body1);
    }
    SET_VECTOR_ELT(ans, 1, PairToVectorList(CDR(body2)));
    UNPROTECT_PTR(body2);
    setAttrib(ans, install("Rd_tag"), header);
    UNPROTECT_PTR(header);
    if (SrcFile)
        setAttrib(ans, R_SrcrefSymbol, makeSrcref(lloc, SrcFile));
#if DEBUGVALS
    Rprintf(" result: %p\n", ans);
#endif
    return ans;
}

static void xxsavevalue(SEXP Rd, YYLTYPE *lloc)
{
    PROTECT(Value = PairToVectorList(CDR(Rd)));
    if (!isNull(Value))
    {
        setAttrib(Value, R_ClassSymbol, mkString("Rd"));
        if (SrcFile)
            setAttrib(Value, R_SrcrefSymbol, makeSrcref(lloc, SrcFile));
    }
    UNPROTECT_PTR(Rd);
}

static SEXP xxtag(SEXP item, int type, YYLTYPE *lloc)
{
    setAttrib(item, install("Rd_tag"), mkString(yytname[YYTRANSLATE(type)]));
    if (SrcFile)
        setAttrib(item, R_SrcrefSymbol, makeSrcref(lloc, SrcFile));
    return item;
}

/*----------------------------------------------------------------------------*/

static int (*ptr_getc)(void);

/* Private pushback, since file ungetc only guarantees one byte.
   We need up to one MBCS-worth and one failed #ifdef or one numeric
   garbage markup match */

#define PUSHBACK_BUFSIZE 30

static int pushback[PUSHBACK_BUFSIZE];
static unsigned int npush = 0;

static int xxgetc(void)
{
    int c;

    if (npush)
        c = pushback[--npush];
    else
        c = ptr_getc();
    if (c == EOF)
        return R_EOF;

    R_ParseContextLast = (R_ParseContextLast + 1) % PARSE_CONTEXT_SIZE;
    R_ParseContext[R_ParseContextLast] = c;

    if (c == '\n')
    {
        xxlineno += 1;
        xxlastlinelen = xxcolno;
        xxcolno = 0;
    }
    else
        xxcolno++;

    return c;
}

static int xxungetc(int c)
{
    if (c == '\n')
    {
        xxlineno -= 1;
        xxcolno = xxlastlinelen; /* FIXME:  could we push back more than one line? */
        xxlastlinelen = 0;
    }
    else
        xxcolno--;

    R_ParseContext[R_ParseContextLast] = '\0';
    /* Mac OS X requires us to keep this non-negative */
    R_ParseContextLast = (R_ParseContextLast + PARSE_CONTEXT_SIZE - 1) % PARSE_CONTEXT_SIZE;
    if (npush >= PUSHBACK_BUFSIZE - 2)
        return EOF;
    pushback[npush++] = c;
    return c;
}

static SEXP makeSrcref(YYLTYPE *lloc, SEXP srcfile)
{
    SEXP val;

    PROTECT(val = allocVector(INTSXP, 4));
    INTEGER(val)[0] = lloc->first_line;
    INTEGER(val)[1] = lloc->first_column;
    INTEGER(val)[2] = lloc->last_line;
    INTEGER(val)[3] = lloc->last_column;
    setAttrib(val, R_SrcfileSymbol, srcfile);
    setAttrib(val, R_ClassSymbol, mkString("srcref"));
    UNPROTECT(1);
    return val;
}

static SEXP mkString2(const char *s, int len)
{
    SEXP t;
    cetype_t enc = CE_NATIVE;

    if (known_to_be_latin1)
        enc = CE_LATIN1;
    else if (known_to_be_utf8)
        enc = CE_UTF8;

    PROTECT(t = allocVector(STRSXP, 1));
    SET_STRING_ELT(t, 0, mkCharLenCE(s, len, enc));
    UNPROTECT(1);
    return t;
}

/* Stretchy List Structures : Lists are created and grown using a special */
/* dotted pair.  The CAR of the list points to the last cons-cell in the */
/* list and the CDR points to the first.  The list can be extracted from */
/* the pair by taking its CDR, while the CAR gives fast access to the end */
/* of the list. */

/* Create a stretchy-list dotted pair */

static SEXP NewList(void)
{
    SEXP s = CONS(R_NilValue, R_NilValue);
    SETCAR(s, s);
    return s;
}

/* Add a new element at the end of a stretchy list */

static SEXP GrowList(SEXP l, SEXP s)
{
    SEXP tmp;
    PROTECT(s);
    tmp = CONS(s, R_NilValue);
    UNPROTECT(1);
    SETCDR(CAR(l), tmp);
    SETCAR(l, tmp);
    return l;
}

/*--------------------------------------------------------------------------*/

/*
 *  Parsing Entry Points:
 *
 *  The Following entry points provide Rd parsing facilities.
 *
 *	SEXP R_ParseRd(Rconnection con, ParseStatus *status, SEXP srcfile)
 *
 */

static SEXP ParseRd(ParseStatus *status, SEXP srcfile)
{
    R_ParseContextLast = 0;
    R_ParseContext[0] = '\0';

    xxlineno = 1;
    xxcolno = 0;

    if (!isNull(srcfile))
        SrcFile = srcfile;
    else
        SrcFile = NULL;

    npush = 0;
    xxmode = LATEXLIKE;
    xxitemType = UNKNOWN;
    xxbraceDepth = 0;
    xxinRString = 0;

    Value = R_NilValue;

    if (yyparse())
        *status = PARSE_ERROR;
    else
        *status = PARSE_OK;

#if DEBUGVALS
    Rprintf("ParseRd result: %p\n", Value);
#endif
    UNPROTECT_PTR(Value);
    return Value;
}

#include "Rconnections.h"
static Rconnection con_parse;

/* need to handle incomplete last line */
static int con_getc(void)
{
    int c;
    static int last = -1000;

    c = Rconn_fgetc(con_parse);
    if (c == EOF && last != '\n')
        c = '\n';
    return (last = c);
}

/* used in source.c */
attribute_hidden SEXP R_ParseRd(Rconnection con, ParseStatus *status, SEXP srcfile)
{
    con_parse = con;
    ptr_getc = con_getc;
    return ParseRd(status, srcfile);
}

/*----------------------------------------------------------------------------
 *
 *  The Lexical Analyzer:
 *
 *  Basic lexical analysis is performed by the following
 *  routines.
 *
 *  The function yylex() scans the input, breaking it into
 *  tokens which are then passed to the parser.
 *
 */

/* Special Symbols */
/* Section and R code headers */

struct
{
    char *name;
    int token;
} static keywords[] = {
    /* These sections contain Latex-like text */

    {"\\author", SECTIONHEADER},
    {"\\concept", SECTIONHEADER},
    {"\\description", SECTIONHEADER},
    {"\\details", SECTIONHEADER},
    {"\\docType", SECTIONHEADER},

    {"\\encoding", SECTIONHEADER},
    {"\\format", SECTIONHEADER},
    {"\\keyword", SECTIONHEADER},
    {"\\name", SECTIONHEADER},
    {"\\note", SECTIONHEADER},

    {"\\references", SECTIONHEADER},
    {"\\section", SECTIONHEADER2},
    {"\\seealso", SECTIONHEADER},
    {"\\source", SECTIONHEADER},
    {"\\title", SECTIONHEADER},

    /* These sections contain R-like text */

    {"\\examples", RSECTIONHEADER},
    {"\\usage", RSECTIONHEADER},

    /* This section contains verbatim text */

    {"\\alias", VSECTIONHEADER},
    {"\\synopsis", VSECTIONHEADER},
    {"\\Rdversion", VSECTIONHEADER},

    /* These macros take no arguments.  One character non-alpha escapes get the
       same token value */

    {"\\cr", ESCAPE},
    {"\\dots", ESCAPE},
    {"\\ldots", ESCAPE},
    {"\\R", ESCAPE},
    {"\\tab", ESCAPE},

    /* These macros take one LaTeX-like argument. */

    {"\\acronym", LATEXMACRO},
    {"\\bold", LATEXMACRO},
    {"\\cite", LATEXMACRO},
    {"\\dfn", LATEXMACRO},
    {"\\dQuote", LATEXMACRO},
    {"\\email", LATEXMACRO},

    {"\\emph", LATEXMACRO},
    {"\\file", LATEXMACRO},
    {"\\linkS4class", LATEXMACRO},
    {"\\pkg", LATEXMACRO},
    {"\\sQuote", LATEXMACRO},

    {"\\strong", LATEXMACRO},

    {"\\var", LATEXMACRO},

    /* These are like SECTIONHEADER/LATEXMACRO, but they change the interpretation of \item */

    {"\\arguments", LISTSECTION},
    {"\\value", LISTSECTION},

    {"\\describe", DESCRIPTION},
    {"\\enumerate", ITEMIZE},
    {"\\itemize", ITEMIZE},

    {"\\item", NOITEM}, /* will change to UNKNOWN, ESCAPE, or LATEXMACRO2 depending on context */

    /* These macros take two LaTeX-like arguments. */

    {"\\enc", LATEXMACRO2},
    {"\\method", LATEXMACRO2},
    {"\\S3method", LATEXMACRO2},
    {"\\S4method", LATEXMACRO2},
    {"\\tabular", LATEXMACRO2},

    /* These macros take one optional bracketed option and always take
       one LaTeX-like argument */

    {"\\link", OPTMACRO},

    /* These markup macros require an R-like text argument */

    {"\\code", RCODEMACRO},
    {"\\dontrun", VERBMACRO}, /* at least for now */
    {"\\dontshow", RCODEMACRO},
    {"\\donttest", RCODEMACRO},
    {"\\testonly", RCODEMACRO},

    /* These macros take one verbatim arg and ignore everything except braces */

    {"\\command", VERBMACRO},
    {"\\env", VERBMACRO},
    {"\\kbd", VERBMACRO},
    {"\\option", VERBMACRO},
    {"\\preformatted", VERBMACRO},

    {"\\samp", VERBMACRO},
    {"\\special", VERBMACRO},
    {"\\url", VERBMACRO},

    /* These ones take one or two verbatim args */

    {"\\eqn", VERBMACRO2},
    {"\\deqn", VERBMACRO2},

    /* We parse IFDEF/IFNDEF as markup, not as a separate preprocessor step */

    {"#ifdef", IFDEF},
    {"#ifndef", IFDEF},
    {"#endif", ENDIF},

    {0, 0}
    /* All other markup macros are rejected. */
};

/* Record the longest # directive here */
#define DIRECTIVE_LEN 7

static int KeywordLookup(const char *s)
{
    int i;
    for (i = 0; keywords[i].name; i++)
    {
        if (strcmp(keywords[i].name, s) == 0)
        {
            return keywords[i].token;
        }
    }
    return UNKNOWN;
}

static void yyerror(char *s)
{
    static const char *const yytname_translations[] = {
    /* the left column are strings coming from bison, the right
       column are translations for users.
       The first YYENGLISH from the right column are English to be translated,
       the rest are to be copied literally.  The #if 0 block below allows xgettext
       to see these.
    */
#define YYENGLISH 16
        "$undefined",     "input",          "SECTIONHEADER",
        "macro",          "RSECTIONHEADER", "macro",
        "VSECTIONHEADER", "macro",          "LISTSECTION",
        "macro",

        "LATEXMACRO",     "macro",          "LATEXMACRO2",
        "macro",          "RCODEMACRO",     "macro",
        "VERBMACRO",      "macro",          "VERBMACRO2",
        "macro",

        "ESCAPE",         "macro",          "ITEMIZE",
        "macro",          "IFDEF",          "conditional",
        "SECTIONHEADER2", "macro",          "OPTMACRO",
        "macro",

        "DESCRIPTION",    "macro",          0};
    static char const yyunexpected[] = "syntax error, unexpected ";
    static char const yyexpecting[] = ", expecting ";
    char *expecting;
#if 0
 /* these are just here to trigger the internationalization */
    _("input"); 	
    _("macro");
    _("conditional");
#endif

    R_ParseError = xxlineno;
    R_ParseErrorFile = SrcFile;

    if (!strncmp(s, yyunexpected, sizeof yyunexpected - 1))
    {
        int i;
        /* Edit the error message */
        expecting = strstr(s + sizeof yyunexpected - 1, yyexpecting);
        if (expecting)
            *expecting = '\0';
        for (i = 0; yytname_translations[i]; i += 2)
        {
            if (!strcmp(s + sizeof yyunexpected - 1, yytname_translations[i]))
            {
                sprintf(R_ParseErrorMsg, _("unexpected %s"),
                        i / 2 < YYENGLISH ? _(yytname_translations[i + 1]) : yytname_translations[i + 1]);
                return;
            }
        }
        sprintf(R_ParseErrorMsg, _("unexpected %s"), s + sizeof yyunexpected - 1);
    }
    else
    {
        strncpy(R_ParseErrorMsg, s, PARSE_ERROR_SIZE - 1);
    }
}

#define TEXT_PUSH(c)                                                                                                   \
    do                                                                                                                 \
    {                                                                                                                  \
        unsigned int nc = bp - stext;                                                                                  \
        if (nc >= nstext - 1)                                                                                          \
        {                                                                                                              \
            char *old = stext;                                                                                         \
            nstext *= 2;                                                                                               \
            stext = malloc(nstext);                                                                                    \
            if (!stext)                                                                                                \
                error(_("unable to allocate buffer for long string at line %d"), xxlineno);                            \
            memmove(stext, old, nc);                                                                                   \
            if (old != st0)                                                                                            \
                free(old);                                                                                             \
            bp = stext + nc;                                                                                           \
        }                                                                                                              \
        *bp++ = (c);                                                                                                   \
    } while (0)

static void setfirstloc(void)
{
    yylloc.first_line = xxlineno;
    yylloc.first_column = xxcolno;
}

static void setlastloc(void)
{
    yylloc.last_line = xxlineno;
    yylloc.last_column = xxcolno;
}

/* Split the input stream into tokens. */
/* This is the lowest of the parsing levels. */

static int token(void)
{
    int c;
    int outsideLiteral = xxmode == LATEXLIKE || xxmode == INOPTION || xxbraceDepth == 0;

    setfirstloc();
    c = xxgetc();

    /* % comments are active everywhere */

    if (c == '%')
        return mkComment(c);

    if (c == '\\')
    {
        int lookahead = xxgetc();
        xxungetc(lookahead);
        if (xxmode == VERBATIM)
        {
            if (lookahead == LBRACE || lookahead == RBRACE)
                return mkVerb(c);
        }
        else
        {
            if (xxinRString && lookahead != 'l')
                return mkCode(c);

            return mkMarkup(c);
        }
    }

    if (xxinRString)
    {
        if (c == R_EOF)
            error(_("Unexpected end of input (in %c quoted string opened at %d:%d)"), xxinRString, xxQuoteLine,
                  xxQuoteCol);
        return mkCode(c);
    }

    if (c == R_EOF)
        return END_OF_INPUT;

    if (c == '#' && xxcolno == 1)
        return mkIfdef(c);

    if (c == LBRACE)
    {
        xxbraceDepth++;
        if (outsideLiteral)
            return c;
    }

    if (c == RBRACE)
    {
        xxbraceDepth--;
        if (outsideLiteral || xxbraceDepth == 0)
            return c;
    }

    if ((c == '[' || c == ']') && xxmode == INOPTION)
        return c;

    switch (xxmode)
    {
    case RLIKE:
        return mkCode(c);
    case INOPTION:
    case LATEXLIKE:
        return mkText(c);
    case VERBATIM:
        return mkVerb(c);
    }

    return ERROR; /* We shouldn't get here. */
}

#define INITBUFSIZE 128

static int mkText(int c)
{
    char st0[INITBUFSIZE];
    unsigned int nstext = INITBUFSIZE;
    char *stext = st0, *bp = st0, lookahead;

    while (1)
    {
        switch (c)
        {
        case '\\':
            lookahead = xxgetc();
            if (lookahead == LBRACE || lookahead == RBRACE || lookahead == '%')
            {
                c = lookahead;
                break;
            }
            xxungetc(lookahead);
            goto stop;
        case ']':
            if (xxmode == INOPTION)
                goto stop;
            break;
        case '%':
        case LBRACE:
        case RBRACE:
        case R_EOF:
            goto stop;
        }
        TEXT_PUSH(c);
        if (c == '\n')
            goto stop;
        c = xxgetc();
    };
stop:
    if (c != '\n')
        xxungetc(c); /* newline causes a break, but we keep it */
    PROTECT(yylval = mkString2(stext, bp - stext));
    if (stext != st0)
        free(stext);
    return TEXT;
}

static int mkComment(int c)
{
    char st0[INITBUFSIZE];
    unsigned int nstext = INITBUFSIZE;
    char *stext = st0, *bp = st0;

    do
        TEXT_PUSH(c);
    while ((c = xxgetc()) != '\n' && c != R_EOF);
    if (c == R_EOF)
        xxungetc(c);
    PROTECT(yylval = mkString2(stext, bp - stext));
    if (stext != st0)
        free(stext);
    return COMMENT;
}

static int mkCode(int c)
{
    char st0[INITBUFSIZE];
    unsigned int nstext = INITBUFSIZE;
    char *stext = st0, *bp = st0;

    /* Avoid double counting initial braces */
    if (c == LBRACE)
        xxbraceDepth--;
    if (c == RBRACE)
        xxbraceDepth++;

    while (1)
    {
        int escaped = 0;
        if (c == '\\')
        {
            int lookahead = xxgetc();
            if (lookahead == '\\' || lookahead == '%')
            {
                c = lookahead;
                escaped = 1;
            }
            else
                xxungetc(lookahead);
        }
        if ((!escaped && c == '%') || c == R_EOF)
            break;
        if (xxinRString)
        {
            /* This stuff is messy, because there are two levels of escaping:
               The Rd escaping and the R code string escaping. */
            if (c == '\\')
            {
                int lookahead = xxgetc();
                if (lookahead == '\\')
                { /* This must be the 3rd backslash */
                    lookahead = xxgetc();
                    if (lookahead == xxinRString || lookahead == '\\')
                    {
                        TEXT_PUSH(c);
                        c = lookahead;
                        escaped = 1;
                    }
                    else
                        xxungetc(lookahead);
                }
                else if (lookahead == xxinRString)
                { /* There could be one or two before this */
                    TEXT_PUSH(c);
                    c = lookahead;
                    escaped = 1;
                }
                else if (!escaped && lookahead == 'l')
                { /* assume \link */
                    xxungetc(lookahead);
                    break;
                }
                else
                    xxungetc(lookahead);
            }
            if (!escaped && c == xxinRString)
                xxinRString = 0;
        }
        else
        {
            if (c == '#')
            {
                do
                {
                    TEXT_PUSH(c);
                    c = xxgetc();
                    if (c == LBRACE)
                        xxbraceDepth++;
                    else if (c == RBRACE)
                        xxbraceDepth--;
                } while (c != '\n' && c != R_EOF && xxbraceDepth > 0);
                if (c == RBRACE)
                    xxbraceDepth++; /* avoid double counting */
            }
            if (c == '\'' || c == '"' || c == '`')
            {
                xxinRString = c;
                xxQuoteLine = xxlineno;
                xxQuoteCol = xxcolno;
            }
            else if (c == '\\' && !escaped)
            {
                int lookahead = xxgetc();
                if (lookahead == LBRACE || lookahead == RBRACE)
                {
                    c = lookahead;
                }
                else if (isalpha(lookahead))
                {
                    xxungetc(lookahead);
                    c = '\\';
                    break;
                }
                else
                {
                    TEXT_PUSH('\\');
                    c = lookahead;
                }
            }
            else if (c == LBRACE)
            {
                xxbraceDepth++;
            }
            else if (c == RBRACE)
            {
                if (xxbraceDepth == 1)
                    break;
                else
                    xxbraceDepth--;
            }
            else if (c == R_EOF)
                break;
        }
        TEXT_PUSH(c);
        if (c == '\n')
            break;
        c = xxgetc();
    }
    if (c != '\n')
        xxungetc(c);
    PROTECT(yylval = mkString2(stext, bp - stext));
    if (stext != st0)
        free(stext);
    return RCODE;
}

static int mkMarkup(int c)
{
    char st0[INITBUFSIZE];
    unsigned int nstext = INITBUFSIZE;
    char *stext = st0, *bp = st0;
    int retval, attempt = 0;

    TEXT_PUSH(c);
    while (isalnum((c = xxgetc())))
        TEXT_PUSH(c);

    while (attempt++ < 2)
    {
        /* character escapes are processed as text, not markup */
        if (bp == stext + 1)
        {
            TEXT_PUSH(c);
            TEXT_PUSH('\0');
            retval = TEXT;
            c = xxgetc();
            break;
        }
        else
        {
            TEXT_PUSH('\0');
            retval = KeywordLookup(stext);
            if (retval == UNKNOWN && attempt == 1)
            {         /* try again, non-digits only */
                bp--; /* pop the \0 */
                while (isdigit(*(bp - 1)))
                {
                    xxungetc(c);
                    c = *(--bp); /* pop the last letter into c */
                }
            }
            else
            {
                if (retval == NOITEM)
                    retval = xxitemType;
                break;
            }
        }
    }
    PROTECT(yylval = mkString2(stext, bp - stext - 1));
    if (stext != st0)
        free(stext);
    xxungetc(c);
    return retval;
}

static int mkIfdef(int c)
{
    char st0[INITBUFSIZE];
    unsigned int nstext = INITBUFSIZE;
    char *stext = st0, *bp = st0;
    int retval;

    TEXT_PUSH(c);
    while (isalpha((c = xxgetc())) && bp - stext <= DIRECTIVE_LEN)
        TEXT_PUSH(c);
    TEXT_PUSH('\0');
    xxungetc(c);
    retval = KeywordLookup(stext);
    PROTECT(yylval = mkString2(stext, bp - stext - 1));

    if (retval == UNKNOWN)
    {
        UNPROTECT(1);
        bp--;
        bp--;
        for (; bp > stext; bp--)
            xxungetc(*bp);
        switch (xxmode)
        {
        case RLIKE:
            retval = mkCode(*bp);
            break;
        case INOPTION:
        case LATEXLIKE:
            retval = mkText(*bp);
            break;
        case VERBATIM:
            retval = mkVerb(*bp);
            break;
        }
    }
    if (stext != st0)
        free(stext);
    return retval;
}

static int mkVerb(int c)
{
    char st0[INITBUFSIZE];
    unsigned int nstext = INITBUFSIZE;
    char *stext = st0, *bp = st0;

    /* Avoid double counting initial braces */
    if (c == LBRACE)
        xxbraceDepth--;
    if (c == RBRACE)
        xxbraceDepth++;

    while (1)
    {
        int escaped = 0;
        if (c == '\\')
        {
            int lookahead = xxgetc();
            if (lookahead == '\\' || lookahead == '%' || lookahead == LBRACE || lookahead == RBRACE)
            {
                c = lookahead;
                escaped = 1;
            }
            else
                xxungetc(lookahead);
        }
        if ((!escaped && c == '%') || c == R_EOF)
            break;
        if (!escaped && c == LBRACE)
            xxbraceDepth++;
        else if (!escaped && c == RBRACE)
        {
            if (xxbraceDepth == 1)
                break;
            else
                xxbraceDepth--;
        }
        else if ((!escaped && c == '%') || c == R_EOF)
            break;
        TEXT_PUSH(c);
        if (c == '\n')
            break;
        c = xxgetc();
    };
    if (c != '\n')
        xxungetc(c);
    PROTECT(yylval = mkString2(stext, bp - stext));
    if (stext != st0)
        free(stext);
    return VERB;
}

static int yylex(void)
{
    int tok = token();

    if (xxDebugTokens)
    {
        Rprintf("%d:%d: %s", yylloc.first_line, yylloc.first_column + 1, yytname[YYTRANSLATE(tok)]);
        if (xxinRString)
            Rprintf("(in %c%c)", xxinRString, xxinRString);
        if (tok > 255 && tok != END_OF_INPUT)
            Rprintf(": %s", CHAR(STRING_ELT(yylval, 0)));
        Rprintf("\n");
    }
    setlastloc();
    return tok;
}

/* "do_parseRd"

 .Internal( parseRd(file, srcfile, encoding, verbose) )
 If there is text then that is read and the other arguments are ignored.
*/

SEXP attribute_hidden do_parseRd(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP s = R_NilValue, source;
    Rconnection con;
    Rboolean wasopen, old_latin1 = known_to_be_latin1, old_utf8 = known_to_be_utf8;
    int ifile;
    const char *encoding;
    ParseStatus status;

#if DEBUGMODE
    yydebug = 1;
#endif

    checkArity(op, args);
    R_ParseError = 0;
    R_ParseErrorMsg[0] = '\0';

    ifile = asInteger(CAR(args));
    args = CDR(args);

    con = getConnection(ifile);
    wasopen = con->isopen;
    source = CAR(args);
    args = CDR(args);
    if (!isString(CAR(args)) || LENGTH(CAR(args)) != 1)
        error(_("invalid '%s' value"), "encoding");
    encoding = CHAR(STRING_ELT(CAR(args), 0)); /* ASCII */
    args = CDR(args);
    known_to_be_latin1 = known_to_be_utf8 = FALSE;
    if (streql(encoding, "latin1"))
        known_to_be_latin1 = TRUE;
    if (streql(encoding, "UTF-8"))
        known_to_be_utf8 = TRUE;
    if (!isLogical(CAR(args)) || LENGTH(CAR(args)) != 1)
        error(_("invalid '%s' value"), "verbose");
    xxDebugTokens = asInteger(CAR(args));

    if (ifile >= 3)
    { /* file != "" */
        if (!wasopen)
        {
            if (!con->open(con))
                error(_("cannot open the connection"));
            if (!con->canread)
            {
                con->close(con);
                error(_("cannot read from this connection"));
            }
        }
        else if (!con->canread)
            error(_("cannot read from this connection"));
        s = R_ParseRd(con, &status, source);
        if (!wasopen)
            con->close(con);
        if (status != PARSE_OK)
            parseError(call, R_ParseError);
    }
    else
        error(_("invalid Rd file"));
    known_to_be_latin1 = old_latin1;
    known_to_be_utf8 = old_utf8;
    return s;
}
