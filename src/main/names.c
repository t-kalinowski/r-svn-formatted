/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997--2002  Robert Gentleman, Ross Ihaka and the
 *                            R Development Core Team
 *  Copyright (C) 2003        The R Foundation
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

#define __R_Names__
#include "Defn.h"
#include "Print.h"
#include "arithmetic.h"

#include <R_ext/RConverters.h>

/* Table of  .Internal(.) and .Primitive(.)  R functions
 * =====     =========	      ==========
 *
 * Each entry is a line with
 *
 *  printname	c-entry	 offset	 eval	arity	  pp-kind   precedence	    rightassoc
 *  ---------	-------	 ------	 ----	-----	  -------   ----------	    ----------
 *2 name	cfun	 code	 eval	arity	  gram.kind gram.precedence gram.rightassoc
 *3 PRIMNAME	PRIMFUN	 PRIMVAL [*]    PRIMARITY PPINFO    PPINFO 	    PPINFO
 *
 * where "2" are the component names of the FUNTAB struct (Defn.h)
 * and	 "3" are the accessor macros. [*]: PRIMPRINT(.) uses the eval component
 *
 * printname:	The function name in R
 *
 * c-entry:	The name of the corresponding C function,
 *		actually declared in ../include/Internal.h .
 *		Convention:
 *		 - all start with "do_",
 *		 - all return SEXP.
 *		 - all have argument list
 *			 (SEXP call, SEXP op, SEXP args, SEXP env)
 *
 * offset:	the 'op' (offset pointer) above; used for C functions
 *		which deal with more than one R function...
 *
 * eval:	= XYZ (three digits)  [New Apr 9/96, before only had "YZ"].
 *		  --- where e.g. '1' means '001'
 *		X=1 says that we should switch R_Visible off
 *		    (the least common situation).
 *		Y=1 says that this is an internal function which must
 *		    be accessed with a	.Internal(.) call, any other value is
 *		    accessible directly and printed in R as ".Primitive(..)".
 *		Z=1 says evaluate arguments before calling and
 *		Z=0 says don't evaluate.
 *
 * arity:	How many arguments are required/allowed;  "-1"	meaning ``any''
 *
 * pp-kind:	Deparsing Info (-> PPkind in ../include/Defn.h )
 *
 * precedence: Operator precedence (-> PPprec in ../include/Defn.h )
 *
 * rightassoc: Right (1) or left (0) associative operator
 *
 */
FUNTAB R_FunTab[] = {

    /* Language Related Constructs */

    /* printname	c-entry		offset	eval	arity	pp-kind	     precedence	rightassoc
     * ---------	-------		------	----	-----	-------      ----------	----------*/
    {"if", do_if, 0, 0, -1, {PP_IF, PREC_FN, 1}},
    {"while", do_while, 0, 0, -1, {PP_WHILE, PREC_FN, 0}},
    {"for", do_for, 0, 0, -1, {PP_FOR, PREC_FN, 0}},
    {"repeat", do_repeat, 0, 0, -1, {PP_REPEAT, PREC_FN, 0}},
    {"break", do_break, CTXT_BREAK, 0, -1, {PP_BREAK, PREC_FN, 0}},
    {"next", do_break, CTXT_NEXT, 0, -1, {PP_NEXT, PREC_FN, 0}},
    {"return", do_return, 0, 0, -1, {PP_RETURN, PREC_FN, 0}},
    {"stop", do_stop, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"warning", do_warning, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
#ifdef NEW_CONDITION_HANDLING
    {".addCondHands", do_addCondHands, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {".resetCondHands", do_resetCondHands, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {".signalCondition", do_signalCondition, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {".dfltStop", do_dfltStop, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {".dfltWarn", do_dfltWarn, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {".addRestart", do_addRestart, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {".getRestart", do_getRestart, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {".invokeRestart", do_invokeRestart, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
#endif
    {".addTryHandlers", do_addTryHandlers, 0, 111, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"geterrmessage", do_geterrmessage, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"restart", do_restart, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"function", do_function, 0, 0, -1, {PP_FUNCTION, PREC_FN, 0}},
    {"as.function.default", do_asfunction, 0, 11, 2, {PP_FUNCTION, PREC_FN, 0}},
    {"<-", do_set, 1, 100, -1, {PP_ASSIGN, PREC_LEFT, 1}},
    {"=", do_set, 3, 100, -1, {PP_ASSIGN, PREC_EQ, 1}},
    {"<<-", do_set, 2, 100, -1, {PP_ASSIGN2, PREC_LEFT, 1}},
    {"{", do_begin, 0, 0, -1, {PP_CURLY, PREC_FN, 0}},
    {"(", do_paren, 0, 1, 1, {PP_PAREN, PREC_FN, 0}},
    {".subset", do_subset_dflt, 1, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {".subset2", do_subset2_dflt, 2, 1, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"[", do_subset, 1, 0, -1, {PP_SUBSET, PREC_SUBSET, 0}},
    {"[[", do_subset2, 2, 0, 2, {PP_SUBSET, PREC_SUBSET, 0}},
    {"$", do_subset3, 3, 0, 2, {PP_DOLLAR, PREC_DOLLAR, 0}},
    /*{"::",	do_NS_get,	?,	?,	2,	{PP_BINARY2, PREC_NS,	  0}},*/
    {"@", do_AT, 0, 0, 2, {PP_DOLLAR, PREC_DOLLAR, 0}},
    {"[<-", do_subassign, 0, 0, 3, {PP_SUBASS, PREC_LEFT, 1}},
    {"[[<-", do_subassign2, 1, 100, 3, {PP_SUBASS, PREC_LEFT, 1}},
    {"$<-", do_subassign3, 1, 0, 3, {PP_SUBASS, PREC_LEFT, 1}},
    {"switch", do_switch, 0, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"browser", do_browser, 0, 100, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"debug", do_debug, 0, 101, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"undebug", do_debug, 1, 101, 1, {PP_FUNCALL, PREC_FN, 0}},
    {".primTrace", do_trace, 0, 101, 1, {PP_FUNCALL, PREC_FN, 0}},
    {".primUntrace", do_trace, 1, 101, 1, {PP_FUNCALL, PREC_FN, 0}},
    {".Internal", do_internal, 0, 0, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"on.exit", do_onexit, 0, 100, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Recall", do_recall, 0, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"delay", do_delay, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    /*{".Alias",	do_alias,	0,	1,	1,	{PP_FUNCALL, PREC_FN,	  0}},*/
    {".Primitive", do_primitive, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"identical", do_ident, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    /* Binary Operators */
    {"+", do_arith, PLUSOP, 1, 2, {PP_BINARY, PREC_SUM, 0}},
    {"-", do_arith, MINUSOP, 1, 2, {PP_BINARY, PREC_SUM, 0}},
    {"*", do_arith, TIMESOP, 1, 2, {PP_BINARY, PREC_PROD, 0}},
    {"/", do_arith, DIVOP, 1, 2, {PP_BINARY2, PREC_PROD, 0}},
    {"^", do_arith, POWOP, 1, 2, {PP_BINARY2, PREC_POWER, 1}},
    {"%%", do_arith, MODOP, 1, 2, {PP_BINARY2, PREC_PERCENT, 0}},
    {"%/%", do_arith, IDIVOP, 1, 2, {PP_BINARY2, PREC_PERCENT, 0}},
    {"%*%", do_matprod, 0, 1, 2, {PP_BINARY, PREC_PERCENT, 0}},
    {"crossprod", do_matprod, 1, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"==", do_relop, EQOP, 1, 2, {PP_BINARY, PREC_COMPARE, 0}},
    {"!=", do_relop, NEOP, 1, 2, {PP_BINARY, PREC_COMPARE, 0}},
    {"<", do_relop, LTOP, 1, 2, {PP_BINARY, PREC_COMPARE, 0}},
    {"<=", do_relop, LEOP, 1, 2, {PP_BINARY, PREC_COMPARE, 0}},
    {">=", do_relop, GEOP, 1, 2, {PP_BINARY, PREC_COMPARE, 0}},
    {">", do_relop, GTOP, 1, 2, {PP_BINARY, PREC_COMPARE, 0}},
    {"&", do_logic, 1, 1, 2, {PP_BINARY, PREC_AND, 0}},
    {"|", do_logic, 2, 1, 2, {PP_BINARY, PREC_OR, 0}},
    {"!", do_logic, 3, 1, 1, {PP_UNARY, PREC_NOT, 0}},
    {"&&", do_logic2, 1, 0, 2, {PP_BINARY, PREC_AND, 0}},
    {"||", do_logic2, 2, 0, 2, {PP_BINARY, PREC_OR, 0}},
    {":", do_seq, 0, 1, 2, {PP_BINARY2, PREC_COLON, 0}},
    {"~", do_tilde, 0, 0, 2, {PP_BINARY, PREC_TILDE, 0}},

    /* Logic Related Functions */

    {"all", do_logic3, 1, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"any", do_logic3, 2, 11, -1, {PP_FUNCALL, PREC_FN, 0}},

    /* Vectors, Matrices and Arrays */

    /* printname	c-entry		offset	eval	arity	pp-kind	     precedence	rightassoc
     * ---------	-------		------	----	-----	-------      ----------	----------*/
    {"vector", do_makevector, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"complex", do_complex, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"matrix", do_matrix, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    /*MM(98/4/22){"array",	do_array,0,	1,	2,	{PP_FUNCALL, PREC_FN,	0}},*/
    {"length", do_length, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"length<-", do_lengthgets, 0, 1, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"row", do_rowscols, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"col", do_rowscols, 2, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"c", /* bind.c:*/ do_c, 0, 0, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"unlist", do_unlist, 0, 10, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"cbind", do_bind, 1, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"rbind", do_bind, 2, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"drop", do_drop, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"oldClass", do_class, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"oldClass<-", do_classgets, 0, 1, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"class", R_do_data_class, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"class<-", R_do_set_class, 0, 1, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"unclass", do_unclass, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"names", do_names, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"names<-", do_namesgets, 0, 11, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"dimnames", do_dimnames, 0, 0, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dimnames<-", do_dimnamesgets, 0, 0, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"all.names", do_allnames, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"dim", do_dim, 0, 0, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dim<-", do_dimgets, 0, 0, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"attributes", do_attributes, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"attributes<-", do_attributesgets, 0, 1, 1, {PP_FUNCALL, PREC_LEFT, 1}},
    {"attr", do_attr, 0, 1, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"attr<-", do_attrgets, 0, 0, 3, {PP_FUNCALL, PREC_LEFT, 1}},
    {"comment", do_comment, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"comment<-", do_commentgets, 0, 11, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"get", do_get, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"exists", do_get, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"assign", do_assign, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"remove", do_remove, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"duplicated", do_duplicated, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"unique", do_duplicated, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"which.min", do_first_min, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"which.max", do_first_max, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"match", do_match, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"pmatch", do_pmatch, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"charmatch", do_charmatch, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"match.call", do_matchcall, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"complete.cases", do_compcases, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"attach", do_attach, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"detach", do_detach, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"search", do_search, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},

    /* Mathematical Functions */
    {"round", do_Math2, 10001, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"atan", do_atan, 10002, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"log", do_log, 10003, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"signif", do_Math2, 10004, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"abs", do_abs, 6, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    /* KH(1999/09/12)-> complex: {"abs", do_math1, 0, 1, 1, {PP_FUNCALL, PREC_FN,	0}}, */
    {"floor", do_math1, 1, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"ceiling", do_math1, 2, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"sqrt", do_math1, 3, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"sign", do_math1, 4, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"trunc", do_math1, 5, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"exp", do_math1, 10, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"expm1", do_math1, 11, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"log1p", do_math1, 12, 11, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"cos", do_math1, 20, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"sin", do_math1, 21, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"tan", do_math1, 22, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"acos", do_math1, 23, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"asin", do_math1, 24, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"cosh", do_math1, 30, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"sinh", do_math1, 31, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"tanh", do_math1, 32, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"acosh", do_math1, 33, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"asinh", do_math1, 34, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"atanh", do_math1, 35, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"lgamma", do_math1, 40, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"gamma", do_math1, 41, 11, 1, {PP_FUNCALL, PREC_FN, 0}},

    /* Polygamma Functions */

    {"digamma", do_math1, 42, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"trigamma", do_math1, 43, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"tetragamma", do_math1, 44, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pentagamma", do_math1, 45, 11, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"gammaCody", do_math1, 46, 11, 1, {PP_FUNCALL, PREC_FN, 0}},

    /* Mathematical Functions of Two Numeric (+ 1-2 int) Variables */

    {"atan2", do_math2, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    {"lbeta", do_math2, 2, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"beta", do_math2, 3, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"lchoose", do_math2, 4, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"choose", do_math2, 5, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dchisq", do_math2, 6, 11, 2 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pchisq", do_math2, 7, 11, 2 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qchisq", do_math2, 8, 11, 2 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dexp", do_math2, 9, 11, 2 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pexp", do_math2, 10, 11, 2 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qexp", do_math2, 11, 11, 2 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dgeom", do_math2, 12, 11, 2 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pgeom", do_math2, 13, 11, 2 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qgeom", do_math2, 14, 11, 2 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dpois", do_math2, 15, 11, 2 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"ppois", do_math2, 16, 11, 2 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qpois", do_math2, 17, 11, 2 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dt", do_math2, 18, 11, 2 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pt", do_math2, 19, 11, 2 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qt", do_math2, 20, 11, 2 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dsignrank", do_math2, 21, 11, 2 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"psignrank", do_math2, 22, 11, 2 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qsignrank", do_math2, 23, 11, 2 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"besselJ", do_math2, 24, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"besselY", do_math2, 25, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    /* Mathematical Functions of a Complex Argument */

    {"Re", do_cmathfuns, 1, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Im", do_cmathfuns, 2, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Mod", do_cmathfuns, 3, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Arg", do_cmathfuns, 4, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Conj", do_cmathfuns, 5, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    /* {"abs",	do_cmathfuns,	6,	1,	1,	{PP_FUNCALL, PREC_FN,	0}},*/

    /* Mathematical Functions of Three Numeric (+ 1-2 int) Variables */

    {"dbeta", do_math3, 1, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pbeta", do_math3, 2, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qbeta", do_math3, 3, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dbinom", do_math3, 4, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pbinom", do_math3, 5, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qbinom", do_math3, 6, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dcauchy", do_math3, 7, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pcauchy", do_math3, 8, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qcauchy", do_math3, 9, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"df", do_math3, 10, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pf", do_math3, 11, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qf", do_math3, 12, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dgamma", do_math3, 13, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pgamma", do_math3, 14, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qgamma", do_math3, 15, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dlnorm", do_math3, 16, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"plnorm", do_math3, 17, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qlnorm", do_math3, 18, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dlogis", do_math3, 19, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"plogis", do_math3, 20, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qlogis", do_math3, 21, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dnbinom", do_math3, 22, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pnbinom", do_math3, 23, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qnbinom", do_math3, 24, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dnorm", do_math3, 25, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pnorm", do_math3, 26, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qnorm", do_math3, 27, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dunif", do_math3, 28, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"punif", do_math3, 29, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qunif", do_math3, 30, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dweibull", do_math3, 31, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pweibull", do_math3, 32, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qweibull", do_math3, 33, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dnchisq", do_math3, 34, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pnchisq", do_math3, 35, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qnchisq", do_math3, 36, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dnt", do_math3, 37, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pnt", do_math3, 38, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qnt", do_math3, 39, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dwilcox", do_math3, 40, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pwilcox", do_math3, 41, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qwilcox", do_math3, 42, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"besselI", do_math3, 43, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"besselK", do_math3, 44, 11, 3, {PP_FUNCALL, PREC_FN, 0}},

    /* Mathematical Functions of Four Numeric (+ 1-2 int) Variables */

    {"dhyper", do_math4, 1, 11, 4 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"phyper", do_math4, 2, 11, 4 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qhyper", do_math4, 3, 11, 4 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dnbeta", do_math4, 4, 11, 4 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pnbeta", do_math4, 5, 11, 4 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qnbeta", do_math4, 6, 11, 4 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dnf", do_math4, 7, 11, 4 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pnf", do_math4, 8, 11, 4 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qnf", do_math4, 9, 11, 4 + 2, {PP_FUNCALL, PREC_FN, 0}},

    {"dtukey", do_math4, 10, 11, 4 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"ptukey", do_math4, 11, 11, 4 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qtukey", do_math4, 12, 11, 4 + 2, {PP_FUNCALL, PREC_FN, 0}},

    /* Random Numbers */

    {"rchisq", do_random1, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"rexp", do_random1, 1, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"rgeom", do_random1, 2, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"rpois", do_random1, 3, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"rt", do_random1, 4, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"rsignrank", do_random1, 5, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    {"rbeta", do_random2, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rbinom", do_random2, 1, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rcauchy", do_random2, 2, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rf", do_random2, 3, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rgamma", do_random2, 4, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rlnorm", do_random2, 5, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rlogis", do_random2, 6, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rnbinom", do_random2, 7, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rnchisq", do_random2, 12, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rnorm", do_random2, 8, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"runif", do_random2, 9, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rweibull", do_random2, 10, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rwilcox", do_random2, 11, 11, 3, {PP_FUNCALL, PREC_FN, 0}},

    {"rhyper", do_random3, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},

    {"rmultinom", do_rmultinom, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"sample", do_sample, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},

    {"RNGkind", do_RNGkind, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"set.seed", do_setseed, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    /* Data Summaries */

    {"sum", do_summary, 0, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    /*MM{"mean",	do_summary,	1,	11,	1,	{PP_FUNCALL, PREC_FN,	0}},*/
    {"min", do_summary, 2, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"max", do_summary, 3, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"prod", do_summary, 4, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"range", do_range, 0, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"cov", do_cov, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"cor", do_cov, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},

    {"cumsum", do_cum, 1, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"cumprod", do_cum, 2, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"cummax", do_cum, 3, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"cummin", do_cum, 4, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    /* Type coercion */

    {"as.character", do_ascharacter, 0, 0, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"as.vector", do_asvector, 0, 10, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"paste", do_paste, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"format", do_format, 0, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"format.info", do_formatinfo, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"cat", do_cat, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"call", do_call, 0, 0, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"do.call", do_docall, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"as.call", do_ascall, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"type.convert", do_typecvt, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"as.environment", do_as_environment, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    /* String Manipulation */

    {"nchar", do_nchar, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"substr", do_substr, 1, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"substrgets", do_substrgets, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"strsplit", do_strsplit, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"abbreviate", do_abbrev, 1, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"make.names", do_makenames, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"grep", do_grep, 1, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"sub", do_gsub, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"gsub", do_gsub, 1, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"regexpr", do_regexpr, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"grep.perl", do_pgrep, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"sub.perl", do_pgsub, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"gsub.perl", do_pgsub, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"regexpr.perl", do_pregexpr, 1, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"agrep", do_agrep, 1, 11, 8, {PP_FUNCALL, PREC_FN, 0}},
    {"tolower", do_tolower, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"toupper", do_toupper, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"chartr", do_chartr, 1, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"sprintf", do_sprintf, 1, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"make.unique", do_makeunique, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    /* Type Checking (typically implemented in ./coerce.c ) */

    {"is.null", do_is, NILSXP, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.logical", do_is, LGLSXP, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.integer", do_is, INTSXP, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.real", do_is, REALSXP, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.double", do_is, REALSXP, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.complex", do_is, CPLXSXP, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.character", do_is, STRSXP, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.symbol", do_is, SYMSXP, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.environment", do_is, ENVSXP, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.list", do_is, VECSXP, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.pairlist", do_is, LISTSXP, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.expression", do_is, EXPRSXP, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"is.object", do_is, 50, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"is.numeric", do_is, 100, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.matrix", do_is, 101, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.array", do_is, 102, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"is.atomic", do_is, 200, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.recursive", do_is, 201, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"is.call", do_is, 300, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.language", do_is, 301, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.function", do_is, 302, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"is.single", do_is, 999, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"is.vector", do_isvector, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"is.na", do_isna, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.nan", do_isnan, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.finite", do_isfinite, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.infinite", do_isinfinite, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    /* Miscellaneous */

    {"proc.time", do_proctime, 0, 1, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"gc.time", do_gctime, 0, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"Version", do_version, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"machine", do_machine, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    /*{"Machine",	do_Machine,	0,	11,	0,	{PP_FUNCALL, PREC_FN,	0}},*/
    {"commandArgs", do_commandArgs, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"int.unzip", do_int_unzip, 0, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
#ifdef Win32
    {"system", do_system, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
#else
    {"system", do_system, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
#endif
#ifdef Win32
    {"unlink", do_unlink, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"help.start", do_helpstart, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"show.help.item", do_helpitem, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"flush.console", do_flushconsole, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"win.version", do_winver, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"saveDevga", do_saveDevga, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"shell.exec", do_shellexec, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dir.create", do_dircreate, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"winDialog", do_windialog, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"winDialogString", do_windialogstring, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"winMenuAdd", do_winmenuadd, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"winMenuDel", do_winmenudel, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"memory.size", do_memsize, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"DLL.version", do_dllversion, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"bringToTop", do_bringtotop, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"select.list", do_selectlist, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"readClipboard", do_readClipboard, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"writeClipboard", do_writeClipboard, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"chooseFiles", do_chooseFiles, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
#endif
#if defined(__APPLE_CC__) && defined(HAVE_AQUA)
    {"wsbrowser", do_wsbrowser, 0, 11, 8, {PP_FUNCALL, PREC_FN, 0}},
    {"pkgbrowser", do_browsepkgs, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"data.manager", do_datamanger, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"package.manager", do_packagemanger, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"flush.console", do_flushconsole, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"hsbrowser", do_hsbrowser, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
#endif
    {"parse", do_parse, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"save", do_save, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"saveToConn", do_saveToConn, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"load", do_load, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"loadFromConn", do_loadFromConn, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"serializeToConn", do_serializeToConn, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"unserializeFromConn", do_unserializeFromConn, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"deparse", do_deparse, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"dput", do_dput, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"dump", do_dump, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"substitute", do_substitute, 0, 0, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"quote", do_quote, 0, 0, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"quit", do_quit, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"interactive", do_interactive, 0, 0, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"readline", do_readln, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"menu", do_menu, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"print.default", do_printdefault, 0, 111, 7, {PP_FUNCALL, PREC_FN, 0}},
    {"prmatrix", do_prmatrix, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"invisible", do_invisible, 0, 101, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"gc", do_gc, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"gcinfo", do_gcinfo, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"gctorture", do_gctorture, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"memory.profile", do_memoryprofile, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"rep", do_rep, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"list", do_makelist, 1, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"split", do_split, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"symbol.C", do_symbol, 0, 1, 1, {PP_FOREIGN, PREC_FN, 0}},
    {"symbol.For", do_symbol, 1, 1, 1, {PP_FOREIGN, PREC_FN, 0}},
    {"is.loaded", do_isloaded, 0, 1, -1, {PP_FOREIGN, PREC_FN, 0}},
    {".C", do_dotCode, 0, 1, -1, {PP_FOREIGN, PREC_FN, 0}},
    {".Fortran", do_dotCode, 1, 1, -1, {PP_FOREIGN, PREC_FN, 0}},
    {".External", do_External, 0, 1, -1, {PP_FOREIGN, PREC_FN, 0}},
    {".Call", do_dotcall, 0, 1, -1, {PP_FOREIGN, PREC_FN, 0}},
    {".External.graphics", do_Externalgr, 0, 1, -1, {PP_FOREIGN, PREC_FN, 0}},
    {".Call.graphics", do_dotcallgr, 0, 1, -1, {PP_FOREIGN, PREC_FN, 0}},
    {"dyn.load", do_dynload, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"dyn.unload", do_dynunload, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"ls", do_ls, 1, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"typeof", do_typeof, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"eval", do_eval, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"eval.with.vis", do_eval, 1, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"expression", do_expression, 1, 0, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.parent", do_sys, 1, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.call", do_sys, 2, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.frame", do_sys, 3, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.nframe", do_sys, 4, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.calls", do_sys, 5, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.frames", do_sys, 6, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.on.exit", do_sys, 7, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.parents", do_sys, 8, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.function", do_sys, 9, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"parent.frame", do_parentframe, 0, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sort", do_sort, 1, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"is.unsorted", do_isunsorted, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"psort", do_psort, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qsort", do_qsort, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"radixsort", do_radixsort, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"order", do_order, 0, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"rank", do_rank, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"missing", do_missing, 1, 0, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"nargs", do_nargs, 1, 0, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"scan", do_scan, 0, 11, 16, {PP_FUNCALL, PREC_FN, 0}},
    {"count.fields", do_countfields, 0, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"readTableHead", do_readtablehead, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"t.default", do_transpose, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"aperm", do_aperm, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"builtins", do_builtins, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"edit", do_edit, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"dataentry", do_dataentry, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"args", do_args, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"formals", do_formals, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"body", do_body, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"bodyCode", do_bodyCode, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"globalenv", do_globalenv, 0, 1, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"environment", do_envir, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"environment<-", do_envirgets, 0, 1, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"reg.finalizer", do_regFinaliz, 0, 1, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"options", do_options, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"sink", do_sink, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"sink.number", do_sinknumber, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"lib.fixup", do_libfixup, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"pos.to.env", do_pos2env, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"lapply", do_lapply, 0, 10, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"apply", do_apply, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"colSums", do_colsum, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"colMeans", do_colsum, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"rowSums", do_colsum, 2, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"rowMeans", do_colsum, 3, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"Rprof", do_Rprof, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"object.size", do_objectsize, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"mem.limits", do_memlimits, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"merge", do_merge, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"capabilities", do_capabilities, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"new.env", do_newenv, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"parent.env", do_parentenv, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"parent.env<-", do_parentenvgets, 0, 11, 2, {PP_FUNCALL, PREC_LEFT, 1}},
#if 0
{"visibleflag", do_visibleflag,	0,	1,	0,	{PP_FUNCALL, PREC_FN,	0}},
#endif

    /* Functions To Interact with the Operating System */

    {"file.show", do_fileshow, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"file.create", do_filecreate, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"file.remove", do_fileremove, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"file.rename", do_filerename, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"file.append", do_fileappend, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"file.symlink", do_filesymlink, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"list.files", do_listfiles, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"file.exists", do_fileexists, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"file.choose", do_filechoose, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"file.info", do_fileinfo, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"file.access", do_fileaccess, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"tempfile", do_tempfile, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"tempdir", do_tempdir, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"R.home", do_Rhome, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"date", do_date, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    /*{"Platform",	do_Platform,	0,	11,	0,	{PP_FUNCALL, PREC_FN,	0}},*/
    {"index.search", do_indexsearch, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"getenv", do_getenv, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"putenv", do_putenv, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"getwd", do_getwd, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"setwd", do_setwd, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"basename", do_basename, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dirname", do_dirname, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.info", do_sysinfo, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.sleep", do_syssleep, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"getlocale", do_getlocale, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"setlocale", do_setlocale, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"localeconv", do_localeconv, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"path.expand", do_pathexpand, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"getpid", do_sysgetpid, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},

    /* Complex Valued Functions */
    {"fft", do_fft, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"mvfft", do_mvfft, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"nextn", do_nextn, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"polyroot", do_polyroot, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},

    /* Device Drivers */

    {"PS", do_PS, 0, 111, 16, {PP_FUNCALL, PREC_FN, 0}},
    {"PicTeX", do_PicTeX, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"XFig", do_XFig, 0, 111, 11, {PP_FUNCALL, PREC_FN, 0}},
    {"PDF", do_PDF, 0, 111, 10, {PP_FUNCALL, PREC_FN, 0}},
#ifdef Win32
    {"devga", do_devga, 0, 111, 14, {PP_FUNCALL, PREC_FN, 0}},
#endif
#ifdef Unix
    {"X11", do_X11, 0, 111, 9, {PP_FUNCALL, PREC_FN, 0}},
    {"gnome", do_Gnome, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"GTK", do_GTK, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"Quartz", do_Quartz, 0, 111, 7, {PP_FUNCALL, PREC_FN, 0}},
#endif

    /* Graphics */

    {"dev.control", do_devcontrol, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dev.copy", do_devcopy, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dev.cur", do_devcur, 0, 111, 0, {PP_FUNCALL, PREC_FN, 0}},
    /*
    {"device",	do_device,	0,	111,	3,	{PP_FUNCALL, PREC_FN,	0}},
    */
    {"dev.next", do_devnext, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dev.off", do_devoff, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dev.prev", do_devprev, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dev.set", do_devset, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"rgb", do_rgb, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"rgb256", do_rgb, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"hsv", do_hsv, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"gray", do_gray, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"colors", do_colors, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"col2rgb", do_col2RGB, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"palette", do_palette, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"plot.new", do_plot_new, 0, 111, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"plot.window", do_plot_window, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"axis", do_axis, 0, 111, 12, {PP_FUNCALL, PREC_FN, 0}},
    {"plot.xy", do_plot_xy, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"text", do_text, 0, 111, 7, {PP_FUNCALL, PREC_FN, 0}},
    {"mtext", do_mtext, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"title", do_title, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"abline", do_abline, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"box", do_box, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rect", do_rect, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"polygon", do_polygon, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"par", do_par, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"readonly.pars", do_readonlypars, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"segments", do_segments, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"arrows", do_arrows, 0, 111, 9, {PP_FUNCALL, PREC_FN, 0}},
    {"layout", do_layout, 0, 111, 10, {PP_FUNCALL, PREC_FN, 0}},
    {"locator", do_locator, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"identify", do_identify, 0, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"strheight", do_strheight, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"strwidth", do_strwidth, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"contour", do_contour, 0, 11, 12, {PP_FUNCALL, PREC_FN, 0}},
    {"contourLines", do_contourLines, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"image", do_image, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"dend", do_dend, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"dend.window", do_dendwindow, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"replay", do_replay, 0, 111, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"erase", do_erase, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    /*{"dotplot",	do_dotplot,	0,	111,	1,	{PP_FUNCALL, PREC_FN,	0}}, */
    {"persp", do_persp, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"filledcontour", do_filledcontour, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    /* {"getDL",	do_getDL,	0,	111,	0,	{PP_FUNCALL, PREC_FN,	0}},
    {"getGPar",	do_getGPar,	0,	111,	0,	{PP_FUNCALL, PREC_FN,	0}}, */
    {"playDL", do_playDL, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"setGPar", do_setGPar, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"getSnapshot", do_getSnapshot, 0, 111, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"playSnapshot", do_playSnapshot, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"symbols", do_symbols, 0, 111, -1, {PP_FUNCALL, PREC_FN, 0}},

    /* Objects */
    {"inherits", do_inherits, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"UseMethod", do_usemethod, 0, 0, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"NextMethod", do_nextmethod, 0, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"standardGeneric", do_standardGeneric, 0, 1, -1, {PP_FUNCALL, PREC_FN, 0}},

    /* Modelling Functionality */

    {"nlm", do_nlm, 0, 11, 11, {PP_FUNCALL, PREC_FN, 0}},
    {"fmin", do_fmin, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"zeroin", do_zeroin, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"optim", do_optim, 0, 11, 7, {PP_FUNCALL, PREC_FN, 0}},
    {"optimhess", do_optimhess, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"terms.formula", do_termsform, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"update.formula", do_updateform, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"model.frame", do_modelframe, 0, 11, 8, {PP_FUNCALL, PREC_FN, 0}},
    {"model.matrix", do_modelmatrix, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    {"D", do_D, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"deriv.default", do_deriv, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},

    /* History manipulation */
    {"loadhistory", do_loadhistory, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"savehistory", do_savehistory, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},

    /* date-time manipulations */
    {"Sys.time", do_systime, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"as.POSIXct", do_asPOSIXct, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"as.POSIXlt", do_asPOSIXlt, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"format.POSIXlt", do_formatPOSIXlt, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"strptime", do_strptime, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

#ifdef BYTECODE
    {"mkCode", do_mkcode, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"bcClose", do_bcclose, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"is.builtin.internal", do_is_builtin_internal, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"disassemble", do_disassemble, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"bcVersion", do_bcversion, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"load.from.file", do_loadfile, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"save.to.file", do_savefile, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"putconst", do_putconst, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
#endif

    /* Connections */
    {"stdin", do_stdin, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"stdout", do_stdout, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"stderr", do_stderr, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"readLines", do_readLines, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"writeLines", do_writelines, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"readBin", do_readbin, 0, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"writeBin", do_writebin, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"readChar", do_readchar, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"writeChar", do_writechar, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"open", do_open, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"isOpen", do_isopen, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"isIncomplete", do_isincomplete, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"isSeekable", do_isseekable, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"close", do_close, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"flush", do_flush, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"file", do_url, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"url", do_url, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"pipe", do_pipe, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"fifo", do_fifo, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"gzfile", do_gzfile, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"unz", do_unz, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"bzfile", do_bzfile, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"seek", do_seek, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"truncate", do_truncate, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pushBack", do_pushback, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"pushBackLength", do_pushbacklength, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"textConnection", do_textconnection, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"socketConnection", do_sockconn, 0, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"sockSelect", do_sockselect, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"getAllConnections", do_getallconnections, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"summary.connection", do_sumconnection, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"download", do_download, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"nsl", do_nsl, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"gzcon", do_gzcon, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},

    {"readDCF", do_readDCF, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    {"getNumRtoCConverters", do_getNumRtoCConverters, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"getRtoCConverterDescriptions", do_getRtoCConverterDescriptions, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"getRtoCConverterStatus", do_getRtoCConverterStatus, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"setToCConverterActiveStatus", do_setToCConverterActiveStatus, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"removeToCConverterActiveStatus", do_setToCConverterActiveStatus, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"lockEnvironment", do_lockEnv, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"environmentIsLocked", do_envIsLocked, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"lockBinding", do_lockBnd, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"unlockBinding", do_lockBnd, 1, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"bindingIsLocked", do_bndIsLocked, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"makeActiveBinding", do_mkActiveBnd, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"bindingIsActive", do_bndIsActive, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"mkUnbound", do_mkUnbound, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"isNamespaceEnv", do_isNSEnv, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"registerNamespace", do_regNS, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"unregisterNamespace", do_unregNS, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"getRegisteredNamespace", do_getRegNS, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"getNamespaceRegistry", do_getNSRegistry, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"importIntoEnv", do_importIntoEnv, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},

    {NULL, NULL, 0, 0, 0, {0, PREC_FN, 0}},
};

SEXP do_primitive(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP name;
    int i;
    checkArity(op, args);
    name = CAR(args);
    if (!isString(name) || length(name) < 1 || STRING_ELT(name, 0) == R_NilValue)
        errorcall(call, "string argument required");
    for (i = 0; R_FunTab[i].name; i++)
        if (strcmp(CHAR(STRING_ELT(name, 0)), R_FunTab[i].name) == 0)
        {
            if ((R_FunTab[i].eval % 100) / 10)
                return mkPRIMSXP(i, R_FunTab[i].eval % 10);
            else
                return mkPRIMSXP(i, R_FunTab[i].eval % 10);
        }
    errorcall(call, "no such primitive function");
    return (R_NilValue); /* -Wall */
}

int StrToInternal(char *s)
{
    int i;
    for (i = 0; R_FunTab[i].name; i++)
        if (strcmp(s, R_FunTab[i].name) == 0)
            return i;
    return 0;
}

static void installFunTab(int i)
{
    if ((R_FunTab[i].eval % 100) / 10)
        SET_INTERNAL(install(R_FunTab[i].name), mkPRIMSXP(i, R_FunTab[i].eval % 10));
    else
        SET_SYMVALUE(install(R_FunTab[i].name), mkPRIMSXP(i, R_FunTab[i].eval % 10));
}

static void SymbolShortcuts()
{
    R_Bracket2Symbol = install("[[");
    R_BracketSymbol = install("[");
    R_BraceSymbol = install("{");
    R_TmpvalSymbol = install("*tmp*");
    R_ClassSymbol = install("class");
    R_DimNamesSymbol = install("dimnames");
    R_DimSymbol = install("dim");
    R_DollarSymbol = install("$");
    R_DotsSymbol = install("...");
    R_DropSymbol = install("drop");
    R_LevelsSymbol = install("levels");
    R_ModeSymbol = install("mode");
    R_NamesSymbol = install("names");
    R_NaRmSymbol = install("na.rm");
    R_RowNamesSymbol = install("row.names");
    R_SeedsSymbol = install(".Random.seed");
    R_LastvalueSymbol = install(".Last.value");
    R_TspSymbol = install("tsp");
    R_CommentSymbol = install("comment");
    R_SourceSymbol = install("source");
    R_DotEnvSymbol = install(".Environment");
    R_RecursiveSymbol = install("recursive");
    R_UseNamesSymbol = install("use.names");
}

extern SEXP framenames; /* from model.c */

/* initialize the symbol table */
void InitNames()
{
    int i;
    /* R_UnboundValue */
    R_UnboundValue = allocSExp(SYMSXP);
    SET_SYMVALUE(R_UnboundValue, R_UnboundValue);
    SET_PRINTNAME(R_UnboundValue, R_NilValue);
    SET_ATTRIB(R_UnboundValue, R_NilValue);
    /* R_MissingArg */
    R_MissingArg = allocSExp(SYMSXP);
    SET_SYMVALUE(R_MissingArg, R_MissingArg);
    SET_PRINTNAME(R_MissingArg, mkChar(""));
    SET_ATTRIB(R_MissingArg, R_NilValue);
    /* R_RestartToken */
    R_RestartToken = allocSExp(SYMSXP);
    SET_SYMVALUE(R_RestartToken, R_RestartToken);
    SET_PRINTNAME(R_RestartToken, mkChar(""));
    SET_ATTRIB(R_RestartToken, R_NilValue);
    /* Parser Structures */
    R_CommentSxp = R_NilValue;
    R_ParseText = R_NilValue;
    /* String constants (CHARSXP values */
    /* Note: changed from mkChar so mkChar can see if it is getting
       "NA" and then return NA_STRING rather than alloc a new CHAR */
    /* NA_STRING */
    NA_STRING = allocString(strlen("NA"));
    strcpy(CHAR(NA_STRING), "NA");
    R_print.na_string = NA_STRING;
    /* R_BlankString */
    R_BlankString = mkChar("");
    /* Initialize the symbol Table */
    if (!(R_SymbolTable = (SEXP *)malloc(HSIZE * sizeof(SEXP))))
        R_Suicide("couldn't allocate memory for symbol table");
    for (i = 0; i < HSIZE; i++)
        R_SymbolTable[i] = R_NilValue;
    /* Set up a set of globals so that a symbol table search can be
       avoided when matching something like dim or dimnames. */
    SymbolShortcuts();
    /*  Builtin Functions */
    for (i = 0; R_FunTab[i].name; i++)
        installFunTab(i);
    /*  Unbound values which are to be preserved through GCs */
    R_PreciousList = R_NilValue;
    framenames = R_NilValue;
#ifdef BYTECODE
    R_initialize_bcode();
#endif
}

/*  install - probe the symbol table */
/*  If "name" is not found, it is installed in the symbol table.
    The symbol corresponding to the string "name" is returned. */

SEXP install(char const *name)
{
    char buf[MAXIDSIZE + 1];
    SEXP sym;
    int i, hashcode;

    if (*name == '\0')
        error("attempt to use zero-length variable name");
    if (strlen(name) > MAXIDSIZE)
        error("symbol print-name too long");
    strcpy(buf, name);
    hashcode = R_Newhashpjw(buf);
    i = hashcode % HSIZE;
    /* Check to see if the symbol is already present;  if it is, return it. */
    for (sym = R_SymbolTable[i]; sym != R_NilValue; sym = CDR(sym))
        if (strcmp(buf, CHAR(PRINTNAME(CAR(sym)))) == 0)
            return (CAR(sym));
    /* Create a new symbol node and link it into the table. */
    sym = mkSYMSXP(mkChar(buf), R_UnboundValue);
    SET_HASHVALUE(PRINTNAME(sym), hashcode);
    SET_HASHASH(PRINTNAME(sym), 1);
    R_SymbolTable[i] = CONS(sym, R_SymbolTable[i]);
    return (sym);
}

/*  do_internal - This is the code for .Internal(). */

SEXP do_internal(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP s, fun;
    int save = R_PPStackTop;
    checkArity(op, args);
    s = CAR(args);
    if (!isPairList(s))
        errorcall(call, "invalid .Internal() argument");
    fun = CAR(s);
    if (!isSymbol(fun))
        errorcall(call, "invalid internal function");
    if (INTERNAL(fun) == R_NilValue)
        errorcall(call, "no internal function \"%s\"", CHAR(PRINTNAME(fun)));
    args = CDR(s);
    if (TYPEOF(INTERNAL(fun)) == BUILTINSXP)
        args = evalList(args, env);
    PROTECT(args);
    R_Visible = 1 - PRIMPRINT(INTERNAL(fun));
    args = PRIMFUN(INTERNAL(fun))(s, INTERNAL(fun), args, env);
    UNPROTECT(1);
    if (save != R_PPStackTop)
    {
        REprintf("stack imbalance in internal %s, %d then %d", PRIMNAME(INTERNAL(fun)), save, R_PPStackTop);
    }
    return (args);
}
#undef __R_Names__
