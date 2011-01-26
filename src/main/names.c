/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997--2010  The R Development Core Team
 *  Copyright (C) 2003, 2004  The R Foundation
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

#define __R_Names__ /* used in Defn.h for extern on R_FunTab */
#include <Defn.h>
#include <Print.h>
#include "arithmetic.h" /* for do_math[1234], do_cmathfuns */

#include <R_ext/RConverters.h>

#include <Rinterface.h>

/* Table of  .Internal(.) and .Primitive(.)  R functions
 * =====     =========	      ==========
 *
 * Each entry is a line with
 *
 *  printname	c-entry	 offset	 eval	arity	  pp-kind   precedence	    rightassoc
 *  ---------	-------	 ------	 ----	-----	  -------   ----------	    ----------
 *2 name	cfun	 code	 eval	arity	  gram.kind gram.precedence gram.rightassoc
 *3 PRIMNAME	PRIMFUN	 PRIMVAL [*]    PRIMARITY PPINFO    PPINFO	    PPINFO
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
 * eval:	= XYZ (three digits) --- where e.g. '1' means '001'
 *		X=1 says that we should force R_Visible off
 *		X=0 says that we should force R_Visible on
 *		X=2 says that we should switch R_Visible on but let the C
 *                  code update this.
 *		Y=1 says that this is an internal function which must
 *		    be accessed with a	.Internal(.) call, any other value is
 *		    accessible directly and printed in R as ".Primitive(..)".
 *		Z=1 says evaluate arguments before calling (BUILTINSXP) and
 *		Z=0 says don't evaluate (SPECIALSXP).
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
attribute_hidden FUNTAB R_FunTab[] = {

    /* printname	c-entry		offset	eval	arity	pp-kind	     precedence	rightassoc
     * ---------	-------		------	----	-----	-------      ----------	----------*/

    /* Language Related Constructs */

    /* Primitives */
    {"if", do_if, 0, 200, -1, {PP_IF, PREC_FN, 1}},
    {"while", do_while, 0, 100, -1, {PP_WHILE, PREC_FN, 0}},
    {"for", do_for, 0, 100, -1, {PP_FOR, PREC_FN, 0}},
    {"repeat", do_repeat, 0, 100, -1, {PP_REPEAT, PREC_FN, 0}},
    {"break", do_break, CTXT_BREAK, 0, -1, {PP_BREAK, PREC_FN, 0}},
    {"next", do_break, CTXT_NEXT, 0, -1, {PP_NEXT, PREC_FN, 0}},
    {"return", do_return, 0, 0, -1, {PP_RETURN, PREC_FN, 0}},
    {"function", do_function, 0, 0, -1, {PP_FUNCTION, PREC_FN, 0}},
    {"<-", do_set, 1, 100, -1, {PP_ASSIGN, PREC_LEFT, 1}},
    {"=", do_set, 3, 100, -1, {PP_ASSIGN, PREC_EQ, 1}},
    {"<<-", do_set, 2, 100, -1, {PP_ASSIGN2, PREC_LEFT, 1}},
    {"{", do_begin, 0, 200, -1, {PP_CURLY, PREC_FN, 0}},
    {"(", do_paren, 0, 1, 1, {PP_PAREN, PREC_FN, 0}},
    {".subset", do_subset_dflt, 1, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {".subset2", do_subset2_dflt, 2, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"[", do_subset, 1, 0, -1, {PP_SUBSET, PREC_SUBSET, 0}},
    {"[[", do_subset2, 2, 0, -1, {PP_SUBSET, PREC_SUBSET, 0}},
    {"$", do_subset3, 3, 0, 2, {PP_DOLLAR, PREC_DOLLAR, 0}},
    {"@", do_AT, 0, 0, 2, {PP_DOLLAR, PREC_DOLLAR, 0}},
    {"[<-", do_subassign, 0, 0, 3, {PP_SUBASS, PREC_LEFT, 1}},
    {"[[<-", do_subassign2, 1, 0, 3, {PP_SUBASS, PREC_LEFT, 1}},
    {"$<-", do_subassign3, 1, 0, 3, {PP_SUBASS, PREC_LEFT, 1}},
    {"switch", do_switch, 0, 200, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"browser", do_browser, 0, 101, 4, {PP_FUNCALL, PREC_FN, 0}},
    {".primTrace", do_trace, 0, 101, 1, {PP_FUNCALL, PREC_FN, 0}},
    {".primUntrace", do_trace, 1, 101, 1, {PP_FUNCALL, PREC_FN, 0}},
    {".Internal", do_internal, 0, 200, 1, {PP_FUNCALL, PREC_FN, 0}},
    {".Primitive", do_primitive, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"call", do_call, 0, 0, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"quote", do_quote, 0, 0, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"substitute", do_substitute, 0, 0, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"missing", do_missing, 1, 0, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"nargs", do_nargs, 1, 1, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"on.exit", do_onexit, 0, 100, 1, {PP_FUNCALL, PREC_FN, 0}},

    /* .Internals */

    {"stop", do_stop, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"warning", do_warning, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"gettext", do_gettext, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"ngettext", do_ngettext, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"bindtextdomain", do_bindtextdomain, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {".addCondHands", do_addCondHands, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {".resetCondHands", do_resetCondHands, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {".signalCondition", do_signalCondition, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {".dfltStop", do_dfltStop, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {".dfltWarn", do_dfltWarn, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {".addRestart", do_addRestart, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {".getRestart", do_getRestart, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {".invokeRestart", do_invokeRestart, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {".addTryHandlers", do_addTryHandlers, 0, 111, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"geterrmessage", do_geterrmessage, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"seterrmessage", do_seterrmessage, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"printDeferredWarnings", do_printDeferredWarnings, 0, 111, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"interruptsSuspended", do_interruptsSuspended, 0, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"restart", do_restart, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"as.function.default", do_asfunction, 0, 11, 2, {PP_FUNCTION, PREC_FN, 0}},
    {"debug", do_debug, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"undebug", do_debug, 1, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"isdebugged", do_debug, 2, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"debugonce", do_debug, 3, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"Recall", do_recall, 0, 210, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"delayedAssign", do_delayed, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"makeLazy", do_makelazy, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"identical", do_identical, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},

    /* Binary Operators, all primitives */
    /* these are group generic and so need to eval args */
    {"+", do_arith, PLUSOP, 1, 2, {PP_BINARY, PREC_SUM, 0}},
    {"-", do_arith, MINUSOP, 1, 2, {PP_BINARY, PREC_SUM, 0}},
    {"*", do_arith, TIMESOP, 1, 2, {PP_BINARY, PREC_PROD, 0}},
    {"/", do_arith, DIVOP, 1, 2, {PP_BINARY2, PREC_PROD, 0}},
    {"^", do_arith, POWOP, 1, 2, {PP_BINARY2, PREC_POWER, 1}},
    {"%%", do_arith, MODOP, 1, 2, {PP_BINARY2, PREC_PERCENT, 0}},
    {"%/%", do_arith, IDIVOP, 1, 2, {PP_BINARY2, PREC_PERCENT, 0}},
    {"%*%", do_matprod, 0, 1, 2, {PP_BINARY, PREC_PERCENT, 0}},

    {"==", do_relop, EQOP, 1, 2, {PP_BINARY, PREC_COMPARE, 0}},
    {"!=", do_relop, NEOP, 1, 2, {PP_BINARY, PREC_COMPARE, 0}},
    {"<", do_relop, LTOP, 1, 2, {PP_BINARY, PREC_COMPARE, 0}},
    {"<=", do_relop, LEOP, 1, 2, {PP_BINARY, PREC_COMPARE, 0}},
    {">=", do_relop, GEOP, 1, 2, {PP_BINARY, PREC_COMPARE, 0}},
    {">", do_relop, GTOP, 1, 2, {PP_BINARY, PREC_COMPARE, 0}},
    {"&", do_logic, 1, 1, 2, {PP_BINARY, PREC_AND, 0}},
    {"|", do_logic, 2, 1, 2, {PP_BINARY, PREC_OR, 0}},
    {"!", do_logic, 3, 1, 1, {PP_UNARY, PREC_NOT, 0}},

    /* specials as conditionally evaluate second arg */
    {"&&", do_logic2, 1, 0, 2, {PP_BINARY, PREC_AND, 0}},
    {"||", do_logic2, 2, 0, 2, {PP_BINARY, PREC_OR, 0}},
    {":", do_colon, 0, 1, 2, {PP_BINARY2, PREC_COLON, 0}},
    /* does not evaluate */
    {"~", do_tilde, 0, 0, 2, {PP_BINARY, PREC_TILDE, 0}},

    /* Logic Related Functions */
    /* these are group generic and so need to eval args */
    {"all", do_logic3, 1, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"any", do_logic3, 2, 1, -1, {PP_FUNCALL, PREC_FN, 0}},

    /* Vectors, Matrices and Arrays */

    /* Primitives */

    {"length", do_length, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"length<-", do_lengthgets, 0, 1, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"c", /* bind.c:*/ do_c, 0, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"oldClass", do_class, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"oldClass<-", do_classgets, 0, 1, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"class", R_do_data_class, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {".cache_class", R_do_data_class, 1, 1, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"class<-", R_do_set_class, 0, 1, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"unclass", do_unclass, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"names", do_names, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"names<-", do_namesgets, 0, 1, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"dimnames", do_dimnames, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dimnames<-", do_dimnamesgets, 0, 1, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"dim", do_dim, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dim<-", do_dimgets, 0, 1, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"attributes", do_attributes, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"attributes<-", do_attributesgets, 0, 1, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"attr", do_attr, 0, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"attr<-", do_attrgets, 0, 1, 3, {PP_FUNCALL, PREC_LEFT, 1}},
    {"levels<-", do_levelsgets, 0, 1, 2, {PP_FUNCALL, PREC_LEFT, 1}},

    /* .Internals */

    {"vector", do_makevector, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"complex", do_complex, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"matrix", do_matrix, 0, 11, 7, {PP_FUNCALL, PREC_FN, 0}},
    {"row", do_rowscols, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"col", do_rowscols, 2, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"unlist", do_unlist, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"cbind", do_bind, 1, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"rbind", do_bind, 2, 10, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"drop", do_drop, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"all.names", do_allnames, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"comment", do_comment, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"comment<-", do_commentgets, 0, 11, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"get", do_get, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"mget", do_mget, 1, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"exists", do_get, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"assign", do_assign, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"list2env", do_list2env, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"remove", do_remove, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"duplicated", do_duplicated, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"unique", do_duplicated, 1, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"anyDuplicated", do_duplicated, 2, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"which", do_which, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"which.min", do_first_min, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pmin", do_pmin, 0, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"pmax", do_pmin, 1, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"which.max", do_first_min, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"match", do_match, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"pmatch", do_pmatch, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"charmatch", do_charmatch, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"match.call", do_matchcall, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"complete.cases", do_compcases, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"crossprod", do_matprod, 1, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"tcrossprod", do_matprod, 2, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    {"attach", do_attach, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"detach", do_detach, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"search", do_search, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},

    /* Mathematical Functions */
    /* primitives: these are group generic and so need to eval args (possibly internally) */
    {"round", do_Math2, 10001, 0, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"signif", do_Math2, 10004, 0, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"log", do_log, 10003, 0, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"log10", do_log1arg, 10, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"log2", do_log1arg, 2, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"abs", do_abs, 6, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"floor", do_math1, 1, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"ceiling", do_math1, 2, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"sqrt", do_math1, 3, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"sign", do_math1, 4, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"trunc", do_trunc, 5, 1, -1, {PP_FUNCALL, PREC_FN, 0}},

    {"exp", do_math1, 10, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"expm1", do_math1, 11, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"log1p", do_math1, 12, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"cos", do_math1, 20, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"sin", do_math1, 21, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"tan", do_math1, 22, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"acos", do_math1, 23, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"asin", do_math1, 24, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"atan", do_math1, 25, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"cosh", do_math1, 30, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"sinh", do_math1, 31, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"tanh", do_math1, 32, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"acosh", do_math1, 33, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"asinh", do_math1, 34, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"atanh", do_math1, 35, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"lgamma", do_math1, 40, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"gamma", do_math1, 41, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"digamma", do_math1, 42, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"trigamma", do_math1, 43, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    /* see "psigamma" below !*/

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

    {"psigamma", do_math2, 26, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    /* Mathematical Functions of a Complex Argument */
    /* these are group generic and so need to eval args */

    {"Re", do_cmathfuns, 1, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Im", do_cmathfuns, 2, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Mod", do_cmathfuns, 3, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Arg", do_cmathfuns, 4, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Conj", do_cmathfuns, 5, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

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

    {"dnbinom_mu", do_math3, 45, 11, 3 + 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pnbinom_mu", do_math3, 46, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qnbinom_mu", do_math3, 47, 11, 3 + 2, {PP_FUNCALL, PREC_FN, 0}},

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
    {"rnbinom_mu", do_random2, 13, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rnchisq", do_random2, 12, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rnorm", do_random2, 8, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"runif", do_random2, 9, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rweibull", do_random2, 10, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rwilcox", do_random2, 11, 11, 3, {PP_FUNCALL, PREC_FN, 0}},

    {"rhyper", do_random3, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},

    {"rmultinom", do_rmultinom, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"sample", do_sample, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},

    {"RNGkind", do_RNGkind, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"set.seed", do_setseed, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},

    /* Data Summaries */
    /* these four are group generic and so need to eval args */
    {"sum", do_summary, 0, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"min", do_summary, 2, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"max", do_summary, 3, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"prod", do_summary, 4, 1, -1, {PP_FUNCALL, PREC_FN, 0}},

    {"mean", do_summary, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"range", do_range, 0, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"cov", do_cov, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"cor", do_cov, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},

    /* Note that the number of arguments in this group only applies
       to the default method */
    {"cumsum", do_cum, 1, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"cumprod", do_cum, 2, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"cummax", do_cum, 3, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"cummin", do_cum, 4, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    /* Type coercion */

    {"as.character", do_ascharacter, 0, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"as.integer", do_ascharacter, 1, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"as.double", do_ascharacter, 2, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"as.complex", do_ascharacter, 3, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"as.logical", do_ascharacter, 4, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"as.raw", do_ascharacter, 5, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"as.call", do_ascall, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"as.environment", do_as_environment, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"storage.mode<-", do_storage_mode, 0, 1, 2, {PP_FUNCALL, PREC_FN, 0}},

    {"as.vector", do_asvector, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"paste", do_paste, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"file.path", do_filepath, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"format", do_format, 0, 11, 8, {PP_FUNCALL, PREC_FN, 0}},
    {"format.info", do_formatinfo, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"cat", do_cat, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"do.call", do_docall, 0, 211, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"type.convert", do_typecvt, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},

    /* String Manipulation */

    {"nchar", do_nchar, 1, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"nzchar", do_nzchar, 1, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"substr", do_substr, 1, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"substr<-", do_substrgets, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"strsplit", do_strsplit, 1, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"abbreviate", do_abbrev, 1, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"make.names", do_makenames, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"grep", do_grep, 0, 11, 8, {PP_FUNCALL, PREC_FN, 0}},
    {"grepl", do_grep, 1, 11, 8, {PP_FUNCALL, PREC_FN, 0}},
    {"sub", do_gsub, 0, 11, 7, {PP_FUNCALL, PREC_FN, 0}},
    {"gsub", do_gsub, 1, 11, 7, {PP_FUNCALL, PREC_FN, 0}},
    {"regexpr", do_regexpr, 0, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"gregexpr", do_regexpr, 1, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"agrep", do_agrep, 1, 11, 9, {PP_FUNCALL, PREC_FN, 0}},
    {"tolower", do_tolower, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"toupper", do_tolower, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"chartr", do_chartr, 1, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"sprintf", do_sprintf, 1, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"make.unique", do_makeunique, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"charToRaw", do_charToRaw, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"rawToChar", do_rawToChar, 1, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"rawShift", do_rawShift, 1, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"intToBits", do_intToBits, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"rawToBits", do_rawToBits, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"packBits", do_packBits, 1, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"utf8ToInt", do_utf8ToInt, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"intToUtf8", do_intToUtf8, 1, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"encodeString", do_encodeString, 1, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"iconv", do_iconv, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"strtrim", do_strtrim, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"strtoi", do_strtoi, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

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
    {"is.raw", do_is, RAWSXP, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

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

    {"is.na", do_isna, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.nan", do_isnan, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.finite", do_isfinite, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"is.infinite", do_isinfinite, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"is.vector", do_isvector, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    /* Miscellaneous */

    /* Primitive */
    {"proc.time", do_proctime, 0, 1, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"gc.time", do_gctime, 0, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
#if 0
{"visibleflag", do_visibleflag,	0,	1,	0,	{PP_FUNCALL, PREC_FN,	0}},
#endif
    {"withVisible", do_withVisible, 1, 10, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"expression", do_expression, 1, 0, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"interactive", do_interactive, 0, 1, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"invisible", do_invisible, 0, 101, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"rep", do_rep, 0, 0, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"rep.int", do_rep_int, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"seq.int", do_seq, 0, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"seq_len", do_seq_len, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"seq_along", do_seq_along, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"list", do_makelist, 1, 1, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"xtfrm", do_xtfrm, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"enc2native", do_enc2, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"enc2utf8", do_enc2, 1, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"emptyenv", do_emptyenv, 0, 1, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"baseenv", do_baseenv, 0, 1, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"globalenv", do_globalenv, 0, 1, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"environment<-", do_envirgets, 0, 1, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"pos.to.env", do_pos2env, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"eapply", do_eapply, 0, 10, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"lapply", do_lapply, 0, 10, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"vapply", do_vapply, 0, 10, 4, {PP_FUNCALL, PREC_FN, 0}},

    {".C", do_dotCode, 0, 1, -1, {PP_FOREIGN, PREC_FN, 0}},
    {".Fortran", do_dotCode, 1, 1, -1, {PP_FOREIGN, PREC_FN, 0}},
    {".External", do_External, 0, 1, -1, {PP_FOREIGN, PREC_FN, 0}},
    {".Call", do_dotcall, 0, 1, -1, {PP_FOREIGN, PREC_FN, 0}},
    {".External.graphics", do_Externalgr, 0, 1, -1, {PP_FOREIGN, PREC_FN, 0}},
    {".Call.graphics", do_dotcallgr, 0, 1, -1, {PP_FOREIGN, PREC_FN, 0}},

    /* .Internal */
    {"Version", do_version, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"machine", do_machine, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"commandArgs", do_commandArgs, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"unzip", do_unzip, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
#ifdef Win32
    {"system", do_system, 0, 211, 5, {PP_FUNCALL, PREC_FN, 0}},
#else
    {"system", do_system, 0, 211, 2, {PP_FUNCALL, PREC_FN, 0}},
#endif
#ifdef Win32
    {"flush.console", do_flushconsole, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"win.version", do_winver, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"shell.exec", do_shellexec, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"winDialog", do_windialog, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"winDialogString", do_windialogstring, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"winMenuNames", do_winmenunames, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"winMenuItems", do_wingetmenuitems, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"winMenuAdd", do_winmenuadd, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"winMenuDel", do_winmenudel, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"memory.size", do_memsize, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"DLL.version", do_dllversion, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"bringToTop", do_bringtotop, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"msgWindow", do_msgwindow, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"select.list", do_selectlist, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"getClipboardFormats", do_getClipboardFormats, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"readClipboard", do_readClipboard, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"writeClipboard", do_writeClipboard, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"chooseFiles", do_chooseFiles, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"chooseDir", do_chooseDir, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"getIdentification", do_getIdentification, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"getWindowHandle", do_getWindowHandle, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"getWindowHandles", do_getWindowHandles, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"getWindowTitle", do_getWindowTitle, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"setWindowTitle", do_setTitle, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"arrangeWindows", do_arrangeWindows, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"setStatusBar", do_setStatusBar, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"shortPathName", do_shortpath, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"loadRconsole", do_loadRconsole, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.which", do_syswhich, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"readRegistry", do_readRegistry, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"winProgressBar", do_winprogressbar, 0, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"closeWinProgressBar", do_closewinprogressbar, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"setWinProgressBar", do_setwinprogressbar, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"useInternet2", do_setInternet2, 0, 211, 1, {PP_FUNCALL, PREC_FN, 0}},
#endif
#if defined(__APPLE_CC__) && defined(HAVE_AQUA)
    {"wsbrowser", do_wsbrowser, 0, 11, 8, {PP_FUNCALL, PREC_FN, 0}},
    {"pkgbrowser", do_browsepkgs, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"data.manager", do_datamanger, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"package.manager", do_packagemanger, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"flush.console", do_flushconsole, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"hsbrowser", do_hsbrowser, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"select.list", do_selectlist, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"aqua.custom.print", do_aqua_custom_print, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
#endif
    {"parse", do_parse, 0, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"parse_Rd", do_parseRd, 0, 11, 7, {PP_FUNCALL, PREC_FN, 0}},
    {"parseLatex", do_parseLatex, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"save", do_save, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"saveToConn", do_saveToConn, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"load", do_load, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"loadFromConn2", do_loadFromConn2, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"serializeToConn", do_serializeToConn, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"unserializeFromConn", do_unserializeFromConn, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"deparse", do_deparse, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"deparseRd", do_deparseRd, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"dput", do_dput, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"dump", do_dump, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"quit", do_quit, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"readline", do_readln, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"menu", do_menu, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"print.default", do_printdefault, 0, 111, 9, {PP_FUNCALL, PREC_FN, 0}},
    {"print.function", do_printfunction, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"prmatrix", do_prmatrix, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"gc", do_gc, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"gcinfo", do_gcinfo, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"gctorture", do_gctorture, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"gctorture2", do_gctorture2, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"memory.profile", do_memoryprofile, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"split", do_split, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"is.loaded", do_isloaded, 0, 11, -1, {PP_FOREIGN, PREC_FN, 0}},
    {"recordGraphics", do_recordGraphics, 0, 211, 3, {PP_FOREIGN, PREC_FN, 0}},
    {"dyn.load", do_dynload, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"dyn.unload", do_dynunload, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"ls", do_ls, 1, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"typeof", do_typeof, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"eval", do_eval, 0, 211, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"eval.with.vis", do_eval, 1, 211, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.parent", do_sys, 1, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.call", do_sys, 2, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.frame", do_sys, 3, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.nframe", do_sys, 4, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.calls", do_sys, 5, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.frames", do_sys, 6, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.on.exit", do_sys, 7, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.parents", do_sys, 8, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sys.function", do_sys, 9, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"browserText", do_sysbrowser, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"browserCondition", do_sysbrowser, 2, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"browserSetDebug", do_sysbrowser, 3, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"parent.frame", do_parentframe, 0, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"sort", do_sort, 1, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"is.unsorted", do_isunsorted, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"psort", do_psort, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"qsort", do_qsort, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"radixsort", do_radixsort, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"order", do_order, 0, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"rank", do_rank, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"scan", do_scan, 0, 11, 18, {PP_FUNCALL, PREC_FN, 0}},
    {"count.fields", do_countfields, 0, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"readTableHead", do_readtablehead, 0, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"t.default", do_transpose, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"aperm", do_aperm, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"builtins", do_builtins, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"edit", do_edit, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"dataentry", do_dataentry, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"dataviewer", do_dataviewer, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"args", do_args, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"formals", do_formals, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"body", do_body, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"bodyCode", do_bodyCode, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"environment", do_envir, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"environmentName", do_envirName, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"env2list", do_env2list, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"reg.finalizer", do_regFinaliz, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"options", do_options, 0, 211, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"sink", do_sink, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"sink.number", do_sinknumber, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"lib.fixup", do_libfixup, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"rapply", do_rapply, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"islistfactor", do_islistfactor, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"colSums", do_colsum, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"colMeans", do_colsum, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"rowSums", do_colsum, 2, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"rowMeans", do_colsum, 3, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"Rprof", do_Rprof, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"Rprofmem", do_Rprofmem, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"tracemem", do_tracemem, 0, 1, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"retracemem", do_retracemem, 0, 201, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"untracemem", do_untracemem, 0, 101, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"object.size", do_objectsize, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"inspect", do_inspect, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"mem.limits", do_memlimits, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"merge", do_merge, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"capabilities", do_capabilities, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"capabilitiesX11", do_capabilitiesX11, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"new.env", do_newenv, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"parent.env", do_parentenv, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"parent.env<-", do_parentenvgets, 0, 11, 2, {PP_FUNCALL, PREC_LEFT, 1}},
    {"l10n_info", do_l10n_info, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"Cstack_info", do_Cstack_info, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"startHTTPD", do_startHTTPD, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"stopHTTPD", do_stopHTTPD, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},

    /* Functions To Interact with the Operating System */

    {"file.show", do_fileshow, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"file.edit", do_fileedit, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"file.create", do_filecreate, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"file.remove", do_fileremove, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"file.rename", do_filerename, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"file.append", do_fileappend, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"codeFiles.append", do_fileappend, 1, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"file.symlink", do_filesymlink, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"file.link", do_filelink, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"file.copy", do_filecopy, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"list.files", do_listfiles, 0, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"file.exists", do_fileexists, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"file.choose", do_filechoose, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"file.info", do_fileinfo, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"file.access", do_fileaccess, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"dir.create", do_dircreate, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"tempfile", do_tempfile, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"tempdir", do_tempdir, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"R.home", do_Rhome, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"date", do_date, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.getenv", do_getenv, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.setenv", do_setenv, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.unsetenv", do_unsetenv, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"getwd", do_getwd, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"setwd", do_setwd, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"basename", do_basename, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dirname", do_dirname, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dirchmod", do_dirchmod, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.chmod", do_syschmod, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.umask", do_sysumask, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.readlink", do_readlink, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.info", do_sysinfo, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.sleep", do_syssleep, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.getlocale", do_getlocale, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.setlocale", do_setlocale, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.localeconv", do_localeconv, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"path.expand", do_pathexpand, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.getpid", do_sysgetpid, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"normalizePath", do_normalizepath, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"Sys.glob", do_glob, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"unlink", do_unlink, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},

    /* Complex Valued Functions */
    {"fft", do_fft, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"mvfft", do_mvfft, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"nextn", do_nextn, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"polyroot", do_polyroot, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},

/* Device Drivers */

#ifdef Unix
    {"X11", do_X11, 0, 111, 16, {PP_FUNCALL, PREC_FN, 0}},
    {"savePlot", do_saveplot, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"cairo", do_cairo, 0, 111, 9, {PP_FUNCALL, PREC_FN, 0}},
#endif

    /* Graphics */

    {"dev.control", do_devcontrol, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dev.displaylist", do_devcontrol, 1, 111, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"dev.copy", do_devcopy, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dev.cur", do_devcur, 0, 111, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"dev.next", do_devnext, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dev.off", do_devoff, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dev.prev", do_devprev, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dev.set", do_devset, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"rgb", do_rgb, 0, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"rgb256", do_rgb, 1, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"rgb2hsv", do_RGB2hsv, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"hsv", do_hsv, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"hcl", do_hcl, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"gray", do_gray, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"colors", do_colors, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"col2rgb", do_col2RGB, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"palette", do_palette, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"plot.new", do_plot_new, 0, 111, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"plot.window", do_plot_window, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"axis", do_axis, 0, 111, 13, {PP_FUNCALL, PREC_FN, 0}},
    {"plot.xy", do_plot_xy, 0, 111, 7, {PP_FUNCALL, PREC_FN, 0}},
    {"text", do_text, 0, 111, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"mtext", do_mtext, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"title", do_title, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"abline", do_abline, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"box", do_box, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rect", do_rect, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"path", do_path, 0, 111, 7, {PP_FUNCALL, PREC_FN, 0}},
    {"raster", do_raster, 0, 111, 7, {PP_FUNCALL, PREC_FN, 0}},
    {"polygon", do_polygon, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"xspline", do_xspline, 0, 111, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"par", do_par, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"segments", do_segments, 0, 111, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"arrows", do_arrows, 0, 111, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"layout", do_layout, 0, 111, 10, {PP_FUNCALL, PREC_FN, 0}},
    {"locator", do_locator, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"identify", do_identify, 0, 211, 8, {PP_FUNCALL, PREC_FN, 0}},
    {"strheight", do_strheight, 0, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"strwidth", do_strwidth, 0, 11, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"contour", do_contour, 0, 11, 12, {PP_FUNCALL, PREC_FN, 0}},
    {"contourLines", do_contourLines, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"image", do_image, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"dend", do_dend, 0, 111, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"dend.window", do_dendwindow, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"erase", do_erase, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"persp", do_persp, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"filledcontour", do_filledcontour, 0, 111, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"getSnapshot", do_getSnapshot, 0, 111, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"playSnapshot", do_playSnapshot, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"symbols", do_symbols, 0, 111, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"getGraphicsEvent", do_getGraphicsEvent, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"getGraphicsEventEnv", do_getGraphicsEventEnv, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"setGraphicsEventEnv", do_setGraphicsEventEnv, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"devAskNewPage", do_devAskNewPage, 0, 211, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"dev.size", do_devsize, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"clip", do_clip, 0, 111, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"grconvertX", do_convertXY, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"grconvertY", do_convertXY, 1, 11, 3, {PP_FUNCALL, PREC_FN, 0}},

    /* Objects */
    {"inherits", do_inherits, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"UseMethod", do_usemethod, 0, 200, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"NextMethod", do_nextmethod, 0, 210, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"standardGeneric", do_standardGeneric, 0, 201, -1, {PP_FUNCALL, PREC_FN, 0}},

    /* Modelling Functionality */

    {"nlm", do_nlm, 0, 11, 11, {PP_FUNCALL, PREC_FN, 0}},
    {"fmin", do_fmin, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"zeroin", do_zeroin, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"zeroin2", do_zeroin2, 0, 11, 7, {PP_FUNCALL, PREC_FN, 0}},
    {"optim", do_optim, 0, 11, 7, {PP_FUNCALL, PREC_FN, 0}},
    {"optimhess", do_optimhess, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"terms.formula", do_termsform, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"update.formula", do_updateform, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"model.frame", do_modelframe, 0, 11, 8, {PP_FUNCALL, PREC_FN, 0}},
    {"model.matrix", do_modelmatrix, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    {"D", do_D, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"deriv.default", do_deriv, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},

    /* History manipulation */
    {"loadhistory", do_loadhistory, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"savehistory", do_savehistory, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"addhistory", do_addhistory, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},

    /* date-time manipulations */
    {"Sys.time", do_systime, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"as.POSIXct", do_asPOSIXct, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"as.POSIXlt", do_asPOSIXlt, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"format.POSIXlt", do_formatPOSIXlt, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"strptime", do_strptime, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"Date2POSIXlt", do_D2POSIXlt, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"POSIXlt2Date", do_POSIXlt2D, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},

#ifdef BYTECODE
    {"mkCode", do_mkcode, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"bcClose", do_bcclose, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"is.builtin.internal", do_is_builtin_internal, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"disassemble", do_disassemble, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"bcVersion", do_bcversion, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"load.from.file", do_loadfile, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"save.to.file", do_savefile, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"putconst", do_putconst, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"enableJIT", do_enablejit, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
#endif

    {"setNumMathThreads", do_setnumthreads, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"setMaxNumMathThreads", do_setmaxnumthreads, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},

    /* Connections */
    {"stdin", do_stdin, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"stdout", do_stdout, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"stderr", do_stderr, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"isatty", do_isatty, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"readLines", do_readLines, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"writeLines", do_writelines, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"readBin", do_readbin, 0, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"writeBin", do_writebin, 0, 211, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"readChar", do_readchar, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"writeChar", do_writechar, 0, 211, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"open", do_open, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"isOpen", do_isopen, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"isIncomplete", do_isincomplete, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"isSeekable", do_isseekable, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"close", do_close, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"flush", do_flush, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"file", do_url, 1, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"url", do_url, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"pipe", do_pipe, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"fifo", do_fifo, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"gzfile", do_gzfile, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"bzfile", do_gzfile, 1, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"xzfile", do_gzfile, 2, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"unz", do_unz, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"seek", do_seek, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"truncate", do_truncate, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pushBack", do_pushback, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"clearPushBack", do_clearpushback, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"pushBackLength", do_pushbacklength, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"rawConnection", do_rawconnection, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"rawConnectionValue", do_rawconvalue, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"textConnection", do_textconnection, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"textConnectionValue", do_textconvalue, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"socketConnection", do_sockconn, 0, 11, 6, {PP_FUNCALL, PREC_FN, 0}},
    {"sockSelect", do_sockselect, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"getConnection", do_getconnection, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"getAllConnections", do_getallconnections, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"summary.connection", do_sumconnection, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"download", do_download, 0, 11, 5, {PP_FUNCALL, PREC_FN, 0}},
    {"nsl", do_nsl, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"gzcon", do_gzcon, 0, 11, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"memCompress", do_memCompress, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"memDecompress", do_memDecompress, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    {"readDCF", do_readDCF, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},

    {"getNumRtoCConverters", do_getNumRtoCConverters, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"getRtoCConverterDescriptions", do_getRtoCConverterDescriptions, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"getRtoCConverterStatus", do_getRtoCConverterStatus, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"setToCConverterActiveStatus", do_setToCConverterActiveStatus, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"removeToCConverterActiveStatus", do_setToCConverterActiveStatus, 1, 11, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"lockEnvironment", do_lockEnv, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"environmentIsLocked", do_envIsLocked, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"lockBinding", do_lockBnd, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"unlockBinding", do_lockBnd, 1, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"bindingIsLocked", do_bndIsLocked, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"makeActiveBinding", do_mkActiveBnd, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"bindingIsActive", do_bndIsActive, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    /* looks like mkUnbound is unused in base R */
    {"mkUnbound", do_mkUnbound, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"isNamespaceEnv", do_isNSEnv, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"registerNamespace", do_regNS, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"unregisterNamespace", do_unregNS, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"getRegisteredNamespace", do_getRegNS, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"getNamespaceRegistry", do_getNSRegistry, 0, 11, 0, {PP_FUNCALL, PREC_FN, 0}},
    {"importIntoEnv", do_importIntoEnv, 0, 11, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"env.profile", do_envprofile, 0, 211, 1, {PP_FUNCALL, PREC_FN, 0}},

    {"write.table", do_writetable, 0, 111, 11, {PP_FUNCALL, PREC_FN, 0}},
    {"Encoding", do_encoding, 0, 11, 1, {PP_FUNCALL, PREC_FN, 0}},
    {"setEncoding", do_setencoding, 0, 11, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"lazyLoadDBfetch", do_lazyLoadDBfetch, 0, 1, 4, {PP_FUNCALL, PREC_FN, 0}},
    {"setTimeLimit", do_setTimeLimit, 0, 111, 3, {PP_FUNCALL, PREC_FN, 0}},
    {"setSessionTimeLimit", do_setSessionTimeLimit, 0, 111, 2, {PP_FUNCALL, PREC_FN, 0}},
    {"icuSetCollate", do_ICUset, 0, 111, -1, {PP_FUNCALL, PREC_FN, 0}},
    {"readRenviron", do_readEnviron, 0, 111, 1, {PP_FUNCALL, PREC_FN, 0}},
    {NULL, NULL, 0, 0, 0, {PP_INVALID, PREC_FN, 0}},
};

SEXP attribute_hidden do_primitive(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP name;
    int i;
    checkArity(op, args);
    name = CAR(args);
    if (!isString(name) || length(name) < 1 || STRING_ELT(name, 0) == R_NilValue)
        errorcall(call, _("string argument required"));
    for (i = 0; R_FunTab[i].name; i++) /* all names are ASCII */
        if (strcmp(CHAR(STRING_ELT(name, 0)), R_FunTab[i].name) == 0)
        {
            if ((R_FunTab[i].eval % 100) / 10)
                return mkPRIMSXP(i, R_FunTab[i].eval % 10);
            else
                return mkPRIMSXP(i, R_FunTab[i].eval % 10);
        }
    errorcall(call, _("no such primitive function"));
    return (R_NilValue); /* -Wall */
}

int StrToInternal(const char *s)
{
    int i;
    for (i = 0; R_FunTab[i].name; i++)
        if (strcmp(s, R_FunTab[i].name) == 0)
            return i;
    return 0;
}

static void installFunTab(int i)
{
    SEXP prim;
    /* prim needs to be protected since install can (and does here) allocate */
    PROTECT(prim = mkPRIMSXP(i, R_FunTab[i].eval % 10));
    if ((R_FunTab[i].eval % 100) / 10)
        SET_INTERNAL(install(R_FunTab[i].name), prim);
    else
        SET_SYMVALUE(install(R_FunTab[i].name), prim);
    UNPROTECT(1);
}

static void SymbolShortcuts(void)
{ /* ../include/Rinternals.h : */
    R_Bracket2Symbol = install("[[");
    R_BracketSymbol = install("[");
    R_BraceSymbol = install("{");
    R_ClassSymbol = install("class");
    R_DeviceSymbol = install(".Device");
    R_DimNamesSymbol = install("dimnames");
    R_DimSymbol = install("dim");
    R_DollarSymbol = install("$");
    R_DotsSymbol = install("...");
    R_DropSymbol = install("drop");
    R_LastvalueSymbol = install(".Last.value");
    R_LevelsSymbol = install("levels");
    R_ModeSymbol = install("mode");
    R_NameSymbol = install("name");
    R_NamesSymbol = install("names");
    R_NaRmSymbol = install("na.rm");
    R_PackageSymbol = install("package");
    R_QuoteSymbol = install("quote");
    R_RowNamesSymbol = install("row.names");
    R_SeedsSymbol = install(".Random.seed");
    R_SourceSymbol = install("source");
    R_TspSymbol = install("tsp");
    /* ../include/Defn.h , i.e. non-public : */
    R_CommentSymbol = install("comment");
    R_DotEnvSymbol = install(".Environment");
    R_ExactSymbol = install("exact");
    R_RecursiveSymbol = install("recursive");
    R_SrcfileSymbol = install("srcfile");
    R_SrcrefSymbol = install("srcref");
    R_WholeSrcrefSymbol = install("wholeSrcref");
    R_TmpvalSymbol = install("*tmp*");
    R_UseNamesSymbol = install("use.names");
    R_DoubleColonSymbol = install("::");
    R_TripleColonSymbol = install(":::");
    R_ConnIdSymbol = install("conn_id");
    R_DevicesSymbol = install(".Devices");

    R_dot_Generic = install(".Generic");
    R_dot_Method = install(".Method");
    R_dot_Methods = install(".Methods");
    R_dot_defined = install(".defined");
    R_dot_target = install(".target");
    R_dot_Group = install(".Group");
    R_dot_Class = install(".Class");
    R_dot_GenericCallEnv = install(".GenericCallEnv");
    R_dot_GenericDefEnv = install(".GenericDefEnv");
}

extern SEXP framenames; /* from model.c */

/* initialize the symbol table */
void InitNames()
{
    int i;

    /* allocate the symbol table */
    if (!(R_SymbolTable = (SEXP *)calloc(HSIZE, sizeof(SEXP))))
        R_Suicide("couldn't allocate memory for symbol table");

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
    /* String constants (CHARSXP values) */
    /* Note: we don't want NA_STRING to be in the CHARSXP cache, so that
       mkChar("NA") is distinct from NA_STRING */
    /* NA_STRING */
    NA_STRING = allocCharsxp(strlen("NA"));
    strcpy(CHAR_RW(NA_STRING), "NA");
    SET_CACHED(NA_STRING); /* Mark it */
    R_print.na_string = NA_STRING;
    /* R_BlankString */
    R_BlankString = mkChar("");
    /* Initialize the symbol Table */
    for (i = 0; i < HSIZE; i++)
        R_SymbolTable[i] = R_NilValue;
    /* Set up a set of globals so that a symbol table search can be
       avoided when matching something like dim or dimnames. */
    SymbolShortcuts();
    /*  Builtin Functions */
    for (i = 0; R_FunTab[i].name; i++)
        installFunTab(i);
    framenames = R_NilValue;
#ifdef BYTECODE
    R_initialize_bcode();
#endif
}

/*  install - probe the symbol table */
/*  If "name" is not found, it is installed in the symbol table.
    The symbol corresponding to the string "name" is returned. */

SEXP install(const char *name)
{
    SEXP sym;
    int i, hashcode;

    if (*name == '\0')
        error(_("attempt to use zero-length variable name"));
    if (strlen(name) > MAXIDSIZE)
        error(_("variable names are limited to %d bytes"), MAXIDSIZE);
    hashcode = R_Newhashpjw(name);
    i = hashcode % HSIZE;
    /* Check to see if the symbol is already present;  if it is, return it. */
    for (sym = R_SymbolTable[i]; sym != R_NilValue; sym = CDR(sym))
        if (strcmp(name, CHAR(PRINTNAME(CAR(sym)))) == 0)
            return (CAR(sym));
    /* Create a new symbol node and link it into the table. */
    sym = mkSYMSXP(mkChar(name), R_UnboundValue);
    SET_HASHVALUE(PRINTNAME(sym), hashcode);
    SET_HASHASH(PRINTNAME(sym), 1);
    R_SymbolTable[i] = CONS(sym, R_SymbolTable[i]);
    return (sym);
}

/*  do_internal - This is the code for .Internal(). */

SEXP attribute_hidden do_internal(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP s, fun, ans;
    int save = R_PPStackTop;
    int flag;
    const void *vmax = vmaxget();

    checkArity(op, args);
    s = CAR(args);
    if (!isPairList(s))
        errorcall(call, _("invalid .Internal() argument"));
    fun = CAR(s);
    if (!isSymbol(fun))
        errorcall(call, _("invalid internal function"));
    if (INTERNAL(fun) == R_NilValue)
        errorcall(call, _("no internal function \"%s\""), CHAR(PRINTNAME(fun)));
    args = CDR(s);
    if (TYPEOF(INTERNAL(fun)) == BUILTINSXP)
        args = evalList(args, env, call, 0);
    PROTECT(args);
    flag = PRIMPRINT(INTERNAL(fun));
    R_Visible = flag != 1;
    ans = PRIMFUN(INTERNAL(fun))(s, INTERNAL(fun), args, env);
    /* This resetting of R_Visible=FALSE  was to fix PR#7397,
       now fixed in GEText */
    if (flag < 2)
        R_Visible = flag != 1;
#ifdef CHECK_VISIBILITY
    if (flag < 2 && flag == R_Visible)
    {
        char *nm = CHAR(PRINTNAME(fun));
        if (strcmp(nm, "eval") && strcmp(nm, "options") && strcmp(nm, "Recall") && strcmp(nm, "do.call") &&
            strcmp(nm, "switch") && strcmp(nm, "recordGraphics") && strcmp(nm, "writeBin") &&
            strcmp(nm, "NextMethod") && strcmp(nm, "eval.with.vis"))
            printf("vis: internal %s\n", nm);
    }
#endif
    UNPROTECT(1);
    check_stack_balance(INTERNAL(fun), save);
    vmaxset(vmax);
    return (ans);
}
#undef __R_Names__
