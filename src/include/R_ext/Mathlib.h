/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2000  The R Development Core Team
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
 *

 * Mathlib.h  should contain ALL headers from R's C code in `src/nmath'
   ---------  such that ``the Math library'' can be used by simply

   ``#include "Mathlib.h" ''

   and nothing else.
*/
#ifndef MATHLIB_H
#define MATHLIB_H

/*-- Mathlib as part of R --  undefine this for standalone : */
#define MATHLIB_IN_R

#include "Rconfig.h"
#include "R_ext/Arith.h"

#ifdef FORTRAN_H
#error __MUST__include "Mathlib.h"  _before_  "Fortran.h"
#endif

#include <errno.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

/* TRUE and FALSE conflict with the Mac --- Fortran.h still defines them... */
#define LTRUE	(1)
#define LFALSE	(0)

/* 30 Decimal-place constants */
/* Computed with bc -l (scale=32; proper round) */

/* SVID & X/Open Constants */
/* Names from Solaris math.h */

#ifndef M_E
#define M_E		2.718281828459045235360287471353	/* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E		1.442695040888963407359924681002	/* log2(e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E	0.434294481903251827651128918917	/* log10(e) */
#endif

#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458	/* ln(2) */
#endif

#ifndef M_LN10
#define M_LN10		2.302585092994045684017991454684	/* ln(10) */
#endif

#ifndef M_PI
#define M_PI		3.141592653589793238462643383280	/* pi */
#endif

#ifndef M_PI_2
#define M_PI_2		1.570796326794896619231321691640	/* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4		0.785398163397448309615660845820	/* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI		0.318309886183790671537767526745	/* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI		0.636619772367581343075535053490	/* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI	1.128379167095512573896158903122	/* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2		1.414213562373095048801688724210	/* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2	0.707106781186547524400844362105	/* 1/sqrt(2) */
#endif

/* Other, R-Specific Constants */

#ifndef M_SQRT_32
#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
#endif

#ifndef M_LOG10_2
#define M_LOG10_2	0.301029995663981195213738894724	/* log10(2) */
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI	1.772453850905516027298167483341	/* sqrt(pi) */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT_2dPI
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#endif


#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI	0.572364942924700087071713675677	/* log(sqrt(pi)) */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi)) */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	/* log(sqrt(pi/2)) */
#endif


#ifdef MATHLIB_IN_R/* Mathlib in R */

#include "Error.h"
# define MATHLIB_ERROR(fmt,x)		error(fmt,x);
# define MATHLIB_WARNING(fmt,x)		warning(fmt,x)
# define MATHLIB_WARNING2(fmt,x,x2)	warning(fmt,x,x2)
# define MATHLIB_WARNING3(fmt,x,x2,x3)	warning(fmt,x,x2,x3)
# define MATHLIB_WARNING4(fmt,x,x2,x3,x4) warning(fmt,x,x2,x3,x4)

#else/* Mathlib standalone */

#include <stdio.h>
# define MATHLIB_ERROR(fmt,x)	{ printf(fmt,x); exit(1) }
# define MATHLIB_WARNING(fmt,x)		printf(fmt,x)
# define MATHLIB_WARNING2(fmt,x,x2)	printf(fmt,x,x2)
# define MATHLIB_WARNING3(fmt,x,x2,x3)	printf(fmt,x,x2,x3)
# define MATHLIB_WARNING4(fmt,x,x2,x3,x4) printf(fmt,x,x2,x3,x4)
#endif

#define ME_NONE		0
/*	no error */
#define ME_DOMAIN	1
/*	argument out of domain */
#define ME_RANGE	2
/*	value out of range */
#define ME_NOCONV	4
/*	process did not converge */
#define ME_PRECISION	8
/*	does not have "full" precision */
#define ME_UNDERFLOW	16
/*	and underflow occured (important for IEEE)*/


#ifdef IEEE_754
extern double m_zero;
extern double m_one;
/* extern double m_tiny; */
#define ML_ERROR(x)	/* nothing */
#define ML_POSINF	(m_one / m_zero)
#define ML_NEGINF	((-m_one) / m_zero)
#define ML_NAN		(m_zero / m_zero)
#define ML_UNDERFLOW	(DBL_MIN * DBL_MIN)
#define ML_VALID(x)	(!ISNAN(x))

#else/*--- NO IEEE: No +/-Inf, NAN,... ---*/
void ml_error(int n);
#define ML_ERROR(x)	ml_error(x)
#define ML_POSINF	DBL_MAX
#define ML_NEGINF	(-DBL_MAX)
#define ML_NAN		(-DBL_MAX)
#define ML_UNDERFLOW	0
#define ML_VALID(x)	(errno == 0)
#endif

#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN); return ML_NAN; }

	/* old-R Compatibility */

#define snorm	norm_rand
#define sunif	unif_rand
#define sexp	exp_rand

	/* Undo SGI Madness */

#ifdef ftrunc
# undef ftrunc
#endif
#ifdef qexp
# undef qexp
#endif
#ifdef qgamma
# undef qgamma
#endif

	/* Name Hiding to Avoid Clashes with Fortran */

#ifdef HIDE_NAMES
# define d1mach	c_d1mach
# define i1mach	c_i1mach
#endif

#ifndef R_NO_REMAP
#define bessel_i	Rf_bessel_i
#define bessel_j	Rf_bessel_j
#define bessel_k	Rf_bessel_k
#define bessel_y	Rf_bessel_y
#define beta		Rf_beta
#define chebyshev_eval	Rf_chebyshev_eval
#define chebyshev_init	Rf_chebyshev_init
#define choose		Rf_choose
#define d1mach		Rf_d1mach
#define dbeta		Rf_dbeta
#define dbinom		Rf_dbinom
#define dcauchy		Rf_dcauchy
#define dchisq		Rf_dchisq
#define dexp		Rf_dexp
#define df		Rf_df
#define dgamma		Rf_dgamma
#define dgeom		Rf_dgeom
#define dhyper		Rf_dhyper
#define digamma		Rf_digamma
#define dlnorm		Rf_dlnorm
#define dlogis		Rf_dlogis
#define dnbeta		Rf_dnbeta
#define dnbinom		Rf_dnbinom
#define dnchisq		Rf_dnchisq
#define dnf		Rf_dnf
#define dnorm4		Rf_dnorm4
#define dnt		Rf_dnt
#define dpois		Rf_dpois
#define dpsifn		Rf_dpsifn
#define dsignrank	Rf_dsignrank
#define dt		Rf_dt
#define dtukey		Rf_dtukey
#define dunif		Rf_dunif
#define dweibull	Rf_dweibull
#define dwilcox		Rf_dwilcox
#define fastchoose	Rf_fastchoose
#define fcube		Rf_fcube
#define fmax2		Rf_fmax2
#define fmin2		Rf_fmin2
#define fmod		Rf_fmod
#define fprec		Rf_fprec
#define fround		Rf_fround
#define fsign		Rf_fsign
#define fsquare		Rf_fsquare
#define ftrunc		Rf_ftrunc
#define gammafn		Rf_gammafn
#define gammalims	Rf_gammalims
#define gamma_cody	Rf_gamma_cody
#define i1mach		Rf_i1mach
#define imax2		Rf_imax2
#define imin2		Rf_imin2
#define I_bessel	Rf_I_bessel
#define J_bessel	Rf_J_bessel
#define K_bessel	Rf_K_bessel
#define lbeta		Rf_lbeta
#define lchoose		Rf_lchoose
#define lfastchoose	Rf_lfastchoose
#define lgammacor	Rf_lgammacor
#define lgammafn	Rf_lgammafn
#define logrelerr	Rf_logrelerr
#define pbeta		Rf_pbeta
#define pbeta_raw	Rf_pbeta_raw
#define pbinom		Rf_pbinom
#define pcauchy		Rf_pcauchy
#define pchisq		Rf_pchisq
#define pentagamma	Rf_pentagamma
#define pexp		Rf_pexp
#define pf		Rf_pf
#define pgamma		Rf_pgamma
#define pgeom		Rf_pgeom
#define phyper		Rf_phyper
#define plnorm		Rf_plnorm
#define plogis		Rf_plogis
#define pnbeta		Rf_pnbeta
#define pnbinom		Rf_pnbinom
#define pnchisq		Rf_pnchisq
#define pnf		Rf_pnf
#define pnorm5		Rf_pnorm5
#define pnorm_both	Rf_pnorm_both
#define pnt		Rf_pnt
#define ppois		Rf_ppois
#define psignrank	Rf_psignrank
#define pt		Rf_pt
#define ptukey		Rf_ptukey
#define punif		Rf_punif
#define pweibull	Rf_pweibull
#define pwilcox		Rf_pwilcox
#define qbeta		Rf_qbeta
#define qbinom		Rf_qbinom
#define qcauchy		Rf_qcauchy
#define qchisq		Rf_qchisq
#define qexp		Rf_qexp
#define qf		Rf_qf
#define qgamma		Rf_qgamma
#define qgeom		Rf_qgeom
#define qhyper		Rf_qhyper
#define qlnorm		Rf_qlnorm
#define qlogis		Rf_qlogis
#define qnbeta		Rf_qnbeta
#define qnbinom		Rf_qnbinom
#define qnchisq		Rf_qnchisq
#define qnf		Rf_qnf
#define qnorm5		Rf_qnorm5
#define qnt		Rf_qnt
#define qpois		Rf_qpois
#define qsignrank	Rf_qsignrank
#define qt		Rf_qt
#define qtukey		Rf_qtukey
#define qunif		Rf_qunif
#define qweibull	Rf_qweibull
#define qwilcox		Rf_qwilcox
#define rbeta		Rf_rbeta
#define rbinom		Rf_rbinom
#define rcauchy		Rf_rcauchy
#define rchisq		Rf_rchisq
#define rexp		Rf_rexp
#define rf		Rf_rf
#define rgamma		Rf_rgamma
#define rgeom		Rf_rgeom
#define rhyper		Rf_rhyper
#define rlnorm		Rf_rlnorm
#define rlogis		Rf_rlogis
#define rnbeta		Rf_rnbeta
#define rnbinom		Rf_rnbinom
#define rnchisq		Rf_rnchisq
#define rnf		Rf_rnf
#define rnorm		Rf_rnorm
#define rnt		Rf_rnt
#define rpois		Rf_rpois
#define rsignrank	Rf_rsignrank
#define rt		Rf_rt
#define rtukey		Rf_rtukey
#define runif		Rf_runif
#define rweibull	Rf_rweibull
#define rwilcox		Rf_rwilcox
#define sign		Rf_sign
#define tetragamma	Rf_tetragamma
#define trigamma	Rf_trigamma
#define Y_bessel	Rf_Y_bessel
#endif

#define	rround	fround
#define	prec	fprec
#undef trunc
#define	trunc	ftrunc


	/* Utilities for `dpq' handling (density/probability/quantile) */

/* give_log in "d";  log_p in "p" & "q" : */
#define give_log log_p
							/* "DEFAULT" */
							/* --------- */
#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)		/* 1 */

#define R_D_Lval(p)	(lower_tail ? p	 : 1 - (p))	/*  p  */
#define R_D_Cval(p)	(lower_tail ? 1 - (p) : p )	/*  1 - p */

#define R_D_val(x)	(log_p	? log(x) : x)		/*  x  in pF(x,..) */
#define R_D_qIv(p)	(log_p	? exp(p) : p)		/*  p  in qF(p,..) */
#define R_D_exp(x)	(log_p	?  x	 : exp(x))	/* exp(x) */
#define R_D_log(p)	(log_p	?  p	 : log(p))	/* log(p) */

#define R_DT_val(x)	R_D_val(R_D_Lval(x))		/*  x  in pF */
#define R_DT_qIv(p)	R_D_qIv(R_D_Lval(p))		/*  p  in qF */
#define R_DT_Cval(x)	R_D_val(R_D_Cval(x))		/*  1 - x */
#define R_DT_CIv(x)	R_D_qIv(R_D_Cval(x))		/*  1 - x */
#define R_DT_exp(x)	R_D_exp(R_D_Lval(x))		/* exp(x) */
#define R_DT_Cexp(x)	R_D_exp(R_D_Cval(x))		/* exp(1 - x) */
#define R_DT_log(p)	R_D_log(R_D_Lval(p))		/* log(p) */
#define R_DT_Clog(p)	R_D_log(R_D_Cval(p))		/* log(1 - p) */

#define R_Q_P01_check(p)			\
    if ((log_p	&& p > 0) ||			\
	(!log_p && (p < 0 || p > 1)) )		\
	ML_ERR_return_NAN


	/* R's versions with !R_FINITE checks */

double R_log(double x);
double R_pow(double x, double y);
double R_pow_di(double, int);

	/* Machine Characteristics */

double	d1mach(int);
double	d1mach_(int*);
int	i1mach(int);
int	i1mach_(int*);

	/* General Support Functions */

int	imax2(int, int);
int	imin2(int, int);
double	fmax2(double, double);
double	fmin2(double, double);
double	fmod(double, double);
double	fprec(double, double);
double	fround(double, double);
double	ftrunc(double);
double	sign(double);
double	fsign(double, double);
double	fsquare(double);
double	fcube(double);

	/* Random Number Generators */

double	norm_rand(void);
double	unif_rand(void);
double	exp_rand(void);

	/* Chebyshev Series */

int	chebyshev_init(double*, int, double);
double	chebyshev_eval(double, double *, int);

	/* Gamma and Related Functions */

double	logrelerr(double); /* = log(1+x) {care for small x} */
void	gammalims(double*, double*);
double	lgammacor(double);
double	gammafn(double);
double	gamma_cody(double);
double	lgammafn(double);
void	dpsifn(double, int, int, int, double*, int*, int*);
double	digamma(double);
double	trigamma(double);
double	tetragamma(double);
double	pentagamma(double);

double	choose(double, double);
double	lchoose(double, double);
double	fastchoose(double, double);
double	lfastchoose(double, double);

	/* Bessel Functions of All Kinds */

double	bessel_i(double, double, double);
double	bessel_j(double, double);
double	bessel_k(double, double, double);
double	bessel_y(double, double);
void	I_bessel(double*, double*, long*, long*, double*, long*);
void	J_bessel(double*, double*, long*,	 double*, long*);
void	K_bessel(double*, double*, long*, long*, double*, long*);
void	Y_bessel(double*, double*, long*,	 double*, long*);

	/* Beta and Related Functions */

double	beta(double, double);
double	lbeta(double, double);

	/* Normal Distribution */

/* These 3 will be re-defined in  ./S_compat.h  (when that is used) : */
#define pnorm pnorm5
#define qnorm qnorm5
#define dnorm dnorm4

double	dnorm(double, double, double, int);
double	pnorm(double, double, double, int, int);
double	qnorm(double, double, double, int, int);
double	rnorm(double, double);
void	pnorm_both(double, double *, double *, int, int);/* both tails */

	/* Uniform Distribution */

double	dunif(double, double, double, int);
double	punif(double, double, double, int, int);
double	qunif(double, double, double, int, int);
double	runif(double, double);

	/* Gamma Distribution */

double	dgamma(double, double, double, int);
double	pgamma(double, double, double, int, int);
double	qgamma(double, double, double, int, int);
double	rgamma(double, double);

	/* Beta Distribution */

double	dbeta(double, double, double, int);
double	pbeta(double, double, double, int, int);
double	qbeta(double, double, double, int, int);
double	rbeta(double, double);
double	pbeta_raw(double, double, double, int);

	/* Lognormal Distribution */

double	dlnorm(double, double, double, int);
double	plnorm(double, double, double, int, int);
double	qlnorm(double, double, double, int, int);
double	rlnorm(double, double);

	/* Chi-squared Distribution */

double	dchisq(double, double, int);
double	pchisq(double, double, int, int);
double	qchisq(double, double, int, int);
double	rchisq(double);

	/* Non-central Chi-squared Distribution */

double	dnchisq(double, double, double, int);
double	pnchisq(double, double, double, int, int);
double	qnchisq(double, double, double, int, int);
double	rnchisq(double, double);

	/* F Distibution */

double	df(double, double, double, int);
double	pf(double, double, double, int, int);
double	qf(double, double, double, int, int);
double	rf(double, double);

	/* Student t Distibution */

double	dt(double, double, int);
double	pt(double, double, int, int);
double	qt(double, double, int, int);
double	rt(double);

	/* Binomial Distribution */

double	dbinom(double, double, double, int);
double	pbinom(double, double, double, int, int);
double	qbinom(double, double, double, int, int);
double	rbinom(double, double);

	/* Cauchy Distribution */

double	dcauchy(double, double, double, int);
double	pcauchy(double, double, double, int, int);
double	qcauchy(double, double, double, int, int);
double	rcauchy(double, double);

	/* Exponential Distribution */

double	dexp(double, double, int);
double	pexp(double, double, int, int);
double	qexp(double, double, int, int);
double	rexp(double);

	/* Geometric Distribution */

double	dgeom(double, double, int);
double	pgeom(double, double, int, int);
double	qgeom(double, double, int, int);
double	rgeom(double);

	/* Hypergeometric Distibution */

double	dhyper(double, double, double, double, int);
double	phyper(double, double, double, double, int, int);
double	qhyper(double, double, double, double, int, int);
double	rhyper(double, double, double);

	/* Negative Binomial Distribution */

double	dnbinom(double, double, double, int);
double	pnbinom(double, double, double, int, int);
double	qnbinom(double, double, double, int, int);
double	rnbinom(double, double);

	/* Poisson Distribution */

double	dpois(double, double, int);
double	ppois(double, double, int, int);
double	qpois(double, double, int, int);
double	rpois(double);

	/* Weibull Distribution */

double	dweibull(double, double, double, int);
double	pweibull(double, double, double, int, int);
double	qweibull(double, double, double, int, int);
double	rweibull(double, double);

	/* Logistic Distribution */

double	dlogis(double, double, double, int);
double	plogis(double, double, double, int, int);
double	qlogis(double, double, double, int, int);
double	rlogis(double, double);

	/* Non-central Beta Distribution */

double	dnbeta(double, double, double, double, int);
double	pnbeta(double, double, double, double, int, int);
double	qnbeta(double, double, double, double, int, int);
double	rnbeta(double, double, double);

	/* Non-central F Distribution */

double	dnf(double, double, double, double, int);
double	pnf(double, double, double, double, int, int);
double	qnf(double, double, double, double, int, int);
double	rnf(double, double, double);

	/* Non-central Student t Distribution */

double	dnt(double, double, double, int);
double	pnt(double, double, double, int, int);
double	qnt(double, double, double, int, int);
double	rnt(double, double);

	/* Studentized Range Distribution */

double	dtukey(double, double, double, double, int);
double	ptukey(double, double, double, double, int, int);
double	qtukey(double, double, double, double, int, int);
double	rtukey(double, double, double);

/* Wilcoxon Rank Sum Distribution */

#define WILCOX_MAX 50
double dwilcox(double, double, double, int);
double pwilcox(double, double, double, int, int);
double qwilcox(double, double, double, int, int);
double rwilcox(double, double);

/* Wilcoxon Signed Rank Distribution */

#define SIGNRANK_MAX 50
double dsignrank(double, double, int);
double psignrank(double, double, int, int);
double qsignrank(double, double, int, int);
double rsignrank(double);

#endif
