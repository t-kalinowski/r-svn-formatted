	/* Utilities for `dpq' handling (density/probability/quantile) */

/* give_log in "d";  log_p in "p" & "q" : */
#define give_log log_p
							/* "DEFAULT" */
							/* --------- */
#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)		/* 1 */

#define R_D_Lval(p)	(lower_tail ? (p) : (1 - (p)))	/*  p  */
#define R_D_Cval(p)	(lower_tail ? (1 - (p)) : (p))	/*  1 - p */

#define R_D_val(x)	(log_p	? log(x) : (x))		/*  x  in pF(x,..) */
#define R_D_qIv(p)	(log_p	? exp(p) : (p))		/*  p  in qF(p,..) */
#define R_D_exp(x)	(log_p	?  (x)	 : exp(x))	/* exp(x) */
#define R_D_log(p)	(log_p	?  (p)	 : log(p))	/* log(p) */

#define R_DT_val(x)	R_D_val(R_D_Lval(x))		/*  x  in pF */
#define R_DT_Cval(x)	R_D_val(R_D_Cval(x))		/*  1 - x */
#define R_DT_qIv(p)	R_D_Lval(R_D_qIv(p))		/*  p  in qF ! */
#define R_DT_CIv(p)	R_D_Cval(R_D_qIv(p))		/*  1 - p in qF */

#define R_DT_exp(x)	R_D_exp(R_D_Lval(x))		/* exp(x) */
#define R_DT_Cexp(x)	R_D_exp(R_D_Cval(x))		/* exp(1 - x) */

#define R_DT_log(p)	(lower_tail ? R_D_log(p) :		\
			 logrelerr(- (log_p ? exp(p) : p)))/* log(p)	in qF */

#define R_DT_Clog(p)	(lower_tail ?				\
			 logrelerr(- (log_p ? exp(p) : p)) :	\
			 R_D_log(p))			/* log(1 - p)	in qF */

#define R_Q_P01_check(p)			\
    if ((log_p	&& p > 0) ||			\
	(!log_p && (p < 0 || p > 1)) )		\
	ML_ERR_return_NAN

