/* Interface routines, callable from R using .Call, for Lapack code */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Lapack.h"

static SEXP modLa_svd(SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, SEXP v, SEXP method)
{
    int *xdims, n, p, lwork, info = 0;
    double *work, *xvals, tmp;
    SEXP val, nm;
    char *meth;

    if (!(isString(jobu) && isString(jobv)))
        error("jobu and jobv must be character objects");
    if (!isString(method))
        error("method must be a character object");
    meth = CHAR(STRING_ELT(method, 0));
#ifndef IEEE_754
    if (strcmp(meth, "dgesdd") == 0)
    {
        warning("method = \"dgesdd\" requires IEEE 754 arithmetic: using \"dgesvd\"");
        meth = "dgesvd";
    }
#endif
    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
    p = xdims[1];
    xvals = (double *)R_alloc(n * p, sizeof(double));
    /* work on a copy of x */
    Memcpy(xvals, REAL(x), (size_t)(n * p));

    if (strcmp(meth, "dgesdd"))
    {
        /* ask for optimal size of work array */
        lwork = -1;
        F77_CALL(dgesvd)
        (CHAR(STRING_ELT(jobu, 0)), CHAR(STRING_ELT(jobv, 0)), &n, &p, xvals, &n, REAL(s), REAL(u),
         INTEGER(getAttrib(u, R_DimSymbol)), REAL(v), INTEGER(getAttrib(v, R_DimSymbol)), &tmp, &lwork, &info);
        lwork = (int)tmp;

        work = (double *)R_alloc(lwork, sizeof(double));
        F77_CALL(dgesvd)
        (CHAR(STRING_ELT(jobu, 0)), CHAR(STRING_ELT(jobv, 0)), &n, &p, xvals, &n, REAL(s), REAL(u),
         INTEGER(getAttrib(u, R_DimSymbol)), REAL(v), INTEGER(getAttrib(v, R_DimSymbol)), work, &lwork, &info);
        if (info != 0)
            error("error code %d from Lapack routine dgesvd", info);
    }
    else
    {
        int ldu = INTEGER(getAttrib(u, R_DimSymbol))[0], ldvt = INTEGER(getAttrib(v, R_DimSymbol))[0];
        int *iwork = (int *)R_alloc(8 * (n < p ? n : p), sizeof(int));

        /* ask for optimal size of work array */
        lwork = -1;

        F77_CALL(dgesdd)
        (CHAR(STRING_ELT(jobu, 0)), &n, &p, xvals, &n, REAL(s), REAL(u), &ldu, REAL(v), &ldvt, &tmp, &lwork, iwork,
         &info);
        lwork = (int)tmp;

        work = (double *)R_alloc(lwork, sizeof(double));
        F77_CALL(dgesdd)
        (CHAR(STRING_ELT(jobu, 0)), &n, &p, xvals, &n, REAL(s), REAL(u), &ldu, REAL(v), &ldvt, work, &lwork, iwork,
         &info);
        if (info != 0)
            error("error code %d from Lapack routine dgesdd", info);
    }

    val = PROTECT(allocVector(VECSXP, 3));
    nm = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nm, 0, mkChar("d"));
    SET_STRING_ELT(nm, 1, mkChar("u"));
    SET_STRING_ELT(nm, 2, mkChar("vt"));
    setAttrib(val, R_NamesSymbol, nm);
    SET_VECTOR_ELT(val, 0, s);
    SET_VECTOR_ELT(val, 1, u);
    SET_VECTOR_ELT(val, 2, v);
    UNPROTECT(2);
    return val;
}

static SEXP modLa_rs(SEXP xin, SEXP only_values, SEXP method)
{
    int *xdims, n, lwork, info = 0, ov;
    char jobv[1], uplo[1], range[1];
    SEXP values, ret, nm, x, z = R_NilValue;
    double *work, *rx, *rvalues, tmp;
    char *meth;

    if (!isString(method))
        error("method must be a character object");
    meth = CHAR(STRING_ELT(method, 0));
#ifndef IEEE_754
    if (strcmp(meth, "dsyevr") == 0)
    {
        warning("method = \"dseyvr\" requires IEEE 754 arithmetic: using \"dsyev\"");
        meth = "dsyev";
    }
#endif
    PROTECT(x = duplicate(xin));
    rx = REAL(x);
    uplo[0] = 'L';
    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
    if (n != xdims[1])
        error("x must be a square numeric matrix");
    ov = asLogical(only_values);
    if (ov == NA_LOGICAL)
        error("invalid `only.values'");
    if (ov)
        jobv[0] = 'N';
    else
        jobv[0] = 'V';

    PROTECT(values = allocVector(REALSXP, n));
    rvalues = REAL(values);
    if (strcmp(meth, "dsyevr"))
    {
        /* ask for optimal size of work array */
        lwork = -1;
        F77_CALL(rsyev)(jobv, uplo, &n, rx, &n, rvalues, &tmp, &lwork, &info);
        lwork = (int)tmp;
        if (lwork < 3 * n - 1)
            lwork = 3 * n - 1; /* Sanity check */
        work = (double *)R_alloc(lwork, sizeof(double));
        F77_CALL(rsyev)(jobv, uplo, &n, rx, &n, rvalues, work, &lwork, &info);
        if (info != 0)
            error("error code %d from Lapack routine dsyev", info);
    }
    else
    {
        int liwork, *iwork, itmp, m;
        double vl, vu, abstol = 0.0;
        int il, iu, *isuppz;

        range[0] = 'A';
        if (!ov)
            PROTECT(z = allocMatrix(REALSXP, n, n));
        isuppz = (int *)R_alloc(2 * n, sizeof(int));
        /* ask for optimal size of work arrays */
        lwork = -1;
        liwork = -1;
        F77_CALL(rsyevr)
        (jobv, range, uplo, &n, rx, &n, &vl, &vu, &il, &iu, &abstol, &m, rvalues, REAL(z), &n, isuppz, &tmp, &lwork,
         &itmp, &liwork, &info);
        lwork = (int)tmp;
        liwork = itmp;

        work = (double *)R_alloc(lwork, sizeof(double));
        iwork = (int *)R_alloc(liwork, sizeof(int));
        F77_CALL(rsyevr)
        (jobv, range, uplo, &n, rx, &n, &vl, &vu, &il, &iu, &abstol, &m, rvalues, REAL(z), &n, isuppz, work, &lwork,
         iwork, &liwork, &info);
        if (info != 0)
            error("error code %d from Lapack routine dsyev", info);
    }

    if (!ov)
    {
        ret = PROTECT(allocVector(VECSXP, 2));
        nm = PROTECT(allocVector(STRSXP, 2));
        SET_STRING_ELT(nm, 1, mkChar("vectors"));
        if (strcmp(meth, "dsyevr"))
        {
            SET_VECTOR_ELT(ret, 1, x);
        }
        else
        {
            SET_VECTOR_ELT(ret, 1, z);
            UNPROTECT_PTR(z);
        }
    }
    else
    {
        ret = PROTECT(allocVector(VECSXP, 1));
        nm = PROTECT(allocVector(STRSXP, 1));
    }
    SET_STRING_ELT(nm, 0, mkChar("values"));
    setAttrib(ret, R_NamesSymbol, nm);
    SET_VECTOR_ELT(ret, 0, values);
    UNPROTECT(4);
    return ret;
}

static SEXP unscramble(const double *imaginary, int n, const double *vecs)
{
    int i, j;
    SEXP s = allocMatrix(CPLXSXP, n, n);

    for (j = 0; j < n; j++)
    {
        if (imaginary[j] != 0)
        {
            int j1 = j + 1;
            for (i = 0; i < n; i++)
            {
                COMPLEX(s)[i + n * j].r = COMPLEX(s)[i + n * j1].r = vecs[i + j * n];
                COMPLEX(s)[i + n * j1].i = -(COMPLEX(s)[i + n * j].i = vecs[i + j1 * n]);
            }
            j = j1;
        }
        else
        {
            for (i = 0; i < n; i++)
            {
                COMPLEX(s)[i + n * j].r = vecs[i + j * n];
                COMPLEX(s)[i + n * j].i = 0.0;
            }
        }
    }
    return s;
}

static SEXP modLa_rg(SEXP x, SEXP only_values)
{
    Rboolean vectors, complexValues;
    int i, n, lwork, info, *xdims, ov;
    double *work, *wR, *wI, *left, *right, *xvals, tmp;
    char jobVL[1], jobVR[1];
    SEXP ret, nm, val;

    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
    if (n != xdims[1])
        error("x must be a square numeric matrix");

    xvals = (double *)R_alloc(n * n, sizeof(double));
    /* work on a copy of x */
    Memcpy(xvals, REAL(x), (size_t)(n * n));
    ov = asLogical(only_values);
    if (ov == NA_LOGICAL)
        error("invalid `only.values'");
    vectors = !ov;
    jobVL[0] = jobVR[0] = 'N';
    left = right = (double *)0;
    if (vectors)
    {
        jobVR[0] = 'V';
        right = (double *)R_alloc(n * n, sizeof(double));
    }
    wR = (double *)R_alloc(n, sizeof(double));
    wI = (double *)R_alloc(n, sizeof(double));
    /* ask for optimal size of work array */
    lwork = -1;
    F77_CALL(dgeev)(jobVL, jobVR, &n, xvals, &n, wR, wI, left, &n, right, &n, &tmp, &lwork, &info);
    lwork = (int)tmp;
    work = (double *)R_alloc(lwork, sizeof(double));
    F77_CALL(dgeev)(jobVL, jobVR, &n, xvals, &n, wR, wI, left, &n, right, &n, work, &lwork, &info);
    if (info != 0)
        error("error code %d from Lapack routine dgeev", info);

    complexValues = FALSE;
    for (i = 0; i < n; i++)
        if (wI[i] != 0.0)
        {
            complexValues = TRUE;
            break;
        }
    ret = PROTECT(allocVector(VECSXP, 2));
    nm = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(nm, 0, mkChar("values"));
    SET_STRING_ELT(nm, 1, mkChar("vectors"));
    setAttrib(ret, R_NamesSymbol, nm);
    SET_VECTOR_ELT(ret, 1, R_NilValue);
    if (complexValues)
    {
        val = allocVector(CPLXSXP, n);
        for (i = 0; i < n; i++)
        {
            COMPLEX(val)[i].r = wR[i];
            COMPLEX(val)[i].i = wI[i];
        }
        SET_VECTOR_ELT(ret, 0, val);

        if (vectors)
            SET_VECTOR_ELT(ret, 1, unscramble(wI, n, right));
    }
    else
    {
        val = allocVector(REALSXP, n);
        for (i = 0; i < n; i++)
            REAL(val)[i] = wR[i];
        SET_VECTOR_ELT(ret, 0, val);
        if (vectors)
        {
            val = allocMatrix(REALSXP, n, n);
            for (i = 0; i < (n * n); i++)
                REAL(val)[i] = right[i];
            SET_VECTOR_ELT(ret, 1, val);
        }
    }
    UNPROTECT(2);
    return ret;
}

static SEXP modLa_zgesv(SEXP A, SEXP B)
{
#ifdef HAVE_DOUBLE_COMPLEX
    int n, p, info, *ipiv, *Adims, *Bdims;
    Rcomplex *avals;

    if (!(isMatrix(A) && isComplex(A)))
        error("A must be a complex matrix");
    if (!(isMatrix(B) && isComplex(B)))
        error("A must be a complex matrix");
    Adims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    Bdims = INTEGER(coerceVector(getAttrib(B, R_DimSymbol), INTSXP));
    n = Adims[0];
    if (n == 0)
        error("A is 0-diml");
    p = Bdims[1];
    if (p == 0)
        error("no rhs in B");
    if (Adims[1] != n)
        error("A (%d x %d) must be square", n, Adims[1]);
    if (Bdims[0] != n)
        error("B (%d x %d) must be square", Bdims[0], p);
    ipiv = (int *)R_alloc(n, sizeof(int));

    avals = (Rcomplex *)R_alloc(n * n, sizeof(Rcomplex));
    /* work on a copy of x */
    Memcpy(avals, COMPLEX(A), (size_t)(n * n));
    F77_CALL(zgesv)(&n, &p, avals, &n, ipiv, COMPLEX(B), &n, &info);
    if (info != 0)
        error("error code %d from Lapack routine zgesv", info);
    return B;
#else
    error("Fortran complex functions are not available on this platform");
    return R_NilValue; /* -Wall */
#endif
}

static SEXP modLa_zgeqp3(SEXP Ain)
{
#ifdef HAVE_DOUBLE_COMPLEX
    int m, n, *Adims, info, lwork;
    Rcomplex *work, tmp;
    double *rwork;
    SEXP val, nm, jpvt, tau, rank, A;

    if (!(isMatrix(Ain) && isComplex(Ain)))
        error("A must be a complex matrix");
    PROTECT(A = duplicate(Ain));
    Adims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    m = Adims[0];
    n = Adims[1];
    rwork = (double *)R_alloc(2 * n, sizeof(double));

    jpvt = PROTECT(allocVector(INTSXP, n));
    tau = PROTECT(allocVector(CPLXSXP, m < n ? m : n));
    lwork = -1;
    F77_CALL(zgeqp3)(&m, &n, COMPLEX(A), &m, INTEGER(jpvt), COMPLEX(tau), &tmp, &lwork, rwork, &info);
    lwork = (int)tmp.r;
    work = (Rcomplex *)R_alloc(lwork, sizeof(Rcomplex));
    F77_CALL(zgeqp3)(&m, &n, COMPLEX(A), &m, INTEGER(jpvt), COMPLEX(tau), work, &lwork, rwork, &info);
    if (info != 0)
        error("error code %d from Lapack routine zqeqp3", info);
    val = PROTECT(allocVector(VECSXP, 4));
    nm = PROTECT(allocVector(STRSXP, 4));
    rank = PROTECT(allocVector(INTSXP, 1));
    INTEGER(rank)[0] = m < n ? m : n;
    SET_STRING_ELT(nm, 0, mkChar("qr"));
    SET_STRING_ELT(nm, 1, mkChar("rank"));
    SET_STRING_ELT(nm, 2, mkChar("qraux"));
    SET_STRING_ELT(nm, 3, mkChar("pivot"));
    setAttrib(val, R_NamesSymbol, nm);
    SET_VECTOR_ELT(val, 0, A);
    SET_VECTOR_ELT(val, 1, rank);
    SET_VECTOR_ELT(val, 2, tau);
    SET_VECTOR_ELT(val, 3, jpvt);
    UNPROTECT(6);
    return val;
#else
    error("Fortran complex functions are not available on this platform");
    return R_NilValue; /* -Wall */
#endif
}

static SEXP modqr_coef_cmplx(SEXP Q, SEXP Bin)
{
#ifdef HAVE_DOUBLE_COMPLEX
    int n, nrhs, lwork, info, k, *Bdims, *Qdims;
    SEXP B, qr = VECTOR_ELT(Q, 0), tau = VECTOR_ELT(Q, 2);
    Rcomplex *work, tmp;

    k = LENGTH(tau);
    if (!(isMatrix(Bin) && isComplex(Bin)))
        error("B must be a complex matrix");

    PROTECT(B = duplicate(Bin));
    Qdims = INTEGER(coerceVector(getAttrib(qr, R_DimSymbol), INTSXP));
    n = Qdims[0];
    Bdims = INTEGER(coerceVector(getAttrib(B, R_DimSymbol), INTSXP));
    if (Bdims[0] != n)
        error("rhs should have %d not %d rows", n, Bdims[0]);
    nrhs = Bdims[1];
    lwork = -1;
    F77_CALL(zunmqr)("L", "C", &n, &nrhs, &k, COMPLEX(qr), &n, COMPLEX(tau), COMPLEX(B), &n, &tmp, &lwork, &info);
    lwork = (int)tmp.r;
    work = (Rcomplex *)R_alloc(lwork, sizeof(Rcomplex));
    F77_CALL(zunmqr)("L", "C", &n, &nrhs, &k, COMPLEX(qr), &n, COMPLEX(tau), COMPLEX(B), &n, work, &lwork, &info);
    if (info != 0)
        error("error code %d from Lapack routine zunmqr", info);
    F77_CALL(ztrtrs)("U", "N", "N", &n, &nrhs, COMPLEX(qr), &n, COMPLEX(B), &n, &info);
    if (info != 0)
        error("error code %d from Lapack routine ztrtrs", info);
    UNPROTECT(1);
    return B;
#else
    error("Fortran complex functions are not available on this platform");
    return R_NilValue; /* -Wall */
#endif
}

static SEXP modqr_qy_cmplx(SEXP Q, SEXP Bin, SEXP trans)
{
#ifdef HAVE_DOUBLE_COMPLEX
    int n, nrhs, lwork, info, k, *Bdims, *Qdims, tr;
    SEXP B, qr = VECTOR_ELT(Q, 0), tau = VECTOR_ELT(Q, 2);
    Rcomplex *work, tmp;

    k = LENGTH(tau);
    if (!(isMatrix(Bin) && isComplex(Bin)))
        error("B must be a complex matrix");
    tr = asLogical(trans);
    if (tr == NA_LOGICAL)
        error("invalid `trans' parameter");

    PROTECT(B = duplicate(Bin));
    Qdims = INTEGER(coerceVector(getAttrib(qr, R_DimSymbol), INTSXP));
    n = Qdims[0];
    Bdims = INTEGER(coerceVector(getAttrib(B, R_DimSymbol), INTSXP));
    if (Bdims[0] != n)
        error("rhs should have %d not %d rows", n, Bdims[0]);
    nrhs = Bdims[1];
    lwork = -1;
    F77_CALL(zunmqr)
    ("L", tr ? "C" : "N", &n, &nrhs, &k, COMPLEX(qr), &n, COMPLEX(tau), COMPLEX(B), &n, &tmp, &lwork, &info);
    lwork = (int)tmp.r;
    work = (Rcomplex *)R_alloc(lwork, sizeof(Rcomplex));
    F77_CALL(zunmqr)
    ("L", tr ? "C" : "N", &n, &nrhs, &k, COMPLEX(qr), &n, COMPLEX(tau), COMPLEX(B), &n, work, &lwork, &info);
    if (info != 0)
        error("error code %d from Lapack routine zunmqr", info);
    UNPROTECT(1);
    return B;
#else
    error("Fortran complex functions are not available on this platform");
    return R_NilValue; /* -Wall */
#endif
}

static SEXP modLa_svd_cmplx(SEXP jobu, SEXP jobv, SEXP xin, SEXP s, SEXP u, SEXP v)
{
#ifdef HAVE_DOUBLE_COMPLEX
    int *xdims, n, p, lwork, info;
    double *rwork;
    Rcomplex *work, tmp;
    SEXP x, val, nm;

    if (!(isString(jobu) && isString(jobv)))
        error("jobu and jobv must be character objects");
    PROTECT(x = duplicate(xin));
    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
    p = xdims[1];
    rwork = (double *)R_alloc(5 * (n < p ? n : p), sizeof(double));
    /* ask for optimal size of work array */
    lwork = -1;
    F77_CALL(zgesvd)
    (CHAR(STRING_ELT(jobu, 0)), CHAR(STRING_ELT(jobv, 0)), &n, &p, COMPLEX(x), &n, REAL(s), COMPLEX(u),
     INTEGER(getAttrib(u, R_DimSymbol)), COMPLEX(v), INTEGER(getAttrib(v, R_DimSymbol)), &tmp, &lwork, rwork, &info);
    lwork = (int)tmp.r;
    work = (Rcomplex *)R_alloc(lwork, sizeof(Rcomplex));
    F77_CALL(zgesvd)
    (CHAR(STRING_ELT(jobu, 0)), CHAR(STRING_ELT(jobv, 0)), &n, &p, COMPLEX(x), &n, REAL(s), COMPLEX(u),
     INTEGER(getAttrib(u, R_DimSymbol)), COMPLEX(v), INTEGER(getAttrib(v, R_DimSymbol)), work, &lwork, rwork, &info);
    if (info != 0)
        error("error code %d from Lapack routine dgesvd", info);
    val = PROTECT(allocVector(VECSXP, 3));
    nm = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nm, 0, mkChar("d"));
    SET_STRING_ELT(nm, 1, mkChar("u"));
    SET_STRING_ELT(nm, 2, mkChar("vt"));
    setAttrib(val, R_NamesSymbol, nm);
    SET_VECTOR_ELT(val, 0, s);
    SET_VECTOR_ELT(val, 1, u);
    SET_VECTOR_ELT(val, 2, v);
    UNPROTECT(3);
    return val;
#else
    error("Fortran complex functions are not available on this platform");
    return R_NilValue; /* -Wall */
#endif
}

static SEXP modLa_rs_cmplx(SEXP xin, SEXP only_values)
{
#ifdef HAVE_DOUBLE_COMPLEX
    int *xdims, n, lwork, info, ov;
    char jobv[1], uplo[1];
    SEXP values, ret, nm, x;
    Rcomplex *work, *rx, tmp;
    double *rwork, *rvalues;

    PROTECT(x = duplicate(xin));
    rx = COMPLEX(x);
    uplo[0] = 'L';
    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
    if (n != xdims[1])
        error("x must be a square numeric matrix");
    ov = asLogical(only_values);
    if (ov == NA_LOGICAL)
        error("invalid `only.values'");
    if (ov)
        jobv[0] = 'N';
    else
        jobv[0] = 'V';

    PROTECT(values = allocVector(REALSXP, n));
    rvalues = REAL(values);
    rwork = (double *)R_alloc((3 * n - 2) > 1 ? 3 * n - 2 : 1, sizeof(double));
    /* ask for optimal size of work array */
    lwork = -1;
    F77_CALL(zheev)(jobv, uplo, &n, rx, &n, rvalues, &tmp, &lwork, rwork, &info);
    lwork = (int)tmp.r;
    work = (Rcomplex *)R_alloc(lwork, sizeof(Rcomplex));
    F77_CALL(zheev)(jobv, uplo, &n, rx, &n, rvalues, work, &lwork, rwork, &info);
    if (info != 0)
        error("error code %d from Lapack routine dsyev", info);
    if (!ov)
    {
        ret = PROTECT(allocVector(VECSXP, 2));
        nm = PROTECT(allocVector(STRSXP, 2));
        SET_STRING_ELT(nm, 1, mkChar("vectors"));
        SET_VECTOR_ELT(ret, 1, x);
    }
    else
    {
        ret = PROTECT(allocVector(VECSXP, 1));
        nm = PROTECT(allocVector(STRSXP, 1));
    }
    SET_STRING_ELT(nm, 0, mkChar("values"));
    setAttrib(ret, R_NamesSymbol, nm);
    SET_VECTOR_ELT(ret, 0, values);
    UNPROTECT(4);
    return ret;
#else
    error("Fortran complex functions are not available on this platform");
    return R_NilValue; /* -Wall */
#endif
}

static SEXP modLa_rg_cmplx(SEXP x, SEXP only_values)
{
#ifdef HAVE_DOUBLE_COMPLEX
    int n, lwork, info, *xdims, ov;
    Rcomplex *work, *left, *right, *xvals, tmp;
    double *rwork;
    char jobVL[1], jobVR[1];
    SEXP ret, nm, values, val = R_NilValue;

    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
    if (n != xdims[1])
        error("x must be a square numeric matrix");

    xvals = (Rcomplex *)R_alloc(n * n, sizeof(Rcomplex));
    /* work on a copy of x */
    Memcpy(xvals, COMPLEX(x), (size_t)(n * n));
    ov = asLogical(only_values);
    if (ov == NA_LOGICAL)
        error("invalid `only.values'");
    jobVL[0] = jobVR[0] = 'N';
    left = right = (Rcomplex *)0;
    if (!ov)
    {
        jobVR[0] = 'V';
        PROTECT(val = allocMatrix(CPLXSXP, n, n));
        right = COMPLEX(val);
    }
    PROTECT(values = allocVector(CPLXSXP, n));
    rwork = (double *)R_alloc(2 * n, sizeof(double));
    /* ask for optimal size of work array */
    lwork = -1;
    F77_CALL(zgeev)(jobVL, jobVR, &n, xvals, &n, COMPLEX(values), left, &n, right, &n, &tmp, &lwork, rwork, &info);
    lwork = (int)tmp.r;
    work = (Rcomplex *)R_alloc(lwork, sizeof(Rcomplex));
    F77_CALL(zgeev)(jobVL, jobVR, &n, xvals, &n, COMPLEX(values), left, &n, right, &n, work, &lwork, rwork, &info);
    if (info != 0)
        error("error code %d from Lapack routine zgeev", info);

    if (!ov)
    {
        ret = PROTECT(allocVector(VECSXP, 2));
        nm = PROTECT(allocVector(STRSXP, 2));
        SET_STRING_ELT(nm, 1, mkChar("vectors"));
        SET_VECTOR_ELT(ret, 1, val);
    }
    else
    {
        ret = PROTECT(allocVector(VECSXP, 1));
        nm = PROTECT(allocVector(STRSXP, 1));
    }
    SET_STRING_ELT(nm, 0, mkChar("values"));
    SET_VECTOR_ELT(ret, 0, values);
    setAttrib(ret, R_NamesSymbol, nm);
    UNPROTECT(ov ? 3 : 4);
    return ret;
#else
    error("Fortran complex functions are not available on this platform");
    return R_NilValue; /* -Wall */
#endif
}

#include "R_ext/Rlapack.h"
#include "R_ext/Rdynload.h"

void R_init_lapack(DllInfo *info)
{
    R_LapackRoutines *tmp;
    tmp = (R_LapackRoutines *)malloc(sizeof(R_LapackRoutines));

    tmp->svd = modLa_svd;
    tmp->rs = modLa_rs;
    tmp->rg = modLa_rg;
    tmp->zgesv = modLa_zgesv;
    tmp->zgeqp3 = modLa_zgeqp3;
    tmp->qr_coef_cmplx = modqr_coef_cmplx;
    tmp->qr_qy_cmplx = modqr_qy_cmplx;
    tmp->svd_cmplx = modLa_svd_cmplx;
    tmp->rs_cmplx = modLa_rs_cmplx;
    tmp->rg_cmplx = modLa_rg_cmplx;

    R_setLapackRoutines(tmp);
}

#ifdef Win32
#define PSIGNAL
#include "psignal.h"
/* force in malloc & free, so ATLAS gets the right ones */
/* also force in signal, although what's using that is unclear */
void lapack_dummy()
{
    char *foo;
    foo = (char *)malloc(1);
    free(foo);
    signal(SIGBREAK, NULL);
}
#endif
