/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997, 1998  Robert Gentleman, Ross Ihaka and the R core team.
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
 */

#include "Defn.h"
#include "Mathlib.h"
#include "Graphics.h"
#include "Print.h"

/*  C o n t o u r   P l o t t i n g  */

typedef struct SEG
{
    struct SEG *next;
    double x0;
    double y0;
    double x1;
    double y1;
} SEG, *SEGP;

static SEGP *ctr_SegDB;

static int ctr_intersect(double z0, double z1, double zc, double *f)
{
    if ((z0 - zc) * (z1 - zc) < 0.0)
    {
        *f = (zc - z0) / (z1 - z0);
        return 1;
    }
    return 0;
}

static SEGP ctr_newseg(double x0, double y0, double x1, double y1, SEGP prev)
{
    SEGP seg = (SEGP)R_alloc(1, sizeof(SEG));
    seg->x0 = x0;
    seg->y0 = y0;
    seg->x1 = x1;
    seg->y1 = y1;
    seg->next = prev;
    return seg;
}

static void ctr_swapseg(SEGP seg)
{
    double x, y;
    x = seg->x0;
    y = seg->y0;
    seg->x0 = seg->x1;
    seg->y0 = seg->y1;
    seg->x1 = x;
    seg->y1 = y;
}

/* Determine the entry direction to the next cell */
/* and update the cell indices */

#ifdef OLD
#define XMATCH(x0, x1) (fabs(x0 - x1) < ctr_xtol)
#define YMATCH(y0, y1) (fabs(y0 - y1) < ctr_ytol)
#else
#define XMATCH(x0, x1) (fabs(x0 - x1) == 0)
#define YMATCH(y0, y1) (fabs(y0 - y1) == 0)
#endif

static double ctr_xtol;
static double ctr_ytol;

static int ctr_segdir(double xend, double yend, double *x, double *y, int *i, int *j, int nx, int ny)
{
    if (YMATCH(yend, y[*j]))
    {
        if (*j == 0)
            return 0;
        *j = *j - 1;
        return 3;
    }
    if (XMATCH(xend, x[*i]))
    {
        if (*i == 0)
            return 0;
        *i = *i - 1;
        return 4;
    }
    if (YMATCH(yend, y[*j + 1]))
    {
        if (*j >= ny - 1)
            return 0;
        *j = *j + 1;
        return 1;
    }
    if (XMATCH(xend, x[*i + 1]))
    {
        if (*i >= nx - 1)
            return 0;
        *i = *i + 1;
        return 2;
    }
    return 0;
}

/* Search seglist for a segment with endpoint (xend, yend). */
/* The cell entry direction is dir, and if tail=1/0 we are */
/* building the tail/head of a contour.	 The matching segment */
/* is pointed to by seg and the updated segment list (with */
/* the matched segment stripped is returned by the funtion. */

static SEGP ctr_segupdate(double xend, double yend, int dir, int tail, SEGP seglist, SEGP *seg)
{
    if (seglist == NULL)
    {
        *seg = NULL;
        return NULL;
    }
    switch (dir)
    {
    case 1:
    case 3:
        if (YMATCH(yend, seglist->y0))
        {
            if (!tail)
                ctr_swapseg(seglist);
            *seg = seglist;
            return seglist->next;
        }
        if (YMATCH(yend, seglist->y1))
        {
            if (tail)
                ctr_swapseg(seglist);
            *seg = seglist;
            return seglist->next;
        }
        break;
    case 2:
    case 4:
        if (XMATCH(xend, seglist->x0))
        {
            if (!tail)
                ctr_swapseg(seglist);
            *seg = seglist;
            return seglist->next;
        }
        if (XMATCH(xend, seglist->x1))
        {
            if (tail)
                ctr_swapseg(seglist);
            *seg = seglist;
            return seglist->next;
        }
        break;
    }
    seglist->next = ctr_segupdate(xend, yend, dir, tail, seglist->next, seg);
    return seglist;
}

static void contour(SEXP x, int nx, SEXP y, int ny, SEXP z, double zc, double atom, DevDesc *dd)
{
    double f, xl, xh, yl, yh, zll, zhl, zlh, zhh, xx[4], yy[4];
    double xend, yend;
    int i, ii, j, jj, k, l, m, nacode, ns, ns2, dir;
    SEGP seglist, seg, s, start, end;
    double *xxx, *yyy;

    for (i = 0; i < nx - 1; i++)
    {
        xl = REAL(x)[i];
        xh = REAL(x)[i + 1];
        for (j = 0; j < ny - 1; j++)
        {
            yl = REAL(y)[j];
            yh = REAL(y)[j + 1];
            k = i + j * nx;
            zll = REAL(z)[k];
            zhl = REAL(z)[k + 1];
            zlh = REAL(z)[k + nx];
            zhh = REAL(z)[k + nx + 1];
            k = 0;

            /* If the value at a corner is */
            /* exactly equal to a contour */
            /* level, change the value at */
            /* corner by a tiny amount. */

            if (zll == zc)
                zll = zll + atom;
            if (zhl == zc)
                zhl = zhl + atom;
            if (zlh == zc)
                zlh = zlh + atom;
            if (zhh == zc)
                zhh = zhh + atom;

            /* Check for intersections with sides */

            nacode = 0;
            if (FINITE(zll))
                nacode += 1;
            if (FINITE(zhl))
                nacode += 2;
            if (FINITE(zlh))
                nacode += 4;
            if (FINITE(zhh))
                nacode += 8;

            switch (nacode)
            {
            case 15:
                if (ctr_intersect(zll, zhl, zc, &f))
                {
                    xx[k] = xl + f * (xh - xl);
                    yy[k] = yl;
                    k++;
                }
                if (ctr_intersect(zll, zlh, zc, &f))
                {
                    yy[k] = yl + f * (yh - yl);
                    xx[k] = xl;
                    k++;
                }
                if (ctr_intersect(zhl, zhh, zc, &f))
                {
                    yy[k] = yl + f * (yh - yl);
                    xx[k] = xh;
                    k++;
                }
                if (ctr_intersect(zlh, zhh, zc, &f))
                {
                    xx[k] = xl + f * (xh - xl);
                    yy[k] = yh;
                    k++;
                }
                break;
            case 14:
                if (ctr_intersect(zhl, zhh, zc, &f))
                {
                    yy[k] = yl + f * (yh - yl);
                    xx[k] = xh;
                    k++;
                }
                if (ctr_intersect(zlh, zhh, zc, &f))
                {
                    xx[k] = xl + f * (xh - xl);
                    yy[k] = yh;
                    k++;
                }
                if (ctr_intersect(zlh, zhl, zc, &f))
                {
                    xx[k] = xl + f * (xh - xl);
                    yy[k] = yh + f * (yl - yh);
                    k++;
                }
                break;
            case 13:
                if (ctr_intersect(zll, zlh, zc, &f))
                {
                    yy[k] = yl + f * (yh - yl);
                    xx[k] = xl;
                    k++;
                }
                if (ctr_intersect(zlh, zhh, zc, &f))
                {
                    xx[k] = xl + f * (xh - xl);
                    yy[k] = yh;
                    k++;
                }
                if (ctr_intersect(zll, zhh, zc, &f))
                {
                    xx[k] = xl + f * (xh - xl);
                    yy[k] = yl + f * (yh - yl);
                    k++;
                }
                break;
            case 11:
                if (ctr_intersect(zhl, zhh, zc, &f))
                {
                    yy[k] = yl + f * (yh - yl);
                    xx[k] = xh;
                    k++;
                }
                if (ctr_intersect(zll, zhl, zc, &f))
                {
                    xx[k] = xl + f * (xh - xl);
                    yy[k] = yl;
                    k++;
                }
                if (ctr_intersect(zll, zhh, zc, &f))
                {
                    xx[k] = xl + f * (xh - xl);
                    yy[k] = yl + f * (yh - yl);
                    k++;
                }
                break;
            case 7:
                if (ctr_intersect(zll, zlh, zc, &f))
                {
                    yy[k] = yl + f * (yh - yl);
                    xx[k] = xl;
                    k++;
                }
                if (ctr_intersect(zll, zhl, zc, &f))
                {
                    xx[k] = xl + f * (xh - xl);
                    yy[k] = yl;
                    k++;
                }
                if (ctr_intersect(zlh, zhl, zc, &f))
                {
                    xx[k] = xl + f * (xh - xl);
                    yy[k] = yh + f * (yl - yh);
                    k++;
                }
                break;
            }

            /* We now have k(=2,4) endpoints */
            /* Decide which to join */

            seglist = NULL;

            if (k > 0)
            {
                if (k == 2)
                {
                    seglist = ctr_newseg(xx[0], yy[0], xx[1], yy[1], seglist);
                }
                else if (k == 4)
                {
                    for (k = 3; k >= 1; k--)
                    {
                        m = k;
                        xl = xx[k];
                        for (l = 0; l < k; l++)
                        {
                            if (xx[l] > xl)
                            {
                                xl = xx[l];
                                m = l;
                            }
                        }
                        if (m != k)
                        {
                            xl = xx[k];
                            yl = yy[k];
                            xx[k] = xx[m];
                            yy[k] = yy[m];
                            xx[m] = xl;
                            yy[m] = yl;
                        }
                    }
                    seglist = ctr_newseg(xx[0], yy[0], xx[1], yy[1], seglist);
                    seglist = ctr_newseg(xx[2], yy[2], xx[3], yy[3], seglist);
                }
                else
                    error("k != 2 or 4\n");
            }
            ctr_SegDB[i + j * nx] = seglist;
        }
    }

    /* The segment database is now assembled. */
    /* Begin following contours. */
    /* 1. Grab a segment */
    /* 2. Follow its tail */
    /* 3. Follow its head */
    /* 4. Draw the contour */

    for (i = 0; i < nx - 1; i++)
        for (j = 0; j < ny - 1; j++)
        {
            while ((seglist = ctr_SegDB[i + j * nx]))
            {
                ii = i;
                jj = j;
                start = end = seglist;
                ctr_SegDB[i + j * nx] = seglist->next;
                xend = seglist->x1;
                yend = seglist->y1;
                while ((dir = ctr_segdir(xend, yend, REAL(x), REAL(y), &ii, &jj, nx, ny)))
                {
                    ctr_SegDB[ii + jj * nx] = ctr_segupdate(xend, yend, dir, 1, ctr_SegDB[ii + jj * nx], &seg);
                    if (!seg)
                        break;
                    end->next = seg;
                    end = seg;
                    xend = end->x1;
                    yend = end->y1;
                }
                ii = i;
                jj = j;
                xend = seglist->x0;
                yend = seglist->y0;
                while ((dir = ctr_segdir(xend, yend, REAL(x), REAL(y), &ii, &jj, nx, ny)))
                {
                    ctr_SegDB[ii + jj * nx] = ctr_segupdate(xend, yend, dir, 0, ctr_SegDB[ii + jj * nx], &seg);
                    if (!seg)
                        break;
                    seg->next = start;
                    start = seg;
                    xend = start->x0;
                    yend = start->y0;
                }
                s = start;
                ns = 0;
                while (s)
                {
                    ns++;
                    s = s->next;
                }

                /* countour midpoint */
                /* use for labelling sometime */

                if (ns > 3)
                    ns2 = ns / 2;
                else
                    ns2 = -1;

                s = start;
                xxx = (double *)C_alloc(ns + 1, sizeof(double));
                yyy = (double *)C_alloc(ns + 1, sizeof(double));
                ns = 0;
                xxx[ns] = s->x0;
                yyy[ns++] = s->y0;
                while (s->next)
                {
                    s = s->next;
                    xxx[ns] = s->x0;
                    yyy[ns++] = s->y0;
                }
                xxx[ns] = s->x1;
                yyy[ns++] = s->y1;
                GMode(1, dd);
                GPolyline(ns, xxx, yyy, USER, dd);
                GMode(0, dd);
                C_free((char *)xxx);
                C_free((char *)yyy);
            }
        }
}

/* .Internal(contour(x,y,z, levels, col, lty)  */

SEXP do_contour(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP oargs, c, x, y, z, col, lty;
    int i, j, nx, ny, nc, ncol, nlty;
    int ltysave, colsave;
    double atom, zmin, zmax;
    char *vmax, *vmax0;
    DevDesc *dd = CurrentDevice();

    GCheckState(dd);

    if (length(args) < 4)
        errorcall(call, "too few arguments\n");

    oargs = args;

    x = CAR(args);
    internalTypeCheck(call, x, REALSXP);
    nx = LENGTH(x);
    args = CDR(args);

    y = CAR(args);
    internalTypeCheck(call, y, REALSXP);
    ny = LENGTH(y);
    args = CDR(args);

    z = CAR(args);
    internalTypeCheck(call, z, REALSXP);
    args = CDR(args);

    c = CAR(args);
    internalTypeCheck(call, c, REALSXP);
    nc = LENGTH(c);
    args = CDR(args);

    PROTECT(col = FixupCol(GetPar("col", args), dd));
    ncol = length(col);

    PROTECT(lty = FixupLty(GetPar("lty", args), dd));
    nlty = length(lty);

    /* col, lwd and lty vectors here --- FIXME: "lwd" ???? */

    if (nx < 2 || ny < 2)
        errorcall(call, "insufficient x or y values\n");

    if (nrows(z) != nx || ncols(z) != ny)
        errorcall(call, "dimension mismatch\n");

    if (nc < 1)
        errorcall(call, "no contour values\n");

    for (i = 0; i < nx; i++)
    {
        if (!FINITE(REAL(x)[i]))
            errorcall(call, "missing x values\n");
        if (i > 0 && REAL(x)[i] < REAL(x)[i - 1])
            errorcall(call, "increasing x values expected\n");
    }

    for (i = 0; i < ny; i++)
    {
        if (!FINITE(REAL(y)[i]))
            errorcall(call, "missing y values\n");
        if (i > 0 && REAL(y)[i] < REAL(y)[i - 1])
            errorcall(call, "increasing y values expected\n");
    }

    ctr_xtol = 1e-3 * fabs(REAL(x)[nx - 1] - REAL(x)[0]);
    ctr_ytol = 1e-3 * fabs(REAL(y)[ny - 1] - REAL(y)[0]);

    for (i = 0; i < nc; i++)
        if (!FINITE(REAL(c)[i]))
            errorcall(call, "illegal NA contour values\n");

    zmin = DBL_MAX;
    zmax = DBL_MIN;
    for (i = 0; i < nx * ny; i++)
        if (FINITE(REAL(z)[i]))
        {
            if (zmax < REAL(z)[i])
                zmax = REAL(z)[i];
            if (zmin > REAL(z)[i])
                zmin = REAL(z)[i];
        }

    if (zmin >= zmax)
    {
        if (zmin == zmax)
            warning("all z values are equal\n");
        else
            warning("all z values are NA\n");
        return R_NilValue;
    }

    /* PREVIOUSLY: atom = DBL_EPSILON * (zmax - zmin); */

    atom = 1e-3 * (zmax - zmin);

    /* Initialize the segment data base */
    /* Note we must be careful about resetting */
    /* the top of the stack, otherwise we run out of */
    /* memory after a sequence of displaylist replays */

    vmax0 = vmaxget();
    ctr_SegDB = (SEGP *)R_alloc(nx * ny, sizeof(SEGP));

    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            ctr_SegDB[i + j * nx] = NULL;

    /* Draw the contours -- note the heap release */

    ltysave = dd->gp.lty;
    colsave = dd->gp.col;
    for (i = 0; i < nc; i++)
    {
        vmax = vmaxget();
        dd->gp.lty = INTEGER(lty)[i % nlty];
        if (dd->gp.lty == NA_INTEGER)
            dd->gp.lty = ltysave;
        dd->gp.col = INTEGER(col)[i % ncol];
        if (dd->gp.col == NA_INTEGER)
            dd->gp.col = colsave;
        contour(x, nx, y, ny, z, REAL(c)[i], atom, dd);
        vmaxset(vmax);
    }
    vmaxset(vmax0);
    dd->gp.lty = ltysave;
    dd->gp.col = colsave;
    UNPROTECT(2);
    /* NOTE: only record operation if no "error"  */
    /* NOTE: on replay, call == R_NilValue */
    if (call != R_NilValue)
        recordGraphicOperation(op, oargs, dd);
    return R_NilValue;
}

/*  I m a g e   R e n d e r i n g  */

SEXP do_image(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP oargs, sx, sy, sz, szlim, sc;
    double *x, *y, *z;
    unsigned *c;
    double xlow, xhigh, ylow, yhigh, zmin, zmax;
    int i, j, nx, ny, nz, ic, nc, colsave, xpdsave;
    DevDesc *dd = CurrentDevice();

    GCheckState(dd);

    checkArity(op, args);
    oargs = args;

    sx = CAR(args);
    internalTypeCheck(call, sx, REALSXP);
    nx = LENGTH(sx);
    args = CDR(args);

    sy = CAR(args);
    internalTypeCheck(call, sy, REALSXP);
    ny = LENGTH(sy);
    args = CDR(args);

    sz = CAR(args);
    internalTypeCheck(call, sz, REALSXP);
    nz = length(sz);
    args = CDR(args);

    szlim = CAR(args);
    internalTypeCheck(call, szlim, REALSXP);
    if (length(szlim) != 2 || !FINITE(REAL(szlim)[0]) || !FINITE(REAL(szlim)[1]) || REAL(szlim)[0] >= REAL(szlim)[1])
        errorcall(call, "invalid z limits\n");
    zmin = REAL(szlim)[0];
    zmax = REAL(szlim)[1];
    args = CDR(args);

    PROTECT(sc = FixupCol(CAR(args), dd));
    nc = length(sc);

    /* Shorthand Pointers */

    x = REAL(sx);
    y = REAL(sy);
    z = REAL(sz);
    c = (unsigned *)INTEGER(sc);

    /* Check of grid coordinates */
    /* We want them to all be finite and in strictly ascending order */

    if (nx < 2 || ny < 2)
        goto badxy;
    if (!FINITE(x[0]))
        goto badxy;
    if (!FINITE(y[0]))
        goto badxy;
    for (i = 1; i < nx; i++)
        if (!FINITE(x[i]) || x[i] <= x[i - 1])
            goto badxy;
    for (j = 1; j < ny; j++)
        if (!FINITE(y[j]) || y[j] <= y[j - 1])
            goto badxy;

    colsave = dd->gp.col;
    xpdsave = dd->gp.xpd;
    dd->gp.xpd = 0;

    GMode(1, dd);

    for (i = 0; i < nx; i++)
    {
        if (i == 0)
            xlow = x[0];
        else
            xlow = 0.5 * (x[i] + x[i - 1]);
        if (i == nx - 1)
            xhigh = x[nx - 1];
        else
            xhigh = 0.5 * (x[i] + x[i + 1]);

        for (j = 0; j < ny; j++)
        {
            if (FINITE(z[i + j * nx]))
            {
                ic = floor((nc - 1) * (z[i + j * nx] - zmin) / (zmax - zmin) + 0.5);
                if (ic >= 0 && ic < nc)
                {
                    if (j == 0)
                        ylow = y[0];
                    else
                        ylow = 0.5 * (y[j] + y[j - 1]);
                    if (j == ny - 1)
                        yhigh = y[ny - 1];
                    else
                        yhigh = 0.5 * (y[j] + y[j + 1]);
                    GRect(xlow, ylow, xhigh, yhigh, USER, c[ic], NA_INTEGER, dd);
                }
            }
        }
    }
    GMode(0, dd);
    dd->gp.col = colsave;
    dd->gp.xpd = xpdsave;
    R_Visible = 0;
    UNPROTECT(1);
    if (call != R_NilValue)
        recordGraphicOperation(op, oargs, dd);
    return R_NilValue;

badxy:
    errorcall(call, "invalid x / y limits\n");
    return R_NilValue; /* never used; to keep -Wall happy */
}

/*  P e r s p e c t i v e   S u r f a c e   P l o t s  */

/* Conversion of degrees to radians */

#define DegToRad(x) (0.01745329251994329576 * x)

/* Definitions of data structures for vectors and */
/* transformations in homogeneous 3d coordinates */

typedef double Vector3d[4];
typedef double Trans3d[4][4];

/* The viewing transformation matrix. */

static SEXP gcall;
static Trans3d VT;

#ifdef NOT_used_currently /*-- out 'def'  (-Wall) --*/
static void MakeVector(double x, double y, double z, Vector3d v)
{
    v[0] = x;
    v[1] = y;
    v[2] = z;
    v[3] = 1;
}
#endif

static void TransVector(Vector3d u, Trans3d T, Vector3d v)
{
    double sum;
    int i, j;

    for (i = 0; i < 4; i++)
    {
        sum = 0;
        for (j = 0; j < 4; j++)
            sum = sum + u[j] * T[j][i];
        v[i] = sum;
    }
}

static void Accumulate(Trans3d T)
{
    Trans3d U;
    double sum;
    int i, j, k;

    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            sum = 0;
            for (k = 0; k < 4; k++)
                sum = sum + VT[i][k] * T[k][j];
            U[i][j] = sum;
        }
    }
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            VT[i][j] = U[i][j];
}

static void SetToIdentity(Trans3d T)
{
    int i, j;
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
            T[i][j] = 0;
        T[i][i] = 1;
    }
}

static void Translate(double x, double y, double z)
{
    Trans3d T;
    SetToIdentity(T);
    T[3][0] = x;
    T[3][1] = y;
    T[3][2] = z;
    Accumulate(T);
}

static void Scale(double x, double y, double z)
{
    Trans3d T;
    SetToIdentity(T);
    T[0][0] = x;
    T[1][1] = y;
    T[2][2] = z;
    Accumulate(T);
}

static void XRotate(double angle)
{
    double c, s;
    Trans3d T;
    SetToIdentity(T);
    c = cos(DegToRad(angle));
    s = sin(DegToRad(angle));
    T[1][1] = c;
    T[2][1] = -s;
    T[2][2] = c;
    T[1][2] = s;
    Accumulate(T);
}

static void YRotate(double angle)
{
    double c, s;
    Trans3d T;
    SetToIdentity(T);
    c = cos(DegToRad(angle));
    s = sin(DegToRad(angle));
    T[0][0] = c;
    T[2][0] = s;
    T[2][2] = c;
    T[0][2] = -s;
    Accumulate(T);
}

static void Perspective(double d)
{
    Trans3d T;

    SetToIdentity(T);
    T[2][3] = -1 / d;
    Accumulate(T);
}

/* Determine the depth ordering of the facets to ensure */
/* that they are drawn in an occlusion compatible order. */

void OrderFacets(double *depth, int *index, int n)
{
    int i, j, h;
    int itmp;

    h = 1;
    do
    {
        h = 3 * h + 1;
    } while (h <= n);

    do
    {
        h = h / 3;
        for (i = h; i < n; i++)
        {
            itmp = index[i];
            j = i;
            while (depth[index[j - h]] < depth[itmp])
            {
                index[j] = index[j - h];
                j = j - h;
                if (j < h)
                    goto next_h;
            }
        next_h:
            index[j] = itmp;
        }
    } while (h != 1);
}

/* For each facet, determine the farthest point from the eye. */
/* Sorting the facets so that these depths are decreasing */
/* yields an occlusion compatible ordering. */
/* Note that we ignore z values when doing this. */

static void DepthOrder(double *z, double *x, double *y, int nx, int ny, double *depth, int *index)
{
    int i, ii, j, jj, nx1, ny1;
    Vector3d u, v;
    double d;
    nx1 = nx - 1;
    ny1 = ny - 1;
    for (i = 0; i < nx1 * ny1; i++)
        index[i] = i;
    for (i = 0; i < nx1; i++)
        for (j = 0; j < ny1; j++)
        {
            d = -DBL_MAX;
            for (ii = 0; ii <= 1; ii++)
                for (jj = 0; jj <= 1; jj++)
                {
                    u[0] = x[i + ii];
                    u[1] = y[j + jj];
                    /* Originally I had the following line here: */
                    /* u[2] = z[i+ii+(j+jj)*nx]; */
                    /* But this leads to artifacts. */
                    /* It has been replaced by the following line: */
                    u[2] = 0;
                    u[3] = 1;
                    if (FINITE(u[0]) && FINITE(u[1]) && FINITE(u[2]))
                    {
                        TransVector(u, VT, v);
                        if (v[3] > d)
                            d = v[3];
                    }
                }
            depth[i + j * nx1] = d;
        }
    OrderFacets(depth, index, nx1 * ny1);
}

static void DrawFacets(double *z, double *x, double *y, int nx, int ny, int *index, int *col, int ncol, int border,
                       double *shade)
{
    double xx[4], yy[4];
    Vector3d u, v;
    int i, j, k, n, nx1, ny1, icol, nv;
    DevDesc *dd;
    dd = CurrentDevice();
    nx1 = nx - 1;
    ny1 = ny - 1;
    n = nx1 * ny1;
    for (k = 0; k < n; k++)
    {
        nv = 0;
        i = index[k] % nx1;
        j = index[k] / nx1;
        icol = (i + j * nx1) % ncol;

        u[0] = x[i];
        u[1] = y[j];
        u[2] = z[i + j * nx];
        u[3] = 1;
        if (FINITE(u[0]) && FINITE(u[1]) && FINITE(u[2]))
        {
            TransVector(u, VT, v);
            xx[nv] = v[0] / v[3];
            yy[nv] = v[1] / v[3];
            nv++;
        }

        u[0] = x[i + 1];
        u[1] = y[j];
        u[2] = z[i + 1 + j * nx];
        u[3] = 1;
        if (FINITE(u[0]) && FINITE(u[1]) && FINITE(u[2]))
        {
            TransVector(u, VT, v);
            xx[nv] = v[0] / v[3];
            yy[nv] = v[1] / v[3];
            nv++;
        }

        u[0] = x[i + 1];
        u[1] = y[j + 1];
        u[2] = z[i + 1 + (j + 1) * nx];
        u[3] = 1;
        if (FINITE(u[0]) && FINITE(u[1]) && FINITE(u[2]))
        {
            TransVector(u, VT, v);
            xx[nv] = v[0] / v[3];
            yy[nv] = v[1] / v[3];
            nv++;
        }

        u[0] = x[i];
        u[1] = y[j + 1];
        u[2] = z[i + (j + 1) * nx];
        u[3] = 1;
        if (FINITE(u[0]) && FINITE(u[1]) && FINITE(u[2]))
        {
            TransVector(u, VT, v);
            xx[nv] = v[0] / v[3];
            yy[nv] = v[1] / v[3];
            nv++;
        }

        if (nv > 2)
        {
            unsigned int newcol = col[icol];
            if (shade)
            {
                unsigned int red, green, blue;
                double shadeval = shade[i + j * (nx - 1)];
                red = shadeval * R_RED(newcol);
                green = shadeval * R_GREEN(newcol);
                blue = shadeval * R_BLUE(newcol);
                newcol = R_RGB(red, green, blue);
            }
            GPolygon(nv, xx, yy, USER, newcol, border, dd);
        }
    }
}

#ifdef NOT_used_currently /*-- out 'def'  (-Wall) --*/
static void CheckRange(double *x, int n, double min, double max)
{
    double xmin, xmax;
    int i;
    xmin = DBL_MAX;
    xmax = -DBL_MAX;
    for (i = 0; i < n; i++)
        if (FINITE(x[i]))
        {
            if (x[i] < xmin)
                xmin = x[i];
            if (x[i] > xmax)
                xmax = x[i];
        }
    if (xmin < min || xmax > max)
        errorcall(gcall, "coordinates outsize specified range\n");
}
#endif

static void PerspWindow(double *xlim, double *ylim, double *zlim, DevDesc *dd)
{
    double pin1, pin2, scale, xdelta, ydelta, xscale, yscale, xadd, yadd;
    double xmax, xmin, ymax, ymin, xx, yy;
    Vector3d u, v;
    int i, j, k;

    xmax = xmin = ymax = ymin = 0;
    u[3] = 1;
    for (i = 0; i < 2; i++)
    {
        u[0] = xlim[i];
        for (j = 0; j < 2; j++)
        {
            u[1] = ylim[j];
            for (k = 0; k < 2; k++)
            {
                u[2] = zlim[k];
                TransVector(u, VT, v);
                xx = v[0] / v[3];
                yy = v[1] / v[3];
                if (xx > xmax)
                    xmax = xx;
                if (xx < xmin)
                    xmin = xx;
                if (yy > ymax)
                    ymax = yy;
                if (yy < ymin)
                    ymin = yy;
            }
        }
    }
    pin1 = GConvertXUnits(1.0, NPC, INCHES, dd);
    pin2 = GConvertYUnits(1.0, NPC, INCHES, dd);
    xdelta = fabs(xmax - xmin);
    ydelta = fabs(ymax - ymin);
    xscale = pin1 / xdelta;
    yscale = pin2 / ydelta;
    scale = (xscale < yscale) ? xscale : yscale;
    xadd = .5 * (pin1 / scale - xdelta);
    yadd = .5 * (pin2 / scale - ydelta);
    GScale(xmin - xadd, xmax + xadd, 1, dd);
    GScale(ymin - yadd, ymax + yadd, 2, dd);
    GMapWin2Fig(dd);
}

static int LimitCheck(double *lim, double *c, double *s)
{
    if (!FINITE(lim[0]) || !FINITE(lim[1]) || lim[0] >= lim[1])
        return 0;
    *s = 0.5 * fabs(lim[1] - lim[0]);
    *c = 0.5 * (lim[1] + lim[0]);
    return 1;
}

/* PerspBox: The following code carries out a visibility test */
/* on the surfaces of the xlim/ylim/zlim box around the plot. */
/* If front = 0, only the faces with their inside toward the */
/* eyepoint are drawn.  If front = 1, only the faces with */
/* their outside toward the eye are drawn.  This lets us carry */
/* out hidden line removal by drawing any faces which will be */
/* obscured before the surface, and those which will not be */
/* obscured after the surface. */

static int Vertex[8][3] = {
    {0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1},
};

static int Face[6][4] = {
    {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 2, 3, 1}, {4, 5, 7, 6}, {0, 4, 6, 2}, {1, 3, 7, 5},
};

static void PerspBox(int front, double *x, double *y, double *z, DevDesc *dd)
{
    Vector3d u0, v0, u1, v1, u2, v2, u3, v3;
    double d[3], e[3];
    int f, i, p0, p1, p2, p3, near;
    for (f = 0; f < 6; f++)
    {
        p0 = Face[f][0];
        p1 = Face[f][1];
        p2 = Face[f][2];
        p3 = Face[f][3];

        u0[0] = x[Vertex[p0][0]];
        u0[1] = y[Vertex[p0][1]];
        u0[2] = z[Vertex[p0][2]];
        u0[3] = 1;
        u1[0] = x[Vertex[p1][0]];
        u1[1] = y[Vertex[p1][1]];
        u1[2] = z[Vertex[p1][2]];
        u1[3] = 1;
        u2[0] = x[Vertex[p2][0]];
        u2[1] = y[Vertex[p2][1]];
        u2[2] = z[Vertex[p2][2]];
        u2[3] = 1;
        u3[0] = x[Vertex[p3][0]];
        u3[1] = y[Vertex[p3][1]];
        u3[2] = z[Vertex[p3][2]];
        u3[3] = 1;

        TransVector(u0, VT, v0);
        TransVector(u1, VT, v1);
        TransVector(u2, VT, v2);
        TransVector(u3, VT, v3);

        /* Visibility test */
        /* Determine whether the surface normal is toward the eye. */

        for (i = 0; i < 3; i++)
        {
            d[i] = v1[i] / v1[3] - v0[i] / v0[3];
            e[i] = v2[i] / v2[3] - v1[i] / v1[3];
        }
        near = (d[0] * e[1] - d[1] * e[0]) < 0;

        if ((front && near) || (!front && !near))
        {
            GLine(v0[0] / v0[3], v0[1] / v0[3], v1[0] / v1[3], v1[1] / v1[3], USER, dd);
            GLine(v1[0] / v1[3], v1[1] / v1[3], v2[0] / v2[3], v2[1] / v2[3], USER, dd);
            GLine(v2[0] / v2[3], v2[1] / v2[3], v3[0] / v3[3], v3[1] / v3[3], USER, dd);
            GLine(v3[0] / v3[3], v3[1] / v3[3], v0[0] / v0[3], v0[1] / v0[3], USER, dd);
        }
    }
}

SEXP do_persp(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP x, y, z, xlim, ylim, zlim;
    SEXP depth, index, originalArgs;
    SEXP col, border, shade;
    double theta, phi, r, d, expand, xc, yc, zc, xs, ys, zs;
    int i, j, scale, ncol;
    DevDesc *dd;

    if (length(args) < 12)
        errorcall(call, "too few parameters\n");
    gcall = call;
    originalArgs = args;

    PROTECT(x = coerceVector(CAR(args), REALSXP));
    if (length(x) < 2)
        errorcall(call, "invalid x argument\n");
    args = CDR(args);

    PROTECT(y = coerceVector(CAR(args), REALSXP));
    if (length(y) < 2)
        errorcall(call, "invalid y argument\n");
    args = CDR(args);

    PROTECT(z = coerceVector(CAR(args), REALSXP));
    if (!isMatrix(z) || nrows(z) != length(x) || ncols(z) != length(y))
        errorcall(call, "invalid z argument\n");
    args = CDR(args);

    PROTECT(xlim = coerceVector(CAR(args), REALSXP));
    if (length(xlim) != 2)
        errorcall(call, "invalid xlim argument\n");
    args = CDR(args);

    PROTECT(ylim = coerceVector(CAR(args), REALSXP));
    if (length(ylim) != 2)
        errorcall(call, "invalid ylim argument\n");
    args = CDR(args);

    PROTECT(zlim = coerceVector(CAR(args), REALSXP));
    if (length(zlim) != 2)
        errorcall(call, "invalid zlim argument\n");
    args = CDR(args);

    /* Checks on x/y/z Limits */

    if (!LimitCheck(REAL(xlim), &xc, &xs))
        errorcall(call, "invalid x limits\n");
    if (!LimitCheck(REAL(ylim), &yc, &ys))
        errorcall(call, "invalid y limits\n");
    if (!LimitCheck(REAL(zlim), &zc, &zs))
        errorcall(call, "invalid z limits\n");

    theta = asReal(CAR(args));
    args = CDR(args);

    phi = asReal(CAR(args));
    args = CDR(args);

    r = asReal(CAR(args));
    args = CDR(args);

    d = asReal(CAR(args));
    args = CDR(args);

    scale = asLogical(CAR(args));
    args = CDR(args);

    expand = asReal(CAR(args));
    args = CDR(args);

    PROTECT(col = FixupCol(CAR(args), dd));
    ncol = LENGTH(col);
    if (ncol < 1)
        errorcall(call, "invalid col specification\n");
    args = CDR(args);

    PROTECT(border = FixupCol(CAR(args), dd));
    if (length(border) < 1)
        errorcall(call, "invalid border specification\n");
    args = CDR(args);

    if (!isNull(CAR(args)))
    {
        PROTECT(shade = coerceVector(CAR(args), REALSXP));
        if (!isMatrix(shade) || nrows(shade) != length(x) - 1 || ncols(shade) != length(y) - 1)
            errorcall(call, "invalid shade argument\n");
    }
    else
        PROTECT(shade = R_NilValue);
    args = CDR(args);

    if (!scale)
    {
        double s;
        s = xs;
        if (s < ys)
            s = ys;
        if (s < zs)
            s = zs;
        xs = s;
        ys = s;
        zs = s;
    }

    /* Parameter Checks */

    if (!FINITE(theta) || !FINITE(phi) || !FINITE(r) || !FINITE(d) || d < 0 || r < 0)
        errorcall(call, "invalid viewing parameters\n");
    if (!FINITE(expand) || expand < 0)
        errorcall(call, "invalid expand value\n");
    if (scale == NA_LOGICAL)
        scale = 0;

    dd = GNewPlot(call != R_NilValue, NA_LOGICAL);
    GSetState(1, dd);
    GSavePars(dd);
    ProcessInlinePars(args, dd);
    if (length(border) > 1)
        dd->gp.fg = INTEGER(border)[0];
    dd->gp.xlog = 0;
    dd->gp.ylog = 0;

    /* Specify the viewing transformation. */

    SetToIdentity(VT);                  /* Initialization */
    Translate(-xc, -yc, -zc);           /* center at the origin */
    Scale(1 / xs, 1 / ys, expand / zs); /* scale extents to [-1,1] */
    XRotate(-90.0);                     /* rotate x-y plane to horizontal */
    YRotate(-theta);                    /* azimuthal rotation */
    XRotate(phi);                       /* elevation rotation */
    Translate(0.0, 0.0, -r - d);        /* translate the eyepoint to the origin */
    Perspective(d);                     /* perspective */

    /* Specify the plotting window. */
    /* Here we map the vertices of the cube */
    /* [xmin,xmax]*[ymin,ymax]*[zmin,zmax] */
    /* to the screen and then chose a window */
    /* which is symmetric about (0,0). */

    PerspWindow(REAL(xlim), REAL(ylim), REAL(zlim), dd);

    /* Compute facet order. */

    PROTECT(depth = allocVector(REALSXP, (nrows(z) - 1) * (ncols(z) - 1)));
    PROTECT(index = allocVector(INTSXP, (nrows(z) - 1) * (ncols(z) - 1)));
    DepthOrder(REAL(z), REAL(x), REAL(y), nrows(z), ncols(z), REAL(depth), INTEGER(index));

    /* Now we order the facets by depth */
    /* and then draw them back to front. */
    /* This is the "painters" algorithm. */

    PerspBox(0, REAL(xlim), REAL(ylim), REAL(zlim), dd);

    DrawFacets(REAL(z), REAL(x), REAL(y), nrows(z), ncols(z), INTEGER(index), INTEGER(col), ncol, INTEGER(border)[0],
               isNull(shade) ? NULL : REAL(shade));

    PerspBox(1, REAL(xlim), REAL(ylim), REAL(zlim), dd);

    GRestorePars(dd);
    UNPROTECT(11);
    if (call != R_NilValue)
        recordGraphicOperation(op, originalArgs, dd);

    PROTECT(x = allocVector(REALSXP, 16));
    PROTECT(y = allocVector(INTSXP, 2));
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
        {
            REAL(x)[i + j * 4] = VT[i][j];
        }
    INTEGER(y)[0] = 4;
    INTEGER(y)[1] = 4;
    setAttrib(x, R_DimSymbol, y);
    UNPROTECT(2);
    return x;
}

/* .Internal(shade(x, y, z, xlim, ylim, zlim,      */
/*                 theta, phi, r, d, scale, ...))  */

SEXP do_shade(SEXP call, SEXP op, SEXP args, SEXP env)
{
    SEXP x, y, z, xlim, ylim, zlim;
    SEXP shading, t;
    int nx, ny;
    int i, j, k;
    double i_red, i_green, i_blue;
    double normalX, normalY, normalZ, sum;
    double NdotL, NdotL_Y, NdotL_Z;
    double xl, xh, yl, yh, topLeft, topRight, bottomLeft, bottomRight;
    double v1x, v1y, v1z, v2x, v2y, v2z;
    double theta, phi, r, d, xc, yc, zc, xs, ys, zs;
    double ambient, diffuse;
    Vector3d light, u;

    PROTECT(x = coerceVector(CAR(args), REALSXP));
    if (length(x) < 2)
        errorcall(call, "invalid x argument\n");
    nx = LENGTH(x);
    args = CDR(args);

    PROTECT(y = coerceVector(CAR(args), REALSXP));
    if (length(y) < 2)
        errorcall(call, "invalid y argument\n");
    ny = LENGTH(y);
    args = CDR(args);

    PROTECT(z = coerceVector(CAR(args), REALSXP));
    if (!isMatrix(z) || nrows(z) != length(x) || ncols(z) != length(y))
        errorcall(call, "invalid z argument\n");
    args = CDR(args);

    PROTECT(xlim = coerceVector(CAR(args), REALSXP));
    if (length(xlim) != 2)
        errorcall(call, "invalid xlim argument\n");
    args = CDR(args);

    PROTECT(ylim = coerceVector(CAR(args), REALSXP));
    if (length(ylim) != 2)
        errorcall(call, "invalid ylim argument\n");
    args = CDR(args);

    PROTECT(zlim = coerceVector(CAR(args), REALSXP));
    if (length(zlim) != 2)
        errorcall(call, "invalid zlim argument\n");
    args = CDR(args);

    /* Checks on x/y/z Limits */

    if (!LimitCheck(REAL(xlim), &xc, &xs))
        errorcall(call, "invalid x limits\n");
    if (!LimitCheck(REAL(ylim), &yc, &ys))
        errorcall(call, "invalid y limits\n");
    if (!LimitCheck(REAL(zlim), &zc, &zs))
        errorcall(call, "invalid z limits\n");

    theta = asReal(CAR(args));
    args = CDR(args);

    phi = asReal(CAR(args));
    args = CDR(args);

    ambient = asReal(CAR(args));
    args = CDR(args);

    diffuse = asReal(CAR(args));
    args = CDR(args);

    /* Calculate light source direction */
    SetToIdentity(VT); /* Initialization */
    XRotate(-90.0);    /* rotate x-y plane to horizontal */
    YRotate(-theta);   /* azimuthal rotation */
    XRotate(phi);      /* elevation rotation */

    u[0] = 0;
    u[1] = 0;
    u[2] = -1;
    u[3] = 1;

    TransVector(u, VT, light);

    PROTECT(shading = allocVector(REALSXP, (nx - 1) * (ny - 1)));
    PROTECT(t = allocVector(INTSXP, 2));
    INTEGER(t)[0] = nx - 1;
    INTEGER(t)[1] = ny - 1;
    setAttrib(shading, R_DimSymbol, t);

    for (i = 0; i < nx - 1; i++)
        for (j = 0; j < ny - 1; j++)
        {
            k = i + j * nx;
            bottomLeft = REAL(z)[k];
            bottomRight = REAL(z)[k + 1];
            topLeft = REAL(z)[k + nx];
            topRight = REAL(z)[k + nx + 1];

            xl = REAL(x)[i];
            xh = REAL(x)[i + 1];
            yl = REAL(y)[j];
            yh = REAL(y)[j + 1];

            /* (xl,yl, bottomLeft) */
            /* (xl,yh, topLeft) */
            /* (xh,yh, topRight) */
            /* (xh,yl, bottomRight) */

            v1x = (xh - xl); // xbr - xtl
            v1y = (yl - yh); // ybr - ytl
            v1z = (bottomRight - topLeft);
            v2x = (xh - xl); // xtr - xbl
            v2y = (yh - yl); // ytr - ybl
            v2z = (topRight - bottomLeft);

            normalX = v1y * v2z - v1z * v2y;
            normalY = v1z * v2x - v1x * v2z;
            normalZ = v1x * v2y - v1y * v2x;

            sum = sqrt(normalX * normalX + normalY * normalY + normalZ * normalZ);

            if (sum == 0)
                sum = 1;
            normalX /= sum;
            normalY /= sum;
            normalZ /= sum;

            NdotL = 0.5 * (normalX * light[0] + normalY * light[1] + normalZ * light[2] + 1);

            /* NdotL_Y = normalY * light[1]; */
            /* NdotL_Z = normalZ * light[2]; */

            /* printf("%f %f %f\n", NdotL_X, NdotL_Y, NdotL_Z); */

            i_red = ambient + diffuse * NdotL;
            /* i_green = ambient + diffuse * NdotL_Y; */
            /* i_blue = ambient + diffuse * NdotL_Z; */

            /* printf("%f %f %f\n", i_red, i_green, i_blue); */
            /* REAL(shading)[k] = (i_red << 16) + (i_green << 8) + i_blue; */
            REAL(shading)[i + j * (nx - 1)] = (i_red);
            /* + i_green + i_blue) / 3.0; */
        }
    UNPROTECT(8);
    return shading;
}
