/* ansari.c
   Compute the exact distribution of the Ansari-Bradley test statistic.
   */

#include <S.h>
#include <Mathlib.h>

#include "ctest.h"

static double ***w;

static void errmsg(char *s)
{
    PROBLEM "%s", s RECOVER(NULL_ENTRY);
}

static void w_init(Sint m, Sint n)
{
    Sint i;

    w = Calloc(m + 1, double **);
    if (!w)
        errmsg("allocation error 1 in `ansari.c'");
    for (i = 0; i <= m; i++)
    {
        w[i] = Calloc(n + 1, double *);
        if (!w[i])
            errmsg("allocation error 2 in `ansari.c'");
    }
}

static void w_free(Sint m, Sint n)
{
    Sint i;
    for (i = m; i >= 0; i--)
    {
        Free(w[i]);
    }
    Free(w);
    w = 0;
}

static double cansari(int k, int m, int n)
{
    int i, l, u;

    l = (m + 1) * (m + 1) / 4;
    u = l + m * n / 2;

    if ((k < l) || (k > u))
        return (0);

    if (w[m][n] == 0)
    {
        w[m][n] = Calloc(u + 1, double);
        if (!w[m][n])
            errmsg("allocation error in cansari()");
        for (i = 0; i <= u; i++)
            w[m][n][i] = -1;
    }

    if (w[m][n][k] < 0)
    {
        if (m == 0)
            w[m][n][k] = (k == 0);
        else if (n == 0)
            w[m][n][k] = (k == l);
        else
            w[m][n][k] = cansari(k, m, n - 1) + cansari(k - (m + n + 1) / 2, m - 1, n);
    }

    return (w[m][n][k]);
}

void dansari(Sint *len, double *x, Sint *m, Sint *n)
{
    Sint i;

    w_init(*m, *n);
    for (i = 0; i < *len; i++)
        if (fabs(x[i] - floor(x[i] + 0.5)) > 1e-7)
        {
            x[i] = 0;
        }
        else
        {
            x[i] = cansari((Sint)x[i], (Sint)*m, (Sint)*n) / choose(*m + *n, *m);
        }
    w_free(*m, *n);
}

void pansari(Sint *len, double *x, Sint *m, Sint *n)
{
    Sint i, j, l, u;
    double p, q;

    w_init(*m, *n);
    l = (*m + 1) * (*m + 1) / 4;
    u = l + *m * *n / 2;
    for (i = 0; i < *len; i++)
    {
        q = floor(x[i] + 1e-7);
        if (q < l)
            x[i] = 0;
        else if (q > u)
            x[i] = 1;
        else
        {
            p = 0;
            for (j = l; j <= q; j++)
            {
                p += cansari((Sint)j, (Sint)*m, (Sint)*n);
            }
            x[i] = p / choose(*m + *n, *m);
        }
    }
    w_free(*m, *n);
}
