#include "embeddedRCall.h"

int callLength(SEXP obj);
int R_embeddedShutdown(Rboolean ask);

main(int argc, char *argv[])
{
    SEXP objs[100];
    int i;
    Rf_initEmbeddedR(sizeof(argv) / sizeof(argv[0]), argv);

    for (i = 0; i < 100; i++)
    {
        objs[i] = allocVector(VECSXP, 1000);
        R_PreserveObject(objs[i]);
        callLength(objs[i]);
    }

    R_embeddedShutdown(FALSE);
}

int callLength(SEXP obj)
{
    SEXP e, val;
    int errorOccurred;
    int len = -1;

    PROTECT(e = allocVector(LANGSXP, 2));
    SETCAR(e, Rf_install("length"));
    SETCAR(CDR(e), obj);

    PROTECT(val = Test_tryEval(e, &errorOccurred));
    len = INTEGER(val)[0];
    UNPROTECT(2);

    return (len);
}

int R_embeddedShutdown(Rboolean ask)
{

    R_dot_Last();
    CleanEd();
    KillAllDevices();
    num_old_gens_to_collect = NUM_OLD_GENERATIONS;
    R_gc();
    return (1);
}
