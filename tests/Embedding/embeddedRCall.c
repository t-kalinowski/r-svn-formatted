#include "embeddedRCall.h"

int eval_R_command(const char *funcName, int argc, char *argv[])
{
    SEXP e;
    SEXP arg;

    int i;
    int errorOccurred;
    init_R(argc, argv);

    PROTECT(arg = allocVector(INTSXP, 10));
    for (i = 0; i < LENGTH(arg); i++)
        INTEGER(arg)[i] = i + 1;

    PROTECT(e = lang2(install((char *)funcName), arg));

    /* Evaluate the call to the R function.
       Ignore the return value.
    */
    R_tryEval(e, R_GlobalEnv, &errorOccurred);

    UNPROTECT(2);
    return (0);
}

extern int Rf_initEmbeddedR(int argc, char *argv[]);

void init_R(int argc, char **argv)
{
    int defaultArgc = 1;
    char *defaultArgv[] = {"Rtest"};

    if (argc == 0 || argv == NULL)
    {
        argc = defaultArgc;
        argv = defaultArgv;
    }
    Rf_initEmbeddedR(argc, argv);
}
