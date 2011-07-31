/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2011   The R Development Core Team.
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

#include <R.h>
#include <Rinternals.h>
#include <stdint.h>

typedef uint_least64_t Uint64;

static Uint64 A1p76[3][3] = {
    {82758667, 1871391091, 4127413238}, {3672831523, 69195019, 1871391091}, {3672091415, 3528743235, 69195019}};

static Uint64 A2p76[3][3] = {
    {1511326704, 3759209742, 1610795712}, {4292754251, 1511326704, 3889917532}, {3859662829, 4292754251, 3708466080}};

static Uint64 A1p127[3][3] = {
    {2427906178, 3580155704, 949770784}, {226153695, 1230515664, 3580155704}, {1988835001, 986791581, 1230515664}};

static Uint64 A2p127[3][3] = {
    {1464411153, 277697599, 1610723613}, {32183930, 1464411153, 1022607788}, {2824425944, 32183930, 2093834863}};

SEXP nextStream(SEXP x)
{
    Uint64 seed[6];
    for (int i = 0; i < 6; i++)
        seed[i] = (unsigned int)INTEGER(x)[i + 1];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            seed[i] += A1p127[i][j] * seed[j];
            seed[j] %= 4294967087;
        }
    for (int i = 3; i < 6; i++)
        for (int j = 3; j < 6; j++)
        {
            seed[i] += A2p127[i][j] * seed[j];
            seed[j] %= 4294944443;
        }
    SEXP ans = allocVector(INTSXP, 7);
    INTEGER(ans)[0] = INTEGER(x)[0];
    for (int i = 0; i < 6; i++)
        INTEGER(ans)[i + 1] = seed[i];
    return ans;
}

SEXP nextSubStream(SEXP x)
{
    Uint64 seed[6];
    for (int i = 0; i < 6; i++)
        seed[i] = (unsigned int)INTEGER(x)[i + 1];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            seed[i] += A1p76[i][j] * seed[j];
            seed[j] %= 4294967087;
        }
    for (int i = 3; i < 6; i++)
        for (int j = 3; j < 6; j++)
        {
            seed[i] += A2p76[i][j] * seed[j];
            seed[j] %= 4294944443;
        }
    SEXP ans = allocVector(INTSXP, 6);
    INTEGER(ans)[0] = INTEGER(x)[0];
    for (int i = 0; i < 6; i++)
        INTEGER(ans)[i + 1] = seed[i];
    return ans;
}

#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
    {"nextStream", (DL_FUNC)&nextStream, 1}, {"nextSubStream", (DL_FUNC)&nextSubStream, 1}, {NULL, NULL, 0}};

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
    __attribute__((visibility("default")))
#endif
    R_init_tools(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
