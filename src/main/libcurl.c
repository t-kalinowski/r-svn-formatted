/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2015 The R Core Team
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

#include <Defn.h>
#include <Internal.h>

#ifdef HAVE_CURL_CURL_H
#include <curl/curl.h>
#endif

SEXP attribute_hidden do_curlVersion(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    checkArity(op, args);
    SEXP ans = PROTECT(allocVector(STRSXP, 1));
#ifdef HAVE_CURL_CURL_H
    curl_version_info_data *d = curl_version_info(CURLVERSION_NOW);
    SET_STRING_ELT(ans, 0, mkChar(d->version));
    setAttrib(ans, install("ssl_version"), mkString(d->ssl_version ? d->ssl_version : "none"));
    setAttrib(ans, install("libssh_version"), mkString(((d->age >= 3) && d->libssh_version) ? d->libssh_version : ""));
    const char *const *p;
    int n, i;
    for (p = d->protocols, n = 0; *p; p++, n++)
        ;
    SEXP protocols = PROTECT(allocVector(STRSXP, n));
    for (p = d->protocols, i = 0; i < n; i++, p++)
        SET_STRING_ELT(protocols, i, mkChar(*p));
    setAttrib(ans, install("protocols"), protocols);
    UNPROTECT(1);
#else
    SET_STRING_ELT(ans, 0, mkChar(""));
#endif
    UNPROTECT(1);
    return ans;
}

static char headers[100][2048];
static int used;

static size_t rcvHeaders(void *buffer, size_t size, size_t nmemb, void *userp)
{
    char *d = (char *)buffer;
    size_t result = size * nmemb, res = result > 2048 ? 2048 : result;
    if (used > 100)
        return result;
    strncpy(headers[used], d, res);
    headers[used][res] = '\0';
    used++;
    return result;
}

static size_t rcvData(void *buffer, size_t size, size_t nmemb, void *userp)
{
    return size * nmemb;
}

SEXP attribute_hidden do_curlGetHeaders(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    checkArity(op, args);
#ifndef HAVE_CURL_CURL_H
    error("curlGetHeaders is not supported on this platform");
    return R_NilValue;
#else
    if (!isString(CAR(args)) || LENGTH(CAR(args)) != 1)
        error("invalid %s argument", "url");
    const char *url = translateChar(STRING_ELT(CAR(args), 0));
    used = 0;
    int redirect = asLogical(CADDR(args));
    if (redirect == NA_LOGICAL)
        error("invalid %s argument", "redirect");

    CURL *hnd = curl_easy_init();
    curl_easy_setopt(hnd, CURLOPT_URL, url);
    curl_easy_setopt(hnd, CURLOPT_NOPROGRESS, 1L);
    curl_easy_setopt(hnd, CURLOPT_NOBODY, 1L);
    curl_easy_setopt(hnd, CURLOPT_HEADER, 1L);
    curl_easy_setopt(hnd, CURLOPT_HEADERFUNCTION, &rcvHeaders);
    curl_easy_setopt(hnd, CURLOPT_WRITEHEADER, &headers);
    curl_easy_setopt(hnd, CURLOPT_WRITEFUNCTION, &rcvData); // unused
    const char *ua = translateChar(STRING_ELT(CADR(args), 0));
    curl_easy_setopt(hnd, CURLOPT_USERAGENT, ua);
    if (redirect)
        curl_easy_setopt(hnd, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(hnd, CURLOPT_MAXREDIRS, 50L);
    char errbuf[CURL_ERROR_SIZE];
    curl_easy_setopt(hnd, CURLOPT_ERRORBUFFER, errbuf);
    CURLcode ret = curl_easy_perform(hnd);
    if (ret != CURLE_OK)
        error("libcurl error code %d\n\t%s\n", ret, errbuf);
    long http_code = 0;
    curl_easy_getinfo(hnd, CURLINFO_RESPONSE_CODE, &http_code);
    curl_easy_cleanup(hnd);

    SEXP ans = PROTECT(allocVector(STRSXP, used));
    for (int i = 0; i < used; i++)
        SET_STRING_ELT(ans, i, mkChar(headers[i]));
    setAttrib(ans, install("status"), ScalarInteger((int)http_code));
    UNPROTECT(1);
    return ans;
#endif
}
