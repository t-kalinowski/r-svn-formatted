link.html.help <- function(verbose=FALSE, lib.loc=.libPaths())
{
    if(!file.exists(file.path(R.home(), "doc", "html", "search")))
       return(invisible(NULL))
    if(verbose) {
        cat("updating HTML package descriptions\n")
        flush.console()
    }
    make.packages.html(lib.loc)
    make.search.html(lib.loc)
    fixup.libraries.URLs(lib.loc)
}

make.packages.html <- function(lib.loc=.libPaths())
{
    f.tg <- file.path(R.home(), "doc/html/packages.html")
    f.hd <- file.path(R.home(), "doc/html/packages-head.html")
    if(!file.create(f.tg)) {
        warning("cannot update HTML package index")
        return(FALSE)
    }
    file.append(f.tg, f.hd)
    out <- file(f.tg, open="a")
    cat('<table width=\"100%\" summary="R Package list">\n',
        file=out)
    rh <- gsub("\\\\", "/", R.home())
    drive <- substring(rh, 1, 2)
    for (lib in lib.loc) {
        pg <- sort(.packages(all.available = TRUE, lib.loc = lib))
        ## use relative indexing for .Library
        if(is.na(pmatch(rh, lib))) {
            libname <- gsub("/", "\\\\", lib)
            lib0 <- if(substring(lib, 2, 2) != ":")
                paste(drive, lib, sep="") else lib
            lib0 <- paste("file:///", lib0, sep="")
        } else {
            lib0 <- "../../library"
            libname <- "the standard library"
        }
        if(length(lib.loc) > 1)
            cat("<p><h3>Packages in ", libname, "</h3>\n",
                sep = "", file = out)
        if(libname != "the standard library")
            cat("<p>Cross-links from this library to other libraries may not work.\n\n",file = out)
        cat("<p>\n<table width=\"100%\">\n", file = out)
        for (i in  pg) {
            title <- package.description(i, fields="Title", lib.loc = lib)[1]
            if (is.na(title)) title <- "-- Title is missing --"
            cat('<tr align="left" valign="top">\n',
                '<td width="25%"><a href="', lib0, '/', i,
                '/html/00Index.html">', i, "</a></td><td>", title,
                "</td></tr>\n", file=out, sep="")
        }
        cat("</table>\n\n", file=out)
    }
    cat("</body></html>\n", file=out)
    close(out)
    invisible(TRUE)
}

make.search.html <- function(lib.loc=.libPaths())
{
    f.tg <- file.path(R.home(), "doc/html/search/index.txt")
    out <- file(f.tg, open = "w")
    if(class(out) == "try-error") {
        warning("cannot update HTML search index")
        return()
    }
    for (lib in lib.loc) {
        rh <- gsub("\\\\", "/", R.home())
        drive <- substring(rh, 1, 2)
        pg <- sort(.packages(all.available = TRUE, lib.loc = lib))
        ## use relative indexing for .Library
        if(is.na(pmatch(rh, lib))) {
            lib0 <- if(substring(lib, 2, 2) != ":") paste(drive, lib, sep="")
            else lib
            lib0 <- paste("URL: file:///", lib0, sep="")
            sed.it <- TRUE
        } else {
            sed.it <- FALSE
        }
        for (i in pg) {
            cfile <- file.path(lib, i, "CONTENTS")
            if(!file.exists(cfile)) next
            contents <- if(sed.it)
                gsub("^URL: ../../../library", lib0, readLines(cfile))
            else readLines(cfile)
            writeLines(contents, out)
        }
    }
    close(out)
}

fixup.package.URLs <- function(pkg, force = FALSE)
{
    top <- paste("file:///", gsub("\\\\", "/", R.home()), sep="")
    fixedfile <- file.path(pkg, "fixedHTMLlinks")
    if(file.exists(fixedfile)) {
        oldtop <- readLines(fixedfile)
        if(!force && (length(oldtop)==1) && top == oldtop) return(TRUE)
        olddoc <- paste(oldtop, "/doc", sep="")
        oldbase <- paste(oldtop, "/library/base", sep="")
    } else {
        olddoc <- "../../../doc"
        oldbase <- "../../base"
    }
    if(!file.create(fixedfile)) return(FALSE)
    cat(top, "\n", sep="", file=fixedfile)
    htmldir <- file.path(pkg, "html")
    if(!file.exists(htmldir)) return(FALSE)
    files <- list.files(htmldir, pattern = "\.html$", full.names = TRUE)
    doc <- paste(top, "/doc", sep="")
    base <- paste(top, "/library/base", sep="")
    for(f in files) {
        page <- readLines(f)
        try(out <- file(f, open = "w"))
        if(class(out) == "try-error") {
            warning("cannot update", f)
            next
        }
        page <- gsub(olddoc, doc, page)
        page <- gsub(oldbase, base, page)
        writeLines(page, out)
        close(out)
    }
    return(TRUE)
}

fixup.libraries.URLs <- function(lib.loc = .libPaths())
{
    for (lib in lib.loc) {
        if(lib == .Library) next
        pg <- sort(.packages(all.available = TRUE, lib.loc = lib))
        for(pkg in pg) fixup.package.URLs(file.path(lib,pkg))
    }
}
