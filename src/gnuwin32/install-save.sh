# needs save lib pkg R_HOME
lib=$2
pkg=$3
R_HOME=$4

case $1 in
    CHECK) if test -r install.R; then R_SAVE_IMAGE=true; else R_SAVE_IMAGE=false; fi;;
    *) R_SAVE_IMAGE=$1;;
esac
export R_SAVE_IMAGE

if ${R_SAVE_IMAGE}; then
    echo "  save image"
    if test -s "R_PROFILE.R"; then true
    else
	echo "options(echo=FALSE)" > R_PROFILE.R
    fi
    R_PROFILE=./R_PROFILE.R
    export R_PROFILE
    (echo " .lib.loc <- c(\"${lib}\", .lib.loc)"; 
      cat ${lib}/${pkg}/R/${pkg};
      echo "rm(.lib.loc)") | ${R_HOME}/bin/Rterm --save --silent \
        || (echo "Execution of package source for ${pkg} failed"; exit 1)
    mv .RData ${lib}/${pkg}/R/all.rda
    mv ${lib}/${pkg}/R/${pkg} ${lib}/${pkg}/R/${pkg}.R
    cat ${R_HOME}/share/R/firstlib.R > ${lib}/${pkg}/R/${pkg}
    ## if install.R is non-empty, arrange to evaluate the R code it
    ## contains after the package source (maybe for some kind of
    ## cleanup).
    if test -s install.R; then
      cat install.R >> ${lib}/${pkg}/R/${pkg}
    fi
fi
