#! /bin/sh

USER_R_HOME=$1/.R; shift
SEARCHINDEX="${USER_R_HOME}/doc/html/search/index.txt"
rm -rf ${USER_R_HOME}

dirs="${USER_R_HOME} ${USER_R_HOME}/doc ${USER_R_HOME}/doc/html ${USER_R_HOME}/doc/html/search ${USER_R_HOME}/library"
for d in ${dirs}; do
    mkdir -p ${d}
done

for f in AUTHORS COPYING THANKS; do
  if test -f ${R_HOME}/${f}; then
    ln -s ${R_HOME}/${f} ${USER_R_HOME}/${f}
  fi
done

if test -d ${R_HOME}/doc/manual; then
  ln -s ${R_HOME}/doc/manual ${USER_R_HOME}/doc/manual
fi

for f in ${R_HOME}/doc/html/*; do
    if test -f $f; then
	ln -s ${f} ${USER_R_HOME}/doc/html
    fi
done
# we are going to recreate this in R code
rm -f ${USER_R_HOME}/doc/html/packages.html

# class files must be copied for mozilla to work
for f in ${R_HOME}/doc/html/search/*.class; do
    if test -f $f; then
	cp ${f} ${USER_R_HOME}/doc/html/search
    fi
done
for f in ${R_HOME}/doc/html/search/*.html; do
    if test -f $f; then
	ln -s ${f} ${USER_R_HOME}/doc/html/search
    fi
done
ln -s ${R_HOME}/doc/html/search/index.txt ${USER_R_HOME}/doc/html/search

rm -f ${SEARCHINDEX}

get_unique () {
  if test -r ${1}; then
    x="1"
    while test -r ${1}.${x}; do
      x=`echo "$x+1" | bc`
    done
    echo ${1}.${x}           
  else
    echo $1
  fi
}
     

for lib in $*; do
    if test -d ${lib}; then
      for pkg in `ls -d ${lib}/* | sed '/CVS$/d; /profile$/d'`; do
	if test -d ${pkg}; then
	    pkgname=`basename ${pkg}`
	    target=`get_unique ${USER_R_HOME}/library/${pkgname}`
	    targetname=`basename ${target}`
	    ln -s ${pkg} ${target}	    

	    if test -r ${pkg}/CONTENTS; then
		cat ${pkg}/CONTENTS | \
		    sed "s/\/library\/${pkgname}\//\/library\/${targetname}\//;"  >> ${SEARCHINDEX}
	    fi
	fi
      done
    fi
done

ln -s ${R_HOME}/doc/html/R.css ${USER_R_HOME}/library

### Local Variables: ***
### mode: sh ***
### sh-indentation: 2 ***
### End: ***
