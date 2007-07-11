#! /bin/sh

revision='$Rev$'
version=`set - ${revision}; echo ${2}`
version="Stangle front-end r${version}

Copyright (C) 2006 The R Core Development Team.
This is free software; see the GNU General Public License version 2
or later for copying conditions.  There is NO warranty."

usage="Usage: R CMD Stangle file

A simple front-end for Stangle()

Options:
  -h, --help		print short help message and exit
  -v, --version		print Sweave version info and exit

Report bugs to <r-bugs@r-project.org>."

case ${1} in
  -h|--help)
     echo "${usage}"; exit 0 ;;
  -v|--version)
     echo "${version}"; exit 0 ;;
esac

R_EXE="${R_HOME}/bin/rterm.exe"
echo "library(\"utils\"); Stangle(\"$1\")" | \
  "${R_EXE}" --no-restore --slave
