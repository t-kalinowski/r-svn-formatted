# -*- Makefile -*-
#
# ${R_HOME}/share/make/tests.mk

ECHO_C =
ECHO_N = -n
ECHO_T =

makevars =
srcdir = .
VPATH = .

test-src = $(test-src-1) $(test-src-auto)
test-out = $(test-src:.R=.Rout)

R = srcdir=$(srcdir) $(R_HOME)/bin/Rterm.exe --vanilla
RDIFF = $(R_HOME)/bin/Rcmd Rdiff.sh

.SUFFIXES:
.SUFFIXES: .R .Rin .Rout

.Rin.R:
	@echo "Creating \`$@'"
	@$(R) < $< > /dev/null

.R.Rout:
	@rm -f $@ $@.fail
	@echo "  Running \`$<'"
	@$(R) R_LIBS="$(R_LIBS)" < $< > $@
	@if test -f $(srcdir)/$@.save; then \
	  mv $@ $@.fail; \
	  echo -n "  Comparing \`$@' to \`$@.save' ..."; \
	  $(RDIFF) $@.fail $(srcdir)/$@.save 0 || exit 1; \
	  mv $@.fail $@; \
	  echo "OK"; \
	fi

all:
	@(out=`echo "$(test-out)" | sed 's/ $$//g'`; \
	  if test -n "$${out}"; then \
	    $(MAKE) -f $(R_HOME)/share/make/tests.mk $(makevars) $${out}; \
	  fi)

clean:
	@rm -f $(test-out) $(test-src-auto) *.fail

