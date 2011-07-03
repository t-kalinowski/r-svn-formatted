## ${R_HOME}/share/make/lazycomp.mk
$(top_builddir)/library/$(pkg)/R/$(pkg).rdb: all.R
	@$(INSTALL_DATA) all.R $(top_builddir)/library/$(pkg)/R/$(pkg)
	@$(ECHO) "byte-compiling package '$(pkg)'"
	@$(ECHO) "tools:::makeLazyLoading(\"$(pkg)\")" | \
	  R_COMPILE_PKGS=1 R_COMPILER_SUPPRESS_ALL=1 \
	  R_DEFAULT_PACKAGES=$(DEFPKGS) LC_ALL=C $(R_EXE) > /dev/null
