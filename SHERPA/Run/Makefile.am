bin_PROGRAMS  =	Sherpa@VERSIONING@

Sherpa@VERSIONING@_SOURCES        =	Main.C 

if HZTOOL_SUPPORT
if CERNLIB_SUPPORT
HZTOOLLIBS    = -L../../AddOns/HZTool -lHZToolSherpa \
		@CONDITIONAL_HZTOOLLIBS@ @CONDITIONAL_CERNLIBLIBS@
else
HZTOOLLIBS    =
endif
else
HZTOOLLIBS    =
endif

Sherpa@VERSIONING@_LDADD  =	-L@SHERPABUILDDIR@/Main -lSherpaMain \
		@CONDITIONAL_ROOTLIBS@ $(HZTOOLLIBS) \
		@CONDITIONAL_HEPMC2LIBS@ @CONDITIONAL_GOLEM95LIBS@

Sherpa@VERSIONING@_LDFLAGS        =	@ldflags@

Sherpa@VERSIONING@_CXXFLAGS   = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@

Sherpa@VERSIONING@_CPPFLAGS   = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@

dist_pkgdata_DATA         = \
			Analysis.dat \
			Fragmentation.dat \
			Integration.dat \
			ISR.dat \
			ME.dat \
			MI.dat \
			Run.dat \
			Shower.dat

bin_SCRIPTS = Sherpa@VERSIONING@-config get_bibtex

install-data-local:
	@echo "Installing Decaydata."; \
	cd "$(srcdir)"; \
	for i in Decaydata/*; do \
	  if test -d "$$i"; then \
	    mkdir -p "$(DESTDIR)$(pkgdatadir)/$$i"; \
	    for j in $$i/*.dat; do \
	      cp "$$j" "$(DESTDIR)$(pkgdatadir)/$$i/"; \
	    done; \
	    if test -d "$$i/PhaseSpaceFunctions" ; then \
	      cp -rf "$$i/PhaseSpaceFunctions" "$(DESTDIR)$(pkgdatadir)/$$i/"; \
	      chmod -R +w "$(DESTDIR)$(pkgdatadir)/$$i/"; \
	      if test "x@SVNINCLUDE@" == "x"; then \
		rm -rf "$(DESTDIR)$(pkgdatadir)/$$i/PhaseSpaceFunctions/.svn"; \
	      fi; \
	    fi; \
	  fi; \
	done; \
	cp Decaydata/*.dat "$(DESTDIR)$(pkgdatadir)/Decaydata/"; \
	cd -;

uninstall-local:
	@echo "Removing Decaydata."; \
	rm -rf "$(DESTDIR)$(pkgdatadir)/Decaydata"; \
	if test "$$(ls -A $(DESTDIR)$(pkgdatadir))" = ""; then \
	rm -rf "$(DESTDIR)$(pkgdatadir)"; \
	fi;

dist-hook:
	cd $(srcdir); \
	for i in Decaydata/*; do \
	  if test -d $$i; then \
	    mkdir -p $(distdir)/$$i; \
	    for j in $$i/*.dat; do \
	      cp $$j $(distdir)/$$i/; \
	    done; \
	    if test -d $$i/PhaseSpaceFunctions; then \
	      cp -rf $$i/PhaseSpaceFunctions $(distdir)/$$i/; \
	      if test "x@SVNINCLUDE@" == "x"; then \
		rm -rf "$(distdir)/$$i/PhaseSpaceFunctions/.svn"; \
	      fi; \
	    fi; \
	  fi; \
	done; \
	cp Decaydata/*.dat $(distdir)/Decaydata/; \
	cd -;
