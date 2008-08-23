bin_PROGRAMS  =	Sherpa

Sherpa_SOURCES        =	Main.C 

Sherpa_LDADD  =	-L@SHERPABUILDDIR@/Main -lSherpaMain \
		@CONDITIONAL_ROOTLIBS@ @CONDITIONAL_LHAPDFLIBS@ \
		@CONDITIONAL_HEPMC2LIBS@

Sherpa_LDFLAGS        =	@ldflags@

AM_CXXFLAGS   =	@SHERPAFLAGS@ @CONDITIONAL_ROOTFLAGS@ ## $(MCATNLOFLAGS)

INCLUDES      = @SHERPAINCS@ @ANALYSISINCS@ @CONDITIONAL_HADRONSINCS@ \
		@APACICINCS@ @CONDITIONAL_AMISICINCS@ @AMEGICINCS@ \
		@AHADICINCS@ @EXTRAXSINCS@ @PHASICINCS@ @MODELINCS@ \
		@CONDITIONAL_PHOTONSINCS@ @PDFINCS@ @BEAMINCS@ @ATOOLSINCS@ \
		@CONDITIONAL_ROOTINCS@

dist_pkgdata_SCRIPTS  =	makelibs grid_makelibs

SUBDIRS       =	CTEQ6Grid MRST99Grid

dist_pkgdata_DATA         = \
			Analysis.dat \
			Fragmentation.dat \
			Integration.dat \
			ISR.dat \
			ME.dat \
			MI.dat \
			Run.dat \
			Shower.dat

EXTRA_DIST    =	@SVNINCLUDE@

install-data-local:
	@echo "Installing Decaydata."; \
	cd $(srcdir); \
	for i in Decaydata/*; do \
	  if test -d $$i; then \
	    mkdir -p $(pkgdatadir)/$$i; \
	    for j in $$i/*.dat; do \
	      cp $$j $(pkgdatadir)/$$i/; \
	    done; \
	    if test -d $$i/PhaseSpaceFunctions; then \
	      cp -rf $$i/PhaseSpaceFunctions $(pkgdatadir)/$$i/; \
	    fi; \
	  fi; \
	done; \
	cp Decaydata/*.dat $(pkgdatadir)/Decaydata/; \
	cd -;

uninstall-local:
	@echo "Removing Decaydata."; \
	rm -rf "$(pkgdatadir)/Decaydata"; \
	if test "$$(ls -A $(pkgdatadir))" = ""; then \
	rm -rf "$(pkgdatadir)"; \
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
	    fi; \
	  fi; \
	done; \
	cp Decaydata/*.dat $(distdir)/Decaydata/; \
	SAMPLERUNS="Tevatron1800 LHC EGamma LEP91"; \
	for i in $$SAMPLERUNS ; do \
	  mkdir -p $(distdir)/$$i; \
	  if test -x $$(which svn) -a -d $$i/.svn; then \
	    for j in $$(svn ls $$i); do \
	      cp $$i/$$j $(distdir)/$$i/; \
	    done; \
          else \
	    for j in $$i/*.dat; do \
	      cp $$j $(distdir)/$$i/; \
	    done; \
	  fi; \
	done; \
	cd -;
