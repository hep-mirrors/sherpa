bin_PROGRAMS  =	Sherpa@VERSIONING@

Sherpa@VERSIONING@_SOURCES        =	Main.C 

Sherpa@VERSIONING@_LDADD  =	-L@SHERPABUILDDIR@/Main -lSherpaMain \
		@CONDITIONAL_ROOTLIBS@ \
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

bin_SCRIPTS = Sherpa@VERSIONING@-config
dist_bin_SCRIPTS = get_bibtex
