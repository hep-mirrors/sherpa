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
if PGS_SUPPORT
PGSLIBS       =	-L../../AddOns/PGS -lPGSSherpa \
		@CONDITIONAL_PGSLIBS@
endif

Sherpa@VERSIONING@_LDADD  =	-L@SHERPABUILDDIR@/Main -lSherpaMain \
		$(HZTOOLLIBS) $(PGSLIBS) @ATOOLSLIBS@

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
