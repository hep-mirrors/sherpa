pkglib_LTLIBRARIES       =	libExtraXS2_2.la

libExtraXS2_2_la_SOURCES      =	\
				XS_QCD.C \
				XS_Drell_Yan.C \
				XS_4F.C \
				Off_Shell_QCD.C

pkginclude_HEADERS        =	\
			XS_QCD.H \
			XS_Drell_Yan.H \
			XS_4F.H \
			Off_Shell_QCD.H

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@EXTRAXSINCS@ @PHASICINCS@ @PDFINCS@ \
		@BEAMINCS@ @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
