pkglib_LTLIBRARIES    =	libExtraXSCOBG.la
libExtraXSCOBG_la_SOURCES     = \
				BG_Amplitude.C \
				N_Gluon_BG.C \
				XS_QCD_COBG.C
pkginclude_HEADERS    =	\
			BG_Amplitude.H \
			N_Gluon_BG.H \
			XS_QCD_COBG.H

EXTRA_DIST    = @SVNINCLUDE@
AM_CXXFLAGS   =	@SHERPAFLAGS@
INCLUDES      = @EXTRAXSINCS@ @PHASICINCS@ @MODELINCS@ \
		@PDFINCS@ @BEAMINCS@ @ATOOLSINCS@
