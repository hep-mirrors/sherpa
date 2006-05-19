pkglib_LTLIBRARIES       =	libExtraXS2_2.la

libExtraXS2_2_la_SOURCES      =	\
				XS_QCD.C \
				XS_Drell_Yan.C \
				XS_4F.C \
				QED_Processes.C \
				QCD_Processes.C

pkginclude_HEADERS        =	\
			XS_QCD.H \
			XS_Drell_Yan.H \
			XS_4F.H \
			QED_Processes.H \
			QCD_Processes.H

EXTRA_DIST    = @SVNINCLUDE@

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@EXTRAXSINCS@ @PHASICINCS@ @PDFINCS@ \
		@BEAMINCS@ @MODELINCS@ @ATOOLSINCS@
