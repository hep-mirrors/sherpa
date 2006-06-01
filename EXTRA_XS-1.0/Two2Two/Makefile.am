pkglib_LTLIBRARIES       =	libExtraXS2_2.la

libExtraXS2_2_la_SOURCES      =	\
				XS_QCD.C \
				XS_Drell_Yan.C \
				XS_4F.C \
				Off_Shell_EW.C \
				Off_Shell_QCD.C \
				QED_Processes.C \
				QCD_Processes.C

pkginclude_HEADERS        =	\
			XS_QCD.H \
			XS_Drell_Yan.H \
			XS_4F.H \
			Off_Shell_EW.H \
			Off_Shell_QCD.H \
			QED_Processes.H \
			QCD_Processes.H

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@EXTRAXSINCS@ @PHASICINCS@ @PDFINCS@ \
		@BEAMINCS@ @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
