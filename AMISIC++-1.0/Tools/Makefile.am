pkglib_LTLIBRARIES       = libAmisicTools.la

libAmisicTools_la_SOURCES     = \
				Amisic_Histogram.C \
				Fixed_Variable_Channel.C \
				Semihard_QCD.C 

pkginclude_HEADERS        = \
			Amisic_Histogram.H \
			Fixed_Variable_Channel.H \
			Semihard_QCD.H 

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@EXTRAXSINCS@ @PHASICINCS@ @PDFINCS@ \
		@BEAMINCS@ @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
