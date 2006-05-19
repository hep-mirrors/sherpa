pkglib_LTLIBRARIES       = libAmisicTools.la

libAmisicTools_la_SOURCES     = \
				Amisic_Histogram.C \
				Fixed_Variable_Channel.C \
				Environment.C \
				Semihard_QCD.C 

pkginclude_HEADERS        = \
			Amisic_Histogram.H \
			Fixed_Variable_Channel.H \
			Environment.H \
			Semihard_QCD.H 

INCLUDES      =	@CONDITIONAL_SHERPAINCS@ @EXTRAXSINCS@ @PHASICINCS@ @PDFINCS@ \
		@BEAMINCS@ @MODELINCS@ @ATOOLSINCS@

AM_CXXFLAGS   =	@SHERPAFLAGS@

EXTRA_DIST    = @SVNINCLUDE@
