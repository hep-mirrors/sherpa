pkglib_LTLIBRARIES       = libAmisic.la

libAmisic_la_SOURCES  = \
			MI_Base.C \
			Amisic.C 
			
pkginclude_HEADERS        = \
			MI_Base.H \
			Amisic.H 

INCLUDES      =	@CONDITIONAL_SHERPAINCS@ @AMISICINCS@ @EXTRAXSINCS@ @PHASICINCS@ \
		@PDFINCS@ @BEAMINCS@ @MODELINCS@ @ATOOLSINCS@

AM_CXXFLAGS   =	@SHERPAFLAGS@

EXTRA_DIST    = @SVNINCLUDE@
