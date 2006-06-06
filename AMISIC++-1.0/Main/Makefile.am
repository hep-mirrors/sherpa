pkglib_LTLIBRARIES       = libAmisic.la

libAmisic_la_SOURCES  = \
			MI_Base.C \
			Amisic.C 
			
pkginclude_HEADERS        = \
			MI_Base.H \
			Amisic.H 

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@CONDITIONAL_SHERPAINCS@ @AMISICINCS@ @EXTRAXSINCS@ @PHASICINCS@ \
		@PDFINCS@ @BEAMINCS@ @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
