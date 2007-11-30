if LHAPDF_SUPPORT
pkglib_LTLIBRARIES       = libLHAPDFSherpa.la 
libLHAPDFSherpa_la_SOURCES  =	\
			LHAPDF_Fortran_Interface.C \
			lhapdf_inter.f 
pkginclude_HEADERS        =	LHAPDF_Fortran_Interface.H
else
LHAPDFEXTRADIST = LHAPDF_Fortran_Interface.H LHAPDF_Fortran_Interface.C \
			lhapdf_inter.f
endif

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@PDFINCS@ @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@ $(LHAPDFEXTRADIST)
