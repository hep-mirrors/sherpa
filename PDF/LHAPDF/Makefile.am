if LHAPDF_SUPPORT
pkglib_LTLIBRARIES       = libLHAPDFSherpa.la 

if LHAPDF_NATIVE_WRAPPER
libLHAPDFSherpa_la_SOURCES  =	\
			LHAPDF_Fortran_Interface.C
pkginclude_HEADERS        =	LHAPDF_Fortran_Interface.H

else
libLHAPDFSherpa_la_SOURCES  =	\
			LHAPDF_Fortran_Interface.C \
			lhapdf_inter.f 
pkginclude_HEADERS        =	LHAPDF_Fortran_Interface.H
endif

else
LHAPDFEXTRADIST = LHAPDF_Fortran_Interface.H LHAPDF_Fortran_Interface.C \
			lhapdf_inter.f
endif

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@PDFINCS@ @ATOOLSINCS@ @CONDITIONAL_LHAPDFINCS@

EXTRA_DIST    = @SVNINCLUDE@ $(LHAPDFEXTRADIST)
