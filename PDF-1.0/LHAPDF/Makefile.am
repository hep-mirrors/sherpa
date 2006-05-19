if LHAPDF_SUPPORT
pkglib_LTLIBRARIES       = libLHAPDF.la 
libLHAPDF_la_SOURCES  =	\
			LHAPDF_Fortran_Interface.C \
			lhapdf_inter.f 
pkginclude_HEADERS        =	LHAPDF_Fortran_Interface.H
endif

EXTRA_DIST    = @SVNINCLUDE@ LHAPDF_Fortran_Interface.C \
		LHAPDF_Fortran_Interface.H lhapdf_inter.f

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@PDFINCS@ @MODELINCS@ @ATOOLSINCS@
