if LHAPDF_SUPPORT
pkglib_LTLIBRARIES       = libLHAPDFSherpa.la 

if LHAPDF_NATIVE_WRAPPER
libLHAPDFSherpa_la_SOURCES  =	\
			LHAPDF_Fortran_Interface.C
LHAPDFEXTRADIST = lhapdf_inter.f
else
libLHAPDFSherpa_la_SOURCES  =	\
			LHAPDF_Fortran_Interface.C \
			lhapdf_inter.f 
endif

localinc_HEADERS      =	LHAPDF_Fortran_Interface.H

else
LHAPDFEXTRADIST = LHAPDF_Fortran_Interface.H LHAPDF_Fortran_Interface.C \
			lhapdf_inter.f
endif

libLHAPDFSherpa_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_LHAPDFINCS@

EXTRA_DIST    = @SVNINCLUDE@ $(LHAPDFEXTRADIST)
