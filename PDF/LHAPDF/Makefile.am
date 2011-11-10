SVNTAG = PDF/LHAPDF
include ../../svn.make

pkglib_LTLIBRARIES       = libLHAPDFSherpa.la 

if LHAPDF_SUPPORT
libLHAPDFSherpa_la_SOURCES  =	\
			SVN_Info.C \
			LHAPDF_Fortran_Interface.C

localinc_HEADERS      =	LHAPDF_Fortran_Interface.H

else
LHAPDFEXTRADIST = SVN_Info.C LHAPDF_Fortran_Interface.H LHAPDF_Fortran_Interface.C
endif

libLHAPDFSherpa_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_LHAPDFINCS@

EXTRA_DIST    = @SVNINCLUDE@ $(LHAPDFEXTRADIST)
