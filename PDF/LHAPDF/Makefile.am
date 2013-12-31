SVNTAG = PDF/LHAPDF
include ../../svn.make

pkglib_LTLIBRARIES       = libLHAPDFSherpa.la 

if LHAPDF_SUPPORT
libLHAPDFSherpa_la_SOURCES    = \
				SVN_Info.C
if LHAPDF_NATIVE_CPP
libLHAPDFSherpa_la_SOURCES += \
				LHAPDF_CPP_Interface.C
else
libLHAPDFSherpa_la_SOURCES += \
				LHAPDF_Fortran_Interface.C
localinc_HEADERS            = \
				LHAPDF_Fortran_Interface.H
endif
else
LHAPDFEXTRADIST = \
		SVN_Info.C \
		LHAPDF_Fortran_Interface.H \
		LHAPDF_Fortran_Interface.C \
		LHAPDF_CPP_Interface.C
NOMD5SUM      =	1
endif

libLHAPDFSherpa_la_LIBADD = @CONDITIONAL_LHAPDFLIBS@

libLHAPDFSherpa_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_LHAPDFINCS@

EXTRA_DIST    = $(LHAPDFEXTRADIST)
