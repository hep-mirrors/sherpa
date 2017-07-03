SVNTAG = PDF/LHAPDF
include ../../svn.make

pkglib_LTLIBRARIES       = libLHAPDFSherpa.la 

libLHAPDFSherpa_la_SOURCES    = \
				SVN_Info.C
if LHAPDF_SUPPORT
if LHAPDF6_SUPPORT
libLHAPDFSherpa_la_SOURCES += \
				LHAPDF_CPP_Interface.C
LHAPDFEXTRADIST = \
		LHAPDF_Fortran_Interface.H \
		LHAPDF_Fortran_Interface.C
else
libLHAPDFSherpa_la_SOURCES += \
				LHAPDF_Fortran_Interface.C
localinc_HEADERS            = \
				LHAPDF_Fortran_Interface.H
LHAPDFEXTRADIST = \
		LHAPDF_CPP_Interface.C
endif
else
libLHAPDFSherpa_la_SOURCES += \
		LHAPDF_Dummy.C

LHAPDFEXTRADIST = \
		LHAPDF_Fortran_Interface.H \
		LHAPDF_Fortran_Interface.C \
		LHAPDF_CPP_Interface.C
endif

libLHAPDFSherpa_la_LIBADD = @CONDITIONAL_LHAPDFLIBS@

libLHAPDFSherpa_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_LHAPDFINCS@

EXTRA_DIST    = $(LHAPDFEXTRADIST)
