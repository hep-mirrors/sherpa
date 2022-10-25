GITTAG = PDF/LHAPDF
include ../../git.make

pkglib_LTLIBRARIES       = libLHAPDFSherpa.la 

libLHAPDFSherpa_la_SOURCES    = \
				Git_Info.C
libLHAPDFSherpa_la_SOURCES += \
				LHAPDF_CPP_Interface.C

libLHAPDFSherpa_la_LIBADD = @CONDITIONAL_LHAPDFLIBS@ @PDFLIBS@ @ATOOLSLIBS@

libLHAPDFSherpa_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_LHAPDFINCS@
