pkglib_LTLIBRARIES       = libMRST.la

libMRST_la_SOURCES   =	\
			mrst99.cc \
			mrst2001lo.f \
			PDF_MRST99.C \
			PDF_MRST01LO.C

pkginclude_HEADERS        = \
			mrst99.h \
			PDF_MRST99.H \
			PDF_MRST01LO.H

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@PDFINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
