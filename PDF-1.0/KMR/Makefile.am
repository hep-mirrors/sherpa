pkglib_LTLIBRARIES       = libDUPDF.la

libDUPDF_la_SOURCES   =	\
			LL_Branching.C \
			LL_Single_Sudakov.C \
			LL_Sudakov.C \
			Doubly_Unintegrated_PDF.C

pkginclude_HEADERS        = \
			LL_Branching.H \
			LL_Single_Sudakov.H \
			LL_Sudakov.H \
			Doubly_Unintegrated_PDF.H

EXTRA_DIST    = @SVNINCLUDE@

AM_CXXFLAGS   =	@SHERPAFLAGS@ @CONDITIONAL_ROOTFLAGS@

INCLUDES      =	@APACICINCS@ @PDFINCS@ @MODELINCS@ @ATOOLSINCS@ @CONDITIONAL_ROOTINCS@
