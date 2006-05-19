if LHAPDF_SUPPORT
else
pkglib_LTLIBRARIES       =	libCTEQ.la
libCTEQ_la_SOURCES     =	\
			Cteq61Pdf.f \
			cteq6_inter.f \
		    	CTEQ6_Fortran_Interface.C 
pkginclude_HEADERS        =	\
			CTEQ6_Fortran_Interface.H 
endif

INCLUDES      =	@APACICINCS@ @PDFINCS@ @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@

AM_CXXFLAGS   =	@SHERPAFLAGS@
