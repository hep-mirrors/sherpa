if LHAPDF_SUPPORT
CTEQEXTRADIST = CTEQ6_Fortran_Interface.H CTEQ6_Fortran_Interface.C \
		cteq6_inter.f Cteq61Pdf.f
else
pkglib_LTLIBRARIES       =	libCTEQ.la
libCTEQ_la_SOURCES     =	\
			Cteq61Pdf.f \
			cteq6_inter.f \
		    	CTEQ6_Fortran_Interface.C 
pkginclude_HEADERS        =	\
			CTEQ6_Fortran_Interface.H 
endif

INCLUDES      =	@PDFINCS@ @ATOOLSINCS@

AM_CXXFLAGS   =	@SHERPAFLAGS@

EXTRA_DIST    = @SVNINCLUDE@ $(CTEQEXTRADIST)
