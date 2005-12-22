SHERPASYS     =	../..

PDFFLAGS      =	# -pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFINCS       =	-I../Main
APACICINCS    =	-I$(SHERPASYS)/APACIC++-2.0/Showers


if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if LHAPDF_SUPPORT
else
lib_LTLIBRARIES       =	libCTEQ.la
libCTEQ_la_SOURCES     =	\
			Cteq61Pdf.f \
			cteq6_inter.f \
		    	CTEQ6_Fortran_Interface.C 
include_HEADERS        =	\
			CTEQ6_Fortran_Interface.H 
endif
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(PDFFLAGS)
INCLUDES      =	$(APACICINCS) $(PDFINCS) $(MODELINCS) $(ATOOLSINCS)

