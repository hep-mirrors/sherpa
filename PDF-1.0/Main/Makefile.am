SHERPASYS     =	../..

PDFFLAGS      =	-pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
BEAMINCS      = -I$(SHERPASYS)/BEAM-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFINCS       =	-I../LHAPDF -I../MRST -I../CTEQ -I../GRV -I../Sudakov \
		-I../KMR -I../Remnant
APACICINCS    =	-I$(SHERPASYS)/APACIC++-2.0/Showers

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       =	libPDF.la
libPDF_la_SOURCES     =	\
		    	PDF_Base.C \
		    	PDF_Electron.C \
			Continued_PDF.C \
		    	PDF_Handler.C \
			ISR_Base.C \
		    	Intact.C \
			Structure_Function.C \
			ISR_Handler.C 
include_HEADERS       =	\
	         	PDF_Base.H \
		 	PDF_Electron.H \
			Continued_PDF.H \
			PDF_Handler.H \
			ISR_Base.H \
			Intact.H \
			Structure_Function.H\
			ISR_Info.H \
		 	ISR_Handler.H 
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(PDFFLAGS)
INCLUDES      =	$(SHERPAINCS) $(APACICINCS) $(PDFINCS) \
		$(MODELINCS) $(BEAMINCS) $(ATOOLSINCS)

