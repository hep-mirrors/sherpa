SHERPASYS     =	../..

PDFFLAGS      =	-pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFINCS       =	-I../LHAPDF -I../Sudakov -I../KMR
APACICINCS    =	-I$(SHERPASYS)/APACIC++-2.0/Showers


if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       =	libPDF.la
libPDF_la_SOURCES     =	\
			Cteq61Pdf.f \
			cteq6_inter.f \
			grvphoton.f \
			mrst99.cc \
		    	PDF_Base.C \
		    	PDF_Electron.C \
		    	CTEQ6_Fortran_Interface.C \
			GRVph_Fortran_Interface.C \
			PDF_MRST99.C \
		    	PDF_Handler.C \
			ISR_Base.C \
		    	Intact.C \
			Structure_Function.C \
			ISR_Handler.C 
noinst_HEADERS        =	\
			mrst99.h \
	         	PDF_Base.H \
		 	PDF_Electron.H \
			CTEQ6_Fortran_Interface.H \
			GRVph_Fortran_Interface.H \
			PDF_MRST99.H \
			PDF_Handler.H \
			ISR_Base.H \
			Intact.H \
			Structure_Function.H\
		 	ISR_Handler.H 
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(PDFFLAGS)
INCLUDES      =	$(APACICINCS) $(PDFINCS) $(MODELINCS) $(ATOOLSINCS)

