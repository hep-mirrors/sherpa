SHERPASYS     =	../..

PDFFLAGS      =	-pedantic -Wall -Winline -Woverloaded-virtual

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFINCS       =	-I../LHAPDF -I../KMR

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       =	libPDF.la
libPDF_la_SOURCES     =	\
			ISR_Handler.C \
		    	Intact.C \
			Structure_Function.C \
		    	PDF_Base.C \
			mrst99.cc \
			PDF_MRST99.C \
		    	PDF_Electron.C \
			Cteq61Pdf.f \
			cteq6_inter.f \
		    	CTEQ6_Fortran_Interface.C \
			GRVph_Fortran_Interface.C \
			grvphoton.f \
		    	PDF_Handler.C 
noinst_HEADERS        =	\
			ISR_Base.H \
		 	ISR_Handler.H \
			Intact.H \
			Structure_Function.H\
	         	PDF_Base.H \
			mrst99.h \
			PDF_MRST99.H \
		 	PDF_Electron.H \
			CTEQ6_Fortran_Interface.H \
			GRVph_Fortran_Interface.H \
			PDF_Handler.H 
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(PDFFLAGS)
INCLUDES      =	$(SHERPAINCS) $(APACICINCS) $(PDFINCS) $(MODELINCS) $(ATOOLSINCS)

