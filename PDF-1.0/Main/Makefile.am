SHERPASYS  = ../..

ATOOLSINC  = $(SHERPASYS)/ATOOLS-2.0
MODELINC   = $(SHERPASYS)/MODEL-1.0

lib_LTLIBRARIES       =	libPDF.la
libPDF_la_SOURCES     =	ISR_Handler.C \
		    	Intact.C \
			Structure_Function.C \
		    	PDF_Base.C \
			mrst99.cc \
			PDF_MRST99.C \
		    	PDF_Electron.C \
			Cteq61Pdf.f \
			cteq6_inter.f \
		    	CTEQ6_Fortran_Interface.C \
		    	PDF_Handler.C 
CXXFLAGS      =	-pedantic -Wall -Winline
INCLUDES      =	-I../LHAPDF \
	   	-I$(MODELINC)/Main \
	   	-I$(ATOOLSINC)/Phys \
	   	-I$(ATOOLSINC)/Org \
	   	-I$(ATOOLSINC)/Math
noinst_HEADERS        =	ISR_Base.H \
		 	ISR_Handler.H \
			Intact.H \
			Structure_Function.H\
	         	PDF_Base.H \
			mrst99.h \
			PDF_MRST99.H \
		 	PDF_Electron.H \
			CTEQ6_Fortran_Interface.H \
			PDF_Handler.H 

