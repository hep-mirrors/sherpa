SHERPASYS     =	../..

SHERPAFLAGS   =	-pedantic -Wall -Winline 

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFDIR	      =	$(SHERPASYS)/PDF-1.0
PDFINCS	      =	-I$(PDFDIR)/Main -I$(PDFDIR)/LHAPDF -I$(PDFDIR)/KMR
SHERPAINCS    =	-I../Single_Events 	

lib_LTLIBRARIES       =	libHerwigTools.la
libHerwigTools_la_SOURCES       =	\
				herwig6500mod.f \
				fiherwig.f \
				hwudat.f \
				Herwig_Interface.C 
noinst_HEADERS        =	\
			herwig6500.inc \
		 	Herwig_Wrapper.H \
			Herwig_Interface.H 

AM_CXXFLAGS   =	$(SHERPAFLAGS)
INCLUDES      =	$(MODELINCS) $(ATOOLSINCS) $(PDFINCS) $(SHERPAINCS)




