SHERPASYS     =	../..

SHERPAFLAGS   =	# -pedantic -Wall -Winline 
DEFS          =

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFDIR	      =	$(SHERPASYS)/PDF-1.0
PDFINCS	      =	-I$(PDFDIR)/Main -I$(PDFDIR)/LHAPDF -I$(PDFDIR)/KMR
SHERPAINCS    =	-I../Single_Events 	

lib_LTLIBRARIES       =	libHerwigTools.la
libHerwigTools_la_SOURCES       =	\
				fiherwig.f \
				fimcatnlo.f \
				mcatnlo_hwlhin.f \
				mcatnlo_str.f \
				herwig6504mod.f \
				MCatNLO_Interface.C \
				Herwig_Interface.C 
noinst_HEADERS        =	\
			MCatNLO_Wrapper.H \
			MCatNLO_Interface.H \
		 	Herwig_Wrapper.H \
			Herwig_Interface.H \
			HERWIG65.INC 

AM_CXXFLAGS   =	$(SHERPAFLAGS)
INCLUDES      =	$(MODELINCS) $(ATOOLSINCS) $(PDFINCS) $(SHERPAINCS)
