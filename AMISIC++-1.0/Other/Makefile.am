SHERPASYS  = ../..

ATOOLSINC  = $(SHERPASYS)/ATOOLS-2.0
PDFINC     = $(SHERPASYS)/PDF-1.0
BEAMINC    = $(SHERPASYS)/BEAM-1.0
PHASICINC  = $(SHERPASYS)/PHASIC++-1.0
MODELINC   = $(SHERPASYS)/MODEL-1.0
EXTRAXSINC = $(SHERPASYS)/EXTRA_XS-1.0

lib_LTLIBRARIES       = libAmisicHelpers.la
libAmisicHelpers_la_SOURCES   =	\
				QCD_Processes_C.C \
				Channel_Interface.C 
CXXFLAGS      =	-pedantic -Wall -Winline
INCLUDES      =	\
		-I../Tools \
		-I$(EXTRAXSINC)/Two2Two \
		-I$(EXTRAXSINC)/Main \
		-I$(PHASICINC)/Main \
		-I$(MODELINC)/Main \
		-I$(PDFINC)/Main \
		-I$(BEAMINC)/Main \
		-I$(ATOOLSINC)/Phys \
		-I$(ATOOLSINC)/Math \
		-I$(ATOOLSINC)/Org 
noinst_HEADERS        =	\
			QCD_Processes_C.H \
			Channel_Interface.H 
