SHERPASYS  = ../..

ATOOLSINC  = $(SHERPASYS)/ATOOLS-2.0
MODELINC   = $(SHERPASYS)/MODEL-1.0
AMISICINC  = $(SHERPASYS)/AMISIC++-1.0

lib_LTLIBRARIES       =	libLundTools.la
libLundTools_la_SOURCES       =	\
				pydata.f \
				pyplst.f \
				Lund_Interface.C 
AM_CXXFLAGS   =	-pedantic -Wall -Winline -Wno-long-long
INCLUDES      =	\
		-I../Single_Events \
		-I$(AMISICINC)/Tools \
		-I$(MODELINC)/Main \
		-I$(ATOOLSINC)/Phys \
		-I$(ATOOLSINC)/Math \
		-I$(ATOOLSINC)/Org 
noinst_HEADERS        =	\
		 	Lund_Wrapper.H \
			Lund_Interface.H 




