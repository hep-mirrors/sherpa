SHERPASYS  = ../..

ATOOLSINC  = $(SHERPASYS)/ATOOLS-2.0
MODELINC   = $(SHERPASYS)/MODEL-1.0
AMISICINC  = $(SHERPASYS)/AMISIC++-1.0

ROOTINC    = $(ROOTSYS)/include
ROOTLIB    = $(ROOTSYS)/lib

CLHEPINC   = $(CLHEPDIR)/include
CLHEPLIB   = $(CLHEPDIR)/lib

lib_LTLIBRARIES       =	libSherpaTools.la
libSherpaTools_la_SOURCES     =	\
				Kinematics.C \
				Output_Handler.C \
			    	HepMC_Interface.C \
				HepEvt_Interface.C \
				fhepevt.f \
			    	Sample_Analysis.C 
AM_CXXFLAGS   =	-pedantic -Wall -Winline -Wno-long-long
INCLUDES      =	\
		-I../Single_Events \
		-I$(AMISICINC)/Analysis \
		-I$(AMISICINC)/Tools \
		-I$(MODELINC)/Main \
		-I$(ATOOLSINC)/Phys \
		-I$(ATOOLSINC)/Math \
		-I$(ATOOLSINC)/Org \
		-I$(ROOTINC) \
		-I$(HEPMCINC) \
		-I$(CLHEPINC)
noinst_HEADERS        =	\
			Kinematics.H \
			Output_Handler.H \
		 	HepMC_Interface.H \
			HepEvt_Interface.H \
			Sample_Analysis.H 




