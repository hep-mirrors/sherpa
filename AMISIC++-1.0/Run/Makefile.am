SHERPASYS  = ../..

ATOOLSINC  = $(SHERPASYS)/ATOOLS-2.0
PDFINC     = $(SHERPASYS)/PDF-1.0
BEAMINC    = $(SHERPASYS)/BEAM-1.0
PHASICINC  = $(SHERPASYS)/PHASIC++-1.0
MODELINC   = $(SHERPASYS)/MODEL-1.0
EXTRAXSINC = $(SHERPASYS)/EXTRA_XS-1.0

ATOOLSLIB  = $(SHERPASYS)/lib
PDFLIB     = $(SHERPASYS)/lib
BEAMLIB    = $(SHERPASYS)/lib
PHASICLIB  = $(SHERPASYS)/lib
MODELLIB   = $(SHERPASYS)/lib
EXTRAXSLIB = $(SHERPASYS)/lib

ROOTINC    = $(ROOTSYS)/include
ROOTLIB    = $(ROOTSYS)/lib

bin_PROGRAMS  =	Amisic
Amisic_SOURCES	      =	\
			Main.C 
Amisic_LDADD  =	\
		@libdir@/libAmisic.la \
		@libdir@/libAmisicHelpers.la \
		@libdir@/libAmisicTools.la \
		-L$(EXTRAXSLIB) -lExtraXS -lExtraXS2_2 \
		-L$(PHASICLIB) -lPhasespace \
		-L$(MODELLIB) -lModelMain -lIsajet -lHdecay \
		-L$(PDFLIB) -lPDF -lLHAPDF \
		-L$(BEAMLIB) -lBeam \
		-L$(ATOOLSLIB) -lToolsPhys -lToolsMath -lToolsOrg \
		-L$(ROOTLIB) -lCint -lCore -lHist -lGraf \
		-lGraf3d -lGpad -lTree -lRint -lPostscript \
		-lMatrix -lPhysics -lGui -ldl \
		-lg2c -lm 
CXXFLAGS      =	-pedantic -Wall -Winline -Wno-long-long
INCLUDES      = \
		-I../Main \
		-I../Tools \
		-I../Other \
		-I$(EXTRAXSINC)/Two2Two \
		-I$(EXTRAXSINC)/Main \
		-I$(PHASICINC)/Main \
		-I$(MODELINC)/Main \
		-I$(PDFINC)/Main \
		-I$(BEAMINC)/Main \
		-I$(ATOOLSINC)/Phys \
		-I$(ATOOLSINC)/Math \
		-I$(ATOOLSINC)/Org \
		-I$(ROOTINC) 
EXTRA_DIST    =	\
		Run.dat \
		Beam.dat \
		ISR.dat \
		Model.dat \
		Integration.dat \
		Particle.dat \
		Hadron.dat \
		Processes.dat \
		XS.dat \
		XSSelector.dat \
		MI.dat 
