SHERPASYS  = ../..

ATOOLSINC  = $(SHERPASYS)/ATOOLS-2.0
PDFINC     = $(SHERPASYS)/PDF-1.0
BEAMINC    = $(SHERPASYS)/BEAM-1.0
PHASICINC  = $(SHERPASYS)/PHASIC++-1.0
MODELINC   = $(SHERPASYS)/MODEL-1.0
EXTRAXSINC = $(SHERPASYS)/EXTRA_XS-1.0
AMISICINC  = $(SHERPASYS)/AMISIC++-1.0
AMEGICINC  = $(SHERPASYS)/AMEGIC++-2.0
APACICINC  = $(SHERPASYS)/APACIC++-2.0

ATOOLSLIB  = @libdir@
PDFLIB     = @libdir@
BEAMLIB    = @libdir@
PHASICLIB  = @libdir@
MODELLIB   = @libdir@
EXTRAXSLIB = @libdir@
AMISICLIB  = @libdir@
AMEGICLIB  = @libdir@
APACICLIB  = @libdir@
SHERPALIB  = @libdir@

ROOTINC    = $(ROOTSYS)/include
ROOTLIB    = $(ROOTSYS)/lib

CLHEPLIB   = $(CLHEPSYS)/lib
CLHEPINC   = $(CLHEPSYS)/include

bin_PROGRAMS  =	Sherpa
Sherpa_SOURCES        =	Main.C 
Sherpa_LDADD  =	\
		-L$(SHERPALIB) -lSherpaMain \
		-lSherpaInitialization -lSherpaSingleEvents \
		-lSherpaPerturbativePhysics -lSherpaSoftPhysics \
		-lLundInterface -lSherpaTools \
		-L$(AMEGICLIB) -lAmegic -lAmplitude -lAmegicPSGen \
		-lZfunctions -lAmplTools -lModel -lProcess -lString \
		-L$(APACICLIB) -lApacicShowers -lApacicMain \
		-L$(EXTRAXSLIB) -lExtraXS -lExtraXS2_2 \
		-L$(PHASICLIB) -lPhasespace \
		-L$(MODELLIB) -lModelMain -lIsajet -lHdecay \
		-L$(PDFLIB) -lPDF -lLHAPDF \
		-L$(BEAMLIB) -lBeam \
		-L$(ATOOLSLIB) -lToolsPhys -lToolsMath -lToolsOrg \
		-lg2c -lm -ldl
#               -L$(AMISICLIB) -lAmisicTools -lAmisicHelpers \
#		-lAmisic -lAmisicAnalysis \
#		-L$(ROOTLIB) -lCint -lCore -lHist -lGraf \
#		-lGraf3d -lGpad -lTree -lRint -lPostscript \
#		-lMatrix -lPhysics -lGui -ldl \
Sherpa_LDFLAGS        =	-rdynamic
INCLUDES      = \
		-I../Main \
		-I../Initialization \
		-I../PerturbativePhysics \
		-I../SoftPhysics \
		-I../SoftPhysics/Lund \
		-I../Single_Events \
		-I../Tools \
		-I$(APACICINC)/Main \
		-I$(APACICINC)/Showers \
		-I$(AMEGICINC)/Amplitude \
		-I$(AMEGICINC)/Amplitude/AmplTools \
		-I$(AMEGICINC)/Amplitude/Zfunctions \
		-I$(AMEGICINC)/Phasespace \
		-I$(AMEGICINC)/Main \
		-I$(AMEGICINC)/Model \
		-I$(AMEGICINC)/String \
		-I$(AMISICINC)/Main \
		-I$(AMISICINC)/Other \
		-I$(AMISICINC)/Tools \
		-I$(AMISICINC)/Tools \
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
SUBDIRS       =	Process
EXTRA_DIST    =	\
		Integration.dat Processes.dat Selector.dat Model.dat Beam.dat \
	     	Particle.dat Hadron.dat Run.dat ISR.dat ME.dat MI.dat simplechain.dat \
		Decays.dat Shower.dat Fragmentation.dat makeclean makelibs Lund.dat
#	     	PDFsets.tar.gz MRST99Grid.tar.gz CTEQ6Grid.tar.gz

