SHERPASYS  = ../..

SHERPAFLAGS   =	-pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
AMEGICDIR     =	$(SHERPASYS)/AMEGIC++-2.0
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0
SHERPAINCS    =	-I../Main -I../Initialization -I../PerturbativePhysics \
		-I../SoftPhysics -I../LundTools -I../Single_Events \
		-I../Tools 
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
EXTRAXSINCS   =	-I$(EXTRAXSDIR)/Two2Two -I$(EXTRAXSDIR)/Main
if AMISIC_SUPPORT
AMISICINCS    =	-I$(SHERPASYS)/AMISIC++-1.0/Main
endif
AMEGICINCS    =	-I$(AMEGICDIR)/Main -I$(AMEGICDIR)/Amplitude/Zfunctions \
		-I$(AMEGICDIR)/Amplitude/AmplTools -I$(AMEGICDIR)/Amplitude \
		-I$(AMEGICDIR)/Phasespace -I$(AMEGICDIR)/Model \
		-I$(AMEGICDIR)/String 
APACICINCS    = -I$(APACICDIR)/Main -I$(APACICDIR)/Showers

SYSLIBS	      =	-lg2c -lm -ldl
ATOOLSLIBS    =	-L@libdir@ -lToolsPhys -lToolsMath -lToolsOrg 
PDFLIBS       =	-L@libdir@ -lPDF -lLHAPDF -lDUPDF
BEAMLIBS      =	-L@libdir@ -lBeam 
PHASICLIBS    =	-L@libdir@ -lPhasespace 
MODELLIBS     =	-L@libdir@ -lModelMain -lIsajet -lHdecay 
EXTRAXSLIBS   =	-L@libdir@ -lExtraXS -lExtraXS2_2 
if AMISIC_SUPPORT
AMISICLIBS    =	-L@libdir@ -lAmisic -lAmisicAnalysis
endif
AMEGICLIBS    =	-L@libdir@ -lAmegic -lAmplitude -lAmegicPSGen \
		-lZfunctions -lAmplTools -lModel -lProcess -lString 
APACICLIBS    =	-L@libdir@ -lApacicShowers -lApacicMain 
SHERPALIBS    =	-L@libdir@ -lSherpaMain -lSherpaInitialization \
		-lSherpaSingleEvents -lSherpaPerturbativePhysics \
		-lSherpaSoftPhysics -lLundTools -lSherpaTools 

if AMISIC_SUPPORT
AMISICFLAGS   =	-DUSING__Amisic
endif

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTLIBS      =	-L$(ROOTSYS)/lib -lCint -lCore -lHist -lGraf \
		-lGraf3d -lGpad -lTree -lRint -lPostscript \
		-lMatrix -lPhysics -lGui -ldl 
endif

if CLHEP_SUPPORT
CLHEPINCS     =	-I$(CLHEPDIR)/include
CLHEPLIBS     =	$(wildcard $(CLHEPDIR)/lib/libCLHEP-g++.*.a)
CLHEPFLAGS    = -DCLHEP_SUPPORT
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

bin_PROGRAMS  =	Sherpa
Sherpa_SOURCES        =	Main.C 
Sherpa_LDADD  =	$(SHERPALIBS) $(AMEGICLIBS) $(APACICLIBS) $(AMISICLIBS) \
		$(EXTRAXSLIBS) $(PHASICLIBS) $(MODELLIBS) $(PDFLIBS) $(BEAMLIBS) \
		$(ATOOLSLIBS) $(SYSLIBS) $(CLHEPLIBS) $(ROOTLIBS)
Sherpa_LDFLAGS        =	-rdynamic
AM_CXXFLAGS   =	$(SHERPAFLAGS) $(AMISICFLAGS) $(CLHEPFLAGS)
INCLUDES      = $(SHERPAINCS) $(APACICINCS) $(AMEGICINCS) $(AMISICINCS) \
		$(EXTRAXSINCS) $(PHASICINCS) $(MODELINCS) $(PDFINCS) \
		$(BEAMINCS) $(ATOOLSINCS) $(CLHEPINCS) $(ROOTINCS)
SUBDIRS       =	Process
EXTRA_DIST    =	\
		$(CVSINCLUDE) \
		Integration.dat Processes.dat Selector.dat Model.dat Beam.dat \
	     	Particle.dat Hadron.dat Run.dat ISR.dat ME.dat MI.dat simplechain.dat \
		Decays.dat Shower.dat Fragmentation.dat makeclean makelibs Lund.dat
#	     	PDFsets.tar.gz MRST99Grid.tar.gz CTEQ6Grid.tar.gz

