SHERPASYS  = ../..

SHERPAFLAGS   =	-pedantic -Wall -Winline 

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ANALYSISDIR   =	$(SHERPASYS)/ANALYSIS-1.0
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
AMEGICDIR     =	$(SHERPASYS)/AMEGIC++-2.0
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0
HADRONSDIR    = $(SHERPASYS)/HADRONS++-0.0
if ADICIC_SUPPORT
ADICICDIR     =	$(SHERPASYS)/ADICIC++-0.0
ADICICFLAGS   = -DUSING__Adicic
endif
if AHADIC_SUPPORT
AHADICDIR     =	$(SHERPASYS)/AHADIC++-0.0
AHADICFLAGS   = -DUSING__Ahadic
endif
if MCATNLO_SUPPORT
MCATNLOLIBS   = -lHerwigTools
MCATNLOFLAGS  = -DUSING__MCatNLO
endif

SHERPAINCS    =	-I../Main -I../Initialization -I../PerturbativePhysics \
		-I../SoftPhysics -I../LundTools -I../Single_Events \
		-I../HerwigTools -I../Tools 
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
ANALYSISINCS  = -I$(ANALYSISDIR)/Main -I$(ANALYSISDIR)/Observables
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
HADRONSINCS   = -I$(HADRONSDIR)/Main
if ADICIC_SUPPORT
ADICICINCS    = -I$(ADICICDIR)/Main -I$(ADICICDIR)/Structure \
                -I$(ADICICDIR)/Handler -I$(ADICICDIR)/Tool -I$(ADICICDIR)/Input
endif
if AHADIC_SUPPORT
AHADICINCS    =	-I$(AHADICDIR)/Main -I$(AHADICDIR)/ClusterFormation \
		-I$(AHADICDIR)/ClusterDecays
endif

SYSLIBS	      =	-lg2c -lm -ldl
ATOOLSLIBS    =	-L@libdir@ -lToolsPhys -lToolsMath -lToolsOrg 
PDFLIBS       =	-L@libdir@ -lPDF -lDUPDF -lSudakov -lMRST -lGRV 



ANALYSISLIBS  = -L@libdir@ -lAnalysis -lObservables
BEAMLIBS      =	-L@libdir@ -lBeam 
PHASICLIBS    =	-L@libdir@ -lPhasespace 
MODELLIBS     =	-L@libdir@ -lModelMain -lIsajet -lHdecay 
EXTRAXSLIBS   =	-L@libdir@ -lExtraXS -lExtraXS2_2 
if AMISIC_SUPPORT
AMISICLIBS    =	-L@libdir@ -lAmisic -lAmisicModel -lAmisicTools
endif
AMEGICLIBS    =	-L@libdir@ -lAmegic -lAmplitude -lAmegicPSGen \
		-lZfunctions -lAmplTools -lModel -lString 
APACICLIBS    =	-L@libdir@ -lApacicShowers -lApacicMain 
HADRONSLIBS   =	-L@libdir@ -lHadronsMain 
if ADICIC_SUPPORT
ADICICLIBS    =	-L@libdir@ -lAdicicMain -lAdicicHdl -lAdicicStruct \
		-lAdicicTool -lAdicicInput
endif
if AHADIC_SUPPORT
AHADICLIBS    =	-L@libdir@ -lAhadicMain -lClusterFormation -lClusterDecays
endif
SHERPALIBS    =	-L@libdir@ -lSherpaMain -lSherpaInitialization \
		-lSherpaSingleEvents -lSherpaPerturbativePhysics \
		-lSherpaSoftPhysics -lLundTools $(MCATNLOLIBS) -lSherpaTools 

if AMISIC_SUPPORT
AMISICFLAGS   =	-DUSING__Amisic
endif

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTLIBS      =	-L$(ROOTSYS)/lib -lCint -lCore -lHist -lGraf \
		-lGraf3d -lGpad -lTree -lRint -lPostscript \
		-lMatrix -lPhysics -lGui -ldl 
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

if CLHEP_SUPPORT
CLHEPINCS     =	-I$(CLHEPDIR)/include
CLHEPLIBS     =	$(wildcard $(CLHEPDIR)/lib/libCLHEP-g++.*.a)
CLHEPFLAGS    = -DCLHEP_SUPPORT
endif

if LHAPDF_SUPPORT
LHAPDFFLAGS    = -DLHAPDF_SUPPORT
LHAPDFLIBS     = -L@libdir@ -lLHAPDF $(LHAPDFDIR)/libLHAPDF.a 
else
CTEQLIBS      = -L@libdir@ -lCTEQ 
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

bin_PROGRAMS  =	Sherpa
Sherpa_SOURCES        =	Main.C 
Sherpa_LDADD  =	$(ANALYSISLIBS) $(SHERPALIBS) $(AMEGICLIBS) $(APACICLIBS) $(ADICICLIBS) \
		$(AMISICLIBS) $(AHADICLIBS) $(EXTRAXSLIBS) $(PHASICLIBS) $(MODELLIBS) \
		$(PDFLIBS) $(LHAPDFLIBS) $(CTEQLIBS) $(BEAMLIBS) $(ATOOLSLIBS) $(HADRONSLIBS)\
		$(SYSLIBS) $(CLHEPLIBS) $(ROOTLIBS) 
Sherpa_LDFLAGS        =	-rdynamic
AM_CXXFLAGS   =	$(SHERPAFLAGS) $(AHADICFLAGS) $(AMISICFLAGS) $(ADICICFLAGS) \
		$(CLHEPFLAGS) $(ROOTFLAGS) $(LHAPDFFLAGS) 
INCLUDES      = $(SHERPAINCS) $(ANALYSISINCS) $(AHADICINCS) $(APACICINCS) $(ADICICINCS) \
		$(AMISICINCS) $(AMEGICINCS) $(EXTRAXSINCS) $(PHASICINCS) $(MODELINCS) \
		$(PDFINCS) $(BEAMINCS) $(ATOOLSINCS) $(CLHEPINCS) $(ROOTINCS) \
		$(HADRONSINCS)
EXTRA_DIST    =	\
		$(CVSINCLUDE) \
		Integration.dat Processes.dat Selector.dat Model.dat Beam.dat \
	     	Particle.dat Hadron.dat Run.dat ISR.dat ME.dat MI.dat \
		Decays.dat Shower.dat Fragmentation.dat makeclean makelibs Lund.dat \
		Cluster.dat
#	     	PDFsets.tar.gz MRST99Grid.tar.gz CTEQ6Grid.tar.gz

