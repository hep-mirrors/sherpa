SHERPASYS  = ../..

SHERPAFLAGS   =	-pedantic -Wall -Winline 

AHADICDIR     =	$(SHERPASYS)/AHADIC++-0.0
AHADICFLAGS   = -DUSING__Ahadic
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ANALYSISDIR   =	$(SHERPASYS)/ANALYSIS-1.0
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
AMEGICDIR     =	$(SHERPASYS)/AMEGIC++-2.0
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0

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
AHADICINCS    =	-I$(AHADICDIR)/Main -I$(AHADICDIR)/Formation \
		-I$(AHADICDIR)/Decays -I$(AHADICDIR)/Tools
AMEGICINCS    =	-I$(AMEGICDIR)/Main -I$(AMEGICDIR)/Amplitude/Zfunctions \
		-I$(AMEGICDIR)/Amplitude/AmplTools -I$(AMEGICDIR)/Amplitude \
		-I$(AMEGICDIR)/Phasespace -I$(AMEGICDIR)/Model \
		-I$(AMEGICDIR)/String 
APACICINCS    = -I$(APACICDIR)/Main -I$(APACICDIR)/Showers

if HADRONS_SUPPORT
HADRONSLIBS   =	-L@libdir@ -lHadronsMain -lHadronsMEs -lHadronsPSs
HADRONSINCS   = -I$(HADRONSDIR)/Main
HADRONSDIR    = $(SHERPASYS)/HADRONS++-0.0
endif
if ADICIC_SUPPORT
ADICICLIBS    =	-L@libdir@ -lAdicicMain -lAdicicHdl -lAdicicStruct \
		-lAdicicTool -lAdicicInput
ADICICINCS    = -I$(ADICICDIR)/Main -I$(ADICICDIR)/Structure \
                -I$(ADICICDIR)/Handler -I$(ADICICDIR)/Tool -I$(ADICICDIR)/Input
ADICICDIR     =	$(SHERPASYS)/ADICIC++-0.0
endif
if MCATNLO_SUPPORT
MCATNLOLIBS   = -lHerwigTools
endif
if LHAPDF_SUPPORT
LHAPDFFLAGS    = -DLHAPDF_SUPPORT
LHAPDFLIBS     = -L@libdir@ -lLHAPDF $(LHAPDFDIR)/libLHAPDF.a 
else
CTEQLIBS      = -L@libdir@ -lCTEQ 
endif
if AMISIC_SUPPORT
AMISICINCS    =	-I$(SHERPASYS)/AMISIC++-1.0/Main
AMISICLIBS    =	-L@libdir@ -lAmisic -lAmisicModel -lAmisicTools
endif

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTLIBS      =	-L$(ROOTSYS)/lib -lCint -lCore -lHist -lGraf \
		-lGraf3d -lGpad -lTree -lRint -lPostscript \
		-lMatrix -lPhysics -lGui -ldl 
ROOTFLAGS     = -Wno-long-long
endif

if CLHEP_SUPPORT
CLHEPINCS     =	-I$(CLHEPDIR)/include
CLHEPLIBS     =	$(wildcard $(CLHEPDIR)/lib/libCLHEP-g++.*.a)
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif



SYSLIBS	      =	-lg2c -lm -ldl
AHADICLIBS    =	-L@libdir@ -lAhadicMain -lAhadicFormation \
		-lAhadicDecays -lAhadicTools

SHERPALIBS    =	-L@libdir@ -lSherpaMain -lSherpaInitialization \
		-lSherpaSingleEvents -lSherpaPerturbativePhysics \
		-lSherpaSoftPhysics -lLundTools $(MCATNLOLIBS) -lSherpaTools 
ATOOLSLIBS    =	-L@libdir@ -lToolsPhys -lToolsMath -lToolsOrg 
ANALYSISLIBS  = -L@libdir@ -lAnalysis -lObservables
BEAMLIBS      =	-L@libdir@ -lBeam 
PHASICLIBS    =	-L@libdir@ -lPhasespace 
PDFLIBS       =	-L@libdir@ -lDUPDF -lSudakov -lPDF -lMRST -lGRV -lRemnant
MODELLIBS     =	-L@libdir@ -lModelMain -lIsajet -lHdecay 
EXTRAXSLIBS   =	-L@libdir@ -lExtraXS -lExtraXS2_2 
AMEGICLIBS    =	-L@libdir@ -lAmegic -lAmplitude -lAmegicPSGen \
		-lZfunctions -lAmplTools -lModel -lString 
APACICLIBS    =	-L@libdir@ -lApacicShowers -lApacicMain 


bin_PROGRAMS  =	Cluster
Cluster_SOURCES = Main.C 
Cluster_LDADD  = $(AHADICLIBS) $(PHASICLIBS) $(ATOOLSLIBS) $(SYSLIBS) \
		 $(MODELLIBS) $(PDFLIBS) $(LHAPDFLIBS) $(CTEQLIBS) \
		 $(BEAMLIBS) $(APACICLIBS)
Cluster_LDFLAGS        =	SHERPA_LDFLAGS
AM_CXXFLAGS   =	$(AHADICFLAGS) $(ROOTFLAGS)
INCLUDES      = $(SHERPAINCS) $(ANALYSISINCS) $(HADRONSINCS) $(AHADICINCS) \
		$(APACICINCS) $(ADICICINCS) $(AMISICINCS) $(AMEGICINCS) \
		$(EXTRAXSINCS) $(PHASICINCS) $(MODELINCS) $(PDFINCS) \
		$(BEAMINCS) $(ATOOLSINCS) $(CLHEPINCS) $(ROOTINCS)
EXTRA_DIST    =	\
		$(CVSINCLUDE)

