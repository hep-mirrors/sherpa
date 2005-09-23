SHERPASYS  = ../..

SHERPAFLAGS   =	-pedantic -Wall -Winline 

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
PDFDIR        =	$(SHERPASYS)/PDF-1.0
PDFINCS       =	-I$(PDFDIR)/Main -I$(PDFDIR)/Remnant
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
EXTRAXSINCS   =	-I$(EXTRAXSDIR)/Two2Two -I$(EXTRAXSDIR)/Main
AMEGICINCS    =	-I$(AMEGICDIR)/Main -I$(AMEGICDIR)/Amplitude/Zfunctions \
		-I$(AMEGICDIR)/Amplitude/AmplTools -I$(AMEGICDIR)/Amplitude \
		-I$(AMEGICDIR)/Phasespace -I$(AMEGICDIR)/Model \
		-I$(AMEGICDIR)/String 
APACICINCS    = -I$(APACICDIR)/Main -I$(APACICDIR)/Showers

SYSLIBS	      =	-lg2c -lm -ldl
ATOOLSLIBS    =	-L@libdir@ -lToolsPhys -lToolsMath -lToolsOrg 
PDFLIBS       =	-L@libdir@ -lPDF -lDUPDF -lSudakov -lMRST -lGRV -lRemnant
ANALYSISLIBS  = -L@libdir@ -lAnalysis -lObservables
BEAMLIBS      =	-L@libdir@ -lBeam 
PHASICLIBS    =	-L@libdir@ -lPhasespace 
MODELLIBS     =	-L@libdir@ -lModelMain -lIsajet -lHdecay 
EXTRAXSLIBS   =	-L@libdir@ -lExtraXS -lExtraXS2_2 
AMEGICLIBS    =	-L@libdir@ -lAmegic -lAmplitude -lAmegicPSGen \
		-lZfunctions -lModel -lString 
APACICLIBS    =	-L@libdir@ -lApacicShowers -lApacicMain 
SHERPALIBS    =	-L@libdir@ -lSherpaMain -lSherpaInitialization \
		-lSherpaSingleEvents -lSherpaPerturbativePhysics \
		-lSherpaSoftPhysics -lLundTools $(MCATNLOLIBS) -lSherpaTools 

if MCATNLO_SUPPORT
MCATNLOLIBS   = -lHerwigTools
MCATNLOFLAGS  = -DUSING__MCatNLO
endif

if HADRONS_SUPPORT
HADRONSDIR    = $(SHERPASYS)/HADRONS++-0.0
HADRONSINCS   = -I$(HADRONSDIR)/Main
HADRONSLIBS   =	-L@libdir@ -lHadronsMain -lHadronsMEs -lHadronsPSs
endif
if AHADIC_SUPPORT
AHADICDIR     =	$(SHERPASYS)/AHADIC++-0.0
AHADICINCS    =	-I$(AHADICDIR)/Main -I$(AHADICDIR)/Formation \
		-I$(AHADICDIR)/Decays -I$(AHADICDIR)/Tools
AHADICLIBS    =	-L@libdir@ -lAhadicMain -lAhadicFormation \
		-lAhadicDecays -lAhadicTools
endif
if AMISIC_SUPPORT
AMISICINCS    =	-I$(SHERPASYS)/AMISIC++-1.0/Main
AMISICLIBS    =	-L@libdir@ -lAmisic -lAmisicModel -lAmisicTools
endif
if ADICIC_SUPPORT
ADICICDIR     =	$(SHERPASYS)/ADICIC++-0.0
ADICICINCS    = -I$(ADICICDIR)/Main -I$(ADICICDIR)/Structure \
                -I$(ADICICDIR)/Handler -I$(ADICICDIR)/Tool -I$(ADICICDIR)/Input
ADICICLIBS    =	-L@libdir@ -lAdicicMain -lAdicicHdl -lAdicicStruct \
		-lAdicicTool -lAdicicInput
endif
if CSS_SUPPORT
CSSDIR        = $(SHERPASYS)/CS_SHOWER++-0.0
CSSINCS       = -I$(CSSDIR)/Tools -I$(CSSDIR)/Showers -I$(CSSDIR)/Main
CSSLIBS       = -L@libdir@ -lCSSTools -lCSShowers -lCSMain
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

if LHAPDF_SUPPORT
LHAPDFLIBS    = -L@libdir@ -lLHAPDF $(LHAPDFDIR)/libLHAPDF.a 
else
CTEQLIBS      = -L@libdir@ -lCTEQ 
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

bin_PROGRAMS  =	Sherpa
Sherpa_SOURCES        =	Main.C 
Sherpa_LDADD  =	$(ANALYSISLIBS) $(SHERPALIBS) $(AMEGICLIBS) $(APACICLIBS) \
		$(ADICICLIBS) $(CSSLIBS) $(AMISICLIBS) $(AHADICLIBS) \
		$(EXTRAXSLIBS) $(PHASICLIBS) $(MODELLIBS) $(BEAMLIBS) \
		$(ATOOLSLIBS) $(HADRONSLIBS) $(PDFLIBS) $(LHAPDFLIBS) \
		$(CTEQLIBS) $(SYSLIBS) $(CLHEPLIBS) $(ROOTLIBS)
Sherpa_LDFLAGS        =	SHERPA_LDFLAGS
AM_CXXFLAGS   =	$(SHERPAFLAGS) $(ROOTFLAGS) 
INCLUDES      = $(SHERPAINCS) $(ANALYSISINCS) $(HADRONSINCS) $(AHADICINCS) \
		$(APACICINCS) $(ADICICINCS) $(AMISICINCS) $(AMEGICINCS) \
		$(EXTRAXSINCS) $(PHASICINCS) $(CSSINCS) $(MODELINCS) \
		$(PDFINCS) $(BEAMINCS) $(ATOOLSINCS) $(CLHEPINCS) $(ROOTINCS)
dist_bin_SCRIPTS      =	makelibs
EXTRA_DIST    =	\
		$(CVSINCLUDE) \
		Integration.dat Processes.dat Selector.dat Model.dat Beam.dat \
	     	Particle.dat Hadron.dat Run.dat ISR.dat ME.dat MI.dat \
		Decays.dat Shower.dat Fragmentation.dat \
		Lund.dat Cluster.dat

#	     	PDFsets.tar.gz MRST99Grid.tar.gz CTEQ6Grid.tar.gz
dist-hook:
	sed -e "s/SHERPA_LDFLAGS/SHERPA_LDFLAGS/1" \
	  < Makefile.am > $(distdir)/Makefile.am
