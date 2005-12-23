SHERPASYS  = ../..

SHERPAFLAGS   =	# -pedantic -Wall -Winline 
DEFS          =

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
ATOOLSLIBS    =	-lToolsPhys -lToolsMath -lToolsOrg 
PDFLIBS       =	-lPDF -lDUPDF -lSudakov -lMRST -lGRV -lRemnant
ANALYSISLIBS  = -lAnalysis -lObservables
BEAMLIBS      =	-lBeam 
PHASICLIBS    =	-lPhasespace 
MODELLIBS     =	-lModelMain
EXTRAXSLIBS   =	-lExtraXS -lExtraXS2_2
AMEGICLIBS    =	-lAmegic -lAmplitude -lAmegicPSGen \
		-lZfunctions -lModel -lString 
APACICLIBS    =	-lApacicShowers -lApacicMain 
SHERPALIBS    =	-lSherpaMain -lSherpaInitialization \
		-lSherpaSingleEvents -lSherpaPerturbativePhysics \
		-lSherpaSoftPhysics -lLundTools -lSherpaTools 

if HADRONS_SUPPORT
HADRONSDIR    = $(SHERPASYS)/HADRONS++-0.0
HADRONSINCS   = -I$(HADRONSDIR)/Main
HADRONSLIBS   =	-lHadronsMain -lHadronsMEs -lHadronsPSs
endif
if AMISIC_SUPPORT
AMISICINCS    =	-I$(SHERPASYS)/AMISIC++-1.0/Main
AMISICLIBS    =	-lAmisic -lAmisicModel -lAmisicTools
endif

if ROOT_SUPPORT
ROOTINCS      =	-I@ROOTINCDIR@
ROOTLDADD     =	-L@ROOTLIBDIR@ @ROOTGLIBS@ @ROOTAUXLIBS@
ROOTFLAGS     = -Wno-long-long
endif

if CLHEP_SUPPORT
CLHEPINCS     =	-I$(CLHEPDIR)/include
CLHEPLIBS     =	$(wildcard $(CLHEPDIR)/lib/libCLHEP-g++.*.a)
endif

if LHAPDF_SUPPORT
LHAPDFLIBS    = -lLHAPDF $(LHAPDFDIR)/libLHAPDF.a 
else
CTEQLIBS      = -lCTEQ 
endif

if ISAJET_SUPPORT
ISAJETLIB     = -L@libdir@ -lIsajet
endif

if HDECAY_SUPPORT
HDECAYLIB     = -L@libdir@ -lHdecay
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

bin_PROGRAMS  =	Sherpa
Sherpa_SOURCES        =	Main.C 
Sherpa_LDADD  =	-L@libdir@ $(ANALYSISLIBS) $(SHERPALIBS) $(AMEGICLIBS) \
		$(APACICLIBS) $(AMISICLIBS) \
		$(EXTRAXSLIBS) $(PHASICLIBS) $(MODELLIBS) \
		$(BEAMLIBS) $(ATOOLSLIBS) $(HADRONSLIBS) $(PDFLIBS) \
		$(LHAPDFLIBS) $(CTEQLIBS) $(SYSLIBS) $(CLHEPLIBS) \
		$(ROOTLDADD) $(ISAJETLIB) $(HDECAYLIB)
Sherpa_LDFLAGS        =	-dynamic -flat_namespace
AM_CXXFLAGS   =	$(SHERPAFLAGS) $(ROOTFLAGS) 
INCLUDES      = $(SHERPAINCS) $(ANALYSISINCS) $(HADRONSINCS) \
		$(APACICINCS) $(AMISICINCS) $(AMEGICINCS) \
		$(EXTRAXSINCS) $(PHASICINCS) $(MODELINCS) \
		$(PDFINCS) $(BEAMINCS) $(ATOOLSINCS) $(CLHEPINCS) $(ROOTINCS)
dist_bin_SCRIPTS      =	makelibs
nobase_data_DATA      = \
			CTEQ6Grid/cteq6l1.tbl \
			CTEQ6Grid/cteq6l.tbl \
			CTEQ6Grid/cteq6m.tbl \
			CTEQ6Grid/cteq6d.tbl \
			MRST99Grid/cor01.dat \
			MRST99Grid/cor02.dat \
			MRST99Grid/cor03.dat \
			MRST99Grid/cor04.dat \
			MRST99Grid/cor05.dat \
			MRST99Grid/cor06.dat \
			MRST99Grid/cor07.dat \
			MRST99Grid/cor08.dat \
			MRST99Grid/cor09.dat \
			MRST99Grid/cor10.dat \
			MRST99Grid/cor11.dat \
			MRST99Grid/cor12.dat \
			MRST01Grid/lo2002.dat
EXTRA_DIST    =	\
		$(CVSINCLUDE) \
		Integration.dat Processes.dat Selector.dat Model.dat Beam.dat \
	     	Particle.dat Hadron.dat Run.dat ISR.dat ME.dat MI.dat \
		Decays.dat Shower.dat Fragmentation.dat \
		Lund.dat Cluster.dat

#	     	PDFsets.tar.gz MRST99Grid.tar.gz CTEQ6Grid.tar.gz

dist-hook:
	sed -e "s/-dynamic -flat_namespace/SHERPA_LDFLAGS/1" \
	  < Makefile.am > $(distdir)/Makefile.am
