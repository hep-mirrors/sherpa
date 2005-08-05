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
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main -I$(SHERPASYS)/PDF-1.0/CTEQ
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
PDFLIBS       =	-L@libdir@ -lPDF -lMRST -lGRV -lCTEQ -lDUPDF -lRemnant -lSudakov
#-lDUPDF


#ANALYSISLIBS  = -L@libdir@ -lAnalysis -lObservables
BEAMLIBS      =	-L@libdir@ -lBeam 
PHASICLIBS    =	-L@libdir@ -lPhasespace 
MODELLIBS     =	-L@libdir@ -lModelMain -lIsajet -lHdecay 
EXTRAXSLIBS   =	-L@libdir@ -lExtraXSns -lExtraXS2_2 
#AMEGICLIBS    =	-L@libdir@ -lAmegic -lAmplitude -lAmegicPSGen \
#		-lZfunctions -lAmplTools -lModel -lString 

APACICLIBS    =	-L@libdir@ -lApacicShowersNS -lApacicMain 
#SHERPALIBS    =	-L@libdir@ -lSherpaMain -lSherpaInitialization \
#		-lSherpaSingleEvents -lSherpaPerturbativePhysics \
#		-lSherpaSoftPhysics -lLundTools $(MCATNLOLIBS) -lSherpaTools 
#

#

CTEQLIBS      = -L@libdir@ -lCTEQ 

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

bin_PROGRAMS  =	ApacicLEP ApacicTevatron
ApacicLEP_SOURCES        =	MainLEP.C 
ApacicLEP_LDADD  =	$(ANALYSISLIBS) $(SHERPALIBS) $(AMEGICLIBS) $(APACICLIBS) \
		$(EXTRAXSLIBS) \
		$(PHASICLIBS) $(MODELLIBS) $(PDFLIBS)  \
		$(BEAMLIBS) $(ATOOLSLIBS) \
		$(SYSLIBS) $(ROOTLIBS)
ApacicLEP_LDFLAGS        =	-rdynamic

ApacicTevatron_SOURCES = MainTevatron.C 
ApacicTevatron_LDADD  =	$(ANALYSISLIBS) $(SHERPALIBS) $(AMEGICLIBS) $(APACICLIBS) \
		$(ADICICLIBS) $(AMISICLIBS) $(AHADICLIBS) $(EXTRAXSLIBS) \
		$(PHASICLIBS) $(MODELLIBS) $(PDFLIBS) \
		$(CTEQLIBS) $(BEAMLIBS) $(ATOOLSLIBS) \
		$(SYSLIBS) $(CLHEPLIBS) $(ROOTLIBS)
ApacicTevatron_LDFLAGS = -rdynamic

AM_CXXFLAGS   =	$(SHERPAFLAGS) $(AHADICFLAGS) $(AMISICFLAGS) $(ADICICFLAGS) \
		$(CLHEPFLAGS) $(ROOTFLAGS) $(LHAPDFFLAGS) $(HADRONFLAGS)
INCLUDES      = $(ANALYSISINCS)  \
		$(APACICINCS)  \
		$(EXTRAXSINCS) $(PHASICINCS) $(MODELINCS) $(PDFINCS) \
		$(BEAMINCS) $(ATOOLSINCS) $(ROOTINCS)
