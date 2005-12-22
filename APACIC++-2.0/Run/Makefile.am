SHERPASYS  = ../..

SHERPAFLAGS   =	# -pedantic -Wall -Winline 

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
PDFLIBS       =	-L@libdir@ -lPDF -lMRST -lGRV -lCTEQ -lDUPDF \
		-lRemnant -lSudakov

if ROOT_SUPPORT
ROOTINCS      =	-I@ROOTINCDIR@
ROOTLDADD     =	-L@ROOTLIBDIR@ @ROOTGLIBS@ @ROOTAUXLIBS@
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

BEAMLIBS      =	-L@libdir@ -lBeam 
PHASICLIBS    =	-L@libdir@ -lPhasespace 
if ISAJET_SUPPORT
ISAJETLIB     = -L@libdir@ -lIsajet
endif
if HDECAY_SUPPORT
HDECAYLIB     = -L@libdir@ -lHdecay
endif
MODELLIBS     =	-L@libdir@ -lModelMain $(HDECAYLIB) $(ISAJETLIB)
EXTRAXSLIBS   =	-L@libdir@ -lExtraXS -lExtraXS2_2

APACICLIBS    =	-L@libdir@ -lApacicShowers -lApacicMain 

CTEQLIBS      = -L@libdir@ -lCTEQ 

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

bin_PROGRAMS  =	ApacicLEP ApacicTevatron
ApacicLEP_SOURCES        =	MainLEP.C 
ApacicLEP_LDADD  =	$(ANALYSISLIBS) $(SHERPALIBS) $(AMEGICLIBS) $(APACICLIBS) \
		$(EXTRAXSLIBS) \
		$(PHASICLIBS) $(MODELLIBS) $(PDFLIBS)  \
		$(BEAMLIBS) $(ATOOLSLIBS) $(ROOTLDADD) \
		$(SYSLIBS)
ApacicLEP_LDFLAGS     =	-rdynamic

ApacicTevatron_SOURCES = MainTevatron.C 
ApacicTevatron_LDADD  =	$(ANALYSISLIBS) $(SHERPALIBS) $(AMEGICLIBS) $(APACICLIBS) \
		$(ADICICLIBS) $(AMISICLIBS) $(AHADICLIBS) $(EXTRAXSLIBS) \
		$(PHASICLIBS) $(MODELLIBS) $(PDFLIBS) \
		$(CTEQLIBS) $(BEAMLIBS) $(ATOOLSLIBS) $(ROOTLDADD) \
		$(SYSLIBS) $(CLHEPLIBS)
ApacicTevatron_LDFLAGS        =	-rdynamic

AM_CXXFLAGS   =	$(SHERPAFLAGS) $(AHADICFLAGS) $(AMISICFLAGS) $(ADICICFLAGS) \
		$(CLHEPFLAGS) $(LHAPDFFLAGS) $(HADRONFLAGS) $(ROOTFLAGS)
INCLUDES      = $(ANALYSISINCS)  \
		$(APACICINCS)  \
		$(EXTRAXSINCS) $(PHASICINCS) $(MODELINCS) $(PDFINCS) \
		$(BEAMINCS) $(ATOOLSINCS) $(ROOTINCS)
EXTRA_DIST    = $(CVSINCLUDE)
dist-hook:
	sed -e "s/-rdynamic/SHERPA_LDFLAGS/1" \
	  < Makefile.am > $(distdir)/Makefile.am
