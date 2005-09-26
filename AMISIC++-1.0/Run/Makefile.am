SHERPASYS  = ../..

AMISICFLAGS   =	-pedantic -Wall -Winline 

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0

AMISICINCS    =	-I../Main  
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
EXTRAXSINCS   =	-I$(EXTRAXSDIR)/Two2Two -I$(EXTRAXSDIR)/Main
APACICINCS    = -I$(APACICDIR)/Main -I$(APACICDIR)/Showers

SYSLIBS	      =	-lg2c -lm -ldl
ATOOLSLIBS    =	-L@libdir@ -lToolsPhys -lToolsMath -lToolsOrg 
PDFLIBS       =	-L@libdir@ -lPDF -lDUPDF -lSudakov -lLHAPDF -lCTEQ -lMRST -lGRV
BEAMLIBS      =	-L@libdir@ -lBeam 
PHASICLIBS    =	-L@libdir@ -lPhasespace -lFoam
if ISAJET_SUPPORT
ISAJETLIB     = -L@libdir@ -lIsajet
endif
if HDECAY_SUPPORT
HDECAYLIB     = -L@libdir@ -lHdecay
endif
MODELLIBS     =	-L@libdir@ -lModelMain $(HDECAYLIB) $(ISAJETLIB)

EXTRAXSLIBS   =	-L@libdir@ -lExtraXS -lExtraXS2_2 
AMISICLIBS    =	-L@libdir@ -lAmisic -lAmisicModel -lAmisicTools
APACICLIBS    =	-L@libdir@ -lApacicShowers -lApacicMain 

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTLIBS      =	-L$(ROOTSYS)/lib -lCint -lCore -lHist -lGraf \
		-lGraf3d -lGpad -lTree -lRint -lPostscript \
		-lMatrix -lPhysics -lGui -ldl 
ROOTFLAGS     = -Wno-long-long
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if SHERPA_SUPPORT
SOURCE	      = Main.C 
else
bin_PROGRAMS  =	Amisic
Amisic_SOURCES        =	Main.C 
Amisic_LDADD  =	$(AMISICLIBS) $(APACICLIBS) $(EXTRAXSLIBS) \
		$(PHASICLIBS) $(MODELLIBS) $(PDFLIBS) \
		$(BEAMLIBS) $(ATOOLSLIBS) $(SYSLIBS) $(CLHEPLIBS) $(ROOTLIBS) 
Amisic_LDFLAGS        =	-rdynamic
AM_CXXFLAGS   =	$(AMISICFLAGS) $(ROOTFLAGS)
INCLUDES      = $(AMISICINCS) $(APACICINCS) $(EXTRAXSINCS) \
		$(PHASICINCS) $(MODELINCS) $(PDFINCS) \
		$(BEAMINCS) $(ATOOLSINCS) $(CLHEPINCS) $(ROOTINCS) 
endif

EXTRA_DIST    =	\
		$(CVSINCLUDE) $(SOURCE) \
		Particle.dat Hadron.dat Integration.dat Model.dat Beam.dat \
		ISR.dat MI.dat simplechain.dat Run.dat 
dist-hook:
	if test -f Makefile_am ; then \
	  cp Makefile_am $(distdir)/Makefile.am; \
	fi
