SHERPASYS     =	../..

SHERPAFLAGS   =	-pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
AMEGICDIR     =	$(SHERPASYS)/AMEGIC++-2.0
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0
SHERPAINCS    =	-I../PerturbativePhysics -I../SoftPhysics \
		-I../SoftPhysics/Lund  -I../Tools 
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
EXTRAXSINCS   =	-I$(EXTRAXSDIR)/Two2Two -I$(EXTRAXSDIR)/Main
AMISICINCS    =	-I$(SHERPASYS)/AMISIC++-1.0/Main
AMEGICINCS    =	-I$(AMEGICDIR)/Main -I$(AMEGICDIR)/Amplitude/Zfunctions \
		-I$(AMEGICDIR)/Amplitude/AmplTools -I$(AMEGICDIR)/Amplitude \
		-I$(AMEGICDIR)/Phasespace -I$(AMEGICDIR)/Model \
		-I$(AMEGICDIR)/String 
APACICINCS    = -I$(APACICDIR)/Main -I$(APACICDIR)/Showers

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       =	libSherpaInitialization.la
libSherpaInitialization_la_SOURCES    =	Initialization_Handler.C 
noinst_HEADERS        =	Initialization_Handler.H 
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(SHERPAFLAGS)
INCLUDES      = $(SHERPAINCS) $(APACICINCS) $(AMEGICINCS) $(AMISICINCS) \
		$(EXTRAXSINCS) $(PHASICINCS) $(MODELINCS) $(PDFINCS) \
		$(BEAMINCS) $(ATOOLSINCS)



