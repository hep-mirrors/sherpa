SHERPASYS     =	../..

SHERPAFLAGS   =	# -pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
AMEGICDIR     =	$(SHERPASYS)/AMEGIC++-2.0
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0
SHERPAINCS    =	-I../Initialization -I../PerturbativePhysics \
		-I../SoftPhysics -I../SoftPhysics/Lund -I../Single_Events \
		-I../LundTools -I../Tools 
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

if CLHEP_SUPPORT
CLHEPINCS     =	-I$(CLHEPDIR)/include
CLHEPFLAGS    = -DCLHEP_SUPPORT
endif

lib_LTLIBRARIES       =	libSherpaMain.la
libSherpaMain_la_SOURCES      =	Sherpa.C 
noinst_HEADERS        =	Sherpa.H 
AM_CXXFLAGS   =	$(SHERPAFLAGS) $(CLHEPFLAGS)
INCLUDES      = $(SHERPAINCS) $(APACICINCS) $(AMEGICINCS) $(AMISICINCS) \
		$(EXTRAXSINCS) $(PHASICINCS) $(MODELINCS) $(PDFINCS) \
		$(BEAMINCS) $(ATOOLSINCS) $(CLHEPINCS)




