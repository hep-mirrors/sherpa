SHERPASYS     =	../..

SHERPAFLAGS   =	# -pedantic -Wall -Winline
DEFS          =

ANALYSISDIR   =	$(SHERPASYS)/ANALYSIS-1.0
AMEGICDIR     =	$(SHERPASYS)/AMEGIC++-2.0
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
PDFDIR        =	$(SHERPASYS)/PDF-1.0
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0

ANALYSISINCS  = -I$(ANALYSISDIR)/Main -I$(ANALYSISDIR)/Observables
AMEGICINCS    =	-I$(AMEGICDIR)/Main -I$(AMEGICDIR)/Amplitude/Zfunctions \
		-I$(AMEGICDIR)/Amplitude/AmplTools \
		-I$(AMEGICDIR)/Amplitude -I$(AMEGICDIR)/Phasespace \
		-I$(AMEGICDIR)/Model -I$(AMEGICDIR)/String 
APACICINCS    = -I$(APACICDIR)/Main -I$(APACICDIR)/Showers
EXTRAXSINCS   =	-I$(EXTRAXSDIR)/Two2Two -I$(EXTRAXSDIR)/Main
SHERPAINCS    =	-I../PerturbativePhysics -I../SoftPhysics \
		-I../LundTools -I../HerwigTools -I../Tools 
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFINCS       =	-I$(PDFDIR)/Main -I$(PDFDIR)/Remnant -I$(PDFDIR)/Sudakov
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org

if HADRONS_SUPPORT
HADRONSDIR    =	$(SHERPASYS)/HADRONS++-0.0
HADRONSINCS   = -I$(HADRONSDIR)/Main -I$(HADRONSDIR)/ME_Library \
		-I$(HADRONSDIR)/PS_Library
endif
if AMISIC_SUPPORT
AMISICDIR     =	$(SHERPASYS)/AMISIC++-1.0
AMISICINCS    =	-I$(SHERPASYS)/AMISIC++-1.0/Main
endif

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTFLAGS     = -Wno-long-long
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES	      =	libSherpaSingleEvents.la
libSherpaSingleEvents_la_SOURCES      =	\
					Event_Handler.C \
					Event_Phase_Handler.C \
				   	Signal_Processes.C \
					Multiple_Interactions.C \
					Jet_Evolution.C \
					Hard_Decays.C \
					Hadronization.C \
					Hadron_Decays.C \
					Analysis_Phase.C \
					EvtReadin_Phase.C \
					MC_Interface.C

include_HEADERS        =	\
			Event_Handler.H \
			Event_Phase_Handler.H \
		 	Signal_Processes.H \
			Multiple_Interactions.H \
			Jet_Evolution.H \
			Hard_Decays.H \
		 	Hadronization.H \
			Hadron_Decays.H \
			Analysis_Phase.H \
			EvtReadin_Phase.H \
			MC_Interface.H

EXTRA_DIST    = $(CVSINCLUDE) 
AM_CXXFLAGS   =	$(SHERPAFLAGS) $(ROOTFLAGS)
INCLUDES      = $(SHERPAINCS) $(ANALYSISINCS) $(APACICINCS) \
		$(AMEGICINCS) $(AMISICINCS) \
		$(EXTRAXSINCS) $(PHASICINCS) $(MODELINCS) \
		$(PDFINCS) $(BEAMINCS) $(HADRONSINCS) $(ATOOLSINCS) \
		$(ROOTINCS)




