SHERPASYS     =	../..

SHERPAFLAGS   =	-pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ANALYSISDIR   =	$(SHERPASYS)/ANALYSIS-1.0
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
AMEGICDIR     =	$(SHERPASYS)/AMEGIC++-2.0
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0
ADICICDIR     =	$(SHERPASYS)/ADICIC++-0.0
SHERPAINCS    =	-I../PerturbativePhysics -I../SoftPhysics \
		-I../LundTools -I../HerwigTools -I../Tools 
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
ANALYSISINCS  = -I$(ANALYSISDIR)/Main -I$(ANALYSISDIR)/Observables
PDFDIR        =	-I$(SHERPASYS)/PDF-1.0
PDFINCS       = $(PDFDIR)/Main $(PDFDIR)/Sudakov
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
ADICICINCS    = -I$(ADICICDIR)/Main -I$(ADICICDIR)/Structure \
		-I$(ADICICDIR)/Handler -I$(ADICICDIR)/Tool

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
					Analysis_Phase.C \
					EvtReadin_Phase.C \
					MC_Interface.C

noinst_HEADERS        =	\
			Event_Handler.H \
			Event_Phase_Handler.H \
		 	Signal_Processes.H \
			Multiple_Interactions.H \
			Jet_Evolution.H \
			Hard_Decays.H \
		 	Hadronization.H \
			Analysis_Phase.H \
			EvtReadin_Phase.H \
			MC_Interface.H

EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(SHERPAFLAGS)
INCLUDES      = $(SHERPAINCS) $(APACICINCS) $(ADICICINCS) $(AMEGICINCS) \
		$(AMISICINCS) $(EXTRAXSINCS) $(PHASICINCS) $(MODELINCS) \
		$(PDFINCS) $(BEAMINCS) $(ATOOLSINCS) $(ANALYSISINCS) 




