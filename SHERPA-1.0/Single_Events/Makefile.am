SHERPASYS     =	../..

SHERPAFLAGS   =	-pedantic -Wall -Winline

ANALYSISDIR   =	$(SHERPASYS)/ANALYSIS-1.0
ANALYSISINCS  = -I$(ANALYSISDIR)/Main -I$(ANALYSISDIR)/Observables
AMEGICDIR     =	$(SHERPASYS)/AMEGIC++-2.0
AMEGICINCS    =	-I$(AMEGICDIR)/Main -I$(AMEGICDIR)/Amplitude/Zfunctions \
		-I$(AMEGICDIR)/Amplitude/AmplTools -I$(AMEGICDIR)/Amplitude \
		-I$(AMEGICDIR)/Phasespace -I$(AMEGICDIR)/Model \
		-I$(AMEGICDIR)/String 
ADICICDIR     =	$(SHERPASYS)/ADICIC++-0.0
ADICICINCS    = -I$(ADICICDIR)/Main -I$(ADICICDIR)/Structure \
                -I$(ADICICDIR)/Handler -I$(ADICICDIR)/Tool -I$(ADICICDIR)/Input
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0
APACICINCS    = -I$(APACICDIR)/Main -I$(APACICDIR)/Showers
AMISICINCS    =	-I$(SHERPASYS)/AMISIC++-1.0/Main
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
EXTRAXSINCS   =	-I$(EXTRAXSDIR)/Two2Two -I$(EXTRAXSDIR)/Main
SHERPAINCS    =	-I../PerturbativePhysics -I../SoftPhysics \
		-I../LundTools -I../HerwigTools -I../Tools 
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFDIR        =	-I$(SHERPASYS)/PDF-1.0
PDFINCS       = $(PDFDIR)/Main $(PDFDIR)/Sudakov
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org

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
INCLUDES      = $(SHERPAINCS) $(ANALYSISINCS) $(APACICINCS) $(ADICICINCS) \
		$(AMEGICINCS) $(AMISICINCS) $(EXTRAXSINCS) $(PHASICINCS) \
		$(MODELINCS) $(PDFINCS) $(BEAMINCS) $(ATOOLSINCS) 




