SHERPASYS     =	../..

SHERPAFLAGS   =	-pedantic -Wall -Winline

ANALYSISDIR   =	$(SHERPASYS)/ANALYSIS-1.0
ANALYSISINCS  = -I$(ANALYSISDIR)/Main -I$(ANALYSISDIR)/Observables
AMEGICDIR     =	$(SHERPASYS)/AMEGIC++-2.0
AMEGICINCS    =	-I$(AMEGICDIR)/Main -I$(AMEGICDIR)/Amplitude/Zfunctions \
		-I$(AMEGICDIR)/Amplitude/AmplTools -I$(AMEGICDIR)/Amplitude \
		-I$(AMEGICDIR)/Phasespace -I$(AMEGICDIR)/Model \
		-I$(AMEGICDIR)/String 
if ADICIC_SUPPORT
ADICICDIR     = $(SHERPASYS)/ADICIC++-0.0
ADICICINCS    = -I$(ADICICDIR)/Main -I$(ADICICDIR)/Structure \
                -I$(ADICICDIR)/Handler -I$(ADICICDIR)/Tool -I$(ADICICDIR)/Input
ADICICFLAGS   = -DUSING__Adicic
endif
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0
APACICINCS    = -I$(APACICDIR)/Main -I$(APACICDIR)/Showers
AMISICINCS    =	-I$(SHERPASYS)/AMISIC++-1.0/Main
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
EXTRAXSINCS   =	-I$(EXTRAXSDIR)/Two2Two -I$(EXTRAXSDIR)/Main
SHERPAINCS    =	-I../PerturbativePhysics -I../SoftPhysics \
		-I../LundTools -I../HerwigTools -I../Tools 
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFDIR        =	$(SHERPASYS)/PDF-1.0
PDFINCS       =	-I$(PDFDIR)/Main -I$(PDFDIR)/Remnant -I$(PDFDIR)/Sudakov
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
AHADICDIR     =	$(SHERPASYS)/AHADIC++-0.0
AHADICINCS    =	-I$(AHADICDIR)/Main -I$(AHADICDIR)/ClusterFormation \
		-I$(AHADICDIR)/ClusterDecays
HADRONSDIR     = $(SHERPASYS)/HADRONS++-0.0
HADRONSINCS    = -I$(HADRONSDIR)/Main -I$(HADRONSDIR)/PS_Library -I$(HADRONSDIR)/ME_Library
	-I$(HADRONSDIR)/ME_Library -I$(HADRONSDIR)/PS_Library \
HADRONSFLAGS   = -DUSING__Hadrons


if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if MCATNLO_SUPPORT
MCATNLOFLAGS  = -DUSING__MCatNLO
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

noinst_HEADERS        =	\
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
AM_CXXFLAGS   =	$(SHERPAFLAGS) $(ADICICFLAGS) $(MCATNLOFLAGS) $(ROOTFLAGS)
INCLUDES      = $(SHERPAINCS) $(ANALYSISINCS) $(APACICINCS) $(ADICICINCS) \
		$(AMEGICINCS) $(AMISICINCS) $(AHADICINCS) $(EXTRAXSINCS) \
		$(PHASICINCS) $(MODELINCS) $(PDFINCS) $(BEAMINCS) $(HADRONSINCS) \
		$(ATOOLSINCS) $(ROOTINCS)




