SHERPASYS  = ../..

ATOOLSINC  = $(SHERPASYS)/ATOOLS-2.0
PDFINC     = $(SHERPASYS)/PDF-1.0
BEAMINC    = $(SHERPASYS)/BEAM-1.0
PHASICINC  = $(SHERPASYS)/PHASIC++-1.0
MODELINC   = $(SHERPASYS)/MODEL-1.0
EXTRAXSINC = $(SHERPASYS)/EXTRA_XS-1.0
AMEGICINC  = $(SHERPASYS)/AMEGIC++-2.0
APACICINC  = $(SHERPASYS)/APACIC++-2.0
AMISICINC  = $(SHERPASYS)/AMISIC++-1.0

lib_LTLIBRARIES	      =	libSherpaSingleEvents.la
libSherpaSingleEvents_la_SOURCES      =	\
					Event_Handler.C \
					Event_Phase_Handler.C \
				   	Signal_Processes.C \
					Multiple_Interactions.C \
					Jet_Evolution.C \
					Hard_Decays.C \
					Hadronization.C
#CXXFLAGS      =	-pedantic -Wall -Winline
INCLUDES      =	\
		-I../../SHERPA-1.0/Tools \
		-I../../SHERPA-1.0/PerturbativePhysics \
		-I../../SHERPA-1.0/SoftPhysics \
		-I../../SHERPA-1.0/SoftPhysics/Lund \
		-I$(AMEGICINC)/Main \
		-I$(AMEGICINC)/Amplitude \
		-I$(AMEGICINC)/Phasespace \
		-I$(AMEGICINC)/Amplitude/AmplTools \
		-I$(AMEGICINC)/Amplitude/Zfunctions \
		-I$(AMEGICINC)/Phasespace \
		-I$(AMEGICINC)/String \
		-I$(AMEGICINC)/Model \
		-I$(APACICINC)/Main \
		-I$(APACICINC)/Showers \
		-I$(AMISICINC)/Main \
		-I$(AMISICINC)/Other \
		-I$(AMISICINC)/Tools \
		-I$(EXTRAXSINC)/Main \
		-I$(EXTRAXSINC)/Two2Two \
		-I$(MODELINC)/Main \
		-I$(PHASICINC)/Main \
		-I$(BEAMINC)/Main \
		-I$(PDFINC)/LHAPDF \
		-I$(PDFINC)/Main \
		-I$(ATOOLSINC)/Phys \
		-I$(ATOOLSINC)/Math \
		-I$(ATOOLSINC)/Org 
noinst_HEADERS        =	\
			Event_Handler.H \
			Event_Phase_Handler.H \
		 	Signal_Processes.H \
			Multiple_Interactions.H \
			Jet_Evolution.H \
			Hard_Decays.H \
		 	Hadronization.H




