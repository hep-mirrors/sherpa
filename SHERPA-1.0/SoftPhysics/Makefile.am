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

SUBDIRS       =	Lund
lib_LTLIBRARIES       =	libSherpaSoftPhysics.la
libSherpaSoftPhysics_la_SOURCES	      =	\
					Primordial_KPerp.C \
					Remnant_Base.C \
					Hadron_Remnant.C \
					Electron_Remnant.C \
					Photon_Remnant.C \
					Beam_Remnant_Handler.C \
				   	Fragmentation_Handler.C \
					Hadron_Decay_Handler.C
AM_CXXFLAGS   =	-pedantic -Wall -Winline
INCLUDES      =	\
		-I../../SHERPA-1.0/Tools \
		-I../../SHERPA-1.0/SoftPhysics/Lund \
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
			Primordial_KPerp.H \
			Primordial_KPerp.IC \
			Remnant_Base.H \
			Hadron_Remnant.H \
			Electron_Remnant.H \
			Photon_Remnant.H \
			Beam_Remnant_Handler.H \
		 	Fragmentation_Handler.H \
			Hadron_Decay_Handler.H
