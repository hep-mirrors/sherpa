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

lib_LTLIBRARIES       =	libSherpaPerturbativePhysics.la
libSherpaPerturbativePhysics_la_SOURCES       =	\
						Matrix_Element_Handler.C \
						MI_Handler.C \
					   	Shower_Handler.C \
						Hard_Decay_Handler.C \
						Interface_Tools.C \
					   	Cluster_Partons.C \
					   	SimpleXS_Apacic_Interface.C \
						Amegic_Apacic_Interface.C \
						Combine_Table.C \
						NLL_Sudakov.C
AM_CXXFLAGS   =	-pedantic -Wall -Winline
INCLUDES      =	\
		-I../../SHERPA-1.0/Tools \
		-I../../SHERPA-1.0/LundTools \
		-I../../SHERPA-1.0/Single_Events \
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
			Matrix_Element_Handler.H \
			MI_Handler.H \
		 	Shower_Handler.H \
			Hard_Decay_Handler.H \
		 	Perturbative_Interface.H \
		 	Interface_Tools.H \
			Cluster_Partons.H \
			SimpleXS_Apacic_Interface.H \
			Amegic_Apacic_Interface.H \
			Combine_Table.H \
			NLL_Sudakov.H



