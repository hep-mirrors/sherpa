SHERPASYS     =	../..

SHERPAFLAGS   =	-pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
AMEGICDIR     =	$(SHERPASYS)/AMEGIC++-2.0
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0
SHERPAINCS    =	-I../SingleEvents -I../Tools 
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
						NLL_Sudakov.C \
						NLL_Branching_Probability_Base.C \
						NLL_Branching_Probabilities.C \
						NLL_Single_Sudakov.C \
						NLL_Combined_Sudakov.C \
						NLL_JetRate.C
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
			NLL_Sudakov.H \
			NLL_Branching_Probabilities.H \
			NLL_Branching_Probability_Base.H \
			NLL_Sudakov_Base.H \
			NLL_Single_Sudakov.H \
			NLL_Combined_Sudakov.H \
			NLL_JetRate.H
AM_CXXFLAGS   =	$(SHERPAFLAGS)
INCLUDES      = $(SHERPAINCS) $(APACICINCS) $(AMEGICINCS) $(AMISICINCS) \
		$(EXTRAXSINCS) $(PHASICINCS) $(MODELINCS) $(PDFINCS) \
		$(BEAMINCS) $(ATOOLSINCS)


