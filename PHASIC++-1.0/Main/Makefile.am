SHERPASYS     =	../..

PHASICFLAGS   =	-pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       =	libPhasespace.la
libPhasespace_la_SOURCES      =	\
				Integrable_Base.C \
			   	Channel_Basics.C \
				Channel_Elements.C \
				Single_Channel.C \
				Simple_Pole.C \
			   	Resonance.C \
				Threshold.C \
				Leading_Log.C \
				LBS_Compton_Peak.C \
			   	Leading_Log_Z.C \
			   	LL_KPerp.C \
			   	LDL_KPerp.C \
				FSR_Channel.C \
				Rambo.C \
				RamboKK.C \
				Sarge.C \
				Vegas.C ISR_Vegas.C \
				Channel_Interface.C \
			   	Multi_Channel.C \
				Phase_Space_Handler.C \
				Phase_Space_Integrator.C 
noinst_HEADERS        =	\
		 	Channel_Basics.H \
			Channel_Elements.H \
			Single_Channel.H \
		 	Simple_Pole.H \
		 	Resonance.H \
			Threshold.H \
			Leading_Log.H \
			LBS_Compton_Peak.H \
		 	Leading_Log_Z.H \
		 	LL_KPerp.H \
		 	LDL_KPerp.H \
			FSR_Channel.H \
			Channel_Interface.H \
			Rambo.H \
			RamboKK.H \
			Sarge.H \
			Vegas.H ISR_Vegas.H \
			Multi_Channel.H \
		 	Integrable_Base.H \
			Phase_Space_Handler.H \
			Phase_Space_Integrator.H 
AM_CXXFLAGS   =	$(PHASICFLAGS)
EXTRA_DIST    = $(CVSINCLUDE)
INCLUDES      =	$(BEAMINCS) $(PDFINCS) $(MODELINCS) $(ATOOLSINCS)
