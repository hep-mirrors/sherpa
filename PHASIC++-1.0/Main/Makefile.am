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
				Info_Key.C \
				Integration_Info.C \
			   	Channel_Basics.C \
				Channel_Elements.C \
				Single_Channel.C \
			 	ISR_Channel.C \
				FSR_Channel.C \
				Rambo.C \
				RamboKK.C \
				Sarge.C \
				Channel_Interface.C \
			   	Multi_Channel.C \
				Phase_Space_Handler.C \
				Phase_Space_Integrator.C 
#				Simple_Pole.C \
#			   	Resonance.C \
#				Threshold.C \
#				Leading_Log.C \
#				LBS_Compton_Peak.C \
#			   	LL_KPerp.C \
noinst_HEADERS        =	\
			Info_Key.H \
			Info_Key.inl.H \
			Integration_Info.H \
			Integration_Info.inl.H \
		 	Channel_Basics.H \
			Channel_Elements.H \
			Single_Channel.H \
		 	ISR_Channel.H \
			FSR_Channel.H \
			Channel_Interface.H \
			Rambo.H \
			RamboKK.H \
			Sarge.H \
			Multi_Channel.H \
		 	Integrable_Base.H \
			Phase_Space_Handler.H \
			Phase_Space_Integrator.H 
#		 	Simple_Pole.H \
#		 	Resonance.H \
#			Threshold.H \
#			Leading_Log.H \
#			LBS_Compton_Peak.H \
#		 	LL_KPerp.H \
AM_CXXFLAGS   =	$(PHASICFLAGS)
EXTRA_DIST    = $(CVSINCLUDE)
INCLUDES      =	$(BEAMINCS) $(PDFINCS) $(MODELINCS) $(ATOOLSINCS)
