SHERPASYS     =	../..

PHASICFLAGS   =	# -pedantic -Wall -Winline
DEFS          =

MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org

if CVS_INCLUDE
CVSINCLUDE    =	.svn
endif

if ROOT_SUPPORT
ROOTINCS      =	-I@ROOTINCDIR@
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

pkglib_LTLIBRARIES       =	libPhasespace.la
libPhasespace_la_SOURCES      =	\
				Regulator_Base.C \
				Regulator.C \
				Integrable_Base.C \
			   	Channel_Basics.C \
				Channel_Elements.C \
				Channel_Elements_KK.C \
				Single_Channel.C \
				Simple_Pole.C \
			   	Resonance.C \
				Threshold.C \
				Leading_Log.C \
				LBS_Compton_Peak.C \
			   	Leading_Log_Z.C \
			   	LL_KPerp.C \
				FSR_Channel.C \
				Rambo.C \
				RamboKK.C \
				Sarge.C \
				Vegas.C ISR_Vegas.C \
				Channel_Interface.C \
			   	Multi_Channel.C \
				PI_Interface.C \
				Phase_Space_Handler.C \
				Phase_Space_Integrator.C 
pkginclude_HEADERS        =	\
			Regulator_Base.H \
			Regulator.H \
		 	Channel_Basics.H \
			Channel_Elements.H \
			Channel_Elements_KK.H \
			Single_Channel.H \
		 	Simple_Pole.H \
		 	Resonance.H \
			Threshold.H \
			Leading_Log.H \
			LBS_Compton_Peak.H \
		 	Leading_Log_Z.H \
		 	LL_KPerp.H \
			FSR_Channel.H \
			Channel_Interface.H \
			Rambo.H \
			RamboKK.H \
			Sarge.H \
			Vegas.H ISR_Vegas.H \
			Multi_Channel.H \
		 	Integrable_Base.H \
			PI_Interface.H \
			Phase_Space_Handler.H \
			Phase_Space_Integrator.H 
AM_CXXFLAGS   =	$(PHASICFLAGS) $(ROOTFLAGS)
INCLUDES      =	$(PHASICINCS) $(MODELINCS) $(PDFINCS) $(BEAMINCS) \
		$(ATOOLSINCS) $(ROOTINCS)
EXTRA_DIST    = $(CVSINCLUDE)
