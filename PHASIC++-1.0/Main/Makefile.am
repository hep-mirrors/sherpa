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
				Integration_Info.C \
			   	Channel_Basics.C \
				Channel_Elements.C \
				Single_Channel.C \
			   	Multi_Channel.C \
			   	ISR_Channel.C \
				FSR_Channel.C \
				Rambo.C \
				RamboKK.C \
				Sarge.C \
				Phase_Space_Handler.C \
				Phase_Space_Integrator.C 
noinst_HEADERS        =	\
		 	Integrable_Base.H \
			Integration_Info.H \
			Integration_Info.IC \
		 	Channel_Basics.H \
			Channel_Elements.H \
			Single_Channel.H \
			Multi_Channel.H \
		 	ISR_Channel.H \
			FSR_Channel.H \
			Rambo.H \
			RamboKK.H \
			Sarge.H \
			Phase_Space_Handler.H \
			Phase_Space_Integrator.H 
AM_CXXFLAGS   =	$(PHASICFLAGS)
EXTRA_DIST    = $(CVSINCLUDE)
INCLUDES      =	$(BEAMINCS) $(PDFINCS) $(MODELINCS) $(ATOOLSINCS)
